//  grid.h
//  Created by zkm on 2024.12.13.
//  Spatial grid utilities for neighbor search.

#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include "matrixtrans.h"
#include <vector>
#include <unordered_map>
#include <tuple>


#include "multibody.h"
#include "particle.h"


using namespace std;

using GridIndex = vector3i;

// Hash for GridIndex.
namespace std {
    template <>
    struct hash<GridIndex> {
        using result_type = std::size_t;

        result_type operator()(const GridIndex& index) const {
            size_t h1 = hash<int>{}(index[0]);
            size_t h2 = hash<int>{}(index[1]);
            size_t h3 = hash<int>{}(index[2]);
            h1 ^= h2 + 0x9e3779b9 + (h1 << 6) + (h1 >> 2);
            h1 ^= h3 + 0x9e3779b9 + (h1 << 6) + (h1 >> 2);
            return h1;
        }
    };
    // ? equal_to
    template<>
    struct equal_to<GridIndex> {
        bool operator()(const GridIndex& a, const GridIndex& b) const {
            return (a[0] == b[0]) && (a[1] == b[1]) && (a[2] == b[2]);
        }
    };
}


// 
class GridCell {
public:
    int n = 0; // cell id (legacy)

    vector3i id;
    
    int particleCount = 0;
    std::vector<int> particleIds;

    int triCount = 0;
    std::vector<int> bodyIds; // body ids
    std::vector<int> triIds;

    //int numNehb = 0;//?
    //int Nehb[26]; // neighbor cells (unused)

    GridCell() {
        particleIds.reserve(9);// ?9
        bodyIds.reserve(10);
        triIds.reserve(10);// ?10
    }

    void addParticle(int id) {
        particleIds.emplace_back(id);
        particleCount++;
    }

    void addTri(int bodyId, int triId) {
        bodyIds.emplace_back(bodyId);
        triIds.emplace_back(triId);
        triCount++;
    }
    void clear() {
        particleIds.clear();
        bodyIds.clear();
        triIds.clear();
        //particleIds.shrink_to_fit();
        //bodyIds.shrink_to_fit();
        //triIds.shrink_to_fit();
    }
};


typedef struct contactPair {
    int idx, flag;// particle index (dynamic)
    int body;
    double deltP[3][3];//deltS,deltR,deltT;

    struct contactPair* next;
}Cp;

bool InitList(Cp* part);
void BackInsertList(Cp* part, int id, double P[][3], int row = 3);

void CreatTree(PARTICLE &pt, BODYSET &bodyset, unordered_map<GridIndex, GridCell> &grid);
void CreatTree(PARTICLE& pt, unordered_map<GridIndex, GridCell>& grid);

void checkTree(PARTICLE& pt, BODYSET& bodyset, unordered_map<GridIndex, GridCell>& grid);
