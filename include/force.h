//  force.h
//  Created by zkm on 2024.12.15.
//
// ?

#pragma once

#include "multibody.h"
#include "particle.h"
#include "wall.h"
#include "grid.h"
#include "state.h"

struct CONTROL;
struct ODE;
struct MECH;

#include "matrixtrans.h" 

void ZeroForceStack(double** F, int n);

class FORCE {
public:
    int ptNum = 0;            // 
    int bodyNum = 0;          // 

     // Fc contact force; Ff field force; Fd drag force; Fw wall force; FR resultant force;FI STL force, EphemList STL;
    double** Fc; //[ptNum][6]??
    double** Ff; //[ptNum][6]
    double** Fd; //[ptNum][6] drag force
    double** Fw; //[ptNum][6] wall force
    double** FR; //[ptNum][6] resultant force

    double*** FI;    //[bodyNum][ptNum][6] body -> particle force
    double*** FI_pt; //[bodyNum][ptNum][6] particle -> body force

    double** poly_F; //[bodyNum][6] body resultant force

    double*** ForceBody; //[bodyNum][ptNum][19] per-body per-particle force terms
    // Detailed per-contact forces (legacy).
    int bodyContactCount = 0;
    double bodyContactMaxDepth = 0.0;
    double bodyContactFz = 0.0;
    
    // Initialize storage.
    FORCE(int nb = 0, int nbody = 0) : ptNum(nb), bodyNum(nbody) {
        Fc = new double* [nb];
        for (int j = 0; j < nb; j++) Fc[j] = new double[6];//??
        Ff = new double* [nb];
        for (int j = 0; j < nb; j++) Ff[j] = new double[6];//
        Fd = new double* [nb];
        for (int j = 0; j < nb; j++) Fd[j] = new double[6];
        Fw = new double* [nb];
        for (int j = 0; j < nb; j++) Fw[j] = new double[6];
        FR = new double* [nb];
        for (int j = 0; j < nb; j++) FR[j] = new double[6];//

        FI_pt = new double** [nbody];
        for (int ib = 0; ib < nbody; ib++) {
            FI_pt[ib] = new double* [nb];
            for (int j = 0; j < nb; j++) FI_pt[ib][j] = new double[6];
        }
        FI = new double** [nbody];
        for (int ib = 0; ib < nbody; ib++) {
            FI[ib] = new double* [nb];
            for (int j = 0; j < nb; j++) FI[ib][j] = new double[6];
        }
        ForceBody = new double** [nbody];
        for (int ib = 0; ib < nbody; ib++) {
            ForceBody[ib] = new double* [nb];
            for (int j = 0; j < nb; j++) ForceBody[ib][j] = new double[19];
        }

        poly_F = new double* [nbody];
        for (int i = 0; i < nbody; i++) poly_F[i] = new double[6];
    }

    void ForceCaculate(double time, PARTICLE& pt, BODYSET& bodyset, CONTROL& CT,
        MECH& mechPP, MECH& mechPT, MECH& mechPW, ODE& ode,
        Cp** part, Cp** poly_part, unordered_map<GridIndex, GridCell>& grid, VHALF& vhalf_pt, VHALF& vhalf_bd);

    void ForceCaculate(double time, PARTICLE& pt, CONTROL& CT,
        MECH& mechPP, MECH& mechPW, ODE& ode,
        Cp** part, Cp** poly_part, unordered_map<GridIndex, GridCell>& grid);

    void ParticlesForceOutput(int i, int nb, int nbody);
    void ParticlesForceOutput(int i, int nb);

    void BodysetForceOutput(int i, int nb, int nbody);

    void setZero()
    {
        {
            ZeroForceStack(Fw, ptNum);
            ZeroForceStack(Ff, ptNum);
            ZeroForceStack(Fc, ptNum);
            ZeroForceStack(Fd, ptNum);
            for (int ib = 0; ib < bodyNum; ib++)
                ZeroForceStack(FI[ib], ptNum);
            ZeroForceStack(FR, ptNum);
            for (int ib = 0; ib < bodyNum; ib++)
                ZeroForceStack(FI_pt[ib], ptNum);
        }

    for (int ib = 0; ib < bodyNum; ib++)
    {
        for (int i = 0; i < ptNum; i++) {
            for (int j = 0; j < 19; j++) {
                ForceBody[ib][i][j] = 0.0;
            }
        }
    }
    for (int ib = 0; ib < bodyNum; ib++) {
        for (int j = 0; j < 6; j++) {
            poly_F[ib][j] = 0.0;
        }
    }
    }

    ~FORCE() {
        for (int j = 0; j < ptNum; j++) delete[]Fc[j];
        delete[]Fc;
        for (int j = 0; j < ptNum; j++) delete[]Ff[j];
        delete[]Ff;
        for (int j = 0; j < ptNum; j++) delete[]Fd[j];
        delete[]Fd;
        for (int j = 0; j < ptNum; j++) delete[]Fw[j];
        delete[]Fw;
        for (int j = 0; j < ptNum; j++) delete[]FR[j];
        delete[]FR;
        for (int ib = 0; ib < bodyNum; ib++)
        {
            for (int j = 0; j < ptNum; j++)delete[] FI[ib][j];
            delete[]FI[ib];
        }
        delete[]FI;
        for (int ib = 0; ib < bodyNum; ib++)
        {
            for (int j = 0; j < ptNum; j++)delete[] FI_pt[ib][j];
            delete[]FI_pt[ib];
        }
        delete[]FI_pt;
        for (int i = 0; i < bodyNum; i++)delete[]poly_F[i];
        delete[] poly_F;
        for (int ib = 0; ib < bodyNum; ib++)
        {
            for (int j = 0; j < ptNum; j++)delete[] ForceBody[ib][j];
            delete[]ForceBody[ib];
        }
        delete[]ForceBody;
    }
};




void Gravity(int id1, int id2, PARTICLE& pt, vector3d &Force1);//
int WallForce(int id1, int id2, WALL& wall, PARTICLE& pt, MECH& mech, vector3d& Force, vector3d& Torque, ODE ode,Cp* part);
int ContactForce(int id1, int id2, PARTICLE& pt, MECH& mech, vector3d &Force1, vector3d &Torque1, vector3d &Torque2, ODE ode, Cp* part);

void Coefficient(int id1, int id2, PARTICLE& pt, MECH& mech);

int Polyhedron_Force(int id1, int id2, PARTICLE& pt, MECH& mech, vector3d& Force1, vector3d& Torque1, vector3d& Torque2,
double* ForceBodyTemp, ODE ode, Cp* part, BODYSET& bodyset, int idbody,
VHALF& vhalf_pt, VHALF& vhalf_bd);

//void ForceCaculate(double time, PARTICLE& pt, BODYSET& bodyset, CONTROL& CT, MECH& mech, ODE& ode, double** Fc, double** Ff, double** Fd, double** Fw, double** FR,
//double*** FI, double*** FI_pt, double*** ForceBody, double** poly_F, Cp** part, Cp** poly_part, unordered_map<GridIndex, GridCell>& grid,
//VHALF vhalf_pt, VHALF vhalf_bd);
//
//void ForceCaculate(double time, PARTICLE& pt, CONTROL& CT, MECH& mech, ODE& ode, double** Fc, double** Ff, double** Fd, double** Fw, double** FR,
//Cp** part, unordered_map<GridIndex, GridCell>& grid);

bool isPointInsideTriangle(const vector3d& A, const vector3d& B, const vector3d& C, const vector3d& Q);
vector3d closestPointOnSegment(const vector3d& A, const vector3d& B, const vector3d& P);
vector3d computeClosestPoint(const vector3d& A, const vector3d& B, const vector3d& C, const vector3d& n, const vector3d& P);
double checkContact(const vector3d& A, const vector3d& B, const vector3d& C, const vector3d& n, const vector3d& P, double radius, vector3d& closestPoint);

vector3d FrictionModelOfAdams(const vector3d& A, const MECH& mech);
