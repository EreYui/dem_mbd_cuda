#pragma once
#include "matrixtrans.h"



class ACC {
public:
    int Num;
    vector3d* acc;
    vector3d* domg;
        explicit ACC(int size = 0) : Num(size) {
        try {
            if (size < 0) throw std::invalid_argument("Negative size");
            
            if (size > 0) {
                acc = new vector3d[Num];
                domg = new vector3d[Num];
            }
        }
        catch (...) {
            clearResources();
            throw;
        }
    }
    ACC(const ACC&) = delete;
    ACC& operator=(const ACC&) = delete;

    // Move constructor
    ACC(ACC&& other) noexcept : Num(0), acc(nullptr), domg(nullptr) {
        *this = std::move(other);
    }

    ACC& operator=(ACC&& other) noexcept {
        if (this != &other) {
            clearResources();
            std::swap(Num, other.Num);
            std::swap(acc, other.acc);
            std::swap(domg, other.domg);
        }
        return *this;
    }

    ~ACC() {
        clearResources();
    }

private:
    void clearResources() noexcept {
        delete[] acc;
        delete[] domg;

        acc = nullptr;
        domg = nullptr;
        Num = 0;
    }
};



class VHALF {
public:
    int Num;
    vector3d* Vel = nullptr;
    vector3d* AngSpd = nullptr;
    vector4d* dQuat = nullptr;
    vector4d* quat = nullptr;

    explicit VHALF(int size = 0) : Num(size) {
        try {
            if (size < 0) throw std::invalid_argument("Negative size");
            
            if (size > 0) {
                AngSpd = new vector3d[Num];
                Vel    = new vector3d[Num];
                dQuat  = new vector4d[Num];
                quat   = new vector4d[Num];
            }
        }
        catch (...) {
            clearResources();
            throw;
        }
    }

    VHALF(const VHALF&) = delete;

    VHALF& operator=(const VHALF&) = delete;

    // Move constructor
    VHALF(VHALF&& other) noexcept : Num(0), Vel(nullptr), AngSpd(nullptr), dQuat(nullptr), quat(nullptr) {
        *this = std::move(other);
    }

    VHALF& operator=(VHALF&& other) noexcept {
        if (this != &other) {
            clearResources();
            std::swap(Num, other.Num);
            std::swap(Vel, other.Vel);
            std::swap(AngSpd, other.AngSpd);
            std::swap(dQuat, other.dQuat);
            std::swap(quat, other.quat);
        }
        return *this;
    }

    ~VHALF() {
        clearResources();
    }

private:
    void clearResources() noexcept {
        delete[] AngSpd;
        delete[] Vel;
        delete[] dQuat;
        delete[] quat;

        AngSpd = nullptr;
        Vel = nullptr;
        dQuat = nullptr;
        quat = nullptr;
        Num = 0;
    }
};