//  matrixtrans.cpp
//  Created by zkm on 2024.12.15.


#include"matrixtrans.h"

matrix3d Quat2IvDCM(vector4d q) // inverse DCM from quaternion
{
    matrix3d m;
    double qq[10];
    //
    qq[0] = q[0] * q[0];
    qq[1] = q[1] * q[1];
    qq[2] = q[2] * q[2];
    qq[3] = q[3] * q[3];
    qq[4] = 2.0 * q[0] * q[1];
    qq[5] = 2.0 * q[0] * q[2];
    qq[6] = 2.0 * q[0] * q[3];
    qq[7] = 2.0 * q[1] * q[2];
    qq[8] = 2.0 * q[1] * q[3];
    qq[9] = 2.0 * q[2] * q[3];
    //
    m[0][0] = qq[0] + qq[1] - qq[2] - qq[3];
    m[0][1] = qq[7] - qq[6];
    m[0][2] = qq[8] + qq[5];
    m[1][0] = qq[7] + qq[6];
    m[1][1] = qq[0] - qq[1] + qq[2] - qq[3];
    m[1][2] = qq[9] - qq[4];
    m[2][0] = qq[8] - qq[5];
    m[2][1] = qq[9] + qq[4];
    m[2][2] = qq[0] - qq[1] - qq[2] + qq[3];

    return m;
}

matrix3d Quat2DCM(vector4d q) // DCM from quaternion
{
    matrix3d m;
    double qq[10];
    //
    qq[0] = q[0] * q[0];
    qq[1] = q[1] * q[1];
    qq[2] = q[2] * q[2];
    qq[3] = q[3] * q[3];
    qq[4] = 2.0 * q[0] * q[1];
    qq[5] = 2.0 * q[0] * q[2];
    qq[6] = 2.0 * q[0] * q[3];
    qq[7] = 2.0 * q[1] * q[2];
    qq[8] = 2.0 * q[1] * q[3];
    qq[9] = 2.0 * q[2] * q[3];
    //
    m[0][0] = qq[0] + qq[1] - qq[2] - qq[3];
    m[0][1] = qq[7] + qq[6];
    m[0][2] = qq[8] - qq[5];
    m[1][0] = qq[7] - qq[6];
    m[1][1] = qq[0] - qq[1] + qq[2] - qq[3];
    m[1][2] = qq[9] + qq[4];
    m[2][0] = qq[8] + qq[5];
    m[2][1] = qq[9] - qq[4];
    m[2][2] = qq[0] - qq[1] - qq[2] + qq[3];

    return m;
}

void QuatNorm(vector4d& q)
{
    double mag = sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
    //assert(mag > 0.0);
    for (int i = 0; i < 4; i++) q[i] = q[i] / mag;
}

long long IJK2N(int I, int J, int K, int D0Sz, int D1Sz, int D2Sz)
{
    return (long long)(I * D1Sz + J) * D2Sz + K;//using "long long int" type in Linux 
}

void N2IJK(int N, int IJK[3], int D0Sz, int D1Sz, int D2Sz)
{
    int tmpI;

    IJK[2] = N % D2Sz;
    tmpI = (N - IJK[2]) / D2Sz;
    IJK[1] = tmpI % D1Sz;
    IJK[0] = (tmpI - IJK[1]) / D1Sz;
}

vector4d gama(vector4d orien, vector3d omg) // quaternion derivative from angular velocity
{
    vector4d dorien;

    dorien[0] = 0.5 * (-orien[1] * omg[0] - orien[2] * omg[1] - orien[3] * omg[2]);
    dorien[1] = 0.5 * (orien[0] * omg[0] + orien[2] * omg[2] - orien[3] * omg[1]);
    dorien[2] = 0.5 * (orien[0] * omg[1] - orien[1] * omg[2] + orien[3] * omg[0]);
    dorien[3] = 0.5 * (orien[0] * omg[2] + orien[1] * omg[1] - orien[2] * omg[0]);
    return dorien;
}
