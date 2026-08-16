//  dyneq.cpp

#include"matrixtrans.h"
#include"multibody.h"
#include"dyneq.h"
#include"particle.h"


void DynEqnParticle(int id, PARTICLE& pt, double* FR, vector3d &acc, vector3d &domg)
{
    //
    vector3d F, T;
    matrix3d DCM;
    vector4d orien = pt.Quat[id];
    vector3d omg = pt.AngSpd[id];

    F = vector3d(FR[0], FR[1], FR[2]);
    T = vector3d(FR[3], FR[4], FR[5]);

    acc = F * (1.0 / pt.Mass[id]);
/*
    dorien[0] = 0.5 * (-orien[1] * omg[0] - orien[2] * omg[1] - orien[3] * omg[2]);
    dorien[1] = 0.5 * (orien[0] * omg[0] + orien[2] * omg[2] - orien[3] * omg[1]);
    dorien[2] = 0.5 * (orien[0] * omg[1] - orien[1] * omg[2] + orien[3] * omg[0]);
    dorien[3] = 0.5 * (orien[0] * omg[2] + orien[1] * omg[1] - orien[2] * omg[0]);
*/
    DCM = Quat2DCM(orien);
    T = DCM * T; // MC_id transformed to body-fixed frame

    domg = T * (1.0 / pt.Inertia[id]);
}

void DynMotionBody(double dt, BODYSET& bodyset)//?//
{
    for (int i = 0; i < bodyset.Num; i++)
    {
        bodyset.body[i].MassCenter = bodyset.body[i].MassCenter + bodyset.body[i].Vel * dt;
    }
}

void DynMotionBody(double stepSize, int step, BODYSET& bodyset, double** oldF, double** F)
{
    int n = bodyset.Num;
    double time = step * stepSize;
    double** ddu = new double* [n];
    for (int i = 0; i < n; i++)  ddu[i] = new double[6];

    double dt = stepSize / 2;
    DynEqnBodySet(time, bodyset, F, ddu); // Evaluate dynamics to get accelerations.

    vector3d Acc, AngularAcc;
    //Matrix IvDCM,DCM;
    vector3d tmpV, tmpP, tmpAcc, tmpAngularAcc, tmpAngularV, omg;
    vector4d dorien;
    // Update position, velocity, and orientation per body.
    for (int i = 0; i < n; i++)
    {
        Acc = vector3d(ddu[i][0], ddu[i][1], ddu[i][2]);
        tmpAcc = Acc * dt; // Half-step linear acceleration.
        tmpV = bodyset.body[i].Vel + tmpAcc; // Half-step velocity.
        tmpP = tmpV * stepSize; // Position delta.
        bodyset.body[i].MassCenter = bodyset.body[i].MassCenter + tmpP; // Position update.
        bodyset.body[i].Vel = tmpV + tmpAcc; // Velocity update.


                QuatNorm(bodyset.body[i].orien);
        matrix3d IvDCM = Quat2IvDCM(bodyset.body[i].orien);
        matrix3d DCM = Quat2DCM(bodyset.body[i].orien);

        AngularAcc = vector3d(ddu[i][3], ddu[i][4], ddu[i][5]); // Angular acceleration.
        tmpAngularAcc = AngularAcc * dt; // Half-step angular acceleration.
        tmpAngularV = bodyset.body[i].AngularVel + tmpAngularAcc; // Half-step angular velocity.

        omg = tmpAngularV;
        dorien[0] = 0.5 * (-bodyset.body[i].orien[1] * omg[0] - bodyset.body[i].orien[2] * omg[1] - bodyset.body[i].orien[3] * omg[2]);
        dorien[1] = 0.5 * (bodyset.body[i].orien[0] * omg[0] + bodyset.body[i].orien[2] * omg[2] - bodyset.body[i].orien[3] * omg[1]);
        dorien[2] = 0.5 * (bodyset.body[i].orien[0] * omg[1] - bodyset.body[i].orien[1] * omg[2] + bodyset.body[i].orien[3] * omg[0]);
        dorien[3] = 0.5 * (bodyset.body[i].orien[0] * omg[2] + bodyset.body[i].orien[1] * omg[1] - bodyset.body[i].orien[2] * omg[0]);

        bodyset.body[i].orien[0] += stepSize * dorien[0];
        bodyset.body[i].orien[1] += stepSize * dorien[1];
        bodyset.body[i].orien[2] += stepSize * dorien[2];
        bodyset.body[i].orien[3] += stepSize * dorien[3];
        QuatNorm(bodyset.body[i].orien);

        bodyset.body[i].AngularVel = tmpAngularV + tmpAngularAcc;
    }

    for (int j = 0; j < n; j++) delete[]ddu[j];
    delete[]ddu;
}

// Rigid body dynamics evaluation.
void DynEqnBodySet(double time, BODYSET& bodyset, double** F, double** ddu)
{
    int n = bodyset.Num;
    BODY x;
    vector3d T;
    matrix3d DCM;
    for (int i = 0; i < n; i++)
    {
        x = bodyset.body[i];
        if (x.state != BODY_DYNAMIC) {
            for (int j = 0; j < 6; ++j) ddu[i][j] = 0.0;
            continue;
        }
        for (int j = 0; j < 3; j++)
        {
            ddu[i][j] = F[i][j] / x.Mass;
        }
        T[0] = F[i][3]; T[1] = F[i][4]; T[2] = F[i][5];
        DCM = Quat2DCM(x.orien);
        T = DCM * T;

        ddu[i][3] = (T[0] - (x.I[2][2] - x.I[1][1]) * x.AngularVel[2] * x.AngularVel[1]) / x.I[0][0];
        ddu[i][4] = (T[1] - (x.I[0][0] - x.I[2][2]) * x.AngularVel[0] * x.AngularVel[2]) / x.I[1][1];
        ddu[i][5] = (T[2] - (x.I[1][1] - x.I[0][0]) * x.AngularVel[1] * x.AngularVel[0]) / x.I[2][2];
    }
// for (int i = 0; i < n; i++)
// {
// for (int j = 0; j < 6; j++)
// {
// ddu[i][j] = 0;
// }
// }
}


void DynEqnBodySet(double time, BODYSET& bodyset, double** F, double** ddu, int step)
{

}

