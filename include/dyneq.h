//  dyneq.h
//  Created by zkm on 2024.12.25.
//

#pragma once

#include"multibody.h"
#include"particle.h"

#include "matrixtrans.h" 

void DynEqnParticle(int id, PARTICLE& pt, double* FR, vector3d& acc, vector3d& domg);

//void DynMotionBody(double dt, BODYSET& bodyset);
//void DynMotionBody(double stepSize, int step, BODYSET& bodyset, double** oldF, double** F);

void DynEqnBodySet(double time, BODYSET& bodyset, double** F, double** ddu);
void DynEqnBodySet(double time, BODYSET& bodyset, double** F, double** ddu, int step);
