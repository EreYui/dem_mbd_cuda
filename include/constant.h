//
//  
//
//  Created by Yang Yu on 14-11-11.
//
//

#ifndef _CONSTANT_H_
#define _CONSTANT_H_

// input & output const settings
const int WidthOutput = 8; // The output number width of file name
const int NumParameter = 15; // The number of input parameters from PAR file
const int WidthInt = 22; // The output width of integer data
const int WidthDouble = 22; // The output width of double data
const int PrecDouble = 12; // The output precession of double data

// solver constants
//const int EqnDimAsteroid = 7; // The dimension of binary equations
//const int EqnDimDebris = 6; // The dimension of debris equations
const double PrecM2E = 1.0e-15; // The iterational precession of Mean Anomaly to Eccentric Anomaly
//const int RebCrNum = 10; // The criteria number of rebounce, i.e., the number of timesteps for the lowest hop on the surface, below which
						 // the velocity is identified as Zero

// solar system physical constants
const double PI = 3.14159265359; // The value of PI
const double G = 6.67384e-11; // The gravitational constant, unit: m^3/kg/s^2
const double SM = 1.9891e30; // The solar mass, unit: kg
const double AU = 1.49597871e11; // The astronomical unit, unit: m
const double SRF = 1.361e3; // The solar radiation flux at 1 AU, unit: kg/s^3
const double LS = 2.988e8; // The light speed value, unit: m/s

// tree params
const int StackLmtLen = 1000000;
const int BranchLmtLen = 10;
const int NghbrUp[27][3] =
{
	{ 0,  0,  0, },//0
	{  1,  0,  0, },//1
	{ -1,  1,  0, },//2
	{  0,  1,  0, },//3
	{  1,  1,  0, },//4
	{ -1, -1,  1, },//5
	{  0, -1,  1, },//6
	{  1, -1,  1, },//7
	{ -1,  0,  1, },//8
	{  0,  0,  1, },//9
	{  1,  0,  1, },//10
	{ -1,  1,  1, },//11
	{  0,  1,  1, },//12
	{  1,  1,  1, },//13
	{ 1,  -1,  0, },//14
	{ -1,  0,  0, },//15
	{ 0,  -1,  0, },//16
	{ -1,  -1,  0, },//17
	{ -1, -1,  -1, },//18
	{ 0, -1,  -1, },//19
	{ 1, -1,  -1, },//20
	{ -1,  0,  -1, },//21
	{ 0,  0,  -1, },//22
	{ 1,  0,  -1, },//23
	{ -1,  1,  -1, },//24
	{ 0,  1,  -1, },//25
	{ 1,  1,  -1, }//26
};


//


#endif

