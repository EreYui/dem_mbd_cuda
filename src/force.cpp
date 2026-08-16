//  force.cpp
// 计算万有引力、颗粒之间作用力、颗粒与墙、颗粒与刚体

#include "force.h"
#include "contact_geometry.h"
#include "multibody.h"
#include "constant.h"
#include "grid.h"
#include "simulation.h"
#include <cassert>
#include <vector>
#include <array>
#include "matrixtrans.h"
//只有颗粒
void FORCE::ForceCaculate(double time, PARTICLE& pt, CONTROL& control,
	MECH& mechPP, MECH& mechPW, ODE& ode, Cp** part, Cp** poly_part,
	unordered_map<GridIndex, GridCell>& grid)
{
	setZero();

	vector3d Force1, Torque1, Torque2;
	for (int j = 0; j < pt.Num; j++)
	{
		if (control.Universal_gravitation_flag)//这个是万有引力
		{
			for (int k = 0; k < pt.Num; k++)
			{
				if (j == k)
					continue;
				Gravity(j, k, pt, Force1);
				Ff[j][0] = Ff[j][0] + Force1[0];
				Ff[j][1] = Ff[j][1] + Force1[1];
				Ff[j][2] = Ff[j][2] + Force1[2];
			}
		}
		if (control.Grav_flag)//这个是重力，或者说小行星对颗粒的引力
		{
			Ff[j][0] = Ff[j][0] + control.g[0] * pt.Mass[j];
			Ff[j][1] = Ff[j][1] + control.g[1] * pt.Mass[j];
			Ff[j][2] = Ff[j][2] + control.g[2] * pt.Mass[j];
		}
		if (control.Wall_flag)//墙壁给颗粒的力
		{
			for (int w = 0; w < control.wall.number; w++)
			{
				if (WallForce(j, w, control.wall, pt, mechPW, Force1, Torque1, ode, poly_part[j]))
				{
					for (int q = 0; q < 3; q++)
					{
						Fw[j][q] = Fw[j][q] + Force1[q];
						Fw[j][q + 3] = Fw[j][q + 3] + Torque1[q];
					}
				}
			}
		}
		for (int k = 0; k < 27; k++)//颗粒间的接触力
		{
			GridIndex index;
			index[0] = pt.Idx[j][0] + NghbrUp[k][0];
			index[1] = pt.Idx[j][1] + NghbrUp[k][1];
			index[2] = pt.Idx[j][2] + NghbrUp[k][2];

			if (grid.find(index) == grid.end())//网格不存在则跳过
				continue;
			auto it = grid.find(index);
			GridCell& cell = it->second;
			for (int s = 0; s < cell.particleIds.size(); s++) //遍历该网格中的所有颗粒
			{
				int ptid = cell.particleIds[s];//当前颗粒的编号
				if (ptid == j)
					continue;
				if (ContactForce(j, ptid, pt, mechPP, Force1, Torque1, Torque2, ode, part[j]))
				{
					Fc[j][0] = Fc[j][0] + Force1[0];
					Fc[j][1] = Fc[j][1] + Force1[1];
					Fc[j][2] = Fc[j][2] + Force1[2];
					Fc[j][3] = Fc[j][3] + Torque1[0];
					Fc[j][4] = Fc[j][4] + Torque1[1];
					Fc[j][5] = Fc[j][5] + Torque1[2];
				}
			}
		}//k<27
	}//颗粒数

	// calculate principal vectors and moments FR in World Frame (inertia)
	for (int j = 0; j < pt.Num; j++)
		for (int k = 0; k < 6; k++)
			FR[j][k] = Fc[j][k] + Ff[j][k] + Fd[j][k] + Fw[j][k];//计算每个颗粒受到的合力、和力矩
}


//有颗粒和刚体系
void FORCE::ForceCaculate(double time, PARTICLE& pt, BODYSET& bodyset, CONTROL& control,
	MECH& mechPP, MECH& mechPT, MECH& mechPW, ODE& ode,
		Cp** part, Cp** poly_part, unordered_map<GridIndex, GridCell>& grid, VHALF &vhalf_pt, VHALF& vhalf_bd)
{

	setZero();

	vector3d Force1, Torque1, Torque2;
	double ForceBodyTemp[19];
	for (int j = 0; j < pt.Num; j++)
	{
		if (control.Universal_gravitation_flag)//这个是万有引力
		{
			for (int k = 0; k < pt.Num; k++)
			{
				if (j == k)
					continue;
				Gravity(j, k, pt, Force1);
				Ff[j][0] = Ff[j][0] + Force1[0];
				Ff[j][1] = Ff[j][1] + Force1[1];
				Ff[j][2] = Ff[j][2] + Force1[2];
			}
		}
		if (control.Grav_flag)//这个是重力，或者说小行星对颗粒的引力
		{
			Ff[j][0] = Ff[j][0] + control.g[0] * pt.Mass[j];
			Ff[j][1] = Ff[j][1] + control.g[1] * pt.Mass[j];
			Ff[j][2] = Ff[j][2] + control.g[2] * pt.Mass[j];
		}
		if (control.Wall_flag)//墙壁给颗粒的力
		{
			for (int w = 0; w < control.wall.number; w++)
			{
				if (WallForce(j, w, control.wall, pt, mechPW, Force1, Torque1, ode, poly_part[j]))
				{
					for (int q = 0; q < 3; q++)
					{
						Fw[j][q] = Fw[j][q] + Force1[q];
						Fw[j][q + 3] = Fw[j][q + 3] + Torque1[q];
					}
				}
			}
		}

		//calculating contact force
		std::vector<std::array<int, 2>> triRecord; // 动态 n x 2 数组,记录已经计算的三角面
		for (int k = 0; k < 27; k++)//
		{
			GridIndex index;
			index[0] = pt.Idx[j][0] + NghbrUp[k][0];
			index[1] = pt.Idx[j][1] + NghbrUp[k][1];
			index[2] = pt.Idx[j][2] + NghbrUp[k][2];

			if (grid.find(index) == grid.end())//网格不存在则跳过
				continue;
			auto it = grid.find(index);
			GridCell& cell = it->second;
			for (int s = 0; s < cell.particleIds.size(); s++) //遍历该网格中的所有颗粒
			{
				int ptid = cell.particleIds[s];//当前颗粒的编号
				if (ptid == j)
					continue;
				if (ContactForce(j, ptid, pt, mechPP, Force1, Torque1, Torque2, ode, part[j]))
				{
					
					Fc[j][0] = Fc[j][0] + Force1[0];
					Fc[j][1] = Fc[j][1] + Force1[1];
					Fc[j][2] = Fc[j][2] + Force1[2];
					Fc[j][3] = Fc[j][3] + Torque1[0];
					Fc[j][4] = Fc[j][4] + Torque1[1];
					Fc[j][5] = Fc[j][5] + Torque1[2];
				}
			}
			if (cell.triIds.size() != cell.bodyIds.size()) {
				std::cout << "Body record count is not equal to triangle record count" << std::endl;
			}
			for (int s = 0; s < cell.triIds.size(); s++) //遍历该网格中的所有三角面
			{
				int triid = cell.triIds[s];
				int bodyid = cell.bodyIds[s]; 
				// 同一个三角面会同时存在于好几个网格中，需要预防重复计算
				// 判断该三角面的力是否已计算
				bool tri_flag = 0;//尚未计算则记为0
				for (int p = 0; p < triRecord.size(); p++)
				{
					if (triRecord[p][0] == bodyid && triRecord[p][1] == triid)
					{
						tri_flag = 1;
						break;//中止循环
					}
				}
				if (tri_flag == 0)//若还未计算，则进行计算，并记录该三角面
				{
					triRecord.push_back({bodyid,triid});
					if (Polyhedron_Force(j, triid, pt, mechPT, Force1, Torque1, Torque2, ForceBodyTemp, ode, poly_part[j], bodyset,bodyid,vhalf_pt, vhalf_bd))//表明相交
					{
						//printf("Fc = %8.4f,%8.4f,%8.4f,  Tc = %8.4f,%8.4f,%8.4f, pt = %d\n", Force1[0], Force1[1], Force1[2], Torque1[0], Torque1[1], Torque1[2], j);
						FI_pt[bodyid][j][0] += Force1[0];//颗粒受到的来自刚体的力
						FI_pt[bodyid][j][1] += Force1[1];//力的三个分量沿惯性系三个轴的方向
						FI_pt[bodyid][j][2] += Force1[2];
						FI_pt[bodyid][j][3] += Torque1[0];//颗粒受到的来自刚体的力矩
						FI_pt[bodyid][j][4] += Torque1[1];//力矩的三个分量沿惯性系的三个轴的方向
						FI_pt[bodyid][j][5] += Torque1[2];

						FI[bodyid][j][0] -= Force1[0];//刚体受到的来自颗粒的力
						FI[bodyid][j][1] -= Force1[1];//力的三个分量沿惯性系三个轴的方向
						FI[bodyid][j][2] -= Force1[2];
						FI[bodyid][j][3] += Torque2[0];//刚体受到的颗粒的作用力产生的力矩
						FI[bodyid][j][4] += Torque2[1];//力矩的三个分量沿惯性系的三个轴的方向
						FI[bodyid][j][5] += Torque2[2];

						for (int k_force = 0; k_force < 19; k_force++) {
							ForceBody[bodyid][j][k_force] = ForceBodyTemp[k_force];
						}
					}
				}
			}
		}//k<27
	}//颗粒数
	
	 //calculate principal vectors and moments FR in World Frame (inertia)---------------------------------------------------------------------------------------
	for (int ib = 0; ib < bodyset.Num; ib++)//计算每个刚体受到的合力、合力矩
	{
		for (int k = 0; k < 6; k++)
		{
			double tmp = 0;
			for (int j = 0; j < pt.Num; j++) {
				tmp = tmp + FI[ib][j][k];
			}
			poly_F[ib][k] = tmp;
		}
	}


	for (int j = 0; j < pt.Num; j++)//计算每个颗粒受到的合力、合力矩
	{
		for (int k = 0; k < 6; k++)
		{
			double tmp = 0;
			for (int ib = 0; ib < bodyset.Num; ib++) {
				tmp = tmp + FI_pt[ib][j][k];
			}
			FR[j][k] = Fc[j][k] + Ff[j][k] + Fd[j][k] + Fw[j][k] + tmp;
		}
	}

	if (control.Grav_flag)
	{
		for (int ib = 0; ib < bodyset.Num; ib++)
		{
			poly_F[ib][0] = poly_F[ib][0] + control.g[0] * bodyset.body[ib].Mass;//加入刚体受到的重力
			poly_F[ib][1] = poly_F[ib][1] + control.g[1] * bodyset.body[ib].Mass;
			poly_F[ib][2] = poly_F[ib][2] + control.g[2] * bodyset.body[ib].Mass;
		}
	}


}

void Gravity(int id1, int id2, PARTICLE& pt, vector3d &Force1)//万有引力
{
    double tmpd;
	vector3d r1, r2, r12;

    r1 = pt.Pos[id1];
    r2 = pt.Pos[id2];
    r12 = r1 - r2;

    tmpd = r12.norm();
    tmpd = G * pt.Mass[id1] * pt.Mass[id2] / (tmpd * tmpd * tmpd);//

    Force1 = r12 * tmpd;
}

int WallForce(int id1, int id2, WALL& wall, PARTICLE& pt, MECH& Mech, vector3d& Force, vector3d& Torque, ODE ode,Cp* part)
{
	int Cflag = 0, flag = 0;
	vector3d tmpF1, tmpF2;
	double tmpd1, tmpd2;

	double l_i, l_j;
	double delta;// 嵌入深度
	vector3d Vr, MT, M_FSj;
	
	vector3d pos1, vel1, omg1,pos2, vel2, omg2;
	vector3d velS, omgT, omgR;
	double dt = ode.StepSize;
	vector3d normal = - wall.N[id2];

	vector3d HR, HS, HT;//分别是滚转、切向、扭转的相对位移*K
	vector3d deltR, deltS, deltT;//分别是dt时间内滚转、切向、扭转的相对位移
	vector3d FN_elastic, FN_damping,FN_total;
	vector3d FrictionForce,FrictionForce_static,FrictionForce_static_elastic;
	vector3d M_FSi;//由摩擦力产生的对球心的力矩
	vector3d TwistTorque, TwistTorque_static,TwistTorque_static_elastic;
	vector3d RollTorque, RollTorque_static, RollTorque_static_elastic;

	MECH mech = Mech;

	Cp* pre = part;
	Cp* head = part;

	// 检索已有接触对中是否存在与id2墙面组成的接触对，若存在则更新part指针指向该接触对；若不存在则part指针指向链表末尾，为新接触对做准备
	while (part->next != nullptr)
	{
		if ((part->next->idx == id2) && (part->next->body == -2))//part->next->body==-2表明这个接触对是颗粒与墙的接触对
		{
			Cflag = 1;// 1 sphere-face contact;2 sphere-edge contact;3 sphere-vertex; 0 no contact
			part = part->next;
			break;
		}
		part = part->next;
	}

	//get the information of the particle
	pos1 = pt.Pos[id1];
	vel1 = pt.Vel[id1];
	//QuatSet(orien1, x.Quat[id1][0], x.Quat[id1][1], x.Quat[id1][2], x.Quat[id1][3]);
	omg1 = pt.AngSpd[id1];
	vector4d orien1 = pt.Quat[id1];
	matrix3d IvDCM1 = Quat2IvDCM(orien1);
	omg1 = IvDCM1 * omg1;

	pos2 = wall.Orig[id2];
	vel2 = vector3d(0, 0, 0);//default speed is 0 
	omg2 = vector3d(0, 0, 0);



	l_j = - (pos1 - pos2).dot(normal);
	delta = pt.Radius[id1] - l_j;
	if (delta > 0)
	{
		flag = 1;//函数返回值
		//读取当前接触对的信息；若接触对不存在，则新建接触对---------------------------------------------------------------------
		if (Cflag == 0)//
		{
			Cp* Cpair = (Cp*)malloc(sizeof(Cp));
			assert(Cpair != nullptr);
			for (int i = 0; i < 3; i++)
			{
				Cpair->deltP[0][i] = 0.0;//roll
				Cpair->deltP[1][i] = 0.0;//tangential
				Cpair->deltP[2][i] = 0.0;//twist
				HR[i] = Cpair->deltP[0][i];
				HS[i] = Cpair->deltP[1][i];
				HT[i] = Cpair->deltP[2][i];
			}
			Cpair->idx = id2;
			Cpair->flag = 1;
			Cpair->body = -2;//-2表明这个接触对是颗粒与墙的接触对
			Cpair->next = nullptr;
			part->next = Cpair;
			part = part->next;
		}
		else
		{
			for (int i = 0; i < 3; i++)
			{
				HR[i] = part->deltP[0][i];
				HS[i] = part->deltP[1][i];
				HT[i] = part->deltP[2][i];
			}
		}

		// Calculate derived contact parameters from the particle-wall set.
		mech.kS = 2 * mech.kN / 7;
		double mu = pt.Mass[id1];//mi*mj/(mi+mj)
		double lnepsN = log(mech.epsN);//initial
		double tmp = mech.kN * mu / (PI * PI + lnepsN * lnepsN);
		mech.cN = -2 * lnepsN * sqrt(tmp);
		double lnepsS = log(mech.epsS);//initial
		tmp = mech.kS * mu / (PI * PI + lnepsS * lnepsS);
		mech.cS = -2 * lnepsS * sqrt(tmp);  //没有推导，随便给的，后面需要重新推导表达式
		double R = pt.Radius[id1];
		tmp = (mech.beta * R) * (mech.beta * R);
		mech.kT = 2 * mech.kS * tmp;
		mech.cT = 2 * mech.cS * tmp;
		mech.kR = mech.kN * tmp;
		mech.cR = mech.cN * tmp;

		//计算接触点的相对速度、相对位移等信息---------------------------------------
		l_i = R - delta;
		Vr = vel2 - vel1;//Vr
		//calculate tangential velocity (V-U+Wxr-[(V-U)*n]n)
		velS = -normal * Vr.dot(normal) + Vr + normal.cross(omg1) * l_i;

		//calculate twist angle velocity
		vector3d RelativeAngularVelocity = omg2 - omg1;
		omgT = normal * RelativeAngularVelocity.dot(normal);//扭转角速度

		//calculate roll angle velocity
		tmp = 0.5 / l_i;//unit:m^(-1);// *(l_j - l_i) / (l_i * l_j);//l_j is infinite on the board
		omgR = RelativeAngularVelocity - omgT +tmp * normal.cross(velS);//滚转角速度

		//integral displacement; 
		//calculate tangential contact force
		deltS = velS * dt;
		//calculate twist contact torque
		deltT = omgT * dt;
		//calculate roll contact torque
		deltR = omgR * dt;

		// 计算法向力--------------------------------------------------------------
		//FN_elastic
		FN_elastic = - normal * mech.kN * delta;
		//FN_damping
		FN_damping =  normal * mech.cN * Vr.dot(normal);
		FN_total = FN_elastic + FN_damping;

		tmp = mech.beta * R * (FN_total.norm());
		double MT_max = mech.mu * mech.mu_T * tmp;
		double MR_max = mech.mu_R * tmp;

		// 计算摩擦力--------------------------------------------------------------
		FrictionForce_static = velS * mech.cS + deltS * mech.kS + HS;
		FrictionForce_static_elastic = deltS * mech.kS + HS;
		tmpd1 = FrictionForce_static_elastic.norm();
		tmpd2 = mech.mu * FN_total.norm();
		if (tmpd1 > tmpd2)// force alone exceeds Ft_max
		{
			double dS = deltS.norm();//
			if (dS > 1.0e-14)
			{
				FrictionForce = tmpd2 * deltS.normalized();
				for (int i = 0; i < 3; i++)
				{
					part->deltP[1][i] = FrictionForce[i];//store force instead of relative displancement
				}
			}
			else
			{
				//FrictionForce.setZero();
				FrictionForce = HS;
				for (int i = 0; i < 3; i++)
				{
				part->deltP[1][i] = FrictionForce[i];//store force instead of relative displancement
				}
			}
			//cout << "摩擦力大小为：" << FrictionForce.norm() << endl;
		}
		else
		{
			FrictionForce = FrictionForce_static;
			for (int i = 0; i < 3; i++)
			{
				part->deltP[1][i] = FrictionForce_static_elastic[i];//store force instead of relative displancement
			}
			//cout << "静摩擦，摩擦力大小为：" << FrictionForce.norm() << endl;
		}
		M_FSi = l_i * normal.cross(FrictionForce);//颗粒受到的由摩擦力产生的对球心的力矩
		//-------------------------------------------------------------------------------------
			TwistTorque_static = omgT * mech.cT + deltT * mech.kT + HT;
			TwistTorque_static_elastic = deltT * mech.kT + HT;
			tmpd1 = TwistTorque_static_elastic.norm();
			if (tmpd1 > MT_max)
			{
				double dT = deltT.norm();
				if (dT > 1.0e-14)
				{
					TwistTorque = MT_max* deltT.normalized();
					for (int i = 0; i < 3; i++)
					{
						part->deltP[2][i] = TwistTorque[i];//store force instead of relative displancement
					}
				}
				else
				{
					TwistTorque.setZero();
					//TwistTorque = HT;
					for (int i = 0; i < 3; i++)
					{
						part->deltP[2][i] = TwistTorque[i];//store force instead of relative displancement
					}
				}
			}
			else
			{
				TwistTorque = TwistTorque_static;
				for (int i = 0; i < 3; i++)
				{
					part->deltP[2][i] = TwistTorque_static_elastic[i];//store force instead of relative displancement
				}
			}
		//---------------------------------------------------------------------------------------------------
			RollTorque_static = omgR * mech.cR + deltR * mech.kR + HR;
			RollTorque_static_elastic = deltR * mech.kR + HR;
			tmpd1 = RollTorque_static.norm();
			if (tmpd1 > MR_max) {
				double dR = deltR.norm();
				if (dR > 1.0e-14)
				{
					RollTorque = MR_max * deltR.normalized();
					for (int i = 0; i < 3; i++)		part->deltP[0][i] = RollTorque[i];//store force instead of relative displancement
				}
				else
				{
					RollTorque.setZero();
					//RollTorque = HR;
					for (int i = 0; i < 3; i++)		part->deltP[0][i] = RollTorque[i];//store force instead of relative displancement
				}
			}
			else {
				RollTorque = RollTorque_static;
				for (int i = 0; i < 3; i++)		part->deltP[0][i] = RollTorque_static_elastic[i];//store force instead of relative displancement
			}
		Force = FN_total + FrictionForce;
		Torque = M_FSi + TwistTorque + RollTorque;
	}
	//若接触对曾经存在但现在不再接触，则删除该接触对，并将part指针指向链表头
	if (Cflag == 1 && flag == 0)//delete the contact pair
	{
		Cp* tmp = head;
		tmp = tmp->next;
		while (tmp)
		{
			if ((tmp->idx == id2) && (tmp->body == -2))
			{
				for (int i = 0; i < 3; i++)
				{
					tmp->deltP[i][0] = 0;
					tmp->deltP[i][1] = 0;
					tmp->deltP[i][2] = 0;
				}
				pre->next = tmp->next;
				free(tmp);
				tmp = nullptr;
				return flag;
			}
			tmp = tmp->next;
			pre = pre->next;
		}
	}
	return flag;
}

//颗粒之间的接触力
int ContactForce(int id1, int id2, PARTICLE& pt, MECH& Mech, vector3d& Force1, vector3d& Torque1, vector3d& Torque2, ODE ode, Cp* part)
{
	//颗粒的速度和位置是在惯性系下定义的，但是颗粒的角速度是在颗粒的本体系下定义的
	int flag = 0, Cflag = 0;
	double tmpd1, tmpd2, l_i, l_j, distance, distanceSq, Aeff, R;
	vector3d deltR, deltS, deltT, HR, HS, HT;
	vector3d tmpF1, tmpF2, Vr, tmpT, tmpR, MT, MR, M_FSi, M_FSj;

	MECH mech = Mech;
	double depth;//嵌入深度
	double dt = ode.StepSize;

	vector3d velS, omgT, omgR;//接触点相对切向速度、相对角速度在扭转方向的分量、相对角速度在滚转方向的分量

	vector3d pos1, vel1, omg1, pos2, vel2, omg2;//位置、速度、角速度；
	vector4d orien1, orien2;//四元数
	matrix3d IvDCM1, IvDCM2;

	Cp* pre = part;
	Cp* head = part;
	while (part->next != nullptr)//在链表中寻找是否存在与id2颗粒组成的接触对，若没有找到，这里不做处理
	{
		if (part->next->idx == id2)
		{
			Cflag = 1;// The contact pair have been built
			part = part->next;
			break;//提前退出循环
		}
		part = part->next;
	}

	pos1 = pt.Pos[id1];
	pos2 = pt.Pos[id2];

	distanceSq = (pos2 - pos1).squaredNorm();
	double SqDepth = (pt.Radius[id1] + pt.Radius[id2]) * (pt.Radius[id1] + pt.Radius[id2]) - distanceSq;

	if (SqDepth > 0)
	{

		flag = 1;//函数返回值，为1则会将接触力叠加
		distance = sqrt(distanceSq);
		depth = pt.Radius[id1] + pt.Radius[id2] - distance;//嵌入深度

		vel1 = pt.Vel[id1];		orien1 = pt.Quat[id1];		omg1 = pt.AngSpd[id1];
		vel2 = pt.Vel[id2];		orien2 = pt.Quat[id2];		omg2 = pt.AngSpd[id2];
		vector3d dirVector = (pos2 - pos1).normalized();//n,由颗粒1指向颗粒2的单位向量
		QuatNorm(orien1);
		QuatNorm(orien2);
		IvDCM1 = Quat2IvDCM(orien1);
		IvDCM2 = Quat2IvDCM(orien2);
		omg2 = IvDCM2 * omg2;//角速度的三个分量转化到惯性系下
		omg1 = IvDCM1 * omg1;

		Coefficient(id1, id2, pt, mech);//计算力学参数

		//法向接触力--------------------------------------------------------------------------------------------------------------------------------
		vector3d ContaceForce_elastic = -mech.kN * depth * dirVector;//FN_elastic
		//cout << mech.kN << endl;
		Vr = vel2 - vel1;//Vr相对速度
		vector3d ContaceForce_damping = mech.cN * (Vr.dot(dirVector)) * dirVector;//FN_damping
		Force1 = ContaceForce_elastic + ContaceForce_damping;

		//calculate cohesion force粘附力-----------------------------------------------------------------------------------------------------------
		R = pt.Radius[id1] * pt.Radius[id2] / (pt.Radius[id1] + pt.Radius[id2]);
		Aeff = 4 * (mech.beta * R) * (mech.beta * R);
		Force1 = Force1 + dirVector * mech.c * Aeff;
		//范德华粘附力对毫米量级以上大颗粒的作用非常微弱,world.par中将系数设置为了0

		double tmp = mech.beta * R * (Force1.norm());
		double MT_max = mech.mu * mech.mu_T * tmp;
		double MR_max = mech.mu_R * tmp;

		//contact histroy; 
		if (Cflag == 0)//若没有与id2颗粒组成的接触对，则建立一个
		{
			Cp* Cpair = (Cp*)malloc(sizeof(Cp));
			if (Cpair == nullptr) {
				// 明确处理内存分配失败
				fprintf(stderr, "Memory allocation failed for contact pair");
			}
			for (int i = 0; i < 3; i++)
			{
				Cpair->deltP[0][i] = 0.0;//roll
				Cpair->deltP[1][i] = 0.0;//tangential
				Cpair->deltP[2][i] = 0.0;//twist
				HR[i] = Cpair->deltP[0][i];
				HS[i] = Cpair->deltP[1][i];
				HT[i] = Cpair->deltP[2][i];
			}
			Cpair->idx = id2;
			Cpair->flag = 1;
			Cpair->next = nullptr;
			part->next = Cpair;
			part = part->next;
		}
		else
		{
			for (int i = 0; i < 3; i++)
			{
				HR[i] = part->deltP[0][i];
				HS[i] = part->deltP[1][i];
				HT[i] = part->deltP[2][i];
			}
		}

		//
		l_i = (pt.Radius[id1] * pt.Radius[id1] - pt.Radius[id2] * pt.Radius[id2] + distanceSq) / (2 * distance);
		l_j = distance - l_i;
		//calculate tangential velocity 两颗粒于接触点处的切向相对速度
		velS = Vr - dirVector * Vr.dot(dirVector) + l_i * dirVector.cross(omg1) + l_j * dirVector.cross(omg2);
		deltS = velS * dt;

		//calculate tangential contact force切向力---------------------------------------------------------------------------------------------------------
		vector3d FrictionForce;
		vector3d FrictionForce_static = velS * mech.cS + deltS * mech.kS + HS;//tengential force
		vector3d FrictionForce_static_elastic = deltS * mech.kS + HS;//仅由弹性变形产生的切向力
		//cout << "切向力为：" << tmpF2.transpose() << endl;
		tmpd1 = FrictionForce_static_elastic.norm();//静摩擦力的大小
		tmpd2 = mech.mu * (Force1.norm());//动摩擦力的大小；此处force1的方向就是沿着法向
		if (tmpd1 > tmpd2)// slipping
		{
			double dS = deltS.norm();//dt步长内接触点相对位移的增量
			if (dS > 1.0e-14)
			{
				FrictionForce = tmpd2 * deltS.normalized();
				for (int i = 0; i < 3; i++)
				{
					part->deltP[1][i] = FrictionForce[i];//store force instead of relative displancement
				}
			}
			else
			{
				//FrictionForce.setZero();
				FrictionForce = HS;
				for (int i = 0; i < 3; i++)
				{
					part->deltP[1][i] = FrictionForce[i];//store force instead of relative displancement
				}
			}
		}
		else
		{
			FrictionForce = FrictionForce_static;
			for (int i = 0; i < 3; i++)
			{
				part->deltP[1][i] = FrictionForce_static_elastic[i];//store force instead of relative displancement
			}
		}
		Force1 = Force1 + FrictionForce;
		M_FSi = l_i * dirVector.cross(FrictionForce);//由切向力产生的力矩

		//calculate twist angle velocity相对扭转角速度------------------------------------------------------------------------------------------
		vector3d RelativeAngularVelocity = omg2 - omg1;
		omgT = dirVector * (RelativeAngularVelocity.dot(dirVector));

		//calculate twist contact torque
		deltT = omgT * dt;

		vector3d TwistTorque;
		vector3d TwistTorque_static = omgT * mech.cT + deltT * mech.kT + HT;
		vector3d TwistTorque_static_elastic = deltT * mech.kT + HT;
		tmpd1 = TwistTorque_static_elastic.norm();
		if (tmpd1 > MT_max)
		{
			double dT = deltT.norm();//dt时间步长内，扭转角的变化量
			if (dT > 1.0e-14)
			{
				TwistTorque = MT_max * deltT.normalized();
				for (int i = 0; i < 3; i++)
				{
					part->deltP[2][i] = TwistTorque[i];//store force instead of relative displancement
				}
			}
			else
			{
				//TwistTorque.setZero();
				TwistTorque = HT;
				for (int i = 0; i < 3; i++)
				{
					part->deltP[2][i] = TwistTorque[i];//store force instead of relative displancement
				}
			}
		}
		else{
			TwistTorque = TwistTorque_static;
			for (int i = 0; i < 3; i++)		part->deltP[2][i] = TwistTorque_static_elastic[i];//store force instead of relative displancement

		}
		//calculate roll angle velocity相对滚转角速度------------------------------------------------------------------------------------------
		tmp = 0.5 * (l_j - l_i) / (l_i * l_j);
		omgR = RelativeAngularVelocity - omgT;// +tmp * dirVector.cross(velS);

		//calculate roll contact torque
		deltR = omgR * dt;

		vector3d RollTorque;
		vector3d RollTorque_static = omgR * mech.cR + deltR * mech.kR + HR;
		vector3d RollTorque_static_elastic = deltR * mech.kR + HR;
		tmpd1 = RollTorque_static_elastic.norm();
		if (tmpd1 > MR_max){
			double dR = deltR.norm();
			if (dR > 1.0e-14)
			{
				RollTorque = MR_max * deltR.normalized();
				for (int i = 0; i < 3; i++)		part->deltP[0][i] = RollTorque[i];//store force instead of relative displancement
			}
			else
			{
				//RollTorque.setZero();
				RollTorque = HR;
				for (int i = 0; i < 3; i++)		part->deltP[0][i] = RollTorque[i];//store force instead of relative displancement
			}
		}
		else{
			RollTorque = RollTorque_static;
			for (int i = 0; i < 3; i++)		part->deltP[0][i] = RollTorque_static_elastic[i];//store force instead of relative displancement
		}

		//resultant torque
		Torque1 = TwistTorque + RollTorque + M_FSi;
	}

	if (Cflag == 1 && flag == 0)//delete the contact pair若接触对已存在并且在本次计算中没有接触，则删除该接触对
	{
		Cp* tmp = head;
		tmp = tmp->next;
		while (tmp)
		{
			if (tmp->idx == id2)
			{
				pre->next = tmp->next;
				free(tmp);
				tmp = nullptr;
				return flag;
			}
			tmp = tmp->next;
			pre = pre->next;
		}
	}
	return flag;
}

void Coefficient(int id1, int id2, PARTICLE& pt, MECH& mech)
{
	int count = 0, maxC = 0;
	double mu = 0, delt = 0;
	//
	mu = pt.Mass[id1] * pt.Mass[id2] / (pt.Mass[id1] + pt.Mass[id2]);

	mech.kS = 2 * mech.kN / 7;//initial 

	double lnepsN = log(mech.epsN);//initial
	double tmp = mech.kN * mu / (PI * PI + lnepsN * lnepsN);
	mech.cN = -2 * lnepsN * sqrt(tmp);

	double lnepsS = log(mech.epsS);//initial
	tmp = mech.kS * mu / (PI * PI + lnepsS * lnepsS);
	mech.cS = -2 * lnepsS * sqrt(tmp);  //没有推导，随便给的，后面需要重新推导表达式

	double R = pt.Radius[id1] * pt.Radius[id2] / (pt.Radius[id1] + pt.Radius[id2]);
	tmp = (mech.beta * R) * (mech.beta * R);
	mech.kT = 2 * mech.kS * tmp;
	mech.cT = 2 * mech.cS * tmp;
	mech.kR = mech.kN * tmp;
	mech.cR = mech.cN * tmp;
}

//颗粒与三角面之间的接触力
int Polyhedron_Force(int id1, int id2, PARTICLE& pt, MECH& Mech, vector3d& Force1, vector3d& Torque1, vector3d& Torque2,
	double* ForceBodyTemp, ODE ode, Cp* part, BODYSET& bodyset, int idbody,VHALF& vhalf_pt, VHALF& vhalf_bd)
{

	//cout << "id1 = " << id1 << ", id2 = " << id2 << ", idbody = " << idbody << endl;

	int flag = 0, Cflag = 0;
	double tmpd1, tmpd2, l_i = 0.0, l_j = 0.0, Aeff, R, C_S;
	vector3d  deltR, deltS, deltT, HR, HS, HT;
	vector3d tmpF1, tmpF2,tmpR, M_FSi, M_FSj;
	vector3d pos1, vel1, omg1, pos2, vel2, omg2, velS, omgT, omgR;

	MECH mech = Mech;
	// {
	// 	mech.mu = 0.13;
	// 	mech.epsN = 0.88;
	// 	mech.epsS = 0.88;
	// 	mech.mu_T = 0.0843;
	// 	mech.mu_R = 0.0843;
	// }

	BODY intruder = bodyset.body[idbody];

	for(int i = 0; i < 19; i++) {
		ForceBodyTemp[i] = 0.0;
	}

	vector4d orien1, orien2;
	matrix3d IvDCM1, IvDCM2;
	orien1 = pt.Quat[id1];
	IvDCM1 = Quat2IvDCM(orien1);
	orien2 = intruder.orien;
	IvDCM2 = bodyset.IvDCM[idbody];

	vector3d v2 = intruder.Shape[id2].v2;
	vector3d v3 = intruder.Shape[id2].v3;
	vector3d v4 = intruder.Shape[id2].v4;
	vector3d normal = intruder.Shape[id2].normal;//三角面外法线方向，已单位化

	vector3d MassCenter = intruder.MassCenter;//刚体的质心位置
	v2 = MassCenter + IvDCM2 * v2;//将其转换到惯性系下
	v3 = MassCenter + IvDCM2 * v3;//将其转换到惯性系下
	v4 = MassCenter + IvDCM2 * v4;//将其转换到惯性系下
	normal = IvDCM2 * normal;//将其转换到惯性系下

	pos1 =  pt.Pos[id1];
	vector3d tmpV3d = pos1 - v2;//某一顶点到球心的矢量

	double distance = normal.dot(tmpV3d);// abs();//球心到三角面的垂直距离
	double depth = pt.Radius[id1] - distance;

	//if (distance < 0)//  return 0;
	//{
	//	distance = -distance;
	//	normal = -normal;
	//}

	Cp* pre = part;
	Cp* head = part;

	while (part->next != nullptr)
	{
		if ((part->next->idx == id2) && (part->next->body == idbody))
		{
			Cflag = 1;// 1 sphere-face contact;2 sphere-edge contact;3 sphere-vertex; 0 no contact
			part = part->next;
			break;
		}
		part = part->next;
	}

	if (depth > 0)//有可能相交
	{
		const double contactArea = contact_geometry::sphereTriangleIntersectionArea(
			pt.Pos[id1], pt.Radius[id1], v2, v3, v4);
		//std::cout << "计算接触面积耗时: " << (end - start) * 1000 << " ms" << std::endl;

		if (contactArea > 0)
		{
			flag = 1;

			//重新计算depth
			//vector3d closest;
			//double depth = checkContact(v2, v3, v4, normal, pt.Pos[id1], pt.Radius[id1], closest);
			//tmpR = closest - MassCenter;//刚体质心指向接触点的矢量//设定接触点为距离球心的最近点


			vector3d pos_circle = pos1 - normal * distance;//接触圆圆心坐标
			tmpR = pos_circle - MassCenter;
			
			double sq_circle = pt.Radius[id1] * pt.Radius[id1] - distance * distance;//接触圆的半径平方

			C_S = contactArea / (PI * sq_circle);//接触区域面积(改为与接触圆面积的比值方便后面乘）
			if (C_S > 1)	C_S = 1.0;

			normal = -normal;//法线方向确定

			vel1 = pt.Vel[id1];
			omg1 = pt.AngSpd[id1];
			omg1 = IvDCM1 * omg1;

			omg2 = intruder.AngularVel;//刚体角速度
			omg2 = IvDCM2 * omg2;
			vel2 = omg2.cross(tmpR) + intruder.Vel;//刚体在接触点处的速度（这个速度已经包括转动速度了）

			mech.kS = 2 * mech.kN / 7;
			double mu = pt.Mass[id1];//mi*mj/(mi+mj)
			double lnepsN = log(mech.epsN);//initial
			double tmp = mech.kN * mu / (PI * PI + lnepsN * lnepsN);
			mech.cN = -2 * lnepsN * sqrt(tmp);
			double lnepsS = log(mech.epsS);//initial
			tmp = mech.kS * mu / (PI * PI + lnepsS * lnepsS);
			mech.cS = -2 * lnepsS * sqrt(tmp);  //没有推导，随便给的，后面需要重新推导表达式
			R = pt.Radius[id1];
			tmp = (mech.beta * R) * (mech.beta * R);
			mech.kT = 2 * mech.kS * tmp;
			mech.cT = 2 * mech.cS * tmp;
			mech.kR = mech.kN * tmp;
			mech.cR = mech.cN * tmp;
			double dt = ode.StepSize;
			Aeff = 4 * tmp;

			//cout << "mech.cN = " << mech.cN << endl;
			//cout << "mech.kR = " << mech.kR << endl;
			//cout << "mech.cR = " << mech.cR << endl;
			//cout << "mech.kS = " << mech.kS << endl;
			//cout << "mech.cS = " << mech.cS << endl;

			l_i = distance;
			vector3d Vr = vel2 - vel1;//Vr

			//法向接触力--------------------------------------------------------------------------------------------------------
			vector3d ContaceForce_elastic = -normal * mech.kN * depth;//FN_elastic
			vector3d ContaceForce_damping = normal * mech.cN * Vr.dot(normal);//FN_damping
			vector3d ContaceForce = ContaceForce_elastic + ContaceForce_damping;
			Force1 = ContaceForce;

			//calculate cohesion force----------------------------------------------------------------------------------------------
			Force1 = Force1 + normal * mech.c * Aeff;
			
			tmp = mech.beta * R * (Force1.norm());
			double MT_max = mech.mu * mech.mu_T * tmp;
			double MR_max = mech.mu_R * tmp;

			//calculate tangential velocity (V-U+Wxr-[(V-U)*n]n)
			velS = -normal * Vr.dot(normal) + Vr + normal.cross(omg1) * l_i;

			//calculate twist angle velocity
			vector3d RelativeAngularVelocity = omg2 - omg1;
			omgT = normal * RelativeAngularVelocity.dot(normal);//扭转角速度

			//calculate roll angle velocity
			tmp = 0.5 / l_i;//unit:m^(-1);// *(l_j - l_i) / (l_i * l_j);//l_j is infinite on the board
			omgR = RelativeAngularVelocity - omgT +tmp * normal.cross(velS);//滚转角速度

			//integral displacement; 
			//calculate tangential contact force
			deltS = velS * dt;
			//calculate twist contact torque
			deltT = omgT * dt;
			//calculate roll contact torque
			deltR = omgR * dt;

			//读取当前接触对的信息；若接触对不存在，则新建接触对---------------------------------------------------------------------
			if (Cflag == 0)//
			{
				Cp* Cpair = (Cp*)malloc(sizeof(Cp));
				assert(Cpair != nullptr);
				for (int i = 0; i < 3; i++)
				{
					Cpair->deltP[0][i] = 0.0;//roll
					Cpair->deltP[1][i] = 0.0;//tangential
					Cpair->deltP[2][i] = 0.0;//twist
					HR[i] = Cpair->deltP[0][i];
					HS[i] = Cpair->deltP[1][i];
					HT[i] = Cpair->deltP[2][i];
				}
				Cpair->idx = id2;
				Cpair->flag = 1;
				Cpair->body = idbody;
				Cpair->next = nullptr;
				part->next = Cpair;
				part = part->next;
			}
			else
			{
				for (int i = 0; i < 3; i++)
				{
					HR[i] = part->deltP[0][i];
					HS[i] = part->deltP[1][i];
					HT[i] = part->deltP[2][i];
				}
			}

			//---------------------------------------------------------------------------------------------------------------
			vector3d FrictionForce;
			vector3d FrictionForce_static = velS * mech.cS + deltS * mech.kS + HS;
			vector3d FrictionForce_static_elastic = deltS * mech.kS + HS;
			tmpd1 = FrictionForce_static_elastic.norm();
			tmpd2 = mech.mu * Force1.norm();
			if (tmpd1 > tmpd2)// force alone exceeds Ft_max
			{
				double dS = deltS.norm();//
				if (dS > 1.0e-14)
				{
					FrictionForce = tmpd2 * deltS.normalized();
					for (int i = 0; i < 3; i++)
					{
						part->deltP[1][i] = FrictionForce[i];//store force instead of relative displancement
					}
				}
				else
				{
					//FrictionForce.setZero();
					FrictionForce = HS;
					for (int i = 0; i < 3; i++)
					{
						part->deltP[1][i] = FrictionForce[i];//store force instead of relative displancement
					}
				}
				//cout << "摩擦力大小为：" << FrictionForce.norm() << endl;
			}
			else
			{
				FrictionForce = FrictionForce_static;
				for (int i = 0; i < 3; i++)
				{
					part->deltP[1][i] = FrictionForce_static_elastic[i];//store force instead of relative displancement
				}
				//cout << "静摩擦，摩擦力大小为：" << FrictionForce.norm() << endl;
			}

			Force1 = Force1 + FrictionForce;
			//FrictionForce =  ContaceForce.norm() * FrictionModelOfAdams(velS, mech);

			vector3d ContaceForceTorque, FrictionForceTorque;//刚体受到的法向接触力产生的力矩、摩擦力产生的力矩
			ContaceForceTorque = -tmpR.cross(ContaceForce);
			FrictionForceTorque = -tmpR.cross(FrictionForce);

			M_FSi = l_i * normal.cross(FrictionForce);//颗粒受到的由摩擦力产生的对球心的力矩
			M_FSj = ContaceForceTorque + FrictionForceTorque; //刚体受到的由法向接触力和摩擦力产生的对质心的力矩

			//-------------------------------------------------------------------------------------
			vector3d TwistTorque;
			vector3d TwistTorque_static = omgT * mech.cT + deltT * mech.kT + HT;
			vector3d TwistTorque_static_elastic = deltT * mech.kT + HT;
			tmpd1 = TwistTorque_static_elastic.norm();
			if (tmpd1 > MT_max)
			{
				double dT = deltT.norm();
				if (dT > 1.0e-14)
				{
					TwistTorque = MT_max* deltT.normalized();
					for (int i = 0; i < 3; i++)
					{
						part->deltP[2][i] = TwistTorque[i];//store force instead of relative displancement
					}
				}
				else
				{
					//TwistTorque.setZero();
					TwistTorque = HT;
					for (int i = 0; i < 3; i++)
					{
						part->deltP[2][i] = TwistTorque[i];//store force instead of relative displancement
					}
				}
			}
			else
			{
				TwistTorque = TwistTorque_static;
				for (int i = 0; i < 3; i++)
				{
					part->deltP[2][i] = TwistTorque_static_elastic[i];//store force instead of relative displancement
				}
			}
			//---------------------------------------------------------------------------------------------------

			vector3d RollTorque;
			vector3d RollTorque_static = omgR * mech.cR + deltR * mech.kR + HR;
			vector3d RollTorque_static_elastic = deltR * mech.kR + HR;
			tmpd1 = RollTorque_static_elastic.norm();
			if (tmpd1 > MR_max) {
				double dR = deltR.norm();
				if (dR > 1.0e-14)
				{
					RollTorque = MR_max * deltR.normalized();
					for (int i = 0; i < 3; i++)		part->deltP[0][i] = RollTorque[i];//store force instead of relative displancement
				}
				else
				{
					//RollTorque.setZero();
					RollTorque = HR;
					for (int i = 0; i < 3; i++)		part->deltP[0][i] = RollTorque[i];//store force instead of relative displancement
				}
			}
			else {
				RollTorque = RollTorque_static;
				for (int i = 0; i < 3; i++)		part->deltP[0][i] = RollTorque_static_elastic[i];//store force instead of relative displancement
			}


			Force1 = Force1 * C_S;

			//resultant torque
			Torque1 = TwistTorque + RollTorque;
			Torque2 = -Torque1;
			
			Torque1 = Torque1 + M_FSi;
			Torque2 = Torque2 + M_FSj;

			Torque1 = Torque1 * C_S;
			Torque2 = Torque2 * C_S;

			for (int i = 0; i < 3; i++) {
				ForceBodyTemp[i] = -ContaceForce[i];
				ForceBodyTemp[i + 3] = -FrictionForce[i];
				ForceBodyTemp[i + 6] = ContaceForceTorque[i];
				ForceBodyTemp[i + 9] = FrictionForceTorque[i];
				ForceBodyTemp[i + 12] = -TwistTorque[i];
				ForceBodyTemp[i + 15] = -RollTorque[i];
			}
			ForceBodyTemp[18] = depth;
			for (int i = 0; i < 19; i++)	ForceBodyTemp[i] *= C_S;
		}
	}
	if (Cflag == 1 && flag == 0)//delete the contact pair
	{
		Cp* tmp = head;
		tmp = tmp->next;
		while (tmp)
		{
			if ((tmp->idx == id2) && (tmp->body == idbody))
			{
				for (int i = 0; i < 3; i++)
				{
					tmp->deltP[i][0] = 0;
					tmp->deltP[i][1] = 0;
					tmp->deltP[i][2] = 0;
				}
				pre->next = tmp->next;
				free(tmp);
				tmp = nullptr;
				return flag;
			}
			tmp = tmp->next;
			pre = pre->next;
		}
	}
	return flag;
}

// 判断点Q是否在三角形ABC内部（使用重心坐标法）
bool isPointInsideTriangle(const vector3d& A, const vector3d& B, const vector3d& C, const vector3d& Q) {
	vector3d AB = B - A;
	vector3d AC = C - A;
	vector3d AQ = Q - A;

	double ab_dot_ab = AB.dot(AB);
	double ab_dot_ac = AB.dot(AC);
	double ac_dot_ac = AC.dot(AC);
	double aq_dot_ab = AQ.dot(AB);
	double aq_dot_ac = AQ.dot(AC);

	double denominator = ab_dot_ab * ac_dot_ac - ab_dot_ac * ab_dot_ac;

	if (fabs(denominator) < 1e-12) {
		return false; // 处理退化三角形
	}

	double u = (ac_dot_ac * aq_dot_ab - ab_dot_ac * aq_dot_ac) / denominator;
	double v = (ab_dot_ab * aq_dot_ac - ab_dot_ac * aq_dot_ab) / denominator;

	return (u >= 0) && (v >= 0) && (u + v <= 1);
}

// 计算点P到线段AB的最近点
vector3d closestPointOnSegment(const vector3d& A, const vector3d& B, const vector3d& P) {
	vector3d AB = B - A;
	double t = (P - A).dot(AB) / AB.squaredNorm();
	t = max(0.0, min(1.0, t));
	return A + t * AB;
}

// 计算三角形上离球心P最近的点
vector3d computeClosestPoint(const vector3d& A, const vector3d& B, const vector3d& C, const vector3d& n, const vector3d& P) {
	vector3d n_normalized = n.normalized();
	vector3d AP = P - A;
	double distance = AP.dot(n_normalized);
	vector3d Q = P - distance * n_normalized;

	if (isPointInsideTriangle(A, B, C, Q)) {
		return Q;
	}

	vector3d closestAB = closestPointOnSegment(A, B, P);
	vector3d closestBC = closestPointOnSegment(B, C, P);
	vector3d closestCA = closestPointOnSegment(C, A, P);

	double distAB = (closestAB - P).squaredNorm();
	double distBC = (closestBC - P).squaredNorm();
	double distCA = (closestCA - P).squaredNorm();

	if (distAB <= distBC && distAB <= distCA) {
		return closestAB;
	}
	else if (distBC <= distCA) {
		return closestBC;
	}
	else {
		return closestCA;
	}
}

// 判断是否存在接触并返回最近点
double checkContact(const vector3d& A, const vector3d& B, const vector3d& C, const vector3d& n, const vector3d& P, double radius, vector3d& closestPoint) {
	closestPoint = computeClosestPoint(A, B, C, n, P);
	double distance = (closestPoint - P).norm();
	double depth = radius - distance;
	return depth;
}

vector3d FrictionModelOfAdams(const vector3d& v, const MECH& mech)
{
	vector3d F_mu;

	double normV = v.norm();

	double alpha = 10;
	double c = 0.0001;

	double p = mech.mu - (1 - exp(-alpha * normV / c));

	if (normV > 1e-14)
		F_mu = -p * v.normalized();
	else
		F_mu.setZero();

	return F_mu;
}

void ZeroForceStack(double** F, int n)
{
	int i, j;
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < 6; j++)
		{
			F[i][j] = 0.0;
		}
	}
}
