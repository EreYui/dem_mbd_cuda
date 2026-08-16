//  tree.cpp
//  Created by zkm on 2024.12.13.
//

#include"grid.h"
#include"multibody.h"
#include"particle.h"

#include"matrixtrans.h"
#include <array>
#include <unordered_set>

void CreatTree(PARTICLE& pt, BODYSET& bodyset, unordered_map<GridIndex, GridCell> &grid)
{
	//cout << "diaoyong CreatTree" << endl;
	for (auto it = grid.begin(); it != grid.end(); ++it) {
		// it->first Ǽit->second ֵ
		it->second.clear();
	}
	grid.clear();
	std::unordered_map<GridIndex, GridCell> empty_map;
	//grid.swap(empty_map);

	bodyset.updateMatrix();
	
	const int numParticles = pt.Num;
	if (numParticles == 0) return;

	// 1. ռ߽
	vector3d lower = pt.Pos[0];
	vector3d upper = pt.Pos[0];

	for (int i = 1; i < numParticles; ++i) {
		lower = lower.min(pt.Pos[i]);
		upper = upper.max(pt.Pos[i]);//ϵͳı߽
	}

	vector3d v0;
	for (int i = 0; i < bodyset.Num; i++)
	{
		v0 = bodyset.body[i].Shape[0].v1;
		bodyset.BodyLowerLmt[i] = v0;//ʼһںΧڵֵ
		bodyset.BodyUpperLmt[i] = v0;
	}

	for (int ib = 0; ib < bodyset.Num; ib++)
	{
		for (int i = 0; i < bodyset.body[ib].Size; i++)//identification the range of the mesh
		{
			vector3d v1 = bodyset.body[ib].Shape[i].v1;
			vector3d v2 = bodyset.body[ib].Shape[i].v2;
			vector3d v3 = bodyset.body[ib].Shape[i].v3;
			vector3d v4 = bodyset.body[ib].Shape[i].v4;

			v1 = bodyset.IvDCM[ib] * v1;//תϵ
			v2 = bodyset.IvDCM[ib] * v2;
			v3 = bodyset.IvDCM[ib] * v3;
			v4 = bodyset.IvDCM[ib] * v4;

			//ibϵĵi4ڹϵµ
			vector3d ver1 = bodyset.body[ib].MassCenter + v1;
			vector3d ver2 = bodyset.body[ib].MassCenter + v2;
			vector3d ver3 = bodyset.body[ib].MassCenter + v3;
			vector3d ver4 = bodyset.body[ib].MassCenter + v4;

			std::array<vector3d, 4> vertics{ ver1, ver2, ver3, ver4 };
			bodyset.body[ib].eAABB[i] = AABB();//ðΧУȻµİΧлԽԽ
			bodyset.body[ib].eAABB[i].include(vertics);//õ˰Χеֵ

			bodyset.BodyLowerLmt[ib] = bodyset.BodyLowerLmt[ib].min(bodyset.body[ib].eAABB[i].min);
			bodyset.BodyUpperLmt[ib] = bodyset.BodyUpperLmt[ib].max(bodyset.body[ib].eAABB[i].max);
		}
	}

	//ÿı߽
	bodyset.BodySetLowerLmt = bodyset.BodyLowerLmt[0];
	bodyset.BodySetUpperLmt = bodyset.BodyUpperLmt[0];
	for (int i = 1; i < bodyset.Num; i++) {
		bodyset.BodySetLowerLmt = bodyset.BodySetLowerLmt.min(bodyset.BodyLowerLmt[i]);
		bodyset.BodySetUpperLmt = bodyset.BodySetUpperLmt.max(bodyset.BodyUpperLmt[i]);
	}

	//ۺϿϵı߽磬ܵı߽磬չһ뾶
	lower = lower.min(bodyset.BodySetLowerLmt);
	upper = upper.max(bodyset.BodySetUpperLmt);

	const double padding = pt.maxRadius;
	pt.LowerLmt = lower - padding;
	pt.UpperLmt = upper + padding;

	// 2. 
	const double invMeshSize = 1.0 / pt.MeshSize;//ܸǿȺ͸ϵСռ仮Ϊȵ
	pt.MeshNum = ((pt.UpperLmt - pt.LowerLmt) * invMeshSize).ceiled();

	///******************** ֤ 5Ǹ ********************/
	//cout << "Mesh Numbers: " << pt.MeshNum.transpose() << endl;
	//if ((pt.MeshNum < 0).any()) {
	//	cerr << "ERROR: Negative mesh numbers detected" << endl;
	//}

	// 3. // get particle index
	int gridnum = 0;
	for (int i = 0; i < pt.Num; i++)
	{
		pt.Idx[i] = ((pt.Pos[i] - pt.LowerLmt) * invMeshSize).floored() + 1;
		GridIndex index = pt.Idx[i];
		if (grid.find(index) == grid.end()) {
			// 񲻴ڣ򴴽
			grid[index] = GridCell();
			gridnum++;
			grid[index].n = gridnum;
		}
		grid[index].addParticle(i);//indexӿi
	}

	// 3. 
	for (int ib = 0; ib < bodyset.Num; ib++)
	{
		for (int i = 0; i < bodyset.body[ib].Size; i++)//identification the range of the mesh
		{
			for (int j = 0; j < 3; j++)
			{
				bodyset.triIdx[ib][i][j] = int(floor((bodyset.body[ib].eAABB[i].min[j] - pt.LowerLmt[j]) / pt.MeshSize)) + 1;
				bodyset.triIdx[ib][i][j+3] = int(floor((bodyset.body[ib].eAABB[i].max[j] - pt.LowerLmt[j]) / pt.MeshSize)) + 1;
			}
			//cout << "min_x = " << bodyset.triIdx[ib][i][0] << ",max_x = " << bodyset.triIdx[ib][i][3] << endl;
			//cout << "min_y = " << bodyset.triIdx[ib][i][1] << ",max_y = " << bodyset.triIdx[ib][i][4] << endl;
			//cout << "min_z = " << bodyset.triIdx[ib][i][2] << ",max_z = " << bodyset.triIdx[ib][i][5] << endl;
			//cout << "-----------------------------" << endl;
			for (int k = bodyset.triIdx[ib][i][0]; k <= bodyset.triIdx[ib][i][3]; k++)
				for (int l = bodyset.triIdx[ib][i][1]; l <= bodyset.triIdx[ib][i][4]; l++)
					for (int m = bodyset.triIdx[ib][i][2]; m <= bodyset.triIdx[ib][i][5]; m++)
					{
						GridIndex index( k, l, m );
						if (grid.find(index) == grid.end()) {
							// 񲻴ڣ򴴽
							grid[index] = GridCell();
							gridnum++;
							grid[index].n = gridnum;
						}
						grid[index].addTri(ib,i);//indexi
					}
		}
	}

	//checkTree(pt, bodyset, grid);
	//cout << "tree finish..." << endl;
}

void CreatTree(PARTICLE& pt, unordered_map<GridIndex, GridCell>& grid)
{
	grid.clear();
	std::unordered_map<GridIndex, GridCell> empty_map;
	grid.swap(empty_map);  // ͷ grid ڴ

	/******************** ֤ 1Ƿ񳹵 ********************/
	cout << "Grid size after clear: " << grid.size() << " (should be 0)" << endl;
	if (!grid.empty()) {
		cerr << "ERROR: Grid was not properly cleared!" << endl;
	}

	//cout << grid.size() << endl;
	// get mesh limits
	const int numParticles = pt.Num;
	if (numParticles == 0) return;

	// 1. ռ߽
	vector3d lower = pt.Pos[0];
	vector3d upper = pt.Pos[0];

	for (int i = 1; i < numParticles; ++i) {
		lower = lower.min(pt.Pos[i]);
		upper = upper.max(pt.Pos[i]);
	}

	const double padding = pt.maxRadius; // 뾶
	pt.LowerLmt = lower - padding;
	pt.UpperLmt = upper + padding;

	/******************** ֤ 2߽ ********************/
	cout << "Particle Lower Limit: " << pt.LowerLmt<< endl;
	cout << "Particle Upper Limit: " << pt.UpperLmt << endl;
	if ((pt.LowerLmt.any_ge(pt.UpperLmt))) {
		cerr << "ERROR: Invalid particle boundaries (Lower >= Upper)" << endl;
	}

	// 2. 
	const double invMeshSize = 1.0 / pt.MeshSize;//ܸǿȺ͸ϵСռ仮Ϊȵ
	pt.MeshNum = ((pt.UpperLmt - pt.LowerLmt) * invMeshSize).ceiled();

	for (int i = 0; i < pt.Num; i++)
	{
		pt.Idx[i] = ((pt.Pos[i] - pt.LowerLmt) * invMeshSize).floored() + 1;
		for (int j = 0; j < 3; j++)
		{
			pt.Idx[i][j] = int(floor((pt.Pos[i][j] - pt.LowerLmt[j]) / pt.MeshSize)) + 1;
		}
		GridIndex index = pt.Idx[i];
		if (grid.find(index) == grid.end()) {
			// 񲻴ڣ򴴽
			grid[index] = GridCell();
		}
		grid[index].addParticle(i);//indexӿi
	}
}


bool InitList(Cp* part)
{
	part = (Cp*)malloc(sizeof(Cp));
	if (part == nullptr)
	{
		return false;
	}
	part->next = nullptr;
	part->flag = 0;
	return true;
}

void BackInsertList(Cp* part, int id, double P[][3], int row)
{
	Cp* End;
	End = part;
	while (part)
	{
		Cp* p = (Cp*)malloc(sizeof(Cp));
		p->idx = id;
		for (int i = 1; i < row; i++)
		{
			for (int j = 1; j < row; j++)
				p->deltP[i][j] = P[i][j];
		}
		End->next = p;
		End = p;
	}
	End->next = nullptr;

}

void checkTree(PARTICLE& pt, BODYSET& bodyset, unordered_map<GridIndex, GridCell>& grid)
{
	int BodyNum = 0, PtNum = 0;
	bool* ex = new bool[pt.Num];
	for (int i = 0; i < pt.Num; i++)	ex[i] = 1;
	bool* ex2 = new bool[pt.Num];
	for (int i = 0; i < bodyset.body[0].Size; i++)	ex2[i] = 1;

	for (auto it = grid.begin(); it != grid.end(); ++it) {
		// it->first Ǽit->second ֵ
		GridCell cell = it->second;
		GridIndex a = it->first;
		int num1 = cell.particleIds.size();
		int num2 = cell.triIds.size();

		if (num1 == 0&&num2==0)
			cout << " " << a<< "   fd"<<endl;

		PtNum += num1;
		BodyNum += num2;

		for (int i = 0; i < num1; i++)
		{
			ex[cell.particleIds[i]] = 1- ex[cell.particleIds[i]];
		}

		for (int i = 0; i < num2; i++)
		{
			ex2[cell.triIds[i]] = 0;
		}

	}
	for (int i = 0; i < pt.Num; i++)
	{
		if (ex[i])	cout << "" << i << "δ¼" << endl;
	}
	for (int i = 0; i < bodyset.body[0].Size; i++)
	{
		if (ex2[i])	cout << "" << i << "δ¼" << endl;
	}
}
