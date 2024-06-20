/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Q4.h"

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

//	Constructor
CQ4::CQ4()
{
	intnum = 0;
	NEN_ = 4;	// Each element has 2 nodes
	nodes_ = new CNode*[NEN_];
    
    ND_ = 8;
    LocationMatrix_ = new unsigned int[ND_];

	ElementMaterial_ = nullptr;

	for (int i = 0; i < 8; i++)
	{
		Nixy[i] = 0;
	}

}

//	Desconstructor
CQ4::~CQ4()
{
}

//	Read element data from stream Input
bool CQ4::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList)
{
	unsigned int MSet, INT;	// Material property set number
	unsigned int N1, N2, N3, N4;	// Left node number and right node number

	Input >> N1 >> N2 >> N3 >> N4 >> MSet>> INT;
    ElementMaterial_ = dynamic_cast<Q4Material*>(MaterialSets) + MSet - 1;
	nodes_[0] = &NodeList[N1 - 1];
	nodes_[1] = &NodeList[N2 - 1];
	nodes_[2] = &NodeList[N3 - 1];
	nodes_[3] = &NodeList[N4 - 1];
	intnum = INT;

	return true;
}

//	Write element data to stream
void CQ4::Write(COutputter& output)
{
	output << setw(11) << nodes_[0]->NodeNumber
		   << setw(9) << nodes_[1]->NodeNumber
		   << setw(9) << nodes_[2]->NodeNumber
		   << setw(9) << nodes_[3]->NodeNumber
		   << setw(12) << ElementMaterial_->nset << endl;
}

//	Calculate element stiffness matrix 
//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element
void CQ4::ElementStiffness(double* Matrix)
{
	clear(Matrix, SizeOfStiffnessMatrix());
	double x1 = nodes_[0]->XYZ[0];
	double y1 = nodes_[0]->XYZ[1];
	double x2 = nodes_[1]->XYZ[0];
	double y2 = nodes_[1]->XYZ[1];
	double x3 = nodes_[2]->XYZ[0];
	double y3 = nodes_[2]->XYZ[1];
	double x4 = nodes_[3]->XYZ[0];
	double y4 = nodes_[3]->XYZ[1];
	// det(jacobian) 中心点
	double JE00 = ((x2 + x3 - x1 - x4) * (y4 + y3 - y1 - y2) - (x4 + x3 - x1 - x2) * (y2 + y3 - y1 - y4)) / 16;
	double JEinv00[4];
	JEinv00[0] = 0.25 * (y4 + y3 - y2 - y1) / JE00;//jacobian逆矩阵左上
	JEinv00[1] = 0.25 * (y1 + y4 - y2 - y3) / JE00;//jacobian逆矩阵右上
	JEinv00[2] = 0.25 * (x1 + x2 - x3 - x4) / JE00;//jacobian逆矩阵左下
	JEinv00[3] = 0.25 * (x2 + x3 - x1 - x4) / JE00;//jacobian逆矩阵右下
	
	//形函数导数N1x,N2x,N3x,N4x,N1y,N2y,N3y,N4y
	Nixy[0] = -JEinv00[0] - JEinv00[1];
	Nixy[1] = JEinv00[0] - JEinv00[1];
	Nixy[2] = -Nixy[0];
	Nixy[3] = -Nixy[1];
	Nixy[4] = -JEinv00[2] - JEinv00[3];
	Nixy[5] = JEinv00[2] - JEinv00[3];
	Nixy[6] = -Nixy[4];
	Nixy[7] = -Nixy[5];
	for (int i = 0; i < 8; i++)
		Nixy[i] = Nixy[i] / 4;
	
	Q4Material* material_ = dynamic_cast<Q4Material*>(ElementMaterial_);
	double P1 = material_->Poisson;
	double P2 = (1 - P1) / 2;
	double P3 = (1 + P1) / 2;

	if (intnum == 1)//完全积分
	{
		for (int j = 0; j < 36; j++) //初始化
			Matrix[j] = 0;
		//2*2 Gauss积分点
		double ksi[4] = { -0.5773502692, 0.5773502692, -0.5773502692, 0.5773502692 };
		double eta[4] = { 0.5773502692, 0.5773502692, -0.5773502692, -0.5773502692 };
		for (int i = 0; i < 4; i++)
		{
			// 求解Jacobian矩阵行列式和逆
			double JE[4];
			double JEinv[4];
			JE[0] = (x2 - x1) * (1 - eta[i]) - (x4 - x3) * (1 + eta[i]);
			JE[1] = (y2 - y1) * (1 - eta[i]) - (y4 - y3) * (1 + eta[i]);
			JE[2] = (x4 - x1) * (1 - ksi[i]) + (x3 - x2) * (1 + ksi[i]);
			JE[3] = (y4 - y1) * (1 - ksi[i]) + (y3 - y2) * (1 + ksi[i]);
			double detJE = (JE[0] * JE[3] - JE[1] * JE[2]) / 16;
			JEinv[0] = 0.25 * JE[3] / detJE;
			JEinv[1] = -0.25 * JE[1] / detJE;
			JEinv[2] = -0.25 * JE[2] / detJE;
			JEinv[3] = 0.25 * JE[0] / detJE;
			// Ni类似于Nixy,但不是全局的
			double Ni[8];
			Ni[0] = (JEinv[0] * (eta[i] - 1) + JEinv[1] * (ksi[i] - 1)) / 4;
			Ni[1] = (JEinv[0] * (1 - eta[i]) + JEinv[1] * (-ksi[i] - 1)) / 4;
			Ni[2] = (JEinv[0] * (eta[i] + 1) + JEinv[1] * (ksi[i] + 1)) / 4;
			Ni[3] = (JEinv[0] * (-eta[i] - 1) + JEinv[1] * (1 - ksi[i])) / 4;
			Ni[4] = (JEinv[2] * (eta[i] - 1) + JEinv[3] * (ksi[i] - 1)) / 4;
			Ni[5] = (JEinv[2] * (1 - eta[i]) + JEinv[3] * (-ksi[i] - 1)) / 4;
			Ni[6] = (JEinv[2] * (eta[i] + 1) + JEinv[3] * (ksi[i] + 1)) / 4;
			Ni[7] = (JEinv[2] * (-eta[i] - 1) + JEinv[3] * (1 - ksi[i])) / 4;

			//Gauss 积分前面的系数
			double correctionfactor = detJE * material_->thickness * material_->E / (1 - P1 * P1);//完全积分系数为1

			Matrix[0] += (Ni[0] * Ni[0] + P2 * Ni[4] * Ni[4]) * correctionfactor;
			Matrix[1] += (Ni[4] * Ni[4] + P2 * Ni[0] * Ni[0]) * correctionfactor;
			Matrix[2] += (P3 * Ni[0] * Ni[4]) * correctionfactor;
			Matrix[3] += (Ni[1] * Ni[1] + P2 * Ni[5] * Ni[5]) * correctionfactor;
			Matrix[4] += (P1 * Ni[1] * Ni[4] + P2 * Ni[0] * Ni[5]) * correctionfactor;
			Matrix[5] += (Ni[0] * Ni[1] + P2 * Ni[4] * Ni[5]) * correctionfactor;
			Matrix[6] += (Ni[5] * Ni[5] + P2 * Ni[1] * Ni[1]) * correctionfactor;
			Matrix[7] += (P3 * Ni[1] * Ni[5]) * correctionfactor;
			Matrix[8] += (Ni[4] * Ni[5] + P2 * Ni[0] * Ni[1]) * correctionfactor;
			Matrix[9] += (P1 * Ni[0] * Ni[5] + P2 * Ni[4] * Ni[1]) * correctionfactor;
			Matrix[10] += (Ni[2] * Ni[2] + P2 * Ni[6] * Ni[6]) * correctionfactor;
			Matrix[11] += (P1 * Ni[2] * Ni[5] + P2 * Ni[1] * Ni[6]) * correctionfactor;
			Matrix[12] += (Ni[1] * Ni[2] + P2 * Ni[5] * Ni[6]) * correctionfactor;
			Matrix[13] += (P1 * Ni[4] * Ni[2] + P2 * Ni[0] * Ni[6]) * correctionfactor;
			Matrix[14] += (Ni[0] * Ni[2] + P2 * Ni[4] * Ni[6]) * correctionfactor;
			Matrix[15] += (Ni[6] * Ni[6] + P2 * Ni[2] * Ni[2]) * correctionfactor;
			Matrix[16] += (P3 * Ni[2] * Ni[6]) * correctionfactor;
			Matrix[17] += (Ni[5] * Ni[6] + P2 * Ni[1] * Ni[2]) * correctionfactor;
			Matrix[18] += (P1 * Ni[1] * Ni[6] + P2 * Ni[5] * Ni[2]) * correctionfactor;
			Matrix[19] += (Ni[4] * Ni[6] + P2 * Ni[2] * Ni[0]) * correctionfactor;
			Matrix[20] += (P1 * Ni[0] * Ni[6] + P2 * Ni[4] * Ni[2]) * correctionfactor;
			Matrix[21] += (Ni[3] * Ni[3] + P2 * Ni[7] * Ni[7]) * correctionfactor;
			Matrix[22] += (P1 * Ni[3] * Ni[6] + P2 * Ni[2] * Ni[7]) * correctionfactor;
			Matrix[23] += (Ni[2] * Ni[3] + P2 * Ni[6] * Ni[7]) * correctionfactor;
			Matrix[24] += (P1 * Ni[3] * Ni[5] + P2 * Ni[1] * Ni[7]) * correctionfactor;
			Matrix[25] += (Ni[1] * Ni[3] + P2 * Ni[5] * Ni[7]) * correctionfactor;
			Matrix[26] += (P1 * Ni[4] * Ni[3] + P2 * Ni[0] * Ni[7]) * correctionfactor;
			Matrix[27] += (Ni[0] * Ni[3] + P2 * Ni[4] * Ni[7]) * correctionfactor;
			Matrix[28] += (Ni[7] * Ni[7] + P2 * Ni[3] * Ni[3]) * correctionfactor;
			Matrix[29] += (P3 * Ni[3] * Ni[7]) * correctionfactor;
			Matrix[30] += (Ni[6] * Ni[7] + P2 * Ni[2] * Ni[3]) * correctionfactor;
			Matrix[31] += (P1 * Ni[2] * Ni[7] + P2 * Ni[6] * Ni[3]) * correctionfactor;
			Matrix[32] += (Ni[5] * Ni[7] + P2 * Ni[1] * Ni[3]) * correctionfactor;
			Matrix[33] += (P1 * Ni[1] * Ni[7] + P2 * Ni[3] * Ni[5]) * correctionfactor;
			Matrix[34] += (Ni[4] * Ni[7] + P2 * Ni[3] * Ni[0]) * correctionfactor;
			Matrix[35] += (P1 * Ni[0] * Ni[7] + P2 * Ni[4] * Ni[3]) * correctionfactor;
		}
	}
	else //减缩积分
	{
		if (intnum != 0)
			std::cerr << "Element intnum not correct, Default intnum = 0, Calculation continue" << std::endl;
				
		double correctionfactor = 4 * JE00 * material_->thickness * material_->E / (1 - P1 * P1);//减缩积分系数为4

		Matrix[0] = Nixy[0] * Nixy[0] + P2 * Nixy[4] * Nixy[4];
		Matrix[1] = Nixy[4] * Nixy[4] + P2 * Nixy[0] * Nixy[0];
		Matrix[2] = P3 * Nixy[0] * Nixy[4];
		Matrix[3] = Nixy[1] * Nixy[1] + P2 * Nixy[5] * Nixy[5];
		Matrix[4] = P1 * Nixy[1] * Nixy[4] + P2 * Nixy[0] * Nixy[5];
		Matrix[5] = Nixy[0] * Nixy[1] + P2 * Nixy[4] * Nixy[5];
		Matrix[6] = Nixy[5] * Nixy[5] + P2 * Nixy[1] * Nixy[1];
		Matrix[7] = P3 * Nixy[1] * Nixy[5];
		Matrix[8] = Nixy[4] * Nixy[5] + P2 * Nixy[0] * Nixy[1];
		Matrix[9] = P1 * Nixy[0] * Nixy[5] + P2 * Nixy[4] * Nixy[1];
		Matrix[10] = Nixy[2] * Nixy[2] + P2 * Nixy[6] * Nixy[6];
		Matrix[11] = P1 * Nixy[2] * Nixy[5] + P2 * Nixy[1] * Nixy[6];
		Matrix[12] = Nixy[1] * Nixy[2] + P2 * Nixy[5] * Nixy[6];
		Matrix[13] = P1 * Nixy[4] * Nixy[2] + P2 * Nixy[0] * Nixy[6];
		Matrix[14] = Nixy[0] * Nixy[2] + P2 * Nixy[4] * Nixy[6];
		Matrix[15] = Nixy[6] * Nixy[6] + P2 * Nixy[2] * Nixy[2];
		Matrix[16] = P3 * Nixy[2] * Nixy[6];
		Matrix[17] = Nixy[5] * Nixy[6] + P2 * Nixy[1] * Nixy[2];
		Matrix[18] = P1 * Nixy[1] * Nixy[6] + P2 * Nixy[5] * Nixy[2];
		Matrix[19] = Nixy[4] * Nixy[6] + P2 * Nixy[2] * Nixy[0];
		Matrix[20] = P1 * Nixy[0] * Nixy[6] + P2 * Nixy[4] * Nixy[2];
		Matrix[21] = Nixy[3] * Nixy[3] + P2 * Nixy[7] * Nixy[7];
		Matrix[22] = P1 * Nixy[3] * Nixy[6] + P2 * Nixy[2] * Nixy[7];
		Matrix[23] = Nixy[2] * Nixy[3] + P2 * Nixy[6] * Nixy[7];
		Matrix[24] = P1 * Nixy[3] * Nixy[5] + P2 * Nixy[1] * Nixy[7];
		Matrix[25] = Nixy[1] * Nixy[3] + P2 * Nixy[5] * Nixy[7];
		Matrix[26] = P1 * Nixy[4] * Nixy[3] + P2 * Nixy[0] * Nixy[7];
		Matrix[27] = Nixy[0] * Nixy[3] + P2 * Nixy[4] * Nixy[7];
		Matrix[28] = Nixy[7] * Nixy[7] + P2 * Nixy[3] * Nixy[3];
		Matrix[29] = P3 * Nixy[3] * Nixy[7];
		Matrix[30] = Nixy[6] * Nixy[7] + P2 * Nixy[2] * Nixy[3];
		Matrix[31] = P1 * Nixy[2] * Nixy[7] + P2 * Nixy[6] * Nixy[3];
		Matrix[32] = Nixy[5] * Nixy[7] + P2 * Nixy[1] * Nixy[3];
		Matrix[33] = P1 * Nixy[1] * Nixy[7] + P2 * Nixy[3] * Nixy[5];
		Matrix[34] = Nixy[4] * Nixy[7] + P2 * Nixy[3] * Nixy[0];
		Matrix[35] = P1 * Nixy[0] * Nixy[7] + P2 * Nixy[4] * Nixy[3];

		for (int i = 0; i < 36; i++)
			Matrix[i] = Matrix[i] * correctionfactor;
	}
}

//	Calculate element stress 
void CQ4::ElementStress(double* stress, double* Displacement)
{
	Q4Material* material_ = dynamic_cast<Q4Material*>(ElementMaterial_);	// Pointer to material of the element
	// Q4单元4个结点的xy位移
	double d[8];
	for (int i = 0; i < 8; i++)
	{
		if (LocationMatrix_[i] == 0)
			d[i] = 0.0;
		else
			d[i] = Displacement[LocationMatrix_[i] - 1];//-1参考CDomain::CalculateEquationNumber
	}
	double P1 = material_->Poisson;
	//最佳应力点应力sigma_x,sigma_y,gamma_xy
	stress[0] = material_->E / (1 - P1 * P1) * (Nixy[0] * d[0] + P1 * Nixy[4] * d[1] + Nixy[1] * d[2] + P1 * Nixy[5] * d[3] + Nixy[2] * d[4] + P1 * Nixy[6] * d[5] + Nixy[3] * d[6] + P1 * Nixy[7] * d[7]);
	stress[1] = material_->E / (1 - P1 * P1) * (P1 * Nixy[0] * d[0] + Nixy[4] * d[1] + P1 * Nixy[1] * d[2] + Nixy[5] * d[3] + P1 * Nixy[2] * d[4] + Nixy[6] * d[5] + P1 * Nixy[3] * d[6] + Nixy[7] * d[7]);
	stress[2] = material_->E / (2 + 2 * P1) * (Nixy[4] * d[0] + Nixy[0] * d[1] + Nixy[5] * d[2] + Nixy[1] * d[3] + Nixy[6] * d[4] + Nixy[2] * d[5] + Nixy[7] * d[6] + Nixy[3] * d[7]);
	//cout << "Displacement[" << i << "]=" << Displacement[i] << "\n" << endl;
}
