/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#pragma once

#include "Element.h"

using namespace std;

//! Q4 element class
class CQ4 : public CElement
{
public:
//0 = reduced intergral 1= complete intergral
	unsigned int intnum;

//	超收敛点（单元中心处）N1x,N2x,N3x,N4x,N1y,N2y,N3y,N4y
	double Nixy[8];

//!	Constructor
	CQ4();

//!	Desconstructor
	~CQ4();

//! 重写
	void GenerateLocationMatrix()override {
		unsigned int i = 0;
		for (unsigned int N = 0; N < NEN_; N++)
			for (unsigned int D = 0; D < 2; D++)
				LocationMatrix_[i++] = nodes_[N]->bcode[D];
	}


//!	Read element data from stream Input
	virtual bool Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList);//&表示引用

//!	Write element data to stream
	virtual void Write(COutputter& output);

//!	Calculate element stiffness matrix
	virtual void ElementStiffness(double* Matrix);

//!	Calculate element stress
	virtual void ElementStress(double* stress, double* Displacement);
};
