#pragma once
/************************************************************************************************/
/*Linear Solver for equations of type Ax = B with Multigrid method*/
/*	which makes use of correction steps that compute the error on a*/
/*	coarser grid.Thus ,lower frequency components of the error reappear*/
/*	as higher ones and allow for an efficient attenuation with standart*/
/*	iterative methods, such as SOR method implemented here.*/
/*	FAS srategy is used so advanced Multigrid strategies such as*/ 
/*	V-cycles and W-cycles may be applied in a hierarchical way.*/
/*	Maximal number of levels and iterations may be choosed.*/
/****************************************************************************************************/
#include "SparseMat.h"
#include "toolsKit.h"
#define DEFECT_FACTOR 10
#define V_CYCLE 1
#define W_CYCLE 2
class LinearSolver
{
public:	
	LinearSolver(int iterations,int levels,vector<float>*(*solver)(SparseMat<float> *A, vector<float> *X,vector<float> *B, float w, int numOfIterations));
	//LinearSolver(1,2,sparseMatSor);
	LinearSolver(void);
	~LinearSolver(void);
	/*cycleType = V_CYCLE or W_CYCLE*/
	/*solves equation of type Ax = B*/
	vector<float> * multigrid(SparseMat<float>& A, vector<float>& B,float w,vector<float>& X,	int preIterations,int postIterations,int cycleType);
	/*SOR function*/
	static vector<float> * sparseMatSor(SparseMat<float>* A, vector<float>* x,vector<float>* B, float w, int numOfIterations);
	

private:
	//maximal number of multigrid iterations
	int maxIterations;
	//maximal number of multigrid levels
	int maxLevels;
	vector<float>*(*solver)(SparseMat<float> *A, vector<float> *X,vector<float> *B, float w, int numOfIterations) ;
	void multigridIteration(SparseMat<float>& A, vector<float>& B,float w,vector<float>& X,	int preIterations,int postIterations,int cycleType);
};

