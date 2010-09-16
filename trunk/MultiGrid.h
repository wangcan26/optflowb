#ifndef MULTIGRID_H_
#define MULTIGRID_H_

#include "FAS.h"
#define MAX_ITERATIONS 200
#define DEFECT_FACTOR 10
#define V_CYCLE 1
#define W_CYCLE 2
//tested on w=1, cycleType = V_CYCLE, preIterations=1 , postIterations=1
bool multigrid(vector<float>* (*solver)(SparseMat<float> A, vector<float> X,vector<float> B, float w, int numOfIterations),
				SparseMat<float>& A,const vector<float>& B,float w,vector<float>& X,
				int preIterations,int postIterations,float tol,int cycleType);
int getMULTIGRIDIterations();
#endif /*MULTIGRID_H_*/
