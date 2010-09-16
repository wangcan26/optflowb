#ifndef FAS_H_
#define FAS_H_
#include "SparseMat.h"
#include "Utils.h"
void restrictionX(const vector<float>& Xh,vector<float>& X2h);
void restrictionV(const vector<float>& Xh,vector<float>& X2h);
void restrictionM(SparseMat<float>& Ah,SparseMat<float>& A2h);
//build fine grid Xh of X2h using , result is saved in Xh
void prolongation(const vector<float>& X2h,vector<float>& Xh);
void prolongationD(vector<float>& X2h,vector<float>& Xh);
void coarseGridOperator(SparseMat<float>& Ah,const vector<float>& X2h,vector<float>& B2h);

#endif /*FAS_H_*/
