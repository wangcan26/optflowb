#ifndef UTILS_H_
#define UTILS_H_

void printVector(vector<float> vec);
float norm(const vector<float>& vec);
void multByVec(SparseMat<float>& mat,const vector<float>& vec,vector<float>& result);
void residual(SparseMat<float>& A,const vector<float>& B,const vector<float>& X,vector<float>& residualVec);
float residualNorm(SparseMat<float>& A,const vector<float>& B,vector<float>& X);

bool achievedTolerance(SparseMat<float>& A,const vector<float>& B,const vector<float>& X,float tol);
//sign = 0 add, sign = 1 sub
void sumVectors(const vector<float>& vec1,const vector<float>& vec2,vector<float>& vecResult,int sign);
#endif /*UTILS_H_*/
