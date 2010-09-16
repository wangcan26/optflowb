#ifndef SOR_H_
#define SOR_H_
#include "SparseMat.h"
#include "Utils.h"
#include <math.h>

vector<float>* MY_SOR(SparseMat<float> A, vector<float> X,vector<float> B, float w, int numOfIterations);


#endif /*SOR_H_*/
