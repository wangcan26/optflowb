#pragma once
#include <cv.h>
#include <highgui.h>
#include <stdio.h>
#include <stdarg.h>
#include "SparseMat.h"
#include <vector>

class SparseToolKit{
public:
	static void printSparseMat(CvSparseMat* mat);
	/* return new Sparse Matrix with the elemets below the main diagonal of mat
		does not include the main diagonal*/
	static CvSparseMat * lowerTriangle(CvSparseMat * mat);
	/* return new Sparse Matrix with the elemets above the main diagonal of mat
		does not include the main diagonal*/
	static CvSparseMat * upperTriangle(CvSparseMat * mat);


	static vector<float>* SOR(SparseMat<float> A, vector<float> x,vector<float> B, double w, int numOfIterations); 

	};