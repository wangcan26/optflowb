#pragma once

#include "cv.h"
#include "coarse2FineCompute.h"
#include "toolsKit.h"
#include "SparseToolKit.h"

class flowMatrix
{
public:
	flowMatrix(void);
	virtual ~flowMatrix(void);
	/*
static void computePdfSum(IplImage* pdfSum,  IplImage* psidashFS1, IplImage* psidashFS2,
				   IplImage* fs1_3222,IplImage* fs1_122ht22,IplImage* fs2_2232,IplImage* fs2_22122wt);

	static void computePsidashFS_brox(IplImage* iterU,IplImage* iterV,int width,int height,int channels,flowUV* UV,double _ERROR_CONST);


	static void psiDash(IplImage* src,IplImage* ans,  double _ERROR_CONST);

	static void computeVectBComponents(IplImage* pdfaltSumXX,IplImage* fs1_3222,IplImage* fs1_122ht22,
							IplImage* fs2_2232,IplImage* fs2_22122wt,IplImage* UorV_Org);*/

	static vector<float> * constructMatrix(IplImage* Ix, IplImage* Iy, 	IplImage* Iz,						
								  flowUV* UV, IplImage* du, IplImage* dv,vector<float>* dUdV,
								    SparseMat<float> *  A,vector<float> * B,
								  double alpha,
								  double _ERROR_CONST);

	
};

