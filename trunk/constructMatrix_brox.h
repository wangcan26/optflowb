#pragma once
#include "cv.h"
#include "coarse2FineCompute.h"
#include "toolsKit.h"
#include "SparseToolKit.h"




class constructMatrix_brox
{
public:
	constructMatrix_brox(void);
	static vector<float> * constructMatrix_b(IplImage* Ikx,
								  IplImage* Iky,
								  IplImage* Ikz,
								  IplImage* Ixx,
								  IplImage* Ixy,
								  IplImage* Iyy,
								  IplImage* Ixz,
								  IplImage* Iyz,							
								  flowUV* UV,
								  IplImage* du,
								  IplImage* dv,
								  SparseMat<float> * A,
								  vector<float> * B,
								  vector<float>* dUdV,
								  double gamma,
								  double alpha,
								  double _ERROR_CONST,
								  int nInnerFPIterations);

	static void computePsidashFS_brox(IplImage* iterU,IplImage* iterV,int width,int height,int channels,flowUV* UV,double _ERROR_CONST);



	virtual ~constructMatrix_brox(void);
private:
	
};
