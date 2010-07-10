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
								  IplImage* psidashFS1,
								  IplImage* psidashFS2,
								  IplImage* u,
								  IplImage* v,
								  IplImage* du,
								  IplImage* dv,
								  double gamma,
								  double alpha,
								  double _ERROR_CONST,
								  int nInnerFPIterations);
	virtual ~constructMatrix_brox(void);
private:
	
};
