#pragma once
#include "cv.h"
#include "coarse2FineCompute.h"
#include "toolsKit.h"



class constructMatrix_brox
{
public:
	constructMatrix_brox(void);
	static void constructMatrix_b(IplImage* Ikx,IplImage* Iky,IplImage* Ikz,IplImage* Ixx,IplImage* Ixy,IplImage* Iyy,IplImage* Ixz,
								  IplImage* Iyz,IplImage* psidash,IplImage* psidashFS1,IplImage* psidashFS2,IplImage*  u,IplImage*  v,double gamma,int _ERROR_CONST );
	virtual ~constructMatrix_brox(void);
private:
	
};
