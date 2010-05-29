#pragma once


#include "highgui.h" 
#include "cv.h"
#include "GaussPyramid.h"
#include "toolsKit.h"
#include "flowUV.h"
#include "IplImageIterator.h"
class coarse2FineCompute
{
public:	
	coarse2FineCompute(int imageDepth,double error);
	void Coarse2FineFlow(IplImage* vx, 
						 IplImage* vy, 
						 IplImage &warpI2,
						 const IplImage &Im1, 
						 const IplImage &Im2, 
						 double alpha,
						 double gamma,
						 double ratio, 
						 int minWidth,
						 int nOuterFPIterations, 
						 int nInnerFPIterations, 
						 int nCGIterations);
	virtual ~coarse2FineCompute(void);
private:
	double _ERROR_CONST;
	int _imageDepth;
	IplImage* LaplaceCompute(IplImage* input,IplImage* input2);
	IplImage* coarse2FineCompute::createWarp(IplImage*WarpImage2,IplImage* img1,IplImage* img2,IplImage* vx,IplImage* vy);
	
	void coarse2FineCompute::constructMatrix_brox(IplImage* Ikx,IplImage* Iky,IplImage* Ikz,IplImage* Ixx,IplImage* Ixy,IplImage* Iyy,IplImage* Ixz,
												  IplImage* Iyz,IplImage* psidash,IplImage* psidashFS1,IplImage* psidashFS2,IplImage*  u,IplImage*  v,double gamma );
	flowUV* SmoothFlowPDE(const IplImage* Im1, 
						   const IplImage* Im2, 
						   IplImage* warpIm2, 
						   IplImage* du, 
						   IplImage* dv, 
						   double alpha,
						   double gamma,
						   int nOuterFPIterations, 
						   int nInnerFPIterations, 
						   int nCGIterations);

	void computePsidashFS_brox(IplImage* iterU,IplImage* iterV,int width,int height,int channels,flowUV* UV);

	//temp
	void opt_flow_lk();

};
