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
	coarse2FineCompute(int imageDepth);
	void Coarse2FineFlow(IplImage* vx, 
						 IplImage* vy, 
						 IplImage &warpI2,
						 const IplImage &Im1, 
						 const IplImage &Im2, 
						 double alpha, 
						 double ratio, 
						 int minWidth,
						 int nOuterFPIterations, 
						 int nInnerFPIterations, 
						 int nCGIterations);
	virtual ~coarse2FineCompute(void);
private:
	int _imageDepth;
	IplImage* LaplaceCompute(IplImage* input,IplImage* input2);
	
	void SmoothFlowPDE(const IplImage* Im1, 
					   const IplImage* Im2, 
					   IplImage* warpIm2, 
					   IplImage* u, 
					   IplImage* v, 
					   double alpha, 
					   int nOuterFPIterations, 
					   int nInnerFPIterations, 
					   int nCGIterations);

	flowUV* SmoothFlowPDE2(const IplImage* Im1, 
									   const IplImage* Im2, 
									   IplImage* warpIm2, 
									   IplImage* du, 
									   IplImage* dv, 
									   double alpha, 
									   int nOuterFPIterations, 
									   int nInnerFPIterations, 
									   int nCGIterations);

	void computePsidashFS_brox(IplImage* iterU,IplImage* iterV,int width,int height,int channels,flowUV* UV);

	//temp
	void opt_flow_lk();

};
