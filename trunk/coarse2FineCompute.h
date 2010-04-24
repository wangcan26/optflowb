#pragma once


#include "highgui.h" 
#include "cv.h"
#include "GaussPyramid.h"
#include "toolsKit.h"

class coarse2FineCompute
{


public:
	coarse2FineCompute(void);
public:
	virtual ~coarse2FineCompute(void);
	void coarse2FineCompute::Coarse2FineFlow(IplImage* vx, IplImage* vy, IplImage &warpI2,const IplImage &Im1, const IplImage &Im2, double alpha, double ratio, int minWidth, 
																	 int nOuterFPIterations, int nInnerFPIterations, int nCGIterations);
	IplImage* coarse2FineCompute::LaplaceCompute(IplImage* input,IplImage* input2);
	void coarse2FineCompute::SmoothFlowPDE(const IplImage* Im1, const IplImage* Im2, IplImage* warpIm2, IplImage* u, IplImage* v, 
																    double alpha, int nOuterFPIterations, int nInnerFPIterations, int nCGIterations);


	//temp
	void coarse2FineCompute::opt_flow_lk();
};
