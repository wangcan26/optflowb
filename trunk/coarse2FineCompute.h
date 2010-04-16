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
};
