#pragma once


#include "highgui.h" 
#include "cv.h"
#include "GaussPyramid.h"
#include "flowUV.h"
#include "IplImageIterator.h"
#include "constructMatrix_brox.h"
#include <math.h>
#include "VectorToolsKit.h"
#include "SparseMat.h"
#include "LinearSolver.h"
#include "toolsKit.h"
typedef toolsKit::vectorTools vtools;


class coarse2FineCompute
{
public:
	bool useMediaFiltering;
	double _ERROR_CONST;
	int _imageDepth;
	coarse2FineCompute(int imageDepth,double error,bool useMedianFiltering=false);
	flowUV* Coarse2FineFlow( const IplImage* Im1, 
							 const IplImage* Im2, 
							 double alpha,
							 double gamma,
							 double ratio, 
							 int minWidth,
							 int nOuterFPIterations, 
							 int nInnerFPIterations,
							  CvMat* velx, CvMat* vely);
						
	virtual ~coarse2FineCompute(void);
	static IplImage** meshgrid(int cols, int rows);
private:
	

	IplImage* LaplaceCompute(IplImage* input,IplImage* input2);
	IplImage* createWarp(IplImage*WarpImage2,IplImage* img1,IplImage* img2,IplImage* vx,IplImage* vy);
	IplImage * warpLayer(IplImage* I, IplImage* u, IplImage * v);
	IplImage* RGBwarp(IplImage* img, IplImage* u, IplImage* v);
	
	void SmoothFlowPDE( IplImage* Im1, 
						   IplImage* Im2, 					
						   double alpha,
						   double gamma,
						   int nOuterFPIterations, 
						   int nInnerFPIterations,
						   flowUV* UV);
	void computePsidashFS_brox(IplImage* iterU,IplImage* iterV,int width,int height,int channels,flowUV* UV);

};
