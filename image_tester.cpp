#include <iostream>
#include <stdio.h>
#include "highgui.h"
#include "GaussPyramid.h"
#include "coarse2FineCompute.h"
#include "SparseToolKit.h"
#include "optical_flow_demo.h"
#include <ctime>
#include "improvements.h"
#include "FlowUtils.h"
using namespace std;

int main (int argc,char** argv) 
{ 
	double error_const=0.001;
	bool useTextureDecomposition=false;
	float textureAlpha=0.95f;
	//IPL_DEPTH_32F IPL_DEPTH_8U
	coarse2FineCompute coarse2fComp(IPL_DEPTH_32F,error_const);
	double ratio=0.75;
	int minWidth=3;
	int outerIter=3;
	int innerIter=50;
	double alpha = 3 ; // Global smoothness variable.
	double gamma = 0 ; // Global weight for derivatives.


	std::clock_t start;
	double diff;

	IplImage* img1= cvLoadImage(argv[1],CV_LOAD_IMAGE_COLOR);//zero is for grayscale  CV_LOAD_IMAGE_GRAYSCALE
	IplImage* img2= cvLoadImage(argv[2],CV_LOAD_IMAGE_COLOR); //1 is for color CV_LOAD_IMAGE_COLOR


	IplImage *img1_32 = cvCreateImage(cvSize(img1->width, img1->height), coarse2fComp._imageDepth, img1->nChannels);
	IplImage *img2_32 = cvCreateImage(cvSize(img1->width, img1->height), coarse2fComp._imageDepth, img1->nChannels);

	cvZero(img1_32);
	cvZero(img2_32);
	//upscale from char to float
	cvConvertScale(img1, img1_32, 1.0/255);
	cvConvertScale(img2, img2_32, 1.0/255);


	//cvSmooth(img1_32,img1_32,CV_GAUSSIAN);
	//cvSmooth(img2_32,img2_32,CV_GAUSSIAN);

	IplImage *img1_32g = cvCreateImage(cvSize(img1->width, img1->height), coarse2fComp._imageDepth, 1);
	IplImage *img2_32g = cvCreateImage(cvSize(img1->width, img1->height), coarse2fComp._imageDepth, 1);

	cvCvtColor( img1_32, img1_32g, CV_BGR2GRAY );
	cvCvtColor( img2_32, img2_32g, CV_BGR2GRAY );

	// Preprocess the gray images.
	if (useTextureDecomposition)
	{
		IplImage *texture1=cvCreateImage(cvSize(img1_32g->width,img1_32g->height),IPL_DEPTH_32F,1);
		IplImage *texture2=cvCreateImage(cvSize(img2_32g->width,img2_32g->height),IPL_DEPTH_32F,1);
		structure_texture_decomposition_rof(img1_32g,img1_32g,texture1,texture2,NULL,NULL,1.0f/8.0f,100,textureAlpha);
		cvCopy(texture1,img1_32g);
		cvCopy(texture2,img2_32g);
	}

	//cvNormalize(img1_32g,img1_32g,127,0,CV_MINMAX); //CV_MINMAX
	//cvNormalize(img2_32g,img2_32g,127,0,CV_MINMAX); 
	/*int t=GetTickCount();
	img1_32g=medianFilter(img1_32g,5);
	cFlowUtils::cvShowManyImages("img1,img2 ",1,img1_32g);	
	cout << (double)(GetTickCount()-t)/1000.0 << endl;
	cvWaitKey(0);
	return 0;*/
	
	//read GT
	CvMat* velx;
	CvMat* vely;
	try {
		
		char *filename = "Media\\flow10.flo";
		cFlowUtils::ReadFlowFile(filename,&velx,&vely);
				
	}
  catch (CError &err) {
	fprintf(stderr, err.message);
	fprintf(stderr, "\n");
	exit(1);
    }

	
	cFlowUtils::DrawFlow(vely,velx);
	


	//test wrap


	start = std::clock();
	flowUV* UV=coarse2fComp.Coarse2FineFlow(img1_32, img2_32, 
											alpha,gamma,
											ratio,minWidth, 
											outerIter,innerIter,
											velx,vely);
	diff = ( std::clock() - start ) / (double)CLOCKS_PER_SEC;
	std::cout<<"BROX pyramid alg. took "<< diff <<" secs"<<endl;

	 cvReleaseImage(&img1_32g);
	 cvReleaseImage(&img2_32g);
	 cvReleaseImage(&img1_32);
	 cvReleaseImage(&img2_32);
	 cvReleaseImage(&img1);
	 cvReleaseImage(&img2);

	 cFlowUtils::DrawFlow(UV->getU(),UV->getV());


	cout<<"fin"<<endl;
} 


