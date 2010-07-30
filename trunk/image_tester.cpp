#include <iostream>
#include <stdio.h>
#include "highgui.h"
#include "GaussPyramid.h"
#include "coarse2FineCompute.h"
#include "SparseToolKit.h"
#include "optical_flow_demo.h"
#include <ctime>

using namespace std;

int main (int argc,char** argv) 
{ 
	double error_const=0.001;
	//IPL_DEPTH_32F IPL_DEPTH_8U
	coarse2FineCompute coarse2fComp(IPL_DEPTH_32F,error_const);
	double ratio=0.75;
	int minWidth=10;
	int outerIter=3;
	int innerIter=100;
	double alpha = 2 ; // Global smoothness variable.
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


	//cvNormalize(img1_32g,img1_32g,127,0,CV_MINMAX); //CV_MINMAX
	//cvNormalize(img2_32g,img2_32g,127,0,CV_MINMAX); 

	
	const IplImage *img1_33_file=NULL;
	const IplImage *img2_33_file=NULL;
	img1_33_file=toolsKit::IplFromFile("c:\\a\\urban1r.txt");
//	img1_33_file=toolsKit::IplFromFile("c:\\a\\1_15_15.txt");
//	img1_33_file=toolsKit::IplFromFile("c:\\a\\Urban3_1s.txt");	
	//img2_33_file=toolsKit::IplFromFile("c:\\a\\2_15_15.txt");
	img2_33_file=toolsKit::IplFromFile("c:\\a\\urban2r.txt");
//	img2_33_file=toolsKit::IplFromFile("c:\\a\\Urban3_2s.txt");
		
	toolsKit::cvShowManyImages("img1,img2 ",2,img1_32,img2_32);	
	cvWaitKey(1);
	
	start = std::clock();
	flowUV* UV=coarse2fComp.Coarse2FineFlow(img1_32, img2_32, 
											alpha,gamma,
											ratio,minWidth, 
											outerIter,innerIter);
	diff = ( std::clock() - start ) / (double)CLOCKS_PER_SEC;
	std::cout<<"BROX pyramid alg. took "<< diff <<" secs"<<endl;

	 cvReleaseImage(&img1_32g);
	 cvReleaseImage(&img2_32g);
	 cvReleaseImage(&img1_32);
	 cvReleaseImage(&img2_32);
	 cvReleaseImage(&img1);
	 cvReleaseImage(&img2);

	 toolsKit::drawFlow(UV->getU(),UV->getV(),0);


	cout<<"fin"<<endl;
} 


