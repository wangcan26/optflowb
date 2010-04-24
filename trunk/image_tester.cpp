#include <iostream>
#include <stdio.h>
#include "highgui.h"
#include "GaussPyramid.h"
#include "coarse2FineCompute.h"


using namespace std;
int main (int argc,char** argv) 
{ 

  GaussPyramid GPyramid1;
  coarse2FineCompute coarse2fComp;
  double ratio=0.75;
  int minWidth=30;
	
	
	//IplImage* img= cvLoadImage(NULL); 
	// cvNamedWindow("TEST",  CV_WINDOW_AUTOSIZE); 
	// cvShowImage("TEST",img); 
	// cvWaitKey(0); 
	// cvReleaseImage(&img); 
	// cvDestroyWindow("TEST"); 
	
	 const IplImage* img1= cvLoadImage(argv[1],CV_LOAD_IMAGE_GRAYSCALE);//zero is for grayscale  CV_LOAD_IMAGE_GRAYSCALE
	 const IplImage* img2= cvLoadImage(argv[2],CV_LOAD_IMAGE_GRAYSCALE); //1 is for color CV_LOAD_IMAGE_COLOR

	/* cvNamedWindow("TEST",CV_WINDOW_AUTOSIZE); 
	 cvShowImage("TEST",img1); 
	 cvNamedWindow("TEST2",CV_WINDOW_AUTOSIZE); 
	 cvShowImage("TEST2",img2); */
	// toolsKit::cvShowManyImages("Image",2, img1,img2);
	 cvWaitKey(0);

	 IplImage* vx= NULL;
	 IplImage* vy= NULL;
	 IplImage* warpI2= NULL;

	 //GPyramid1.ConstructPyramid(*img,ratio,minWidth);

	 coarse2fComp.Coarse2FineFlow(  vx, 
									vy, 
									*warpI2,
									*img1, 
									*img2, 
									0, 
									ratio, 
									minWidth, 
									3, 
									500, 
									0);



} 