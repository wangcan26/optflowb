#include <iostream>
#include <stdio.h>
#include "highgui.h"
#include "GaussPyramid.h"
#include "coarse2FineCompute.h"


using namespace std;
int main (int argc,char** argv) 
{ 
  char c;
  char *mystring;
  GaussPyramid GPyramid1;
  coarse2FineCompute coarse2fComp;
  double ratio=0.75;
	int minWidth=30;
	
	//std::cout<<"aaaa"<<std::cout;
//		puts ("Enter text. Include a dot ('.') in a sentence to exit:");

	if(argv[1]==NULL){		
	//IplImage* img= cvLoadImage(NULL); 
	// cvNamedWindow("TEST",  CV_WINDOW_AUTOSIZE); 
	// cvShowImage("TEST",img); 
	// cvWaitKey(0); 
	// cvReleaseImage(&img); 
	// cvDestroyWindow("TEST"); 
	}
	else{	
	 const IplImage* img1= cvLoadImage(argv[1]); 
	 const IplImage* img2= cvLoadImage(argv[2]); 

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
									0, 
									0, 
									0);

	 //cvNamedWindow("TEST",  CV_WINDOW_AUTOSIZE); 
	 //cvShowImage("TEST",img); 
	 cvWaitKey(0); 
	// cvReleaseImage(&img); 
	 //cvDestroyWindow("TEST"); 
	}
//cin >> mystring;

	/* do {
    c=getchar();
    putchar (c);
  } while (c != '.');
  return 0;*/

} 