#include <iostream>
#include <stdio.h>
#include "highgui.h"
#include "GaussPyramid.h"
#include "coarse2FineCompute.h"
#include "SparseToolKit.h"

using namespace std;


void PrintMat(CvMat *A)
{
	int i, j;
	for (i = 0; i < A->rows; i++)
	{
		printf("\n");
		switch (CV_MAT_DEPTH(A->type))
		{
		case CV_32F:
		case CV_64F:
			for (j = 0; j < A->cols; j++)
				printf ("%8.3f ", (float)cvGetReal2D(A, i, j));
			break;
		case CV_8U:
		case CV_16U:
			for(j = 0; j < A->cols; j++)
				printf ("%6d",(int)cvGetReal2D(A, i, j));
			break;
		default:
			break;
		}
	}
	printf("\n");
}



int main (int argc,char** argv) 
{ 

	GaussPyramid GPyramid1;
	double error_const=0.001;
	//IPL_DEPTH_32F IPL_DEPTH_8U
	coarse2FineCompute coarse2fComp(IPL_DEPTH_32F,error_const);
	double ratio=0.75;
	int minWidth=30;
	int outerIter=3;
	int innerIter=100;
	double alpha = 2 ; // Global smoothness variable.
	double gamma = 0 ; // Global weight for derivatives.

	//IplImage* img= cvLoadImage(NULL); 
	// cvNamedWindow("TEST",  CV_WINDOW_AUTOSIZE); 
	// cvShowImage("TEST",img); 
	// cvWaitKey(0); 
	// cvReleaseImage(&img); 
	// cvDestroyWindow("TEST"); 

	 IplImage* img1= cvLoadImage(argv[1],CV_LOAD_IMAGE_COLOR);//zero is for grayscale  CV_LOAD_IMAGE_GRAYSCALE
	 IplImage* img2= cvLoadImage(argv[2],CV_LOAD_IMAGE_COLOR); //1 is for color CV_LOAD_IMAGE_COLOR


	//IplImage* img1g= cvCreateImage(cvSize(img1->width, img1->height),img1->depth, 1);
	//IplImage* img2g= cvCreateImage(cvSize(img1->width, img1->height), img1->depth, 1);
	//IplImage* img1RGB= cvCreateImage(cvSize(img1->width, img1->height),img1->depth, img1->nChannels);
	// IplImage* img2RGB= cvCreateImage(cvSize(img1->width, img1->height), img1->depth, img1->nChannels);

	//cvZero(img1g);
	//cvZero(img2g);

	//cvConvertImage(img1,img1RGB, CV_CVTIMG_SWAP_RB);
	//cvConvertImage(img2,img2RGB, CV_CVTIMG_SWAP_RB);


	//cvCvtColor( img1RGB, img1g, CV_RGB2GRAY );
	//cvCvtColor( img2RGB, img2g, CV_RGB2GRAY );

	

	IplImage *img1_32 = cvCreateImage(cvSize(img1->width, img1->height), coarse2fComp._imageDepth, img1->nChannels);
	IplImage *img2_32 = cvCreateImage(cvSize(img1->width, img1->height), coarse2fComp._imageDepth, img1->nChannels);

	cvZero(img1_32);
	cvZero(img2_32);
	//upscale from char to float
	cvConvertScale(img1, img1_32, 1.0/255);
	cvConvertScale(img2, img2_32, 1.0/255);


	cvSmooth(img1_32,img1_32,CV_GAUSSIAN);
	cvSmooth(img2_32,img2_32,CV_GAUSSIAN);

	IplImage *img1_32g = cvCreateImage(cvSize(img1->width, img1->height), coarse2fComp._imageDepth, 1);
	IplImage *img2_32g = cvCreateImage(cvSize(img1->width, img1->height), coarse2fComp._imageDepth, 1);

	cvCvtColor( img1_32, img1_32g, CV_BGR2GRAY );
	cvCvtColor( img2_32, img2_32g, CV_BGR2GRAY );

	/*cout<<"Im1:before smooth"<<endl;
	toolsKit::IPL_print(img1);
	
	   cvSmooth(img1,img1,CV_GAUSSIAN);
	   cvSmooth(img1,img2,CV_GAUSSIAN);*/


//	cvNormalize(img1_32g,img1_32g,127,0,CV_MINMAX); //CV_MINMAX
//	cvNormalize(img2_32g,img2_32g,127,0,CV_MINMAX); 


	/*cout<<"Im1:after up scale"<<endl;
	toolsKit::IPL_print(img1_32g);
	cout<<"Im2:after up scale"<<endl;
	toolsKit::IPL_print(img2_32g);*/
	

	IplImage* vx= NULL;
	IplImage* vy= NULL;
	//IplImage* warpI2= NULL;

	/*double an[] ={0,0,0,0,0,0,0,0,0};
	double a[] = {1,2,3,4,5,6,7,8,9};
	double b[]={2,2,2,3,3,3,4,5,6};
	CvMat* matans = &cvMat( 3, 3, CV_64FC1, an ); // 64FC1 for double
	CvMat* matans2 = &cvMat( 3, 3, CV_64FC1, an ); // 64FC1 for double
	CvMat* mata = &cvMat( 3, 3, CV_64FC1, a );
	CvMat* matb = &cvMat( 3, 3, CV_64FC1, b );
	cvMul(mata,matb,matans);
	cvMul(matb,mata,matans2);
	PrintMat(mata);
	PrintMat(matb);
	PrintMat(matans);
	PrintMat(matans2);
	cout<<"a"<<endl;*/
	//GPyramid1.ConstructPyramid(*img,ratio,minWidth);

	/*
	IplImage *imga1=cvCreateImage(cvSize(3,3), IPL_DEPTH_32F,1);
	IplImage *imga2=cvCreateImage(cvSize(3,3), IPL_DEPTH_32F,1);
	IplImage *imgans=cvCreateImage(cvSize(3,3), IPL_DEPTH_32F,1);
	((float*)imga1->imageData)[0]=1;
	((float*)imga1->imageData)[1]=2;
	((float*)imga1->imageData)[2]=3;
	((float*)imga1->imageData)[3]=4;
	((float*)imga1->imageData)[4]=5;
	((float*)imga1->imageData)[5]=6;
	((float*)imga1->imageData)[6]=7;
	((float*)imga1->imageData)[7]=8;
	((float*)imga1->imageData)[8]=9;

	//for top bottom
	//((float*)imga2->imageData)[0]=1;
	//((float*)imga2->imageData)[1]=1;
	//((float*)imga2->imageData)[2]=1;
	//((float*)imga2->imageData)[3]=2;
	//((float*)imga2->imageData)[4]=2;
	//((float*)imga2->imageData)[5]=2;
	//((float*)imga2->imageData)[6]=3;
	//((float*)imga2->imageData)[7]=3;
	//((float*)imga2->imageData)[8]=3;
	//for left right
	((float*)imga2->imageData)[0]=1;
	((float*)imga2->imageData)[1]=2;
	((float*)imga2->imageData)[2]=3;
	((float*)imga2->imageData)[3]=1;
	((float*)imga2->imageData)[4]=2;
	((float*)imga2->imageData)[5]=3;
	((float*)imga2->imageData)[6]=1;
	((float*)imga2->imageData)[7]=2;
	((float*)imga2->imageData)[8]=3;
	toolsKit::IPL_print(imga1);
	toolsKit::IPL_print(imga2);
	// cvZero(imgans);
	// toolsKit::IPL_add_left(imga1,imga2,imga1);
	// toolsKit::IPL_print(imga1);
	toolsKit::IPL_add_right(imga1,imga2,imga1);

	toolsKit::IPL_print(imga1);
	toolsKit::IPL_print(imga1);
	*/

	toolsKit::cvShowManyImages("img1,img2 color",2,img1_32,img2_32);
	cvWaitKey(0);
	toolsKit::cvShowManyImages("img1,img2",2,img1_32g,img2_32g);
	cvWaitKey(0);
	 coarse2fComp.Coarse2FineFlow(vx,vy, 
								  *img1_32g, *img2_32g, 
								  alpha,gamma,
								  ratio,minWidth, 
								  outerIter,innerIter);

	cout<<"fin"<<endl;
} 


