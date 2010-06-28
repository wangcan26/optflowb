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
  double ratio=0.80;
  int minWidth=18;
  double alpha = 2 ; // Global smoothness variable.
  double gamma = 0 ; // Global weight for derivatives.
	
	//IplImage* img= cvLoadImage(NULL); 
	// cvNamedWindow("TEST",  CV_WINDOW_AUTOSIZE); 
	// cvShowImage("TEST",img); 
	// cvWaitKey(0); 
	// cvReleaseImage(&img); 
	// cvDestroyWindow("TEST"); 
	
	 const IplImage* img1= cvLoadImage(argv[1],CV_LOAD_IMAGE_GRAYSCALE);//zero is for grayscale  CV_LOAD_IMAGE_GRAYSCALE
	 const IplImage* img2= cvLoadImage(argv[2],CV_LOAD_IMAGE_GRAYSCALE); //1 is for color CV_LOAD_IMAGE_COLOR
	
	 IplImage *img1_32 = cvCreateImage(cvSize(img1->width, img1->height), coarse2fComp._imageDepth, img1->nChannels);
	 IplImage *img2_32 = cvCreateImage(cvSize(img2->width, img2->height), coarse2fComp._imageDepth, img2->nChannels);

	cvConvertScale(img1, img1_32, 1/255.);
	cvConvertScale(img2, img2_32, 1/255.);
	
	cvNormalize(img1_32,img1_32,0,255,CV_MINMAX); 
	cvNormalize(img2_32,img2_32,0,255,CV_MINMAX); 

	
	// toolsKit::cvShowManyImages("Image",2, img1_32,img1_32);
	 cvWaitKey(0);

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
	 
	 coarse2fComp.Coarse2FineFlow(  vx, 
									vy, 
									*img1_32, 
									*img2_32, 
									alpha, 
									gamma,
									ratio, 
									minWidth, 
									3, 
									1, 
									0);




/*
	double a[] = {1,0,1,
		1,1,1,
		1,0,1};
	CvMat A = cvMat(3,3,CV_64FC1,a);
	printf("A:\n");
	PrintMat(&A);
	int aDim[] = {3,3};
	CvSparseMat * SA = cvCreateSparseMat(2,aDim,CV_32FC1);
	
	for (int i=0; i<A.rows; i++) for (int j=0; j<A.cols;j++){
		CvScalar val = cvGet2D(&A,i,j);
		if (val.val[0]!=0)
			cvSet2D(SA,i,j,val);
		val = cvGet2D(SA,i,j);
		printf("(%d, %d):%f %c",i,j,val.val[0],j==A.cols-1?'\n':','); 
		}
	SparseToolKit::printSparseMat(SA);
	double b[] = {4,6,4};
	CvMat B = cvMat(3,1,CV_64FC1,b);
	printf("B:\n");
	PrintMat(&B);
	int bDim[] = {3,1};
	CvSparseMat * SB = cvCreateSparseMat(2,bDim,CV_32FC1);
	
	for (int i=0; i<B.rows; i++) for (int j=0; j<B.cols;j++){
		CvScalar val = cvGet2D(&B,i,j);
		if (val.val[0] !=0 ) cvSet2D(SB,i,j,val);
		printf("(%d, %d):%f %c",i,j,val.val[0],j==B.cols-1?'\n':','); 
		}
	SparseToolKit::printSparseMat(SB);
	double x[] = {0,0,0};
	int xDim[]={3,1};
	CvSparseMat *SX =NULL;//cvCreateSparseMat(2,xDim,CV_32FC1);
	//for (int i=0; i<B.rows; i++) for (int j=0; j<B.cols;j++){
	//	cvSet2D(SX,i,j,cvScalar(0));
	//}
	//	cvSet2D(SX,0,0,cvScalar(1));
	//cvSet2D(SX,0,1,cvScalar(2));
	//cvSet2D(SX,0,2,cvScalar(3));
	//CvMat X = cvMat(3,1,CV_32FC1,x);
	//	cvSolve(SA, SB, SX);    // solve (Ax=b) for x
	//toolsKit::SparseSor(SA,SB,SX);
	toolsKit::split(SA,SB,cvScalar(0.5));
	//SparseToolKit::printSparseMat(SB);
	cvReleaseSparseMat(&SA);
	cvReleaseSparseMat(&SB);
	cvReleaseSparseMat(&SX);*/
	} 


