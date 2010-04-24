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
	cvReleaseSparseMat(&SX);
	} 
