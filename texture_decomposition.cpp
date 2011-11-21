#include <vector>
#include <algorithm>
#include "highgui.h" 
#include "cv.h"
#include "IplImageIterator.h"
#include "toolsKit.h"
#include <cmath>


void findMinMax(IplImage *img,float &min,float &max)
{
	IplImageIterator<float> it(img);

	min=max=*it;
	it++;
	while (!it)
	{
		if (*it<min)
			min=*it;
		if (*it>max)
			max=*it;
		it++;
	}
}

// Reproject so each vector given by [p0 p1] is of length <=1
void Reproject(IplImage *p0,IplImage *p1)
{
	IplImageIterator<float> it0(p0);	
	IplImageIterator<float> it1(p1);		
	float max=1.0f;
	while (!it0)
	{
		max=std::max(max,sqrt((*it0)*(*it0)+(*it1)*(*it1)));
		it0++;
		it1++;
	}

	IplImageIterator<float> i0(p0);	
	IplImageIterator<float> i1(p1);			
	while (!i0)
	{		
		*i0=*i0/max;
		*i1=*i1/max;
		i0++;
		i1++;
	}

}
// Rescale an image to be in the range [minv maxv]
void rescale(IplImage *img1,IplImage *img2,float minv,float maxv)
{
	float min1,max1,min2,max2,min,max;
	findMinMax(img1,min1,max1);
	findMinMax(img2,min2,max2);
	min=std::min(min1,min2);
	max=std::max(max1,max2);
	IplImageIterator<float> it1(img1);		
	IplImageIterator<float> it2(img2);		
	while (!it1)
	{
		*it1=((*it1-min)/(max-min))*(maxv-minv)+minv;
		*it2=((*it2-min)/(max-min))*(maxv-minv)+minv;
		it1++;
		it2++;
	}
}
// Decompose the input IMAGE into structure and texture parts using the
// Rudin-Osher-Fatemi method. The final output is a linear combination 
// of the decomposed texture and the structure parts.

// in1 - first image (input)
// in2 - second image (input)
// texture1 - texture part for first image (ouput)
// texture2 - texture part for second image (ouput)
// structure1 - structure part for first image or NULL if not needed(ouput)
// structure2 - structure part for second image or NULL if not needed (ouput)
// Example usage:
// Assuming img1_32g and img2_32g grayscale images
//
// IplImage *texture1=cvCreateImage(cvSize(img1_32g->width,img1_32g->height),IPL_DEPTH_32F,1);
// IplImage *texture2=cvCreateImage(cvSize(img2_32g->width,img2_32g->height),IPL_DEPTH_32F,1);
// IplImage *structure1=cvCreateImage(cvSize(img1_32g->width,img1_32g->height),IPL_DEPTH_32F,1);
// IplImage *structure2=cvCreateImage(cvSize(img2_32g->width,img2_32g->height),IPL_DEPTH_32F,1);
// structure_texture_decomposition_rof(img1_32g,img1_32g,texture1,texture2,structure1,structure2,1.0f/8.0f,100,textureAlpha);
void structure_texture_decomposition_rof(IplImage *in1,IplImage *in2,IplImage *texture1,IplImage *texture2,IplImage *structure1,IplImage *structure2,float theta,int nIters,float alp)
{	
	in1=cvCloneImage(in1);
	in2=cvCloneImage(in2);
	IplImage *im[]={in1,in2};
	// Rescale the input image to [-1 1]
	rescale(in1,in2,-1,1);
	// Backup orginal images
	IplImage *o1=cvCloneImage(in1);
	IplImage *o2=cvCloneImage(in2);

	IplImage *p1=cvCreateImage(cvSize(in1->width,in1->height),in1->depth,1);		
	IplImage *p2=cvCreateImage(cvSize(in1->width,in1->height),in1->depth,1);	
	IplImage *tmp=cvCreateImage(cvSize(in1->width,in1->height),in1->depth,1);		
	IplImage *tmp1=cvCreateImage(cvSize(in1->width,in1->height),in1->depth,1);		
	IplImage *tmp2=cvCreateImage(cvSize(in1->width,in1->height),in1->depth,1);	
	IplImage *I_x=cvCreateImage(cvSize(in1->width,in1->height),in1->depth,1);		
	IplImage *I_y=cvCreateImage(cvSize(in1->width,in1->height),in1->depth,1);	
	IplImage *div_p=cvCreateImage(cvSize(in1->width,in1->height),in1->depth,1);	
	IplImage *p[]={p1,p2};
	// stepsize
	float delta = 1.0f/(4.0f*theta);
	double x1[]={0,1,-1};
	double x2[]={1,-1};
	CvMat* ker1 = cvCreateMat(1,3,CV_64FC1);
	CvMat* ker4 = cvCreateMat(3,1,CV_64FC1);
	CvMat* ker2 = cvCreateMat(1,2,CV_64FC1);
	CvMat* ker3 = cvCreateMat(2,1,CV_64FC1);
	CvMat* test = cvCreateMat(1,3,CV_64FC1);

	cvmSet(ker1,0,0,0);
	cvmSet(ker1,0,1,1);
	cvmSet(ker1,0,2,-1);

	cvmSet(test,0,0,1);
	cvmSet(test,0,1,1);
	cvmSet(test,0,2,1);

	cvmSet(ker4,0,0,0);
	cvmSet(ker4,1,0,1);
	cvmSet(ker4,2,0,-1);

	cvmSet(ker2,0,0,1);
	cvmSet(ker2,0,1,-1);

	cvmSet(ker3,0,0,1);
	cvmSet(ker3,1,0,-1);

	for (int i=0;i<2;i++)
	{
		// Initialize dual variable p to be 0
		cvZero(p[0]);
		cvZero(p[1]);
		for (int iter=0;iter<nIters;iter++)
		{
			// Compute divergence
			cvFilter2D(p[0],tmp1,ker1);
			cvFilter2D(p[1],tmp2,ker4);
			cvAdd(tmp1,tmp2,div_p);			

			toolsKit::cvMulScalar(div_p,theta);
			cvCopyImage(im[i],tmp);			
			cvAdd(tmp,div_p,tmp);						
			cvFilter2D(tmp,I_x,ker2);			
			cvFilter2D(tmp,I_y,ker3);


			toolsKit::cvMulScalar(I_x,delta);
			toolsKit::cvMulScalar(I_y,delta);
			// Update dual variable
			cvAdd(p[0],I_x,p[0]);
			cvAdd(p[1],I_y,p[1]);
			// Reproject to |p| <= 1
			Reproject(p[0],p[1]);			
		}
		// compute divergence    
		cvFilter2D(p[0],tmp1,ker1);
		cvFilter2D(p[1],tmp2,ker4);
		cvAdd(tmp1,tmp2,div_p);
		toolsKit::cvMulScalar(div_p,theta);		
		cvAdd(im[i],div_p,im[i]);		
	}
	if (structure1 && structure2)
	{
		cvCopyImage(im[0],structure1);
		cvCopyImage(im[1],structure2);
		rescale(structure1,structure2,0,1);
	}
	toolsKit::cvMulScalar(im[0],alp);		
	toolsKit::cvMulScalar(im[1],alp);	
	cvSub(o1,im[0],texture1);	
	cvSub(o2,im[1],texture2);
	rescale(texture1,texture2,0,1);

}