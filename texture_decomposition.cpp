#include <vector>
#include <algorithm>
#include "highgui.h" 
#include "cv.h"
#include <cmath>
#include "UtilsMat.h"

void findMinMax(cv::Mat img,float &min,float &max)
{
	cv::Mat_<float>::iterator it = img.begin<float>();
	cv::Mat_<float>::iterator itEnd = img.end<float>();

	min = max = *it;
	it++;
	while (it != itEnd)
	{
		if (*it<min)
			min=*it;
		if (*it>max)
			max=*it;
		it++;
	}
}

// Reproject so each vector given by [p0 p1] is of length <=1
void Reproject(cv::Mat p0,cv::Mat p1)
{
	cv::Mat_<float>::iterator it0 = p0.begin<float>();
	cv::Mat_<float>::iterator it1 = p1.begin<float>();
	cv::Mat_<float>::iterator it0End = p0.end<float>();
	float max=1.0f;
	while (it0 != it0End)
	{
		max=std::max(max,sqrt((*it0)*(*it0)+(*it1)*(*it1)));
		it0++;
		it1++;
	}

	cv::Mat_<float>::iterator i0 = p0.begin<float>();
	cv::Mat_<float>::iterator i1 = p1.begin<float>();
	cv::Mat_<float>::iterator i0End = p0.end<float>();
	while (i0 != i0End)
	{		
		*i0=*i0/max;
		*i1=*i1/max;
		i0++;
		i1++;
	}

}
// Rescale an image to be in the range [minv maxv]
void rescale(cv::Mat img1,cv::Mat img2,float minv,float maxv)
{
	float min1,max1,min2,max2,min,max;
	findMinMax(img1,min1,max1);
	findMinMax(img2,min2,max2);
	min=std::min(min1,min2);
	max=std::max(max1,max2);
	cv::Mat_<float>::iterator it1 = img1.begin<float>();
	cv::Mat_<float>::iterator it2 = img2.begin<float>();
	cv::Mat_<float>::iterator it1End = img1.end<float>();
	cv::Mat_<float>::iterator it2End = img2.end<float>();	
	while (it1!=it1End)
	{
		*it1=((*it1-min1)/(max1-min1))*(maxv-minv)+minv;
		*it2=((*it2-min2)/(max2-min2))*(maxv-minv)+minv;
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
void structureTextureDecompositionRof(const cv::Mat& in1,const cv::Mat& in2,cv::Mat& texture1,cv::Mat& texture2, cv::Mat* structure1, cv::Mat* structure2, const float theta, const int nIters, const float alpha)
{	
	cv::Mat new_in1 = in1.clone();
	cv::Mat new_in2 = in2.clone();
	cv::Mat im[] = {new_in1,new_in2};
	 //Rescale the input image to [-1 1]
	rescale(new_in1,new_in2,-1,1);
	// Backup orginal images
	cv::Mat o1 = new_in1.clone();
	cv::Mat o2 = new_in2.clone();

	cv::Mat p1(new_in1.rows, new_in1.cols, new_in1.type());
	cv::Mat p2(new_in1.rows, new_in1.cols, new_in1.type());
	cv::Mat tmp(new_in1.rows, new_in1.cols, new_in1.type());
	cv::Mat tmp1(new_in1.rows, new_in1.cols, new_in1.type());
	cv::Mat tmp2(new_in1.rows, new_in1.cols, new_in1.type());
	cv::Mat I_x(new_in1.rows, new_in1.cols, new_in1.type());
	cv::Mat I_y(new_in1.rows, new_in1.cols, new_in1.type());
	cv::Mat div_p(new_in1.rows, new_in1.cols, new_in1.type());
	cv::Mat p[] = {p1, p2};
	
	// stepsize
	float delta = 1.0f/(4.0f*theta);

	cv::Mat ker1 = (cv::Mat_<double>(1,3) << -1, 1, 0);
	cv::Mat ker2 = (cv::Mat_<double>(1,2) << -1, 1);
	cv::Mat ker3 = (cv::Mat_<double>(2,1) << -1, 1);
	cv::Mat ker4 = (cv::Mat_<double>(3,1) << -1, 1, 0);

	for (int i=0;i<2;i++)
	{
		// Initialize dual variable p to be 0
		p[0].setTo(cv::Scalar(0));
		p[1].setTo(cv::Scalar(0));
		for (int iter=0;iter<nIters;iter++)
		{
			// Compute divergence TODO::check the filter borders
			//cv::Mat bla = cv::copyMakeBorder(p[0],bb,1,1,0,0,

			cv::filter2D(p[0],tmp1,tmp1.depth(),ker1,cv::Point(0,0),cv::BORDER_CONSTANT);
			cv::filter2D(p[1],tmp2,tmp2.depth(),ker4,cv::Point(0,0),cv::BORDER_CONSTANT);
			div_p = tmp1 + tmp2;

			div_p *= theta;
			im[i].copyTo(tmp);
			tmp += div_p;

			cv::filter2D(tmp,I_x,I_x.depth(),ker2,cv::Point(0,0),cv::BORDER_REPLICATE);
			cv::filter2D(tmp,I_y,I_y.depth(),ker3,cv::Point(0,0),cv::BORDER_REPLICATE);

			I_x *= delta;
			I_y *= delta;
			// Update dual variable
			p[0] += I_x;
			p[1] += I_y;
			// Reproject to |p| <= 1
			Reproject(p[0],p[1]);			
		}
		// compute divergence    
		cv::filter2D(p[0],tmp1,tmp1.depth(),ker1,cv::Point(0,0),cv::BORDER_CONSTANT);
		cv::filter2D(p[1],tmp2,tmp2.depth(),ker4,cv::Point(0,0),cv::BORDER_CONSTANT);
		div_p = tmp1 + tmp2;
		div_p *= theta;
		im[i] += div_p;		
	}
	if (structure1 && structure2)
	{
		im[0].copyTo(*structure1);
		im[1].copyTo(*structure2);
		rescale(*structure1,*structure2,0,1);
	}
	im[0] *= alpha;		
	im[1] *= alpha;	
	texture1 = o1 - im[0];	
	texture2 = o2 - im[1];
	rescale(texture1,texture2,0,1);
}

