#include <vector>
#include <algorithm>
#include "highgui.h" 
#include "cv.h"
#include "IplImageIterator.h"
#include "toolsKit.h"
#include <cmath>

IplImage *medianFilter(IplImage *src,int n)
{
	CvScalar s;
	std::vector<double> neighbours(n*n);
	IplImage *ret=cvCloneImage(src);
	int edge=n/2;	
	for (int y=edge;y<ret->height-edge;y++)
	{
		for (int x=edge;x<ret->width-edge;x++)
		{
			neighbours.clear();
			for (int j=-edge;j<=edge;j++)
			{
				for (int i=-edge;i<=edge;i++)
				{					
					s=cvGet2D(src,y+j,x+i);
					neighbours.push_back(s.val[0]);
				}
			}
			std::sort(neighbours.begin(),neighbours.end());
			s.val[0]=neighbours[n*n/2+1];
			cvSet2D(ret,y,x,s);
		}
	}
	return ret;	
}
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
void rescale(IplImage *img,float minv,float maxv)
{
	float min,max;
	findMinMax(img,min,max);
	IplImageIterator<float> it(img);		
	while (!it)
	{
		*it=((*it-min)/(max-min))*(maxv-minv)+minv;
		it++;
	}
}
void rescale2(IplImage *img1,IplImage *img2,float minv,float maxv)
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

void Reproject(IplImage *p0,IplImage *p1)
{
	IplImageIterator<float> it0(p0);	
	IplImageIterator<float> it1(p1);		
	float max=1.0f;
	while (!it0)
	{
		max=std::max(1.0f,sqrt((*it0)*(*it0)+(*it1)*(*it1)));
		*it0=*it0/max;
		*it1=*it1/max;
		it0++;
		it1++;
	}
	
}

void ShowImage(char *title,IplImage *img)
{
	cvNamedWindow( title, 1 );
    cvShowImage( title, img);

    cvWaitKey();
    cvDestroyWindow(title);
}
void ShowImage2(char *title,IplImage *img)
{
	cvNamedWindow( title, 1 );
    cvShowImage( title, img);    
}

void structure_texture_decomposition_rof(IplImage *in1,IplImage *in2,IplImage *texture1,IplImage *texture2,IplImage *structure1,IplImage *structure2,float theta,int nIters,float alp)
{	
	in1=cvCloneImage(in1);
	in2=cvCloneImage(in2);
	IplImage *im[]={in1,in2};
	//rescale(in1,-1,1);
	//rescale(in2,-1,1);
	rescale2(in1,in2,-1,1);
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
		cvZero(p[0]);
		cvZero(p[1]);
		for (int iter=0;iter<nIters;iter++)
		{
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
			cvAdd(p[0],I_x,p[0]);
			cvAdd(p[1],I_y,p[1]);
			Reproject(p[0],p[1]);			
		}
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
		rescale2(structure1,structure2,0,1);
	}
	toolsKit::cvMulScalar(im[0],alp);		
	toolsKit::cvMulScalar(im[1],alp);	
	cvSub(o1,im[0],texture1);	
	cvSub(o2,im[1],texture2);
	rescale2(texture1,texture2,0,1);
	
}




inline double f(double x) 
{ 
        return (x-floor(x)); 


} 


void GetSurrounding(IplImage *img,int px,int py,double a[4][4])
{
	for (int y=-1;y<=2;y++)
	{
		for (int x=-1;x<=2;x++)
		{
			if ((px+x>=0) && (py+y>=0) && (px+x<img->width) && (py+y<img->height))
				a[y+1][x+1]=cvGet2D(img,py+y,px+x).val[0];
			else
				a[y+1][x+1]=0;
		}
	}
}

inline double P(double x) 
{ 
        return (x>0.0)?x:0.0; 

} 


inline double R(double x) 
{ 
        //      return (pow(P(x+2),3) - 4* pow(P(x+1),3) + 6*pow(P(x),3) - 4* pow(P(x-1),3) ) / 6; 
        double p0=P(x+2); 
        double p1=P(x+1); 
        double p2=P(x); 
        double p3=P(x-1); 
        return ( (p0*p0*p0) - 4*(p1*p1*p1) + 6*(p2*p2*p2) - 4*(p3*p3*p3) ) / 6.0; 
}

inline double bicubic( 
        double tx, double ty, 
        double P00, double P01, double P02, double P03, 
        double P10, double P11, double P12, double P13, 
        double P20, double P21, double P22, double P23, 
        double P30, double P31, double P32, double P33 ) 
{ 
        double dx=f(tx); 
        double dy=f(ty); 


        double Rdx[4],Rdy[4]; 
        for(int n=0;n<=3;n++) 
        { 
                Rdx[n]=R(n-1-dx); 
                Rdy[n]=R(n-1-dy); 
        } 


        double s; 


        s =P00*Rdx[0]*Rdy[0]; 
        s+=P01*Rdx[1]*Rdy[0]; 
        s+=P02*Rdx[2]*Rdy[0]; 
        s+=P03*Rdx[3]*Rdy[0]; 


        s+=P10*Rdx[0]*Rdy[1]; 
        s+=P11*Rdx[1]*Rdy[1]; 
        s+=P12*Rdx[2]*Rdy[1]; 
        s+=P13*Rdx[3]*Rdy[1]; 


        s+=P20*Rdx[0]*Rdy[2]; 
        s+=P21*Rdx[1]*Rdy[2]; 
        s+=P22*Rdx[2]*Rdy[2]; 
        s+=P23*Rdx[3]*Rdy[2]; 


        s+=P30*Rdx[0]*Rdy[3]; 
        s+=P31*Rdx[1]*Rdy[3]; 
        s+=P32*Rdx[2]*Rdy[3]; 
        s+=P33*Rdx[3]*Rdy[3]; 


        return s; 


}


void ImageResize(IplImage *src,IplImage *dst,bool interp=true)
{
	double s[4][4];	
	double ratio_x=(double)src->width/(double)dst->width;
	double ratio_y=(double)src->height/(double)dst->height;
	for (unsigned int y=0;y<dst->height;y++)
	{
		for (unsigned int x=0;x<dst->width;x++)
		{
			double src_x=(double)x*ratio_x;
			double src_y=(double)y*ratio_y;
			unsigned int src_ix=(int)floor(src_x);
			unsigned int src_iy=(int)floor(src_y);
			GetSurrounding(src,src_ix,src_iy,s);
			CvScalar sc;
			if (interp)
				sc.val[0]=bicubic(src_x,src_y,s[0][0],s[0][1],s[0][2],s[0][3],s[1][0],s[1][1],s[1][2],s[1][3],s[2][0],s[2][1],s[2][2],s[2][3],s[3][0],s[3][1],s[3][2],s[3][3]);
			else
				sc.val[0]=s[1][1];
			cvSet2D(dst,y,x,sc);
		}
	}
	//rescale(dst,0,1);
}
/*
void GetSurrounding2(IplImage *img,int px,int py,double a[4][4])
{
	for (int i=1;i<=4;i++)
	{
		if (px+i<img->width)
			a[0][i-1]=cvGet2D(img,py,px+i).val[0];
		else
			a[0][i-1]=0;		

		if (px-i>=0)
			a[1][i-1]=cvGet2D(img,py,px-i).val[0];
		else
			a[1][i-1]=0;

		if (py+i<img->height)
			a[2][i-1]=cvGet2D(img,py+i,px).val[0];
		else
			a[2][i-1]=0;

		if (py-i>=0)
			a[3][i-1]=cvGet2D(img,py-i,px).val[0];
		else
			a[3][i-1]=0;
	}
}

double bicubicInterpolate (double p[4][4], double x, double y) {
	double a00 = p[1][1];
	double a01 = -.5*p[1][0] + .5*p[1][2];
	double a02 = p[1][0] - 2.5*p[1][1] + 2*p[1][2] - .5*p[1][3];
	double a03 = -.5*p[1][0] + 1.5*p[1][1] - 1.5*p[1][2] + .5*p[1][3];
	double a10 = -.5*p[0][1] + .5*p[2][1];
	double a11 = .25*p[0][0] - .25*p[0][2] - .25*p[2][0] + .25*p[2][2];
	double a12 = -.5*p[0][0] + 1.25*p[0][1] - p[0][2] + .25*p[0][3] + .5*p[2][0] - 1.25*p[2][1] + p[2][2] - .25*p[2][3];
	double a13 = .25*p[0][0] - .75*p[0][1] + .75*p[0][2] - .25*p[0][3] - .25*p[2][0] + .75*p[2][1] - .75*p[2][2] + .25*p[2][3];
	double a20 = p[0][1] - 2.5*p[1][1] + 2*p[2][1] - .5*p[3][1];
	double a21 = -.5*p[0][0] + .5*p[0][2] + 1.25*p[1][0] - 1.25*p[1][2] - p[2][0] + p[2][2] + .25*p[3][0] - .25*p[3][2];
	double a22 = p[0][0] - 2.5*p[0][1] + 2*p[0][2] - .5*p[0][3] - 2.5*p[1][0] + 6.25*p[1][1] - 5*p[1][2] + 1.25*p[1][3] + 2*p[2][0] - 5*p[2][1] + 4*p[2][2] - p[2][3] - .5*p[3][0] + 1.25*p[3][1] - p[3][2] + .25*p[3][3];
	double a23 = -.5*p[0][0] + 1.5*p[0][1] - 1.5*p[0][2] + .5*p[0][3] + 1.25*p[1][0] - 3.75*p[1][1] + 3.75*p[1][2] - 1.25*p[1][3] - p[2][0] + 3*p[2][1] - 3*p[2][2] + p[2][3] + .25*p[3][0] - .75*p[3][1] + .75*p[3][2] - .25*p[3][3];
	double a30 = -.5*p[0][1] + 1.5*p[1][1] - 1.5*p[2][1] + .5*p[3][1];
	double a31 = .25*p[0][0] - .25*p[0][2] - .75*p[1][0] + .75*p[1][2] + .75*p[2][0] - .75*p[2][2] - .25*p[3][0] + .25*p[3][2];
	double a32 = -.5*p[0][0] + 1.25*p[0][1] - p[0][2] + .25*p[0][3] + 1.5*p[1][0] - 3.75*p[1][1] + 3*p[1][2] - .75*p[1][3] - 1.5*p[2][0] + 3.75*p[2][1] - 3*p[2][2] + .75*p[2][3] + .5*p[3][0] - 1.25*p[3][1] + p[3][2] - .25*p[3][3];
	double a33 = .25*p[0][0] - .75*p[0][1] + .75*p[0][2] - .25*p[0][3] - .75*p[1][0] + 2.25*p[1][1] - 2.25*p[1][2] + .75*p[1][3] + .75*p[2][0] - 2.25*p[2][1] + 2.25*p[2][2] - .75*p[2][3] - .25*p[3][0] + .75*p[3][1] - .75*p[3][2] + .25*p[3][3];

	double x2 = x * x;
	double x3 = x2 * x;
	double y2 = y * y;
	double y3 = y2 * y;

	return a00 + a01 * y + a02 * y2 + a03 * y3 +
	       a10 * x + a11 * x * y + a12 * x * y2 + a13 * x * y3 +
	       a20 * x2 + a21 * x2 * y + a22 * x2 * y2 + a23 * x2 * y3 +
	       a30 * x3 + a31 * x3 * y + a32 * x3 * y2 + a33 * x3 * y3;
}


*/