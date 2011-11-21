#include <vector>
#include <algorithm>
#include <string>
#include <cmath>
#include <limits>
using namespace std;
#include "highgui.h" 
#include "cv.h"
#include "IplImageIterator.h"
#include "toolsKit.h"
#include "FlowUtils.h"

int ncols = 0;
#define MAXCOLS 60
int colorwheel[MAXCOLS][3];
void makecolorwheel();
// Default values for error range
const float cFlowUtils::minValueAngleError=0.0f;
const float cFlowUtils::maxValueAngleError=0.436f;
const float cFlowUtils::minValueEndpointError=0.0f;
const float cFlowUtils::maxValueEndpointError=20.0f;
const float cFlowUtils::minValueOriError=0.0f;
const float cFlowUtils::maxValueOriError=0.436f;
const float cFlowUtils::minValueMagError=0.0f;
const float cFlowUtils::maxValueMagError=20.0f;
cFlowUtils::cFlowUtils() : U(NULL),V(NULL),gtU(NULL),gtV(NULL),eeImage(NULL),aeImage(NULL),orientationImage(NULL),magnitudeImage(NULL),mask(NULL)
{	
}
cFlowUtils::cFlowUtils(IplImage *U,IplImage *V,IplImage *gtU,IplImage *gtV)
{
	SetFlow(U,V);
	SetGroundTruthFlow(gtU,gtV);
	CalculateError(errorAE,0.0f);
	CalculateError(errorEP,0.0f);
	CalculateError(errorMagAndOri,0.0f);
}
cFlowUtils::~cFlowUtils()
{
	if (mask)
	{
		delete [] mask;
		mask=NULL;
	}
	if (U)
	{
		cvReleaseImage(&U);
		U=NULL;
	}
	if (V)
	{
		cvReleaseImage(&V);
		V=NULL;
	}
	if (gtU)
	{
		cvReleaseImage(&gtU);
		gtU=NULL;
	}
	if (gtV)
	{
		cvReleaseImage(&gtV);
		gtV=NULL;
	}
	if (aeImage)
	{
		cvReleaseImage(&aeImage);
		aeImage=NULL;
	}
	if (eeImage)
	{
		cvReleaseImage(&eeImage);
		eeImage=NULL;
	}
	if (orientationImage)
	{
		cvReleaseImage(&orientationImage);
		orientationImage=NULL;
	}
	if (magnitudeImage)
	{
		cvReleaseImage(&magnitudeImage);
		magnitudeImage=NULL;
	}
}
void cFlowUtils::SetFlow(IplImage *U,IplImage *V)
{
	this->U=cvCloneImage(U);
	this->V=cvCloneImage(V);
	CalculateMask();
}
void cFlowUtils::SetFlow(CvMat *U,CvMat *V)
{
	this->U=cvCreateImage(cvSize(U->rows,U->cols),IPL_DEPTH_32F,1);
	this->V=cvCreateImage(cvSize(U->rows,U->cols),IPL_DEPTH_32F,1);
	for (int y=0;y<U->cols;y++)
	{
		for (int x=0;x<U->rows;x++)
		{
			float u=(float)cvmGet(U,x,y);
			float v=(float)cvmGet(V,x,y);
			cvSet2D(this->U,y,x,cvScalar(u));
			cvSet2D(this->V,y,x,cvScalar(v));
		}
	}
	CalculateMask();
}
void cFlowUtils::SetGroundTruthFlow(IplImage *U,IplImage *V)
{
	this->gtU=cvCloneImage(U);
	this->gtV=cvCloneImage(V);
	CalculateMask();
}
void cFlowUtils::SetGroundTruthFlow(CvMat *U,CvMat *V)
{
	this->gtU=cvCreateImage(cvSize(U->rows,U->cols),IPL_DEPTH_32F,1);
	this->gtV=cvCreateImage(cvSize(U->rows,U->cols),IPL_DEPTH_32F,1);
	for (int y=0;y<U->cols;y++)
	{
		for (int x=0;x<U->rows;x++)
		{
			float u=(float)cvmGet(U,x,y);
			float v=(float)cvmGet(V,x,y);
			cvSet2D(this->gtU,y,x,cvScalar(u));
			cvSet2D(this->gtV,y,x,cvScalar(v));
		}
	}
	CalculateMask();
}
bool cFlowUtils::LoadFlow(const string &filename)
{
	bool ret = cFlowUtils::ReadFlowFile(filename,&U,&V,IPL_DEPTH_32F);
	CalculateMask();
	return ret;
}
bool cFlowUtils::LoadGroundTruthFlow(const string &filename)
{
	bool ret = cFlowUtils::ReadFlowFile(filename,&gtU,&gtV,IPL_DEPTH_32F);
	CalculateMask();
	return ret;
}
bool cFlowUtils::ReadFlowFile(const string &filename,IplImage **U,IplImage **V,int depth)
{	
	const float tag_float = 202021.25f;
	FILE *file=fopen(filename.c_str(),"rb");
	if (!file)
	{
		throw CError("ReadFlowFile: could not open %s", filename.c_str());
		return false;
	}
	int width, height;
	float tag;

	fread(&tag,sizeof(float),1,file);
	if (tag != tag_float)
	{
		fclose(file);
		return false;
	}
	fread(&width,sizeof(int),1,file);
	fread(&height,sizeof(int),1,file);
	
	*U=cvCreateImage(cvSize(width,height),depth,1);
	*V=cvCreateImage(cvSize(width,height),depth,1);
	
	for (int y=0;y<height;y++)
	{
		for (int x=0;x<width;x++)
		{
			float u,v;
			fread(&u,sizeof(float),1,file);
			fread(&v,sizeof(float),1,file);
			
			CvScalar su=cvScalar(u);
			CvScalar sv=cvScalar(v);
			cvSet2D(*U,y,x,su);				
			cvSet2D(*V,y,x,sv);			
		}
	}
	fclose(file);
	return true;
}
void cFlowUtils::CalculateMask()
{
	if (U==NULL || V==NULL || gtU==NULL || gtV==NULL)
		return;

	validPixels=0;	
	if (!mask)
		mask=new char[U->width*U->height];
	for (int y=0;y<U->height;y++)
	{
		for (int x=0;x<U->width;x++)
		{
			float u=(float)cvGet2D(U,y,x).val[0];
			float v=(float)cvGet2D(V,y,x).val[0];
			float gtu=(float)cvGet2D(gtU,y,x).val[0];
			float gtv=(float)cvGet2D(gtV,y,x).val[0];			
			if (u>1000 || v>1000 || gtu>1000 || gtv>1000)
				mask[y*U->width+x]=0;
			else
			{
				mask[y*U->width+x]=1;
				validPixels++;
			}
		}
	}
}
bool cFlowUtils::WriteFlowFile(const string &filename,IplImage *U,IplImage *V)
{
	if (!U || !V)
		return false;
	if ((U->width!=V->width) || (U->height!=V->height))
		return false;

	FILE *file=fopen(filename.c_str(),"wb");
	if (!file)
	{
		throw CError("ReadFlowFile: could not open %s", filename.c_str());
		return false;
	}

	fwrite("PIEH",1,4,file);
	int width=U->width;
	int height=U->height;
	fwrite(&width,sizeof(int),1,file);
	fwrite(&height,sizeof(int),1,file);

	for (int y=0;y<U->height;y++)
	{
		for (int x=0;x<U->width;x++)
		{
			float u=(float)cvGet2D(U,y,x).val[0];
			float v=(float)cvGet2D(V,y,x).val[0];
			fwrite(&u,sizeof(float),1,file);
			fwrite(&v,sizeof(float),1,file);
		}
	}
	fclose(file);
	return true;
}
bool cFlowUtils::ReadFlowFile(const string &filename,CvMat **U,CvMat **V)
{	
	const float tag_float = 202021.25f;
	FILE *file=fopen(filename.c_str(),"rb");
	if (!file)
	{
		throw CError("ReadFlowFile: could not open %s", filename.c_str());
		return false;
	}
	int width, height;
	float tag;

	fread(&tag,sizeof(float),1,file);
	if (tag != tag_float)
	{
		fclose(file);
		return false;
	}
	fread(&width,sizeof(int),1,file);
	fread(&height,sizeof(int),1,file);
	
	*U=cvCreateMat(width,height,CV_32FC1);
	*V=cvCreateMat(width,height,CV_32FC1);

	for (int y=0;y<height;y++)
	{
		for (int x=0;x<width;x++)
		{
			float u,v;
			fread(&u,sizeof(float),1,file);
			fread(&v,sizeof(float),1,file);
			cvmSet(*U,x,y,u);
			cvmSet(*V,x,y,v);			
		}
	}
	fclose(file);
	return true;
}
bool cFlowUtils::WriteFlowFile(const string &filename,CvMat *U,CvMat *V)
{
	if (!U || !V)
		return false;
	if ((U->width!=V->width) || (U->height!=V->height))
		return false;

	FILE *file=fopen(filename.c_str(),"wb");
	if (!file)
	{
		throw CError("ReadFlowFile: could not open %s", filename.c_str());
		return false;
	}

	fwrite("PIEH",1,4,file);
	int width=U->width;
	int height=U->height;
	fwrite(&width,sizeof(int),1,file);
	fwrite(&height,sizeof(int),1,file);

	for (int y=0;y<U->height;y++)
	{
		for (int x=0;x<U->width;x++)
		{
			float u=(float)cvmGet(U,x,y);
			float v=(float)cvmGet(V,x,y);			
			fwrite(&u,sizeof(float),1,file);
			fwrite(&v,sizeof(float),1,file);
		}
	}
	fclose(file);
	return true;
}
// Calculates the various errors.
// [minValue maxValue] is the error range for the angle error endpoint error and magnitude error.
// [minValue2 maxValue2] is the error range for the magnitude error.
void cFlowUtils::CalculateError(eErrorType errorType,float minValue,float maxValue,float minValue2,float maxValue2)
{
	const float pi=3.1415926535897932384626433832795f;
	double totalError=0.0;
	double totalError2=0.0;
	if (errorType==errorAE)
	{	
		if (minValue==0.0f && maxValue==0.0f)
		{
			minValue=cFlowUtils::minValueAngleError;
			maxValue=cFlowUtils::maxValueAngleError;
		}
		aeImage=cvCreateImage(cvSize(U->width,U->height),IPL_DEPTH_32F,1);		
		for (int y=0;y<U->height;y++)
		{
			for (int x=0;x<U->width;x++)
			{
				if (mask[y*U->width+x])
				{				
					float u=(float)cvGet2D(U,y,x).val[0];
					float v=(float)cvGet2D(V,y,x).val[0];
					float gtu=(float)cvGet2D(gtU,y,x).val[0];
					float gtv=(float)cvGet2D(gtV,y,x).val[0];
					// get the angle between (u,v,1), (gtu,gtv,1).
					float ae=(float)acos( (1.0f+u*gtu+v*gtv) / (sqrt(1.0f+u*u+v*v)*sqrt(1.0f+gtu*gtu+gtv*gtv)));
					if (ae<minValue)
						ae=minValue;
					else if (ae>maxValue)
						ae=maxValue;
					cvSet2D(aeImage,y,x,cvScalar(ae/pi));			
					totalError+=ae;
				}
			}
		}
		averageAE = totalError/((double)validPixels);
	}
	else if (errorType==errorEP)
	{
		if (minValue==0.0f && maxValue==0.0f)
		{
			minValue=cFlowUtils::minValueEndpointError;
			maxValue=cFlowUtils::maxValueEndpointError;
		}
		eeImage=cvCreateImage(cvSize(U->width,U->height),IPL_DEPTH_32F,1);
		for (int y=0;y<U->height;y++)
		{
			for (int x=0;x<U->width;x++)
			{
				if (mask[y*U->width+x])
				{
					float u=(float)cvGet2D(U,y,x).val[0];
					float v=(float)cvGet2D(V,y,x).val[0];
					float gtu=(float)cvGet2D(gtU,y,x).val[0];
					float gtv=(float)cvGet2D(gtV,y,x).val[0];
					// end-point error in pixels.
					float ee=(float)sqrt( (u-gtu) * (u-gtu) + (v-gtv) * (v-gtv) );
					if (ee<minValue)
						ee=minValue;
					else if (ee>maxValue)
						ee=maxValue;

					cvSet2D(eeImage,y,x,cvScalar(ee));			
					totalError+=ee;
				}
			}
		}
		averageEE = totalError/((double)validPixels);
	}
	else if (errorType==errorMagAndOri)
	{
		if (minValue==0.0f && maxValue==0.0f)
		{
			minValue=cFlowUtils::minValueMagError;
			maxValue=cFlowUtils::maxValueMagError;
		}
		if (minValue2==0.0f && maxValue2==0.0f)
		{
			minValue2=cFlowUtils::minValueOriError;
			maxValue2=cFlowUtils::maxValueOriError;
		}
		orientationImage=cvCreateImage(cvSize(U->width,U->height),IPL_DEPTH_32F,1);
		magnitudeImage=cvCreateImage(cvSize(U->width,U->height),IPL_DEPTH_32F,1);
		for (int y=0;y<U->height;y++)
		{
			for (int x=0;x<U->width;x++)
			{
				if (mask[y*U->width+x])
				{
					float u=(float)cvGet2D(U,y,x).val[0];
					float v=(float)cvGet2D(V,y,x).val[0];
					float gtu=(float)cvGet2D(gtU,y,x).val[0];
					float gtv=(float)cvGet2D(gtV,y,x).val[0];
					float mag1=sqrt(u*u+v*v);
					float mag2=sqrt(gtu*gtu+gtv*gtv);
					float mage=fabs(mag1-mag2);
					float ori=(float)acos( (u*gtu+v*gtv) / (mag1*mag2));
					if (ori != ori)
						continue;
					if (ori<minValue2)
						ori=minValue2;
					else if (ori>maxValue2)
						ori=maxValue2;
					if (mage<minValue)
						mage=minValue;
					else if (mage>maxValue)
						mage=maxValue;
					cvSet2D(orientationImage,y,x,cvScalar(ori));			
					totalError+=ori;					
					cvSet2D(magnitudeImage,y,x,cvScalar(mage));
					totalError2+=mage;
				}
			}
		}
		averageOri = totalError/((double)validPixels);
		averageMag = totalError2/((double)validPixels);

	}
}
double cFlowUtils::GetAverageAE()
{
	return averageAE;
}
double cFlowUtils::GetAverageEP()
{
	return averageEE;
}
double cFlowUtils::GetAverageMag()
{
	return averageMag;
}
double cFlowUtils::GetAverageOri()
{
	return averageOri;
}
void cFlowUtils::DisplayFlowError(eErrorType errorType)
{
	if (errorType==errorAE)
	{
		cvShowManyImages("Angle error",1,aeImage);
	}
	else if (errorType==errorEP)
	{
		cvShowManyImages("End-point error",1,eeImage);
	}
	else if (errorType==errorMagAndOri)
	{
		cvShowManyImages("Magnitude and orientation error",2,magnitudeImage,orientationImage);
	}
}
void cFlowUtils::DisplayFlowAndError(eErrorType errorType,float maxMotion)
{
	IplImage *flow=GetFlowImage(U,V,maxMotion);
	if (errorType==errorAE)
	{
		cvShowManyImages("Angle error",2,flow,aeImage);
	}
	else if (errorType==errorEP)
	{
		cvShowManyImages("End-point error",2,flow,eeImage);
	}
	else if (errorType==errorMagAndOri)
	{
		cvShowManyImages("Magnitude and orientation error",3,flow,magnitudeImage,orientationImage);
	}
}

float FindMaxLength(IplImage *U,IplImage *V)
{
	float maxLength=numeric_limits<float>::lowest();
	for (int y=0;y<U->height;y++)
	{
		for (int x=0;x<U->width;x++)
		{
			float u,v;
			u=(float)cvGet2D(U,y,x).val[0];
			v=(float)cvGet2D(V,y,x).val[0];
			if (fabs(u)>1e9 || fabs(v)>1e9)
				continue;
			float l=sqrt(u*u+v*v);
			if (l>maxLength)
				maxLength=l;
		}
	}
	return maxLength;
}

void cFlowUtils::DisplayGroundTruthFlow(float maxMotion)
{
	IplImage *flow=GetFlowImage(gtU,gtV,maxMotion);
	cvShowManyImages("Ground truth flow",1,flow);	
	cvReleaseImage(&flow);
}
void cFlowUtils::DisplayFlow(float maxMotion)
{
	IplImage *flow=GetFlowImage(U,V,maxMotion);
	cvShowManyImages("Flow",1,flow);		
	cvReleaseImage(&flow);
}

IplImage *cFlowUtils::GetFlowImage(IplImage *U,IplImage *V,float maxMotion)
{
	const float pi=3.1415926535897932384626433832795f;
	IplImage* color_img = cvCreateImage( cvSize(U->width,U->height), IPL_DEPTH_8U, 3 );
	cvZero(color_img);
	float maxLength=FindMaxLength(U,V);

	if (maxMotion>0.0f)
		maxLength=maxMotion;

	if (maxLength==0.0f)
		maxLength=1.0f;

	if (ncols==0)
		makecolorwheel();

	for (int y=0;y<U->height;y++)
	{
		for (int x=0;x<U->width;x++)
		{
			CvScalar out=cvGet2D(color_img,y,x);
			float u,v;
			u=(float)cvGet2D(U,y,x).val[0];
			v=(float)cvGet2D(V,y,x).val[0];			
			if (fabs(u)<1e9 && fabs(v)<1e9)
			{
				u=u/maxLength;
				v=v/maxLength;
				float rad = sqrt(u * u + v * v);
				float a = atan2(-v, -u) / pi;
				float fk = (a + 1.0f) / 2.0f * (ncols-1);
				int k0 = (int)fk;
				int k1 = (k0 + 1) % ncols;
				float f = fk - k0;
				//f = 0; // uncomment to see original color wheel
				for (int b = 0; b < 3; b++) {
					float col0 = colorwheel[k0][b] / 255.0f;
					float col1 = colorwheel[k1][b] / 255.0f;
					float col = (1 - f) * col0 + f * col1;
					if (rad <= 1)
						col = 1 - rad * (1 - col); // increase saturation with radius
					else
					{
						col *= .75; // out of range
						out.val[b]=(int)(255.0 * col);
					}
				}				
			}
			cvSet2D(color_img,y,x,out);
		}
	}
	
		
	return color_img;
}
// Tests error calculation and example usage
void cFlowUtils::Test()
{
	LoadFlow("Media\\testCflowUtil.flo");
	LoadGroundTruthFlow("Media\\testCflowUtilGT.flo");
	CalculateError(cFlowUtils::errorAE);
	CalculateError(cFlowUtils::errorEP);
	CalculateError(cFlowUtils::errorMagAndOri);
	cout << "Average Angle error: " << GetAverageAE() << endl;
	cout << "Average endpoint error: " << GetAverageEP() << endl;
	cout << "Average magnitue error: " << GetAverageMag() << endl;
	cout << "Average orientation error: " << GetAverageOri() << endl;
}
void cFlowUtils::DrawFlow(IplImage* U,IplImage* V)
{
	IplImage *flow=GetFlowImage(U,V);
	cvShowManyImages("flow",1,flow);
	cvWaitKey(0);			
	cvReleaseImage(&flow);
}
void cFlowUtils::DrawFlow(CvMat* U,CvMat* V)
{
	IplImage stub, *u,*v;
	u = cvGetImage(U, &stub);
	v = cvGetImage(V, &stub);
	DrawFlow(u,v);
}
void cFlowUtils::DrawFlow2(IplImage* du,IplImage* u,IplImage* dv,IplImage* v)
{
	IplImage* tempsumU =cvCreateImage( cvSize(u->width,u->height), u->depth, u->nChannels );
	IplImage* tempsumV =cvCreateImage( cvSize(u->width,u->height), u->depth, u->nChannels );

	cvAdd(u,du,tempsumU);
	cvAdd(v,dv,tempsumV);

	DrawFlow(tempsumU,tempsumV);

	cvReleaseImage(&tempsumU);
	cvReleaseImage(&tempsumV);

}
void cFlowUtils::cvShowManyImages(char* title, int nArgs, ...) {

	// img - Used for getting the arguments 
	IplImage *img;
	// DispImage - the image in which input images are to be copied
	IplImage *DispImage	;
	int size;
	int i;
	int m, n;
	int x, y;
	// w - Maximum number of images in a row 
	// h - Maximum number of images in a column 
	int w, h;
	// scale - How much we have to resize the image
	float scale;
	int max;
	// If the number of arguments is lesser than 0 or greater than 12
	// return without displaying 
	if(nArgs <= 0) {
		printf("Number of arguments too small....\n");
		return;
	}
	else if(nArgs > 12) {
		printf("Number of arguments too large....\n");
		return;
	}
	// Determine the size of the image, 
	// and the number of rows/cols 
	// from number of arguments 
	else if (nArgs == 1) {
		w = h = 1;
		size = 300;
	}
	else if (nArgs == 2) {
		w = 2; h = 1;
		size = 300;
	}
	else if (nArgs == 3 || nArgs == 4) {
		w = 2; h = 2;
		size = 300;
	}
	else if (nArgs == 5 || nArgs == 6) {
		w = 3; h = 2;
		size = 200;
	}
	else if (nArgs == 7 || nArgs == 8) {
		w = 4; h = 2;
		size = 200;
	}
	else {
		w = 4; h = 3;
		size = 150;
	}
	va_list args1;
	va_start(args1, nArgs);
	// Create a new 3 channel image
	img = va_arg(args1, IplImage*);
	if(img->nChannels==1)
		DispImage = cvCreateImage( cvSize(100 + size*w, 60 + size*h), img->depth, 1 );
	else
		DispImage = cvCreateImage( cvSize(100 + size*w, 60 + size*h), img->depth, 3 );
	va_list args;
	va_start(args, nArgs);
	// Used to get the arguments passed
	// Loop for nArgs number of arguments
	for (i = 0, m = 20, n = 20; i < nArgs; i++, m += (20 + size)) {
		// Get the Pointer to the IplImage
		img = va_arg(args, IplImage*);
		// Check whether it is NULL or not
		// If it is NULL, release the image, and return
		if(img == 0) {
			printf("Invalid arguments");
			cvReleaseImage(&DispImage);
			return;
		}
		// Find the width and height of the image
		x = img->width;
		y = img->height;
		// Find whether height or width is greater in order to resize the image
		max = (x > y)? x: y;
		// Find the scaling factor to resize the image
		scale = (float) ( (float) max / size );
		// Used to Align the images
		if( i % w == 0 && m!= 20) {
			m = 20;
			n+= 20 + size;
		}
		// Set the image ROI to display the current image
		cvSetImageROI(DispImage, cvRect(m, n, (int)( x/scale ), (int)( y/scale )));
		// Resize the input image and copy the it to the Single Big Image
		cvResize(img, DispImage);
		// Reset the ROI in order to display the next image
		cvResetImageROI(DispImage);
	}
	// Create a new window, and show the Single Big Image
	cvNamedWindow( title, 1 );
	cvShowImage( title, DispImage);
	cvWaitKey(1);
	//cvDestroyWindow(title);
	// End the number of arguments
	va_end(args);
	// Release the Image Memory
	cvReleaseImage(&DispImage);
}

void setcols(int r, int g, int b, int k)
{
	colorwheel[k][0] = r;
	colorwheel[k][1] = g;
	colorwheel[k][2] = b;
}

void makecolorwheel()
{
	// relative lengths of color transitions:
	// these are chosen based on perceptual similarity
	// (e.g. one can distinguish more shades between red and yellow 
	//  than between yellow and green)
	int RY = 15;
	int YG = 6;
	int GC = 4;
	int CB = 11;
	int BM = 13;
	int MR = 6;
	ncols = RY + YG + GC + CB + BM + MR;
	//printf("ncols = %d\n", ncols);
	if (ncols > MAXCOLS)
		exit(1);
	int i;
	int k = 0;
	for (i = 0; i < RY; i++) setcols(255,	   255*i/RY,	 0,	       k++);
	for (i = 0; i < YG; i++) setcols(255-255*i/YG, 255,		 0,	       k++);
	for (i = 0; i < GC; i++) setcols(0,		   255,		 255*i/GC,     k++);
	for (i = 0; i < CB; i++) setcols(0,		   255-255*i/CB, 255,	       k++);
	for (i = 0; i < BM; i++) setcols(255*i/BM,	   0,		 255,	       k++);
	for (i = 0; i < MR; i++) setcols(255,	   0,		 255-255*i/MR, k++);
}