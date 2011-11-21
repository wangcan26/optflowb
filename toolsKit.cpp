#include "toolsKit.h"
#include "optical_flow_demo.h"
toolsKit::toolsKit()
	{
	}

toolsKit::~toolsKit(void)
	{
	}
void toolsKit::IPL_mul_inverse(IplImage* img,int opType){
	if(img->depth==IPL_DEPTH_8U){
		IplImageIterator<unsigned char> it(img);		
		(opType==1?IPL_mul_inverse_loop(it):IPLsqrt_mul2(it));
		}
	else if (img->depth==IPL_DEPTH_32F){
		IplImageIterator<float> it(img);
		(opType==1?IPL_mul_inverse_loop(it):IPLsqrt_mul2(it));
		}
	else if (img->depth==IPL_DEPTH_64F){
		IplImageIterator<double> it(img);
		(opType==1?IPL_mul_inverse_loop(it):IPLsqrt_mul2(it));
		}
	else{
		cout<<"IPL_mul_inverse got unsupported depth"<<endl;
		return;
		}

	}

void toolsKit::IPL_add(IplImage* img,IplImage* img2,IplImage* dest){
	IplImageIterator<float> it(img);
	IplImageIterator<float> it2(img2);
	IplImageIterator<float> it3(dest);
	while (!it) {      
		//	cout<<"bef:"<<it.data;
		*it3= ((float)*it)+((float)*it2); 
		// cout<<"=>"<<it.data<<endl;
		++it;
		++it2;
		++it3;
		}
	}
//srcHorizontal is wt+1
void toolsKit::IPL_mul_different_sizes(IplImage* src,IplImage* srcHorizontal,IplImage* dest){
	int j=0;
	int k=src->width-1;
	int m=0;

	for (j=0,m=0; m < dest->width*dest->height;j++,m++)					
		{	

		float temp1=((float*)src->imageData)[m];
		float temp2=((float*)srcHorizontal->imageData)[j];
		((float*)dest->imageData)[m]=((float*)src->imageData)[m]+((float*)srcHorizontal->imageData)[j];			

		if(k==0){//smaller pic is at row end		
			k=srcHorizontal->width-1;
			j++;			
			}
		else
			k--;			
		}		
	}
//srcVertical is ht+1
void toolsKit::IPL_mul_different_sizes2(IplImage* src,IplImage* srcVertical,IplImage* dest){	
	for (int m=0; m < dest->width*dest->height;m++)					
		{	
		((float*)dest->imageData)[m]=((float*)src->imageData)[m]+((float*)srcVertical->imageData)[m];					
		}		
	}
//wt*(ht+1) +(ht+1)*wt matrixs only
//dest size is x*y
void toolsKit::IPL_add_different_sizes(IplImage* imgHorizonal,IplImage* imgVertical,IplImage* dest){
	int i=0;
	int j=0;
	int m=0;
	int k=imgHorizonal->width-1;
	
	for (i=imgHorizonal->width,j=0,m=0; m < dest->width*dest->height; i++,j++,m++)					
		{	

		//float temp1=((float*)imgHorizonal->imageData)[i];
		//float temp2=((float*)imgVertical->imageData)[j];
		((float*)dest->imageData)[m]=((float*)imgHorizonal->imageData)[i]+((float*)imgVertical->imageData)[j];			

		if(k==0){//smaller pic is at row end		
			k=imgHorizonal->width-1;
			j++;			
			}
		else
			k--;	

		}		
	}
void toolsKit::IPL_add_different_sizes2(IplImage* imgVertical,IplImage* imgVertical2,IplImage* dest){
	int i=1;
	int j=1;
	int k=dest->width-1;
	cvZero(dest);
	for (i = 1,j=1; i < dest->width*dest->height-1; i++,j++)					
		{				
		if(k==0){//smaller pic is at row end		
			k=dest->width-1;
			i++;
			}
		else
			k--;
		if (k!=dest->width)
			((float*)dest->imageData)[j]=((float*)imgVertical->imageData)[i]+((float*)imgVertical2->imageData)[i];			
		}		
	}

void toolsKit::IPL_add_different_sizes3(IplImage* imgVertical,IplImage* imgVertical2,IplImage* dest){
	int i=1;
	int j=1;
	int k=dest->width-1;
	cvZero(dest);	
	for (i = 1,j=1; i < imgVertical->width*dest->height; i++,j++)					
		{				
		if(k==0){//smaller pic is at row end		
			k=imgVertical->width-1;
			j++;
			}
		else
			k--;
		if (k!=dest->width)
			((float*)dest->imageData)[j]=((float*)imgVertical->imageData)[i]+((float*)imgVertical2->imageData)[i];			
		}		
	}

void toolsKit::IPL_sub(IplImage* img,IplImage* img2,IplImage* dest){
	IplImageIterator<float> it(img);
	IplImageIterator<float> it2(img2);
	IplImageIterator<float> it3(dest);
	while (!it) {      
		//		cout<<"bef:"<<*it<<" "<<*it2;
		*it3=(float)( ((float)*it)-((float)*it2)); 
		//	 cout<<"=>"<<*it<<" "<<*it2<<endl;
		++it;
		++it2;
		++it3;
		}
	}
//img2 is the shifted image
void toolsKit::IPL_operate_top(IplImage* img,IplImage* shiftImg2,IplImage* dest,toolsKit::operations operation){
	int i=0;
	int j=0;
	int width=img->width;	
	for (i = width,j=0; i < width*img->height; i++,j++)					
		{	
		switch (operation){					  
					  case ADD : ((float*)dest->imageData)[i]=((float*)img->imageData)[i]+((float*)shiftImg2->imageData)[j];				
						  break;
					  case SUB : ((float*)dest->imageData)[i]=((float*)img->imageData)[i]-((float*)shiftImg2->imageData)[j];				
						  break;
					  case MUL : ((float*)dest->imageData)[i]=((float*)img->imageData)[i]*((float*)shiftImg2->imageData)[j];				
						  break;
			}					
		}
	}
//img2 is the shifted image
void toolsKit::IPL_operate_bottom(IplImage* img,IplImage* shiftImg2,IplImage* dest,toolsKit::operations operation){
	int i=0;
	int j=0;
	int width=img->width;				
	for (i = 0,j=width; i < width*img->height-width; i++,j++)					
		{	
		switch (operation){					  
					  case ADD : ((float*)dest->imageData)[i]=((float*)img->imageData)[i]+((float*)shiftImg2->imageData)[j];
						  break;
					  case SUB : ((float*)dest->imageData)[i]=((float*)img->imageData)[i]-((float*)shiftImg2->imageData)[j];
						  break;
					  case MUL : ((float*)dest->imageData)[i]=((float*)img->imageData)[i]*((float*)shiftImg2->imageData)[j];
						  break;
			}

		}
	}
//img2 is the shifted image
void toolsKit::IPL_operate_left(IplImage* img,IplImage* shiftImg2,IplImage* dest,toolsKit::operations operation){
	int width=img->width;
	if (img->nChannels==3)
		for (int i = 1; i < img->width*img->height*3; i+=3)
			{
			// 
			}
	else				
		for (int i = 1; i < width*img->height; i++)					
			{
			if(i %width  !=0)//do not compute first column
				switch (operation){					  
						  case ADD : ((float*)dest->imageData)[i]=((float*)img->imageData)[i]+((float*)shiftImg2->imageData)[i-1];
							  break;
						  case SUB : ((float*)dest->imageData)[i]=((float*)img->imageData)[i]-((float*)shiftImg2->imageData)[i-1];
							  break;
						  case MUL : ((float*)dest->imageData)[i]=((float*)img->imageData)[i]*((float*)shiftImg2->imageData)[i-1];
							  break;
				}

			}
	}
//img2 is the shifted image
void toolsKit::IPL_operate_right(IplImage* img,IplImage* shiftImg2,IplImage* dest,toolsKit::operations operation){
	int width=img->width;					
	int k=width-1;	
	//for first cell	
	for (int i = 0; i < width*img->height; i++)
		{
		if(k)//do not compute last column
			switch (operation){					  
					  case ADD : ((float*)dest->imageData)[i]=((float*)img->imageData)[i]+((float*)shiftImg2->imageData)[i+1];
						  break;
					  case SUB : ((float*)dest->imageData)[i]=((float*)img->imageData)[i]-((float*)shiftImg2->imageData)[i+1];
						  break;
					  case MUL : ((float*)dest->imageData)[i]=((float*)img->imageData)[i]*((float*)shiftImg2->imageData)[i+1];
						  break;
			}									
		k==0?k=width-1:k--;
		}
	}


//adders
void toolsKit::IPL_add_top(IplImage* img,IplImage* shiftImg2,IplImage* dest){
	IPL_operate_top(img,shiftImg2,dest,ADD);
	}
void toolsKit::IPL_add_bottom(IplImage* img,IplImage* shiftImg2,IplImage* dest){
	IPL_operate_bottom(img,shiftImg2,dest,ADD);
	}
void toolsKit::IPL_add_left(IplImage* img,IplImage* shiftImg2,IplImage* dest){
	IPL_operate_left(img,shiftImg2,dest,ADD);
	}
void toolsKit::IPL_add_right(IplImage* img,IplImage* shiftImg2,IplImage* dest){
	IPL_operate_right(img,shiftImg2,dest,ADD);
	}
//subs
void toolsKit::IPL_sub_top(IplImage* img,IplImage* shiftImg2,IplImage* dest){
	IPL_operate_top(img,shiftImg2,dest,SUB);
	}
void toolsKit::IPL_sub_bottom(IplImage* img,IplImage* shiftImg2,IplImage* dest){
	IPL_operate_bottom(img,shiftImg2,dest,SUB);
	}
void toolsKit::IPL_sub_left(IplImage* img,IplImage* shiftImg2,IplImage* dest){
	IPL_operate_left(img,shiftImg2,dest,SUB);
	}
void toolsKit::IPL_sub_right(IplImage* img,IplImage* shiftImg2,IplImage* dest){
	IPL_operate_right(img,shiftImg2,dest,SUB);
	}
//multipliers
void toolsKit::IPL_mul_top(IplImage* img,IplImage* shiftImg2,IplImage* dest){
	IPL_operate_top(img,shiftImg2,dest,MUL);
	}
void toolsKit::IPL_mul_bottom(IplImage* img,IplImage* shiftImg2,IplImage* dest){
	IPL_operate_bottom(img,shiftImg2,dest,MUL);
	}
void toolsKit::IPL_mul_left(IplImage* img,IplImage* shiftImg2,IplImage* dest){
	IPL_operate_left(img,shiftImg2,dest,MUL);
	}
void toolsKit::IPL_mul_right(IplImage* img,IplImage* shiftImg2,IplImage* dest){
	IPL_operate_right(img,shiftImg2,dest,MUL);
	}

void toolsKit::PrintMat(CvMat *A)
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
void toolsKit::IPL_print(const IplImage *image) {
	int ht= image->height; // number of lines
	int wt= image->width ; // total number of element per line
	//int step= image->widthStep; // effective width
	// get the pointer to the image buffer
	//unsigned char *data= reinterpret_cast<unsigned char *>(image->imageData);
	cout<<"=============================width:"<<image->width <<" height:"<< image->height<<" channels:"<<image->nChannels<<"[B,G,R]============================="<<endl;
	for (int i=0; i<ht; i++) {

		for (int j=0; j<wt; j+= 1) {
			CvScalar s= cvGet2D(image,i,j);
			// process each pixel ---------------------
			if(image->nChannels==3){
				cout<<"["<<s.val[0]<<","<<s.val[1]<<","<<s.val[2]<<"] ";
				}
			if(image->nChannels==2){
				cout<<"["<<s.val[0]<<","<<s.val[1]<<"] ";
				}
			else
				cout<<s.val[0]<<" ";
			// end of pixel processing ----------------
			} // end of line          
		cout<<endl;
		//data+= step;  // next line
		//break;
		}
	cout<<"======================================================================================================================================="<<endl;
	}

IplImage*  toolsKit::transposeImage(IplImage* image) {

	IplImage *rotated = cvCreateImage(cvSize(image->height,image->width), image->depth,image->nChannels);

  CvPoint2D32f center;

  float center_val = (float)((image->width)-1) / 2;
  center.x = center_val;
  center.y = center_val;
  CvMat *mapMatrix = cvCreateMat( 2, 3, CV_32FC1 );

  cv2DRotationMatrix(center, 90, 1.0, mapMatrix);
  cvWarpAffine(image, rotated, mapMatrix, CV_INTER_LINEAR + CV_WARP_FILL_OUTLIERS, cvScalarAll(0));

  cvReleaseMat(&mapMatrix);

  return rotated;
}

IplImage*  toolsKit::transposeImage2(IplImage* image) {

	IplImage *rotated = cvCreateImage(cvSize(image->height,image->width), image->depth,image->nChannels);

  CvPoint2D32f center;

  float center_val = (float)((image->width)-1) / 2;
  center.x = center_val;
  center.y = center_val;
  CvMat *mapMatrix = cvCreateMat( 2, 3, CV_32FC1 );
  
  cv2DRotationMatrix(center, 90, 1.0, mapMatrix);
  cvWarpAffine(image, rotated, mapMatrix, CV_INTER_LINEAR + CV_WARP_FILL_OUTLIERS, cvScalarAll(0));

  cvReleaseMat(&mapMatrix);
  cvFlip(rotated, NULL, 0);
  return rotated;
}
void toolsKit::cvMulScalar(IplImage* img,float scalar){
	IplImageIterator<float> it(img);
	while (!it) {   
		if(*it!=0)
			*it= (float)scalar*(*it); 
		++it;
		}
	}

void toolsKit::increaseImageSize(IplImage* src,IplImage* dst,int select){
	int i=0;
	int j=0;
	int width=src->width;					
	int k=width;	
	//for first cell	
	cvZero(dst);

	if (!select)//increase by one col
		for (i = 0,j=0; i < width*src->height;j++, i++)
			{
			if(k){//do not compute last column			  
				((float*)dst->imageData)[j]=((float*)src->imageData)[i];							
				//cout<<((float*)dst->imageData)[j]<<" ";
				}
			if(k==0){//smaller pic is at row end
				//cout<<endl;
				k=width-1;
				j++;
				}
			else
				k--;
			}
	if (select==2)//pad with zeros
		for (i = 0,j=dst->width+1; i < width*src->height;j++, i++)
			{
			
			if(k){//do not compute last column			  
				((float*)dst->imageData)[j]=((float*)src->imageData)[i];							
				//cout<<((float*)dst->imageData)[j]<<" ";			
				}
			if(k==0){//smaller pic is at row end
				//cout<<endl;
				k=width-1;
				j=j+2;		
				((float*)dst->imageData)[j]=((float*)src->imageData)[i];	
				}
			else
				k--;


			}
			
	else if(select==1)//increase by one row
		for (i = 0,j=0; i < width*src->height; i++,j++)					
			{	
			((float*)dst->imageData)[j]=((float*)src->imageData)[i];		
			}

	}
void toolsKit::cvZeroTop(IplImage* img){
	int i=0;
	int width=img->width;					
	//for first cell
	for (i = 0; i < width; i++)
		{
		if(i<width)//first row
			((float*)img->imageData)[i]=0;
		}
	}

void toolsKit::cvZeroBottom(IplImage* img){
	int i=0;
	int width=img->width;					
	//for first cell
	for (i = width*img->height-width; i < width*img->height; i++)
		{
		//if(i> width*img->height-width-1)
		((float*)img->imageData)[i]=0;
		}
	}
void toolsKit::cvZeroRight(IplImage* img){
	int i=0;
	int width=img->width;					
	int k=width-1;
		
	for (i = 0; i < width*img->height; i++)
		{		
		if(!k)//apply to first column only										  
			((float*)img->imageData)[i]=0;
		k==0?k=width-1:k--;
		}
	}

void toolsKit::cvZeroLeftRight(IplImage* img){
	int i=0;
	int width=img->width;					
	int k=width-1;		
	for (i = 0; i < width*img->height-width+1; )
		{	
		((float*)img->imageData)[i]=0;
		((float*)img->imageData)[i+width-1]=0;
		float a2=((float*)img->imageData)[i+width-1];
		i=i+width;
		}
	}

void toolsKit::cvZeroBottomLeft(IplImage* img){
	int i=0;
	int width=img->width;					
	int k=width-1;
		
	for (i = 0; i < width*img->height; i++)
		{
		if(i> width*img->height-width-1)
			((float*)img->imageData)[i]=0;
		if(!k==width-1)//apply to first column only										  
			((float*)img->imageData)[i]=0;
		k==0?k=width-1:k--;
		}

	}



void toolsKit::cvNormalizeEdges(IplImage* img){
	int i=0;
	int k=0;
	int width=img->width;					
	//for first cell
	k=width-1;
	for (i = 0; i < width*img->height; i++)
		{
		if(i<width)//first row
			((float*)img->imageData)[i]=0;
		if(i> width*img->height-width-1)//last row
			((float*)img->imageData)[i]=0;
		if(!k)//apply to last column only										  
			((float*)img->imageData)[i]=0;
		if(k==width-1)//apply to first column only										  
			((float*)img->imageData)[i]=0;

		k==0?k=width-1:k--;
		}
	}


void toolsKit::seperateDuDv(IplImage* du,IplImage* dv,vector<float> * dUdV){
	int i=0;
	int k=0;
	int mult=0;
	int element=0;
	int imageSize=du->width*du->height;

	for (vector<float>::iterator it = dUdV->begin(); it!= dUdV->end(); it++, i++){
				if (i<dUdV->size()/2)						
					((float*)du->imageData)[element+k] = *it;											
				else{
					if (imageSize==i){
						k=0;
						mult=0;
						element=0;										
						
					}
					((float*)dv->imageData)[element+k] = *it;
						
					}	
	

			//toolsKit::IplToFile(du,"c:\\a\\du_cpp.txt");
			//toolsKit::IplToFile(dv,"c:\\a\\dv_cpp.txt");
				if(mult==du->height-1){
					mult=0;
					element++;
					k=0;
				}
				else{
					mult++;
					k=k+du->width;
				}
	}	
}

//averaging as per the numerics section of the second chapter.for x and y
//{UP=1,DOWN=2,LEFT=3,RIGHT=4};
void toolsKit::shiftImage(IplImage* src,IplImage* emptyImage,int select){
	cvZero(emptyImage);
	if(select==RIGHT)//right 
		toolsKit::IPL_add_right(emptyImage,src,emptyImage);
	else if(select==LEFT)//down 
		toolsKit::IPL_add_left(emptyImage,src,emptyImage);
	else if(select==DOWN)//down 
		toolsKit::IPL_add_bottom(emptyImage,src,emptyImage);
	else if(select==UP)//down 
		toolsKit::IPL_add_top(emptyImage,src,emptyImage);

	}

bool toolsKit::AlmostEqualRelativeOrAbsolute(float A, float B,float maxRelativeError, float maxAbsoluteError)
	{
	if (fabs(A - B) < maxAbsoluteError)
		return true;
	float relativeError;
	if (fabs(B) > fabs(A))
		relativeError = fabs((A - B) / B);
	else
		relativeError = fabs((A - B) / A);
	if (relativeError <= maxRelativeError)
		return true;
	return false;
	}

bool toolsKit::IsNan(float A)
	{
	// A NAN has an exponent of 255 (shifted left 23 positions) and
	// a non-zero mantissa.
	int exp = *(int*)&A & 0x7F800000;
	int mantissa = *(int*)&A & 0x007FFFFF;
	if (exp == 0x7F800000 && mantissa != 0)
		return true;
	return false;
	}

void toolsKit::cvZeroNans(IplImage* img){				
	float nan1 = sqrt(-1.0f);
	for (int i = 0; i < img->width*img->height; i++)					
		{					
		//if(AlmostEqualRelativeOrAbsolute(FLT_MAX,((float*)img->imageData)[i],0.00001,0.00001))
		//	((float*)img->imageData)[i]=0;
		//if(AlmostEqualRelativeOrAbsolute(FLT_MIN,((float*)img->imageData)[i],0.0001,0.0001))
		//	((float*)img->imageData)[i]=0;
		if(IsNan(((float*)img->imageData)[i]))
			((float*)img->imageData)[i]=0;

		}
	}
//( var1 + var2*var3 + var4*var5 )^ 2==>ans
void toolsKit::costumeLineCompute(IplImage* ans,IplImage* var1,IplImage* var2,IplImage* var3,IplImage* var4,IplImage* var5){
	IplImage* temp=cvCreateImage(cvSize( var1->width, var1->height ),var1->depth,var1->nChannels);


	cvMul(var2,var3,ans);
	cvMul(var4,var5,temp);			
	cvAdd(var1,ans,ans);
	cvAdd(ans,temp,temp);
	cvPow(temp,ans,2);
	cvReleaseImage(&temp);
	}

IplImage*  toolsKit::psiDerivative(IplImage* x,double epsilon){	
	cvAddS(x,cvScalarAll(epsilon),x);
	toolsKit::IPL_mul_inverse(x,0);	
	return x;
	}

/*Function///////////////////////////////////////////////////////////////

Name:       cvShowManyImages

Purpose:    This is a function illustrating how to display more than one 
image in a single window using Intel OpenCV

Parameters: char *title: Title of the window to be displayed
int nArgs:   Number of images to be displayed
...:         IplImage*, which contains the images

Language:   C++

The method used is to set the ROIs of a Single Big image and then resizing 
and copying the input images on to the Single Big Image.

This function does not stretch the image... 
It resizes the image without modifying the width/height ratio..

This function can be called like this:

cvShowManyImages("Images", 2, img1, img2);
or
cvShowManyImages("Images", 5, img2, img2, img3, img4, img5);

This function can display upto 12 images in a single window.
It does not check whether the arguments are of type IplImage* or not.
The maximum window size is 700 by 660 pixels.
Does not display anything if the number of arguments is less than
one or greater than 12.

If you pass a pointer that is not IplImage*, Error will occur.
Take care of the number of arguments you pass, and the type of arguments, 
which should be of type IplImage* ONLY.

Idea was from BettySanchi of OpenCV Yahoo! Groups.

If you have trouble compiling and/or executing
this code, I would like to hear about it.

You could try posting on the OpenCV Yahoo! Groups
[url]http://groups.yahoo.com/group/OpenCV/messages/ [/url]


Parameswaran, 
Chennai, India.

cegparamesh[at]gmail[dot]com            

...
///////////////////////////////////////////////////////////////////////*/

/*void toolsKit::cvShowManyImages(char* title, int nArgs, ...) {

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

	*/
vector<float> * toolsKit::IplImageToCoulmnVector(IplImage* img){
	vector<float> * ans = new vector<float>(img->height*img->width);
	IplImageIterator<float> it(img);
	for (int i = 0; i<img->height; i++)	
		for (int j=0; j<img->width; j++){
			(*ans)[i+j*img->height] = *it;
			it++;
			}
		return ans;
	}

void toolsKit::ColumnVectorToIplImage(vector<float>* vCol, IplImage* image){
		int i=0, j=0;
		for (vector<float>::iterator it = vCol->begin(); it != vCol->end(); it++){
			cvSet2D(image, i,j,cvScalar(*it));
			i++;
			if (i==image->height){
				j = (j==image->width?0:j+1);
				i=0;	
				}
			}	
	}

IplImage * toolsKit::IplFromFile(string filename){
	ifstream thefile(filename.c_str(),ifstream::in);
	if (!thefile.good())
		return NULL;
	int height, width, depth;
	thefile>>height;
	thefile>>width;
	thefile>>depth;
	IplImage * ans = cvCreateImage(cvSize(width,height), depth,1);
	float val;
	for (int i = 0; i<height; i++)
		for (int j = 0; j<width; j++){
			thefile>>val;	
			cvSet2D(ans,i,j,cvScalar(val));
			}
		thefile.close();

		return ans;
	}

void toolsKit::IplToFile(IplImage* img, string  filename){
	ofstream thefile(filename.c_str());
	thefile<<img->height;
	thefile<<" ";
	thefile<<img->width;
	thefile<<" ";
	thefile<<img->depth;
	thefile<<"\n";
	for (int i =0; i<img->height; i++){
		for (int j=0; j<img->width; j++){
			thefile<<cvGet2D(img,i,j).val[0];
			if (j!=img->width-1) thefile<<" ";
			}
		thefile<<"\n";
		}
	thefile.close();
	}

void toolsKit::vectorTools::vectorMin(vector<float>* v, float val){
	for (vector<float>::iterator it = v->begin(); it != v->end(); it++)
		*it = (*it<val?*it:val);
	}

void toolsKit::vectorTools::vectorMax(vector<float>* v, float val){
	for (vector<float>::iterator it = v->begin(); it != v->end(); it++)
		*it = (*it>val?*it:val);
	}

vector<float> * toolsKit::vectorTools::vectorFloor(vector<float> * v){
	vector<float> * ans = new vector<float>();
	for (vector<float>::iterator it = v->begin(); it != v->end(); it++){
		ans->push_back((float)(((int)*it)/1));
		}
	return ans;
	}

vector<float> * toolsKit::vectorTools::vectorCeil(vector<float> * v){
	vector<float> * ans = new vector<float>();
	for (vector<float>::iterator it = v->begin(); it != v->end(); it++){
		ans->push_back(ceil(*it));
		}
	return ans;
	}

vector<float> * toolsKit::vectorTools::vectorSub(vector<float> * a, vector<float>* b){
	vector<float>* ans = new vector<float>();
	vector<float>::iterator ait = a->begin();
	vector<float>::iterator bit = b->begin();

	while(ait != a->end() && bit != b->end()){
		ans->push_back(*ait - *bit);
		ait++; bit++;
		}
	return ans;
	}

vector<float>* toolsKit::vectorTools::vectorSub(float val, vector<float>* b){
	vector<float>* ans = new vector<float>();
	for (vector<float>::iterator it = b->begin(); it != b->end(); it++)
		ans->push_back(val-*it);
	return ans;
	}

vector<float>* toolsKit::vectorTools::vectorSub(vector<float>* a, float val){
	vector<float>* ans = new vector<float>();
	for (vector<float>::iterator it = a->begin(); it != a->end(); it++)
		ans->push_back(*it-val);
	return ans;

	}

vector<float>* toolsKit::vectorTools::vectorMul(vector<float>* a, vector<float>* b){
	vector<float>* ans = new vector<float>();
	vector<float>::iterator ait = a->begin();
	vector<float>::iterator bit = b->begin();

	while(ait != a->end() && bit != b->end()){
		ans->push_back((*ait) * (*bit));
		ait++; bit++;
		}
	return ans;
	}

vector<float>* toolsKit::vectorTools::vectorMul(float val, vector<float>* b){
	vector<float>* ans = new vector<float>();
	for (vector<float>::iterator it = b->begin(); it != b->end(); it++)
		ans->push_back((*it) * val);
	return ans;
	}

vector<float>* toolsKit::vectorTools::vectorMul(vector<float>* a, float val){
	return toolsKit::vectorTools::vectorMul(val,a);
	}

vector<float>* toolsKit::vectorTools::vectorAdd(vector<float>* a, vector<float>* b){
	vector<float>* ans = new vector<float>();
	vector<float>::iterator ait = a->begin();
	vector<float>::iterator bit = b->begin();

	while(ait != a->end() && bit != b->end()){
		ans->push_back((*ait) + (*bit));
		ait++; bit++;
		}
	return ans;
	}

vector<float>* toolsKit::vectorTools::vectorAdd(int count, ...){
	if (count==0) return NULL;
	va_list vectors;
	va_start(vectors,count);
	vector<float> * elem = va_arg(vectors,vector<float>*);
	vector<float> * ans = new vector<float>();
	//copy the first vector:
	for (vector<float>::iterator it = elem->begin(); it != elem->end(); it++){
		ans->push_back(*it);
	}
	if(count ==1) return ans;
	//add the rest of the vectors:
	for (int i=1; i< count; i++){
		elem = va_arg(vectors,vector<float>*);
		for (vector<float>::iterator it = ans->begin(), it1 = elem->begin(); it != ans->end() && it1 != elem->end(); it++, it1++){
			*it += *it1;
		}
	}
	va_end(vectors);
	return ans;
}

vector<float>* toolsKit::vectorTools::vectorAdd(vector<float>* a, float val){
	vector<float>* ans = new vector<float>();
	for (vector<float>::iterator it = a->begin(); it != a->end(); it++)
		ans->push_back(*it+val);
	return ans;
	}

vector<float>* toolsKit::vectorTools::vectorAdd(float val, vector<float>* b){
	return toolsKit::vectorTools::vectorAdd(b,val);
	}


vector<float>* toolsKit::vectorTools::elementsFromIpl(IplImage* I, vector<float> * p){
		vector<float>* ans = new vector<float>();
		vector<float>* colVect = toolsKit::IplImageToCoulmnVector(I);
		for(vector<float>::iterator it = p->begin(); it!= p->end(); it++)
			ans->push_back((*colVect)[(*it)-1]);
		return ans;
	}

//A + B * (C - 1)
vector<float> * toolsKit::vectorTools::elementsForIpl(vector<float> * A, int B, vector<float> * C){
		vector<float> * args = new vector<float>();
		args = toolsKit::vectorTools::vectorSub(C,1);
		args = toolsKit::vectorTools::vectorMul(B,args);
		args = toolsKit::vectorTools::vectorAdd(A,args);
		return args;
	}

void toolsKit::vectorTools::vectorToFile(vector<float>* v, string filename){
	ofstream thefile(filename.c_str(),ios::out & ios::trunc);
	for (vector<float>::iterator it = v->begin(); it != v->end();it++)
		thefile<<*it<<endl;

	}	

