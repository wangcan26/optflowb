#include "toolsKit.h"

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
		IplImageIterator<float> it2(img);
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
	//img2 is the shifted image
	void toolsKit::IPL_operate_top(IplImage* img,IplImage* shiftImg2,IplImage* dest,toolsKit::operations operation){
			int i,j;
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
			int i,j;
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
			int i;
			int width=img->width;
			if (img->nChannels==3)
				for (i = 1; i < img->width*img->height*3; i+=3)
				{
				// 
				}
			else				
				for (i = 1; i < width*img->height; i++)					
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
			int i,k;
			int width=img->width;					
			//for first cell
			k=width-1;
			for (i = 0; i < width*img->height; i++)
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
	void toolsKit::IPL_print(const IplImage *image) {
		int nl= image->height; // number of lines
		int nc= image->width * image->nChannels; // total number of element per line
		int step= image->widthStep; // effective width
		// get the pointer to the image buffer
		unsigned char *data= reinterpret_cast<unsigned char *>(image->imageData);
		cout<<"=============================width:"<<image->width <<" height:"<< image->height<<" channels:"<<image->nChannels<<"============================="<<endl;
		for (int i=0; i<nl; i++) {

			for (int j=0; j<nc; j+= image->nChannels) {
				CvScalar sca= cvGet2D(image,i,j);
				// process each pixel ---------------------
				// data[j]= data[j]/div * div + div/2;
				cout<<sca.val[0]<<" ";
				//  data[j+1]= data[j+1]/div * div + div/2;
				// data[j+2]= data[j+2]/div * div + div/2;
				// end of pixel processing ----------------
			} // end of line          
			cout<<endl;
			data+= step;  // next line
			//break;
		}
		cout<<"======================================================================================================================================="<<endl;
	}

	void toolsKit::cvMulScalar(IplImage* img,float scalar){
		IplImageIterator<float> it(img);
		while (!it) {      
			*it= (float)scalar*(*it); 
			++it;
		}
	}

	
	void toolsKit::cvNormalizeEdges(IplImage* img){
			int i,k;
			int width=img->width;					
			//for first cell
			k=width-1;
			for (i = 0; i < width*img->height; i++)
			{
				if(!k)//apply to last column only										  
					 ((float*)img->imageData)[i]=0;
							 				
				k==0?k=width-1:k--;
			}

			/*
			const double x_max = -log(FLT_EPSILON);		
			for (int i = 0; i < img->width*img->height; i++)					
				{	
					float t=img->imageData[i]+img->imageData[i];
					//if(((float*)img->imageData)[i] > x_max); 
					//	((float*)img->imageData)[i]=0;
					if(t==((float*)img->imageData)[i])
						((float*)img->imageData)[i]=0;

				}*/

		
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

	void toolsKit::cvNormalizeEdges2(IplImage* img){				
			float nan1 = sqrt(-1.0f);
			for (int i = 0; i < img->width*img->height; i++)					
				{					
					if(AlmostEqualRelativeOrAbsolute(nan1,((float*)img->imageData)[i],0.00001,0.00001))
						((float*)img->imageData)[i]=0;
					if(AlmostEqualRelativeOrAbsolute(FLT_MAX,((float*)img->imageData)[i],0.00001,0.00001))
						((float*)img->imageData)[i]=0;
					if(AlmostEqualRelativeOrAbsolute(FLT_MIN,((float*)img->imageData)[i],0.0001,0.0001))
						((float*)img->imageData)[i]=0;
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
		//double y=1 / (2 * sqrt( x + epsilon ) ) ;
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

	void toolsKit::cvShowManyImages(char* title, int nArgs, ...) {

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
		cvWaitKey();
		cvDestroyWindow(title);
		// End the number of arguments
		va_end(args);
		// Release the Image Memory
		cvReleaseImage(&DispImage);
	}


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
