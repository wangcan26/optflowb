#include "coarse2FineCompute.h"
coarse2FineCompute::coarse2FineCompute(void)
{
}

coarse2FineCompute::~coarse2FineCompute(void)
{
}



int sobelImplementation() {

/*If you do to find Explanations for some functions i'd suggest you

that you check out the older posts of this blog*/

int height,width,step,channels;
uchar *data,*data2;
int i,j,k;

//CvCapture* capture = cvCaptureFromCAM( CV_CAP_ANY );
//IplImage* frame = cvQueryFrame( capture );
IplImage* src1 = cvLoadImage("c:\\a\\car1.jpg",0);/*This is the image where a mono image is stored*/
IplImage* src2 = cvCreateImage(cvGetSize(src1), 8, 1 );/*This is the image
/*if( !capture ) {
fprintf( stderr, "ERROR: capture is NULL n" );
getchar();
exit(0);
}*/
// Create a window in which the captured images will be presented

cvNamedWindow( "mywindow", CV_WINDOW_AUTOSIZE );
// Show the image captured from the camera in the window and repeat
/*
while( 1 ) {
// Get one frame

IplImage* frame = cvQueryFrame( capture );

if( !frame ) {
fprintf( stderr, "ERROR: frame is null...n" );
getchar();

break;
}*/
height = src1 -> height;
width  = src1 -> width;
step   = src1 -> widthStep;
channels = src1->nChannels;
data = (uchar *)src1->imageData;
data2 = (uchar *)src2->imageData;

/*Just try to figure out why i did not use the third for Loop..?
Well this is not a perfect way to get a monochrome Image of a color Image but this is a ver good Example of how to

play with Images */
for(i=0;i< (height);i++)
	for(j=0;j< width; j++)
		data2[i*(src1->widthStep)+j*(src1->nChannels)]=data[(height-i)*step+j*channels];
/*Here i used height-i because the image is fliped vertically.reflipping it removes the Error
We can use a flipping function from OpenCV but since a small manipulation is working why do it any other way..*/
cvShowImage( "mywindow", src1 );

/* Do not release the frame!
//If ESC key pressed, Key=0x10001B under OpenCV 0.9.7(linux version),

//remove higher bits using AND operator*/
//if( (cvWaitKey(10) && 255) == 27 )) break;

//}

/* Release the capture device housekeeping*/

//cvReleaseCapture( &capture );
//cvDestroyWindow( "mywindow" );
return 0;

}


int test0() 
{
	char* fileAddress="c:\\a\\Dumptruck1.png";
	IplImage* orginalImage = cvLoadImage(fileAddress,0);
	cvNamedWindow("Orginal Image");
	cvShowImage("Orginal Image", orginalImage);

	IplImage* edgeImage =	cvCreateImage(cvGetSize(orginalImage),IPL_DEPTH_16S,1);
	cvSobel(orginalImage,edgeImage,0,1,3);
	cvNamedWindow("Edge Image");
	cvShowImage("Edge Image", edgeImage);

	
	cvWaitKey(0);
	cvReleaseImage(&orginalImage);
	cvReleaseImage(&edgeImage);
	cvDestroyWindow("orginal Image");
	cvDestroyWindow("Edge Image");

	return 0;
}


int testSobel(){
	char* fileAddress="c:\\a\\Dumptruck1.png";
	IplImage* orginalImage = cvLoadImage(fileAddress,0);
	//IplImage* dest_dx,dest_dy,df_dx,df_dy;//: IPL_DEPTH_8U image.

	/*create temp images*/
	IplImage* df_dx = cvCreateImage(cvGetSize(orginalImage),IPL_DEPTH_16S,1);
	IplImage* df_dy = cvCreateImage(cvGetSize(orginalImage),IPL_DEPTH_16S,1);
	
	IplImage* dest_dx=cvCreateImage(cvGetSize(orginalImage),IPL_DEPTH_8U,1); 
	IplImage* dest_dy=cvCreateImage(cvGetSize(orginalImage),IPL_DEPTH_8U,1); 
	/* use sobel to find derivatives */
	cvSobel( orginalImage, df_dx, 1, 0, 3);
	cvSobel( orginalImage, df_dy, 0, 1, 3);

	/* Convert signed to unsigned 8*/
	cvConvertScaleAbs( df_dx , dest_dx, 1, 0);
	cvConvertScaleAbs( df_dy , dest_dy, 1, 0); 
	cvNamedWindow("dx");
	cvShowImage("dx", dest_dx);
	cvNamedWindow("dy");
	cvShowImage("dy", dest_dy);

return 0;
}






IplImage* coarse2FineCompute::LaplaceCompute(IplImage* input,IplImage* input2){
	IplImage* output= cvCreateImage(cvSize(input->width, input->height), input->depth, input->nChannels);;
	cvLaplace( input, input2, 3 );
	
	return output;
}



void coarse2FineCompute::Coarse2FineFlow(IplImage* vx, 
										 IplImage* vy, 
										 IplImage &warpI2,
										 const IplImage &Im1, 
										 const IplImage &Im2, 
										 double alpha, 
										 double ratio, 
										 int minWidth,
										 int nOuterFPIterations, 
										 int nInnerFPIterations, 
										 int nCGIterations)
{
	// first build the pyramid of the two images
	GaussPyramid Pyramid1;
	GaussPyramid Pyramid2;		
	cout<<"Constructing pyramid..."<<endl;
	Pyramid1.ConstructPyramid(Im1,ratio,minWidth);


	Pyramid2.ConstructPyramid(Im2,ratio,minWidth);
	cout<<"done!"<<endl;
	
	// now iterate from the top level to the bottom
	//IplImage* Image1=NULL;
	//IplImage* Image2=NULL;
	IplImage* WarpImage2=NULL;

		//opt_flow_lk();

	testSobel();
	return ;
for(int k=Pyramid1.nlevels()-1;k>=0;k--)
	{		
		cout<<"Pyramid level "<<k<<"-";
		
		int width=Pyramid1.getImageFromPyramid(k)->width;
		int height=Pyramid1.getImageFromPyramid(k)->height;
		int depth=Pyramid1.getImageFromPyramid(k)->depth;
		int nChannels=Pyramid1.getImageFromPyramid(k)->nChannels;
		cout<<"width:"<<width<<"  height:"<<height<<endl;

		//on top level(first iteration)
		if(k==Pyramid1.nlevels()-1)
		{
			cout<<"first iteration:"<<k<<endl;
			vx=cvCreateImage(cvSize(width,height),depth,nChannels);
			vy=cvCreateImage(cvSize(width,height),depth,nChannels);		
			//clone image2 to warpImage2
			WarpImage2 = cvCreateImage(cvSize(Pyramid2.getImageFromPyramid(k)->width,Pyramid2.getImageFromPyramid(k)->height ),Pyramid2.getImageFromPyramid(k)->depth, Pyramid2.getImageFromPyramid(k)->nChannels );
			WarpImage2=  cvCloneImage(Pyramid2.getImageFromPyramid(k));
		}
		else
		{
			IplImage *tempVx = cvCreateImage(cvSize(width, height), vx->depth, vx->nChannels);
			cvResize(vx, tempVx); 
			vx=tempVx;
			//vx.Multiplywith(1/ratio);
			IplImage *tempVy = cvCreateImage(cvSize(width, height), vy->depth, vy->nChannels);
			cvResize(vy, tempVy); 
			vy=tempVy;
			//vy.Multiplywith(1/ratio);

			//warpFL(WarpImage2,Pyramid1.getImageFromPyramid(k),Pyramid2.getImageFromPyramid(k),vx,vy);

			toolsKit::cvShowManyImages("Image",4, Pyramid1.getImageFromPyramid(k),Pyramid2.getImageFromPyramid(k),vx,vy);

			//IplImage* out=LaplaceCompute(Pyramid1.getImageFromPyramid(k),Pyramid2.getImageFromPyramid(k));
			//toolsKit::cvShowManyImages("Image",1, out);
		}						
		SmoothFlowPDE( Pyramid1.getImageFromPyramid(k),Pyramid2.getImageFromPyramid(k),WarpImage2,vx,vy,alpha,nOuterFPIterations,nInnerFPIterations,nCGIterations);	
	}
	//warpFL(WarpImage2,Pyramid1.getImageFromPyramid(k),Pyramid2.getImageFromPyramid(k),vx,vy);
}

//--------------------------------------------------------------------------------------------------------
// function to compute optical flow field using two fixed point iterations
// Input arguments:
// Im1, Im2:				    frame 1 and frame 2
// warpIm2:						the warped frame 2 according to the current flow field u and v
// u,v:							the current flow field, NOTICE that they are also output arguments
//	
//--------------------------------------------------------------------------------------------------------
void coarse2FineCompute::SmoothFlowPDE(const IplImage* Im1, 
									   const IplImage* Im2, 
									   IplImage* warpIm2, 
									   IplImage* u, 
									   IplImage* v, 
									   double alpha, 
									   int nOuterFPIterations, 
									   int nInnerFPIterations, 
									   int nCGIterations)
{
	IplImage* mask,imdx,imdy,imdt;
	int width,height,depth,nChannels,nPixels;
	width=Im1->width;
	height=Im1->height;
	depth=Im1->depth;
	nChannels=Im1->nChannels;
	nPixels=width*height;

	
	IplImage* du=cvCreateImage(cvSize(width,height),depth,NULL);
	IplImage* dv=cvCreateImage(cvSize(width,height),depth,NULL);
	IplImage* uu=cvCreateImage(cvSize(width,height),depth,NULL);
	IplImage* vv=cvCreateImage(cvSize(width,height),depth,NULL);
	IplImage* ux=cvCreateImage(cvSize(width,height),depth,NULL);
	IplImage* uy=cvCreateImage(cvSize(width,height),depth,NULL);
	IplImage* vx=cvCreateImage(cvSize(width,height),depth,NULL);
	IplImage* vy=cvCreateImage(cvSize(width,height),depth,NULL);
	IplImage* Phi_1st=cvCreateImage(cvSize(width,height),depth,NULL);
	IplImage* Psi_1st=cvCreateImage(cvSize(width,height),depth,nChannels);

	IplImage* imdxy,imdx2,imdy2,imdtdx,imdtdy;
	IplImage* ImDxy,ImDx2,ImDy2,ImDtDx,ImDtDy;
	IplImage* A11,A12,A22,b1,b2;
	IplImage* foo1,foo2;

	// variables for conjugate gradient
	//DImage r1,r2,p1,p2,q1,q2;
	//double* rou;
	//rou=new double[nCGIterations];

	//double varepsilon_phi=pow(0.001,2);
	//double varepsilon_psi=pow(0.001,2);

	//--------------------------------------------------------------------------
	// the outer fixed point iteration
	//--------------------------------------------------------------------------
	for(int count=0;count<nOuterFPIterations;count++)
	{
		/*
		// compute the gradient
		getDxs(imdx,imdy,imdt,Im1,warpIm2);

		// generate the mask to set the weight of the pix	els moving outside of the image boundary to be zero
		genInImageMask(mask,vx,vy);

		// set the derivative of the flow field to be zero
		du.reset();
		dv.reset();
	*/
		//--------------------------------------------------------------------------
		// the inner fixed point iteration
		//--------------------------------------------------------------------------
		for(int hh=0;hh<nInnerFPIterations;hh++)
		{
			/*
			// compute the derivatives of the current flow field
			if(hh==0)
			{
				uu.copyData(u);
				vv.copyData(v);
			}
			else
			{
				uu.Add(u,du);
				vv.Add(v,dv);
			}
			uu.dx(ux);
			uu.dy(uy);
			vv.dx(vx);
			vv.dy(vy);

			// compute the weight of phi
			Phi_1st.reset();
			double* phiData=Phi_1st.data();
			double temp;
			const double *uxData,*uyData,*vxData,*vyData;
			uxData=ux.data();
			uyData=uy.data();
			vxData=vx.data();
			vyData=vy.data();
			for(int i=0;i<nPixels;i++)
			{
				temp=uxData[i]*uxData[i]+uyData[i]*uyData[i]+vxData[i]*vxData[i]+vyData[i]*vyData[i];
				phiData[i]=1/(2*sqrt(temp+varepsilon_phi));
			}

			// compute the nonlinear term of psi
			Psi_1st.reset();
			double* psiData=Psi_1st.data();
			const double *imdxData,*imdyData,*imdtData;
			const double *duData,*dvData;
			imdxData=imdx.data();
			imdyData=imdy.data();
			imdtData=imdt.data();
			duData=du.data();
			dvData=dv.data();
		
			double _a  = 10000, _b = 0.1;
			if(nChannels==1)
			{
				for(int i=0;i<nPixels;i++)
				{
					temp=imdtData[i]+imdxData[i]*duData[i]+imdyData[i]*dvData[i];
					//if(temp*temp<0.04)
					psiData[i]=1/(2*sqrt(temp*temp+varepsilon_psi));
					//psiData[i] = _a*_b/(1+_a*temp*temp);
				}
			}
			else
			{
				for(int i=0;i<nPixels;i++)
					for(int k=0;k<nChannels;k++)
					{
						int offset=i*nChannels+k;
						temp=imdtData[offset]+imdxData[offset]*duData[i]+imdyData[offset]*dvData[i];
						//if(temp*temp<0.04)
						psiData[offset]=1/(2*sqrt(temp*temp+varepsilon_psi));
						//psiData[offset] =  _a*_b/(1+_a*temp*temp);
					}
			}

			// prepare the components of the large linear system
			ImDxy.Multiply(Psi_1st,imdx,imdy);
			ImDx2.Multiply(Psi_1st,imdx,imdx);
			ImDy2.Multiply(Psi_1st,imdy,imdy);
			ImDtDx.Multiply(Psi_1st,imdx,imdt);
			ImDtDy.Multiply(Psi_1st,imdy,imdt);

			if(nChannels>1)
			{
				ImDxy.collapse(imdxy);
				ImDx2.collapse(imdx2);
				ImDy2.collapse(imdy2);
				ImDtDx.collapse(imdtdx);
				ImDtDy.collapse(imdtdy);
			}
			else
			{
				imdxy.copyData(ImDxy);
				imdx2.copyData(ImDx2);
				imdy2.copyData(ImDy2);
				imdtdx.copyData(ImDtDx);
				imdtdy.copyData(ImDtDy);
			}

			// filtering
			imdx2.smoothing(A11,3);
			imdxy.smoothing(A12,3);
			imdy2.smoothing(A22,3);

			// add epsilon to A11 and A22
			A11.Add(alpha*0.1);
			A22.Add(alpha*0.1);

			// form b
			imdtdx.smoothing(b1,3);
			imdtdy.smoothing(b2,3);
			// laplacian filtering of the current flow field
		    Laplacian(foo1,u,Phi_1st);
			Laplacian(foo2,v,Phi_1st);
			double *b1Data,*b2Data;
			const double *foo1Data,*foo2Data;
			b1Data=b1.data();
			b2Data=b2.data();
			foo1Data=foo1.data();
			foo2Data=foo2.data();

			for(int i=0;i<nPixels;i++)
			{
				b1Data[i]=-b1Data[i]-alpha*foo1Data[i];
				b2Data[i]=-b2Data[i]-alpha*foo2Data[i];
			}
*/
			//-----------------------------------------------------------------------
			// conjugate gradient algorithm
			//-----------------------------------------------------------------------
			//r1.copyData(b1);
			//r2.copyData(b2);
			//du.reset();
			//dv.reset();
			//-----------------------------------------------------------------------
			// end of conjugate gradient algorithm
			//-----------------------------------------------------------------------
		}// end of inner fixed point iteration
		
		/*
		u.Add(du,1);
		v.Add(dv,1);
		warpFL(warpIm2,Im1,Im2,u,v);
	*/
	}// end of outer fixed point iteration
	
	
	
	
}


void coarse2FineCompute::opt_flow_lk(){
	//int MAX_CORNERS = 500;
		// Load two images and allocate other structures
	IplImage* imgA = cvLoadImage("c:\\a\\Dumptruck1.png", CV_LOAD_IMAGE_GRAYSCALE);
	IplImage* imgB = cvLoadImage("c:\\a\\Dumptruck2.png", CV_LOAD_IMAGE_GRAYSCALE);

	CvSize img_sz = cvGetSize( imgA );
	int win_size = 15;

	IplImage* imgC = cvLoadImage("c:\\a\\Dumptruck_of.png", CV_LOAD_IMAGE_UNCHANGED);

	// Get the features for tracking
	IplImage* eig_image = cvCreateImage( img_sz, IPL_DEPTH_32F, 1 );
	IplImage* tmp_image = cvCreateImage( img_sz, IPL_DEPTH_32F, 1 );

	int corner_count = 500;
	CvPoint2D32f* cornersA = new CvPoint2D32f[500];

	cvGoodFeaturesToTrack( imgA, eig_image, tmp_image, cornersA, &corner_count,
		0.05, 5.0, 0, 3, 0, 0.04 );

	cvFindCornerSubPix( imgA, cornersA, corner_count, cvSize( win_size, win_size ),
		cvSize( -1, -1 ), cvTermCriteria( CV_TERMCRIT_ITER | CV_TERMCRIT_EPS, 20, 0.03 ) );

	// Call Lucas Kanade algorithm
	char features_found[500];
	float feature_errors[500];

	CvSize pyr_sz = cvSize( imgA->width+8, imgB->height/3 );

	IplImage* pyrA = cvCreateImage( pyr_sz, IPL_DEPTH_32F, 1 );
	IplImage* pyrB = cvCreateImage( pyr_sz, IPL_DEPTH_32F, 1 );

	CvPoint2D32f* cornersB = new CvPoint2D32f[500];

	cvCalcOpticalFlowPyrLK( imgA, imgB, pyrA, pyrB, cornersA, cornersB, corner_count, 
		cvSize( win_size, win_size ), 5, features_found, feature_errors,
		 cvTermCriteria( CV_TERMCRIT_ITER | CV_TERMCRIT_EPS, 20, 0.3 ), 0 );

	// Make an image of the results

	for( int i=0; i<500;i++ )
		{
		/*	printf("Error is %f/n", feature_errors[i]);
			continue;
		}*/
		//printf("Got it/n");
		CvPoint p0 = cvPoint( cvRound( cornersA[i].x ), cvRound( cornersA[i].y ) );
		CvPoint p1 = cvPoint( cvRound( cornersB[i].x ), cvRound( cornersB[i].y ) );
		cvLine( imgC, p0, p1, CV_RGB(255,0,0), 2 );
	}

	//cvNamedWindow( "ImageA", 0 );
	//cvNamedWindow( "ImageB", 0 );
	//cvNamedWindow( "LKpyr_OpticalFlow", 0 );

	cvShowImage( "ImageA", imgA );
	cvShowImage( "ImageB", imgB );
	cvShowImage( "LKpyr_OpticalFlow", imgC );

	cvWaitKey(0);
	//return 0;
}



/////////////////////////////////////////////////tests///////////////////////////////////////
