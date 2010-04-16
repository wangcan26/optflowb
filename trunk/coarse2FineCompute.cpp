#include "coarse2FineCompute.h"

coarse2FineCompute::coarse2FineCompute(void)
{
}

coarse2FineCompute::~coarse2FineCompute(void)
{
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
	IplImage* Image1=NULL;
	IplImage* Image2=NULL;
	IplImage* WarpImage2=NULL;

for(int k=Pyramid1.nlevels()-1;k>=0;k--)
	{		
		cout<<"Pyramid level "<<k<<"-";
		
		int width=Pyramid1.getImageFromPyramid(k)->width;
		int height=Pyramid1.getImageFromPyramid(k)->height;
		int depth=Pyramid1.getImageFromPyramid(k)->depth;
		int nChannels=Pyramid1.getImageFromPyramid(k)->nChannels;
		cout<<"width:"<<width<<"  height:"<<height<<endl;
		/*
		im2feature(Image1,Pyramid1.Image(k));
		im2feature(Image2,Pyramid2.Image(k));
*/

		//on top level(first iteration)
		if(k==Pyramid1.nlevels()-1)
		{
			cout<<"first iteration:"<<k<<endl;
			vx=cvCreateImage(cvSize(width,height ),depth,nChannels);
			vy=cvCreateImage(cvSize(width,height ),depth,nChannels);		
			//clone image2 to warpImage2
		//	WarpImage2 = cvCreateImage(cvSize(Image2->width,Image2->height ),Image2->depth, Image2->nChannels );
		//	WarpImage2=  cvCloneImage(Image2);
		}
		else
		{
			
			

			IplImage *tempVx = cvCreateImage(cvSize(width, height), 
											   vx->depth,
											   vx->nChannels);
			cvResize(vx, tempVx); 
			vx=tempVx;

			IplImage *tempVy = cvCreateImage(cvSize(width, height), 
											   vy->depth,
											   vy->nChannels);
			cvResize(vy, tempVy); 
			vy=tempVy;
			//toolsKit::cvShowManyImages("Image", 2, Pyramid1.getImageFromPyramid(k),Pyramid1.getImageFromPyramid(k));
			toolsKit::cvShowManyImages("Image", 2,vx,vy);

			
			//vx.imresize(width,height);
			/*vx.Multiplywith(1/ratio);
			vy.imresize(width,height);
			vy.Multiplywith(1/ratio);			
			warpFL(WarpImage2,Image1,Image2,vx,vy);
			*/
		}						
		 // SmoothFlowPDE(Image1,Image2,WarpImage2,vx,vy,alpha,nOuterFPIterations,nInnerFPIterations,nCGIterations);	
	}
	//warpFL(warpI2,Im1,Im2,vx,vy);
}