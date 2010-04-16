#include "GaussPyramid.h"
#include "math.h"

GaussPyramid::GaussPyramid(void)
{
	ImPyramid=NULL;
}

GaussPyramid::~GaussPyramid(void)
{
	if(ImPyramid!=NULL)
		delete []ImPyramid;
}

void GaussPyramid::ConstructPyramid(const IplImage &image, double ratio, int minWidth)
{
	// the ratio cannot be arbitrary numbers
	if(ratio>0.98 || ratio<0.4)
		ratio=0.75;
	// first decide how many levels
	nLevels=log((double)minWidth/image.width)/log(ratio);
	if(ImPyramid!=NULL)
		delete []ImPyramid;
	ImPyramid=new IplImage*[nLevels];
	

	ImPyramid[0] = cvCreateImage( cvSize( image.width,image.height ),image.depth, image.nChannels );
	ImPyramid[0]=cvCloneImage(&image);
/* create destination image */
	cvNamedWindow("TEST",  CV_WINDOW_AUTOSIZE); 
	//cvShowImage("TEST",ImPyramid[0]); 
	//cvWaitKey(0); 
	cvDestroyWindow("TEST"); 
	
/* copy from source to dest */
//cvCopy( src, dst, NULL );
	//ImPyramid[0].copyData(image);
	
	double baseSigma=(1/ratio-1);
	int n=log(0.25)/log(ratio);
	double nSigma=baseSigma*n;

	cout<<"nlevels:"<<nLevels<<endl;
	for(int i=1;i<nLevels;i++)
	{
		IplImage* foo=cvCreateImage( cvSize( ImPyramid[i-1]->width,ImPyramid[i-1]->height ),ImPyramid[i-1]->depth, ImPyramid[i-1]->nChannels );
		
		/*if(i<=n)
		{*/
			double sigma=baseSigma*i;
			//image.GaussSmoothing(foo,sigma,sigma*3);

			//cvSmooth(&image,&foo,2,3,0,sigma,sigma*3);
			cvSmooth(ImPyramid[i-1],foo,2);
			
			//cvNamedWindow("smooth",  CV_WINDOW_AUTOSIZE); 
			//cvShowImage("smooth",foo); 
			//cvWaitKey(0); 

			ImPyramid[i] = cvCreateImage( cvSize( foo->width*ratio,foo->height*ratio ),foo->depth, foo->nChannels );

			cvResize(foo,ImPyramid[i]);
			//foo.imresize(ImPyramid[i],pow(ratio,i));
			
			//cvNamedWindow("resize",  CV_WINDOW_AUTOSIZE); 
			//cvShowImage("resize",ImPyramid[i]); 
			//cvWaitKey(0); 
			//cvDestroyWindow("resize"); 
			//cvDestroyWindow("smooth"); 

		/*}
		else
		{
		/*	cvSmooth(&ImPyramid[i-n],&foo,2,3,0,nSigma,nSigma*3);
			//ImPyramid[i-n].GaussSmoothing(foo,nSigma,nSigma*3);
			double rate=(double)pow(ratio,i)*image.width/foo->width;

			ImPyramid[i] = cvCreateImage( cvSize( image.width*ratio,image.height*ratio ),image.depth, image.nChannels );
			cvResize(&foo,ImPyramid[i]);
			//foo.imresize(ImPyramid[i],rate);
		}*/
	}
}

/*void GaussPyramid::displayTop(const char *filename)
{
//	ImPyramid[nLevels-1].imwrite(filename);
}*/