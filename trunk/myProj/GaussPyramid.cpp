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

void GaussPyramid::ConstructPyramid(const IplImage* image, double ratio, int minWidth)
{
	
	// the ratio cannot be arbitrary numbers
	//if(ratio>0.99 || ratio<0.4)
		ratio=0.75;
	// first decide how many levels
	//nLevels=log((double)minWidth/image.width)/log(ratio);
	if(ImPyramid!=NULL)
		delete []ImPyramid;


	//minWidth=1;

	ImPyramid=new IplImage*[minWidth];
	
	ImPyramid[minWidth-1]=cvCreateImage( cvSize( image->width,image->height ),image->depth, image->nChannels );
	ImPyramid[minWidth-1]=cvCloneImage(image);
	


	

	cout<<"nlevels:"<<minWidth<<endl;
	for(int i=0;i<minWidth-1;i++)
	{
		IplImage* smoothImage=cvCreateImage( cvSize( image->width,image->height ),image->depth, image->nChannels );
		double altRatio=pow( ratio, minWidth-i);		
		cvSmooth(image,smoothImage,2);
		int reduceWidth=smoothImage->width*altRatio+1;
		int reduceHeight=smoothImage->height*altRatio+1 ;
		ImPyramid[i] = cvCreateImage( cvSize(reduceWidth ,reduceHeight),smoothImage->depth, smoothImage->nChannels );
		cvResize(smoothImage,ImPyramid[i]);
	}

}
