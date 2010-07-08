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
	if(ratio>0.99 || ratio<0.4)
		ratio=0.75;
	// first decide how many levels
	nLevels=log((double)minWidth/image.width)/log(ratio);
	if(ImPyramid!=NULL)
		delete []ImPyramid;


	nLevels=1;

	ImPyramid=new IplImage*[nLevels];
	
	ImPyramid[0] = cvCreateImage( cvSize( image.width,image.height ),image.depth, image.nChannels );
	ImPyramid[0]=cvCloneImage(&image);
	


	/*

	cout<<"nlevels:"<<nLevels<<endl;
	for(int i=1;i<nLevels;i++)
	{
		IplImage* foo=cvCreateImage( cvSize( ImPyramid[i-1]->width,ImPyramid[i-1]->height ),ImPyramid[i-1]->depth, ImPyramid[i-1]->nChannels );
		 double altRatio=pow( ratio, i );		
		cvSmooth(ImPyramid[i-1],foo,2);
		ImPyramid[i] = cvCreateImage( cvSize( foo->width*ratio,foo->height*ratio ),foo->depth, foo->nChannels );
		cvResize(foo,ImPyramid[i]);
	}
*/
}
