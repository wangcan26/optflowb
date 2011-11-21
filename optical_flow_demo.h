#ifndef OPTICAL_FLOW_DEMO
#define OPTICAL_FLOW_DEMO

#include <stdio.h>
#include <cv.h>
#include <highgui.h>
#include <math.h>

#define M_PI 3.1415926535897932384626433832795


// the "official" threshold - if the absolute value of either 
// flow component is greater, it's considered unknown
#define UNKNOWN_FLOW_THRESH 1e9

// value to use to represent unknown flow
#define UNKNOWN_FLOW 1e10


void MotionToColor(CvMat* velx, CvMat* vely, IplImage* color_img, float maxmotion);
void computeColor(float fx, float fy, IplImage* color_img,  int x, int y);
void makecolorwheel();
void setcols(int r, int g, int b, int k);


inline bool unknown_flow(float u, float v) {
    return (fabs(u) >  UNKNOWN_FLOW_THRESH) 
	|| (fabs(v) >  UNKNOWN_FLOW_THRESH) ;
//yair	|| isnan(u) || isnan(v);
}

inline bool unknown_flow(float *f) {
    return unknown_flow(f[0], f[1]);
}


/* This is just an inline that allocates images.  I did this to reduce clutter in the
 * actual computer vision algorithmic code.  Basically it allocates the requested image
 * unless that image is already non-NULL.  It always leaves a non-NULL image as-is even
 * if that image's size, depth, and/or channels are different than the request.
 */
inline static void allocateOnDemand( IplImage **img, CvSize size, int depth, int channels )
{
	if ( *img != NULL )	return;

	*img = cvCreateImage( size, depth, channels );
	if ( *img == NULL )
	{
		fprintf(stderr, "Error: Couldn't allocate image.  Out of memory?\n");
		exit(-1);
	}
}





#endif
