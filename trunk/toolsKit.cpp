#include "toolsKit.h"

toolsKit::toolsKit(void)
{
}

toolsKit::~toolsKit(void)
{
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
    IplImage *DispImage;

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
	  	  DispImage = cvCreateImage( cvSize(100 + size*w, 60 + size*h), 8, 1 );
	  else
		  DispImage = cvCreateImage( cvSize(100 + size*w, 60 + size*h), 8, 3 );

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


void toolsKit::opt_flow_lk(){
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

