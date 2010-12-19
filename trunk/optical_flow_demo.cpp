/* --Sparse Optical Flow Demo Program--
 * Written by Yair Adato
 */


#include "optical_flow_demo.h"

static const double pi = 3.14159265358979323846;


//int main__( int argc, char** argv )
//{
//	  /* data structure for the image */
//    IplImage *img1 = 0;
//	IplImage *img2 = 0;
//   
//    /* check for supplied argument */
//    if( argc < 3 ) {
//        fprintf( stderr, "Usage: loadimg <filename1> <filename2>\n" );
//        return 1;
//    }
//   
//    /* load the image,
//       use CV_LOAD_IMAGE_GRAYSCALE or CV_LOAD_IMAGE_COLOR*/
//    img1 = cvLoadImage(argv[1], CV_LOAD_IMAGE_GRAYSCALE );
//	img2 = cvLoadImage(argv[2], CV_LOAD_IMAGE_GRAYSCALE );
//   
//    /* always check */ 
//    if( 0 == img1  || 0 == img2 ) {
//        fprintf( stderr, "Cannot load file %s!\n", argv[1] );
//        return 1;
//    }
//   
//    /* create a window */
//    cvNamedWindow( "yair", CV_WINDOW_AUTOSIZE );
//   
//    /* display the image */
//    cvShowImage( "yair", img1 );
//   
//    /* wait until user press a key */
// //   cvWaitKey(0);
//
//	/* display the image */
//    cvShowImage( "yair", img2 );
//   
//    /* wait until user press a key */
// //   cvWaitKey(0);
//
//
//	CvSize size = cvSize(11 , 11);
//	
//	CvMat* velx = cvCreateMat(img1->height, img1->width, CV_32FC1 );
//	CvMat* vely = cvCreateMat(img1->height, img1->width, CV_32FC1 );
//
//	cvCalcOpticalFlowLK( img1, img2, size, velx, vely );
//	cvShowImage( "yair", velx );
//	
//
//	IplImage* color_img = cvCreateImage( cvSize(velx->height,velx->width), IPL_DEPTH_8U, 3 );
////	MotionToColor( velx,  vely,  color_img,  0.1f);
//	cvShowImage( "yair", color_img );
//	 cvWaitKey(0);
//
//  try {
//    middlebury::CShape sh(img1->height, img1->width, 2);
//	middlebury::CFloatImage img(sh);
//	img.ClearPixels();
//	
//	for (int y = 0; y < velx->height; y++) {
//		const float* ptr = (const float*)(velx->data.ptr + y * velx->step);
//		for (int x = 0; x < velx->width; x++) {
//			img.Pixel(x, y, 0) =  *ptr++;
//		}
//	}
//
//	for (int y = 0; y < vely->height; y++) {
//		const float* ptr = (const float*)(vely->data.ptr + y * vely->step);
//		for (int x = 0; x < vely->width; x++) {
//			img.Pixel(x, y, 1) =  *ptr++;
//		}
//	}
//	char *filename = "test.flo";
//
//	middlebury::WriteFlowFile(img, filename);
//	middlebury::ReadFlowFile(img, filename);
//    }
//  catch (middlebury::CError &err) {
//	fprintf(stderr, err.message);
//	fprintf(stderr, "\n");
//	exit(1);
//    }
//
//   
//
//
//
//
//
//
//
//
//		CvSize frame_size = cvSize( img1->width, img1->height );
//
//
//
//		IplImage  *eig_image = NULL, *temp_image = NULL, *pyramid1 = NULL, *pyramid2 = NULL;
//
//		/* Preparation: Allocate the necessary storage. */
//		allocateOnDemand( &eig_image, frame_size, IPL_DEPTH_32F, 1 );
//		allocateOnDemand( &temp_image, frame_size, IPL_DEPTH_32F, 1 );
//
//		/* Preparation: This array will contain the features found in frame 1. */
//		CvPoint2D32f frame1_features[400];
//
//		/* Preparation: BEFORE the function call this variable is the array size
//		 * (or the maximum number of features to find).  AFTER the function call
//		 * this variable is the number of features actually found.
//		 */
//		int number_of_features;
//		
//		/* I'm hardcoding this at 400.  But you should make this a #define so that you can
//		 * change the number of features you use for an accuracy/speed tradeoff analysis.
//		 */
//		number_of_features = 400;
//
//		/* Actually run the Shi and Tomasi algorithm!!
//		 * "frame1_1C" is the input image.
//		 * "eig_image" and "temp_image" are just workspace for the algorithm.
//		 * The first ".01" specifies the minimum quality of the features (based on the eigenvalues).
//		 * The second ".01" specifies the minimum Euclidean distance between features.
//		 * "NULL" means use the entire input image.  You could point to a part of the image.
//		 * WHEN THE ALGORITHM RETURNS:
//		 * "frame1_features" will contain the feature points.
//		 * "number_of_features" will be set to a value <= 400 indicating the number of feature points found.
//		 */
//		cvGoodFeaturesToTrack(img1, eig_image, temp_image, frame1_features, &number_of_features, .01, .01, NULL);
//
//		/* Pyramidal Lucas Kanade Optical Flow! */
//
//		/* This array will contain the locations of the points from frame 1 in frame 2. */
//		CvPoint2D32f frame2_features[400];
//
//		/* The i-th element of this array will be non-zero if and only if the i-th feature of
//		 * frame 1 was found in frame 2.
//		 */
//		char optical_flow_found_feature[400];
//
//		/* The i-th element of this array is the error in the optical flow for the i-th feature
//		 * of frame1 as found in frame 2.  If the i-th feature was not found (see the array above)
//		 * I think the i-th entry in this array is undefined.
//		 */
//		float optical_flow_feature_error[400];
//
//		/* This is the window size to use to avoid the aperture problem (see slide "Optical Flow: Overview"). */
//		CvSize optical_flow_window = cvSize(3,3);
//		
//		/* This termination criteria tells the algorithm to stop when it has either done 20 iterations or when
//		 * epsilon is better than .3.  You can play with these parameters for speed vs. accuracy but these values
//		 * work pretty well in many situations.
//		 */
//		CvTermCriteria optical_flow_termination_criteria
//			= cvTermCriteria( CV_TERMCRIT_ITER | CV_TERMCRIT_EPS, 20, .3 );
//
//		/* This is some workspace for the algorithm.
//		 * (The algorithm actually carves the image into pyramids of different resolutions.)
//		 */
//     	allocateOnDemand( &pyramid1, frame_size, IPL_DEPTH_8U, 1 );
//		allocateOnDemand( &pyramid2, frame_size, IPL_DEPTH_8U, 1 );
//
//		/* Actually run Pyramidal Lucas Kanade Optical Flow!!
//		 * "frame1_1C" is the first frame with the known features.
//		 * "frame2_1C" is the second frame where we want to find the first frame's features.
//		 * "pyramid1" and "pyramid2" are workspace for the algorithm.
//		 * "frame1_features" are the features from the first frame.
//		 * "frame2_features" is the (outputted) locations of those features in the second frame.
//		 * "number_of_features" is the number of features in the frame1_features array.
//		 * "optical_flow_window" is the size of the window to use to avoid the aperture problem.
//		 * "5" is the maximum number of pyramids to use.  0 would be just one level.
//		 * "optical_flow_found_feature" is as described above (non-zero iff feature found by the flow).
//		 * "optical_flow_feature_error" is as described above (error in the flow for this feature).
//		 * "optical_flow_termination_criteria" is as described above (how long the algorithm should look).
//		 * "0" means disable enhancements.  (For example, the second array isn't pre-initialized with guesses.)
//		 */
//		cvCalcOpticalFlowPyrLK(img1, img2, pyramid1, pyramid2, frame1_features, frame2_features, number_of_features, optical_flow_window, 5, optical_flow_found_feature, optical_flow_feature_error, optical_flow_termination_criteria, 0 );
//	
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//    /* free memory */
//    cvDestroyWindow( "yair" );
//    cvReleaseImage( &img1 );
//	cvReleaseImage( &img2 );
//	cvReleaseMat(&velx);
//	cvReleaseMat(&vely);
//   
//    return 0;
//
//
//}



