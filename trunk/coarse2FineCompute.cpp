#include "coarse2FineCompute.h"

#include <ctime>

coarse2FineCompute::coarse2FineCompute(int imageDepth,double error)
{
	_imageDepth=imageDepth;
	_ERROR_CONST=error;
}

coarse2FineCompute::~coarse2FineCompute(void)
{
}



IplImage** coarse2FineCompute::meshgrid(int cols, int rows){
	IplImage ** ans = new IplImage*[2];
	ans[0] = cvCreateImage(cvSize(cols,rows),IPL_DEPTH_32F,1);
	ans[1] = cvCreateImage(cvSize(cols,rows),IPL_DEPTH_32F,1);
	int count2=1;
	for(int i =0; i<rows; i++){
		int count1 = 1;
		for(int j =0; j<cols; j++){
			cvSet2D(ans[0],i,j,cvScalar(count1++));
			cvSet2D(ans[1],i,j,cvScalar( count2));
			}
		count2++;
		}


	return ans;
	}



IplImage * coarse2FineCompute::warpLayer(IplImage* I, IplImage* u, IplImage * v){
	int height = I->height;
	int width = I->width;
	int nChannels = I->nChannels;
	IplImage ** XY = meshgrid(width, height);
	IplImage * X = XY[0];
	IplImage * Y = XY[1];

	IplImage * Xu = cvCreateImage(cvSize(X->width,X->height),X->depth,X->nChannels);
	cvAdd(X,u,Xu);
	cvReleaseImage(&X);

	IplImage * Yv = cvCreateImage(cvSize(Y->width,Y->height),Y->depth,Y->nChannels);
	cvAdd(Y,v,Yv);
	cvReleaseImage(&Y);

	vector<float>* XI = toolsKit::IplImageToCoulmnVector(Xu);
	vector<float>* YI = toolsKit::IplImageToCoulmnVector(Yv);

	float eM6 = 0.00247875218; //1E-6
	//XI = max(1, min(sx - 1E-6, XI));
	vtools::vectorMin(XI, height-eM6);
	vtools::vectorMax(XI, 1);

	//XI = max(1, min(sx - 1E-6, XI));
	vtools::vectorMin(YI, width-eM6);
	vtools::vectorMax(YI, 1);
	//toolsKit::vectorTools::vectorToFile(XI, "c:\\a\\XI_cpp.txt");
	//toolsKit::vectorTools::vectorToFile(YI, "c:\\a\\YI_cpp.txt");
	//fXI = floor(XI);
	vector<float>* fXI = vtools::vectorFloor(XI);
	//toolsKit::vectorTools::vectorToFile(fXI, "c:\\a\\fXI_cpp.txt");
	//cXI = ceil(XI);
	vector<float>* cXI = vtools::vectorCeil(XI);
	//toolsKit::vectorTools::vectorToFile(cXI, "c:\\a\\cXI_cpp.txt");
	//fYI = floor(YI);
	vector<float>* fYI = vtools::vectorFloor(YI);
	//toolsKit::vectorTools::vectorToFile(fYI, "c:\\a\\fYI_cpp.txt");
	//cYI = ceil(YI);
	vector<float>* cYI = vtools::vectorCeil(YI);
	//toolsKit::vectorTools::vectorToFile(cYI, "c:\\a\\cYI_cpp.txt");

	//alpha_x = XI - fXI;
	vector<float>* alpha_x = vtools::vectorSub(XI, fXI);
	//toolsKit::vectorTools::vectorToFile(alpha_x, "c:\\a\\alpha_x_cpp.txt");
	//alpha_y = YI - fYI;
	vector<float>* alpha_y = vtools::vectorSub(YI, fYI);
	//	toolsKit::vectorTools::vectorToFile(alpha_y, "c:\\a\\alpha_y_cpp.txt");

	//A1 = (1 - alpha_x) .* (1 - alpha_y) .* I(fYI + sy * (fXI - 1))
	//vector<float>*  A1 = (1 - *alpha_x) * (1 - *alpha_y) * (I<<=*fYI+height*(*fXI-1));
	vector<float> * Iargs = vtools::elementsForIpl(fYI,height,fXI);
	vector<float>*  A1 = vtools::vectorMul(vtools::vectorMul((vtools::vectorSub(1 , alpha_x)) , (vtools::vectorSub(1 ,alpha_y))) , vtools::elementsFromIpl(I,Iargs));
	delete Iargs;
	//A2 = alpha_x .* (1 - alpha_y) .* I(fYI + sy * (cXI - 1))

	Iargs = vtools::elementsForIpl(fYI,height,cXI);
	vector<float>*  A2 = vtools::vectorMul(alpha_x, vtools::vectorMul(vtools::vectorSub(1,alpha_y),(vtools::elementsFromIpl(I,Iargs))));
	delete Iargs;

	//A3 = (1 - alpha_x) .* alpha_y .* I(cYI + sy * (fXI - 1))
	Iargs = vtools::elementsForIpl(cYI,height,fXI);
	vector<float> * A3 = vtools::vectorMul(vtools::vectorSub(1 , alpha_x) , vtools::vectorMul(alpha_y , (vtools::elementsFromIpl(I,Iargs))));
	delete Iargs;

	//A4 = alpha_x .* alpha_y .* I(cYI + sy * (cXI - 1))
	Iargs = vtools::elementsForIpl(cYI,height,cXI);
	vector<float>  *A4 = vtools::vectorMul(alpha_x , vtools::vectorMul(alpha_y , (vtools::elementsFromIpl(I,Iargs))));
	delete Iargs;
	vector<float> * O = vtools::vectorAdd(4,A1 , A2 , A3 , A4);
	delete A1; delete A2; delete A3; delete A4;
	//O = reshape(O, sy, sx); -> from vector to IPLImage
	IplImage* IplO = cvCreateImage(cvSize(width, height), I->depth, nChannels);
	toolsKit::ColumnVectorToIplImage(O,IplO);
	delete O;
	delete XI;
	delete YI;
	delete fXI;
	delete cXI;
	delete fYI;
	delete cYI;
	//delete alpha_x;
	//delete alpha_y;
	cvReleaseImage(&X);
	cvReleaseImage(&Y);
	delete[] XY;
	return IplO;
	}

//// Allocate image planes
//	IplImage* r = cvCreateImage( cvGetSize(src), IPL_DEPTH_8U, 1 );
//	IplImage* g = cvCreateImage( cvGetSize(src), IPL_DEPTH_8U, 1 );
//	IplImage* b = cvCreateImage( cvGetSize(src), IPL_DEPTH_8U, 1 );
//
//	// Split image onto the color planes
//	cvSplit( src, r, g, b, NULL );

IplImage* coarse2FineCompute::RGBwarp(IplImage* I, IplImage* u, IplImage* v){

	if(I->nChannels ==1)
		return warpLayer(I,u,v);
	//==3
	
	IplImage* r = cvCreateImage( cvGetSize(I), I->depth, 1 );
	IplImage* g = cvCreateImage( cvGetSize(I), I->depth, 1 );
	IplImage* b = cvCreateImage( cvGetSize(I), I->depth, 1 );

	cvSplit(I,r,g,b,NULL);
	r = warpLayer(r,u,v);
	g = warpLayer(g,u,v);
	b = warpLayer(b,u,v);
	//toolsKit::cvShowManyImages("rgb?",3,r,g,b);

	IplImage * ans = cvCreateImage(cvGetSize(I),I->depth,I->nChannels);
	cvMerge(r,g,b,0,ans);
	cvReleaseImage( &r );
	cvReleaseImage( &g );
	cvReleaseImage( &b );
	toolsKit::cvShowManyImages("I ans",2,I,ans);
	//cvWaitKey(0);
	return ans;
	
	}



//need to recieve cvSobel with function pointer;
int getDXsCV(const IplImage* src1,IplImage* dest_dx,IplImage* dest_dy){		
	double x[7] =	{0.016667,
						-0.15,
						 0.75,
							0,
						-0.75,
						 0.15,
					-0.016667}	;					
	double y[7] =  {0.016667,-0.15, 0.75,0,-0.75,0.15,-0.016667};

	cvZero(dest_dx);
	cvZero(dest_dy);

    //CvPoint point = cvPoint(1,1);
	//x derivative
	CvMat* weickert = &cvMat(1, 7, CV_64FC1, x ); // 64FC1 for double
	cvFilter2D(src1,dest_dx,weickert);
	
	//y derivative
	weickert = &cvMat( 7, 1, CV_64FC1, y );
	cvFilter2D(src1,dest_dy,weickert);
	
	//old div with sobel
	//cvSobel(src, dest_dx, 1, 0, 1);
	//cvSobel( src, dest_dy, 0, 1, 1);
	//fix to fit matlab
	toolsKit::cvMulScalar(dest_dx,-1);
	toolsKit::cvMulScalar(dest_dy,-1);
	//toolsKit::cvMulScalar(dest_dz,-2);

	return 0;
}

flowUV* coarse2FineCompute::Coarse2FineFlow( const IplImage* Im1, 
											 const IplImage* Im2, 
											 double alpha, 
											 double gamma,
											 double ratio, 
											 int minWidth,
											 int nOuterFPIterations, 
											 int nInnerFPIterations)
										
{
	// first build the pyramid of the two images
	GaussPyramid Pyramid1;	
	GaussPyramid Pyramid2;			
	Pyramid1.SetNlevels(minWidth);
	Pyramid2.SetNlevels(minWidth);
	Pyramid1.ConstructPyramid(Im1,ratio,minWidth);
	Pyramid2.ConstructPyramid(Im2,ratio,minWidth);	
	
	//for timing
	std::clock_t start;
	double diff;

	//clone image2 to warpImage2 (in first iter)
	int firstIter=0;
	IplImage* WarpImage2=cvCreateImage(cvSize(Pyramid2.getImageFromPyramid(firstIter)->width,Pyramid2.getImageFromPyramid(firstIter)->height ),Pyramid2.getImageFromPyramid(firstIter)->depth, Pyramid2.getImageFromPyramid(firstIter)->nChannels );
	WarpImage2=  cvCloneImage(Pyramid2.getImageFromPyramid(firstIter));	

	IplImage* vx1=cvCreateImage(cvSize(Pyramid1.getImageFromPyramid(firstIter)->width,Pyramid1.getImageFromPyramid(firstIter)->height),Pyramid1.getImageFromPyramid(firstIter)->depth,1);
	IplImage* vy1=cvCreateImage(cvSize(Pyramid1.getImageFromPyramid(firstIter)->width,Pyramid1.getImageFromPyramid(firstIter)->height),Pyramid1.getImageFromPyramid(firstIter)->depth,1);
	cvZero(vx1);
	cvZero(vy1);
	//flow result
	flowUV* UV=new flowUV(vx1,vy1);

	// now iterate from the top level to the bottom
	for(int k=0;k<Pyramid1.getNlevels();k++)
	{		
		cout<<"==================================Pyramid level "<<k<<"-";
		int width =Pyramid1.getImageFromPyramid(k)->width;
		int height=Pyramid1.getImageFromPyramid(k)->height;
		int depth =Pyramid1.getImageFromPyramid(k)->depth;
		int nChannels=Pyramid1.getImageFromPyramid(k)->nChannels;		
		cout<<"width:"<<width<<"  height:"<<height<<"============================================"<<endl;

		
		IplImage *pyramid1_gray = cvCreateImage(cvSize(Pyramid1.getImageFromPyramid(k)->width, Pyramid1.getImageFromPyramid(k)->height),Pyramid1.getImageFromPyramid(firstIter)->depth, 1);
		IplImage *pyramid2_gray = cvCreateImage(cvSize(Pyramid1.getImageFromPyramid(k)->width, Pyramid1.getImageFromPyramid(k)->height), Pyramid1.getImageFromPyramid(firstIter)->depth, 1);

		cvCvtColor( Pyramid1.getImageFromPyramid(k), pyramid1_gray, CV_BGR2GRAY );
		cvCvtColor( Pyramid2.getImageFromPyramid(k), pyramid2_gray, CV_BGR2GRAY );
		//on all levels but the first one
		if (k!=0){
			IplImage *tempVx = cvCreateImage(cvSize(width, height), UV->getU()->depth, UV->getU()->nChannels);
			cvResize(UV->getU(), tempVx,CV_INTER_LINEAR); 	
			UV->setU(tempVx);
			IplImage *tempVy = cvCreateImage(cvSize(width, height), UV->getU()->depth, UV->getU()->nChannels);
			cvResize(UV->getV(), tempVy,CV_INTER_LINEAR); 
			UV->setV(tempVy);
			UV->setPsidashFSAns1(cvCreateImage(cvSize( tempVx->width, tempVx->height+1 ),tempVx->depth,tempVx->nChannels));
			UV->setPsidashFSAns2(cvCreateImage(cvSize( tempVx->width+1, tempVx->height ),tempVx->depth,tempVx->nChannels)); 						
			
			WarpImage2 = RGBwarp(Pyramid2.getImageFromPyramid(k),UV->getU(),UV->getV());		
			//cvShowImage("warp:",WarpImage2);
			//cvWaitKey(0);

			cvCvtColor( WarpImage2, pyramid2_gray, CV_BGR2GRAY );
		}	
			
					  
							
		start = std::clock();

		
		SmoothFlowPDE( pyramid1_gray,pyramid2_gray,alpha,gamma,nOuterFPIterations,nInnerFPIterations,UV);
		diff = ( std::clock() - start ) / (double)CLOCKS_PER_SEC;
		std::cout<<"SmoothFlowPDE took: "<< diff <<" secs"<<endl;

	}
	return UV;
	
}



void coarse2FineCompute::SmoothFlowPDE(  IplImage* Im1, 
											IplImage* Im2, 											
											double alpha,
											double gamma,
											int nOuterFPIterations, 
											int nInnerFPIterations,
											flowUV* UV){
		
		//dimentions
		int height=Im1->height;
		int width=Im1->width;
		int channels=Im1->nChannels;		
		//init for the different DX,DY & DT		
		IplImage* Ikx=cvCreateImage(cvSize( width, height ),_imageDepth,channels); 
		IplImage* Iky=cvCreateImage(cvSize( width, height ),_imageDepth,channels);
		IplImage* Ikx2=cvCreateImage(cvSize( width, height ),_imageDepth,channels); 
		IplImage* Iky2=cvCreateImage(cvSize( width, height ),_imageDepth,channels); 
		IplImage* Ikt_Org=cvCreateImage(cvSize( width, height ),_imageDepth,channels); //IKZ
		IplImage* IXt_axis=cvCreateImage(cvSize( width, height ),_imageDepth,channels); 
		IplImage* IYt_ayis=cvCreateImage(cvSize( width, height ),_imageDepth,channels); 
		//the gradient of the gradient
		IplImage* Ixx=cvCreateImage(cvSize( width, height ),_imageDepth,channels); 
		IplImage* Ixy=cvCreateImage(cvSize( width, height ),_imageDepth,channels);
		IplImage* Iyx=cvCreateImage(cvSize( width, height ),_imageDepth,channels); 
		IplImage* Iyy=cvCreateImage(cvSize( width, height ),_imageDepth,channels);	
		//the addition in each iter to u&v
		IplImage* Du=cvCreateImage(cvSize( width, height ),_imageDepth,channels); 
		IplImage* Dv=cvCreateImage(cvSize( width, height ),_imageDepth,channels);
			
		//clear all derivatives
		cvZero(Ikx); cvZero(Iky); cvZero(Ikt_Org); cvZero(Ixx); cvZero(Ixy); 
		cvZero(Iyy); cvZero(IXt_axis); cvZero(IYt_ayis);cvZero(Du);cvZero(Dv);
		//create the different DX of the pictures
	
		getDXsCV(Im1,Ikx,Iky);	
		getDXsCV(Im2,Ikx2,Iky2);		
		//by brox we need to take the gradient of the gradient:
		getDXsCV(Ikx,Ixx,Ixy);
		getDXsCV(Iky,Iyx,Iyy);
		
		//DXT of original images and their x&y gradiants	 			
		cvSub(Im1,Im2,Ikt_Org);
		cvSub(Ikx2,Ikx,IXt_axis);
		cvSub(Iky2,Iky,IYt_ayis);
		//////////////////////////
		//cvZero(Ikx); cvZero(Iky); cvZero(Ikt_Org); cvZero(Ixx); cvZero(Ixy); cvZero(Iyy); cvZero(IXt_axis); cvZero(IYt_ayis);
		////Ikx, Iky, Ikt_Org, Ixx, Ixy, Iyy, IXt_axis, IYt_ayis
		//
		//Ikx=toolsKit::IplFromFile("c:\\a\\Ix.txt");
		//Iky=toolsKit::IplFromFile("c:\\a\\Iy.txt");
		//Ikt_Org=toolsKit::IplFromFile("c:\\a\\Iz.txt");
		//Ixx=toolsKit::IplFromFile("c:\\a\\Ixx.txt");
		//Ixy=toolsKit::IplFromFile("c:\\a\\Ixy.txt");
		//Iyy=toolsKit::IplFromFile("c:\\a\\Iyy.txt");
		//IXt_axis=toolsKit::IplFromFile("c:\\a\\Ixz.txt");
		//IYt_ayis=toolsKit::IplFromFile("c:\\a\\Iyz.txt");
		//////////////////////////
	
		
		//outer fixed point iteration
		vector<float> * dUdV = new vector<float>(Ikx->height*Ikx->width*2);
		SparseMat<float> * A = new SparseMat<float>(2*Ikx->height*Ikx->width);
		vector<float> * B = new vector<float>(A->getM());
		for(int iter=0;iter<nOuterFPIterations;iter++){						
			///construct Matrix and solve it
			dUdV = constructMatrix_brox::constructMatrix_b(Ikx, Iky, Ikt_Org, Ixx, Ixy, Iyy, IXt_axis, IYt_ayis, 
																		   UV,Du,Dv,A,B,dUdV, gamma ,alpha, _ERROR_CONST,nInnerFPIterations);						
			//cout<<"dUdV size = "<<dUdV->size()<<endl;
			//cout<<"Du size is: "<<Du->height<<","<<Du->width<<endl;
			//cout<<"Dv size is: "<<Dv->height<<","<<Dv->width<<endl;
			
			toolsKit::seperateDuDv(Du,Dv,dUdV);				

			//delete dUdV;			
			//erase edges as in matlab				
			toolsKit::cvNormalizeEdges(Du);
			toolsKit::cvNormalizeEdges(Dv);
		
			//print flow
			toolsKit::drawFlow2(UV->getU(),Du,UV->getV(),Dv,1);
		}

		//now add the most accurate du/dv to the flow		
		cvAdd(UV->getU(),Du,UV->getU());
		cvAdd(UV->getV(),Dv,UV->getV());		



	//clean temp vars
		delete A;
		delete dUdV;
	cvReleaseImage( &Ikx ); 
	cvReleaseImage( &Iky ); 
	cvReleaseImage( &Ikx2 ); 
	cvReleaseImage( &Iky2 ); 
	cvReleaseImage( &Ikt_Org ); 
	cvReleaseImage( &IXt_axis ); 
	cvReleaseImage( &IYt_ayis ); 
	cvReleaseImage( &Ixx ); 
	cvReleaseImage( &Ixy ); 
	cvReleaseImage( &Iyx ); 
	cvReleaseImage( &Iyy ); 
	cvReleaseImage( &Du ); 
	cvReleaseImage( &Dv ); 
	UV->releaseAns1and2();			

}

