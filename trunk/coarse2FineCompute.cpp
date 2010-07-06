#include "coarse2FineCompute.h"

#include <ctime>
#include "optical_flow_demo.h"
coarse2FineCompute::coarse2FineCompute(int imageDepth,double error)
{
	_imageDepth=imageDepth;
	_ERROR_CONST=error;
}

coarse2FineCompute::~coarse2FineCompute(void)
{
}



	//template <class T>
	int EnforceRange(const int& x,const int& MaxValue) {return __min(__max(x,0),MaxValue-1);};

//--------------------------------------------------------------------------------------------------
// function to interplate multi-channel image plane for (x,y)
// --------------------------------------------------------------------------------------------------
 void BilinearInterpolate(const IplImage* pImage,int width,int height,int nChannels,double x,double y,IplImage* result,int i,int j)
{
	int xx,yy,u,v,k,offset;
	xx=x;
	yy=y;
	double dx,dy,s;
	dx=__max(__min(x-xx,1),0);
	dy=__max(__min(y-yy,1),0);
	//memset(result,0,sizeof(IplImage*)*nChannels);
	
	for(int m=0;m<=1;m++)
		for(int n=0;n<=1;n++)
		{
			u=EnforceRange(xx+m,width);
			v=EnforceRange(yy+n,height);
			offset=(v*width+u)*nChannels;
			s=fabs(1-m-dx)*fabs(1-n-dy);
			for(k=0;k<nChannels;k++)				
				//result[k]+=pImage[offset+k]*s;
				result->imageData[i*width+j+k]+=(float)fabs(((float)pImage->imageData[offset+k]*(float)s));
		}
		//cout<<"bilinier(i,wid,j):"<<i<<","<<width<<","<<j<<":"<<(float)result->imageData[i*width+j]<<endl;
}
//------------------------------------------------------------------------------------------------------------
// function to warp an image with respect to flow field
// pWarpIm2 has to be allocated before hands
//------------------------------------------------------------------------------------------------------------
//template <class T1,class T2>
int warpImage(IplImage* pWarpIm2,const IplImage* pIm1, const IplImage* pIm2, const IplImage* pVx, const IplImage* pVy)
{
	int ans,i,j,k,t;
	int width,height,nChannels;
	width=pIm2->width;
	height=pIm2->height;
	nChannels=pIm2->nChannels;
	
	for(int i=0;i<height;i++){
		for(int j=0;j<width;j++){

			
			int offset=i*width+j;
			double x,y;
		//	y=i+pVy[offset];		
			y=((float)pVy->imageData[offset])+(float)i;			
		//	x=j+pVx[offset];
			x=((float)pVx->imageData[offset])+(float)j;
			offset*=nChannels;
			//edges only
			
			BilinearInterpolate(pIm2,width,height,nChannels,x,y,pWarpIm2,i,j);					
						
			/*if(j+2>width){
				cout<<"i:"<<i<<" j:"<<j<<" i*j:"<<i*j<<endl;
			}*/
		}
		
	}
		/*	k=pIm1->width-1;
			i=0;
			j=0;
			for (t=0; t < pIm1->width*pIm1->height; t++)					
				{					  
					 //((float*)pWarpIm2->imageData)[i]=((float*)pIm1->imageData)[i]+((float*)pIm1->imageData)[i];	
					 if(k==0){
						k=width-1;
						j=0;
						i=0;
						//cout<<"k:"<<k<<" t:"<<t/width<<endl;
					 }
					 else{	 
						k--;
						j++;
						i++;
						int offset=i*width+j;
						double x,y;					
						y=((float)pVy->imageData[offset])+(float)i;								
						x=((float)pVx->imageData[offset])+(float)j;
						offset*=nChannels;
						BilinearInterpolate(pIm2,width,height,nChannels,x,y,pWarpIm2,i,j);
					 }
				}*/

		return 0;
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

void vectorMin(vector<float>* v, float val){
	for (vector<float>::iterator it = v->begin(); it != v->end(); it++)
		*it = (*it<val?*it:val);
	}

void vectorMax(vector<float>* v, float val){
		for (vector<float>::iterator it = v->begin(); it != v->end(); it++)
			*it = (*it>val?*it:val);
	}

vector<float> * vectorFloor(vector<float> * v){
		vector<float> * ans = new vector<float>();
		for (vector<float>::iterator it = v->begin(); it != v->end(); it++){
			ans->push_back((float)(((int)*it)/1));
			}
		return ans;
	}

vector<float> * vectorCeil(vector<float> * v){
		vector<float> * ans = new vector<float>();
		for (vector<float>::iterator it = v->begin(); it != v->end(); it++){
			ans->push_back(ceil(*it));
			}
		return ans;
	}

vector<float> * vectorSub(vector<float> * a, vector<float>* b){
		vector<float>* ans = new vector<float>();
		vector<float>::iterator ait = a->begin();
		vector<float>::iterator bit = b->begin();
		
		while(ait != a->end() && bit != b->end()){
			ans->push_back(*ait - *bit);
			ait++; bit++;
			}
		return ans;
	}

vector<float>* vectorMul(vector<float>* a, vector<float>* b){
	vector<float>* ans = new vector<float>();
	vector<float>::iterator ait = a->begin();
		vector<float>::iterator bit = b->begin();
		
		while(ait != a->end() && bit != b->end()){
			ans->push_back((*ait) * (*bit));
			ait++; bit++;
			}
	return ans;
	}

vector<float>* vectorAdd(vector<float>* a, vector<float>* b){
	vector<float>* ans = new vector<float>();
	vector<float>::iterator ait = a->begin();
		vector<float>::iterator bit = b->begin();
		
		while(ait != a->end() && bit != b->end()){
			ans->push_back((*ait) + (*bit));
			ait++; bit++;
			}
	return ans;
	}

IplImage** coarse2FineCompute::RGBwarp(IplImage* img, IplImage* u, IplImage* v){

	int height = img->height;
	int width = img->width;
	int nChannels = img->nChannels;
	IplImage ** XY = meshgrid(width, height);
	IplImage * X = XY[0];
	IplImage * Y = XY[1];
	
	IplImage * Xu = cvCreateImage(cvSize(X->width,X->height),X->depth,X->nChannels);
	cvAdd(X,u,Xu);
	cvReleaseImage(&Xu);

	IplImage * Yv = cvCreateImage(cvSize(Y->width,Y->height),Y->depth,Y->nChannels);
	cvAdd(Y,v,Yv);
	cvReleaseImage(&Yv);

	vector<float>* XI = toolsKit::IplImageToCoulmnVector(Xu);
	vector<float>* YI = toolsKit::IplImageToCoulmnVector(Yv);
	float eM6 = 0.00247875218; //1E-6
	//XI = max(1, min(sx - 1E-6, XI));
	vectorMin(XI, height-eM6);
	vectorMax(XI, 1);

	//XI = max(1, min(sx - 1E-6, XI));
	vectorMin(YI, width-eM6);
	vectorMax(YI, 1);
 
	//fXI = floor(XI);
	vector<float>* fXI = vectorFloor(XI);
	//cXI = ceil(XI);
	vector<float>* cXI = vectorCeil(XI);
	//fYI = floor(YI);
	vector<float>* fYI = vectorFloor(YI);
	//cYI = ceil(YI);
	vector<float>* cYI = vectorCeil(YI);

	//alpha_x = XI - fXI;
	vector<float>* alpha_x = vecotrSub(XI, fXI);
	//alpha_y = YI - fYI;
	vector<float>* alpha_y = vecotrSub(YI, fYI);


	/*note:
		A([i1 i2 i3]) is a vector of matrix A elements (as a column vecotr) in index i1 i2 i3...
		example:
			A = [1 2 3; 4 5 6; 7 8 9]
			A(:) = 1 4 7 2 5 8 3 6 9
			A([1 3 5]) = 1 7 5
			A([3 3 3]) = 7 7 7*/

	//O = (1 - alpha_x) .* (1 - alpha_y) .* I(fYI + sy * (fXI - 1)) + ...
	//	  alpha_x .* (1 - alpha_y) .* I(fYI + sy * (cXI - 1)) + ...
	//    (1 - alpha_x) .* alpha_y .* I(cYI + sy * (fXI - 1)) + ...
	//    alpha_x .* alpha_y .* I(cYI + sy * (cXI - 1));



	//O = reshape(O, sy, sx); -> from vector to IPLImage

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
	return NULL;
	}



//need to recieve cvSobel with function pointer;
int getDXsCVSobel(const IplImage* src1,IplImage* dest_dx,IplImage* dest_dy){	
	
	double x[7] =	{0.016667,
						-0.15,
						 0.75,
							0,
						-0.75,
						 0.15,
					-0.016667}	;
						

	double y[7] =  {0.016667,-0.15, 0.75,0,-0.75,0.15,-0.016667};

    CvPoint point = cvPoint(-1,-1);
	//x derivative
	CvMat* weickert = &cvMat(1, 7, CV_64FC1, x ); // 64FC1 for double
	cvFilter2D(src1,dest_dx,weickert,point);
	
	//y derivative
	weickert = &cvMat( 7, 1, CV_64FC1, y );
	cvFilter2D(src1,dest_dy,weickert);
	//z derivative
	//matScharr_3x3x2 = &cvMat( 3, 3, CV_64FC1, z );
	//cvFilter2D(src1,dest_dz,matScharr_3x3x2);
	
	//cvSobel(src, dest_dx, 1, 0, 1);
	//cvSobel( src, dest_dy, 0, 1, 1);
	//fix to fit matlab
	toolsKit::cvMulScalar(dest_dx,-1);
	toolsKit::cvMulScalar(dest_dy,-1);
	//toolsKit::cvMulScalar(dest_dz,-2);

	return 0;
}





IplImage* coarse2FineCompute::createWarp(IplImage* WarpImage2, IplImage* img1,IplImage* img2,IplImage* vx,IplImage* vy){

	//IplImage* WarpImage2 = cvCreateImage(cvSize(img2->width,img2->height ),img2->depth,img2->nChannels );
	int status=warpImage(WarpImage2,img1, img2,vx,vy);
	if (status==-1)
		cout<<"warp corrently support only 1 channel pics"<<endl;
	return WarpImage2;
}

void coarse2FineCompute::Coarse2FineFlow(IplImage* vx, 
										 IplImage* vy, 
										 const IplImage &Im1, 
										 const IplImage &Im2, 
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
	cout<<"Constructing pyramid..."<<endl;
	Pyramid1.ConstructPyramid(Im1,ratio,minWidth);


	Pyramid2.ConstructPyramid(Im2,ratio,minWidth);
	cout<<"done!"<<endl;

	// now iterate from the top level to the bottom
	IplImage* WarpImage2=NULL;
	
	std::clock_t start;
	double diff;

	
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
			cvZero(vx);
			cvZero(vy);
			//clone image2 to warpImage2 (in first iter)
			WarpImage2 = cvCreateImage(cvSize(Pyramid2.getImageFromPyramid(k)->width,Pyramid2.getImageFromPyramid(k)->height ),Pyramid2.getImageFromPyramid(k)->depth, Pyramid2.getImageFromPyramid(k)->nChannels );
			WarpImage2=  cvCloneImage(Pyramid2.getImageFromPyramid(k));			
		}
		else
		{
			IplImage *tempVx = cvCreateImage(cvSize(width, height), vx->depth, vx->nChannels);
			cvResize(vx, tempVx); 
			vx=tempVx;			
			IplImage *tempVy = cvCreateImage(cvSize(width, height), vy->depth, vy->nChannels);
			cvResize(vy, tempVy); 
			vy=tempVy;			
			//create the warp image
			WarpImage2 = cvCreateImage(cvSize(Pyramid2.getImageFromPyramid(k)->width,Pyramid2.getImageFromPyramid(k)->height ),Pyramid2.getImageFromPyramid(k)->depth, Pyramid2.getImageFromPyramid(k)->nChannels );
			cvZero(WarpImage2);
			//WarpImage2=createWarp(WarpImage2,Pyramid1.getImageFromPyramid(k),Pyramid2.getImageFromPyramid(k),vx,vy);
			RGBwarp(WarpImage2,vx,vy);
			
					  
		}
		//cvNormalize(img1_32,img1_32,0,255,CV_MINMAX); 
		IplImage *temp1 = cvCreateImage(cvSize(width, height), WarpImage2->depth, WarpImage2->nChannels);
		IplImage *temp2 = cvCreateImage(cvSize(width, height), WarpImage2->depth, WarpImage2->nChannels);
		IplImage *temp3 = cvCreateImage(cvSize(width, height), WarpImage2->depth, WarpImage2->nChannels);
		
		//cvNormalize(WarpImage2,temp1,0,1,CV_MINMAX);
		//cvNormalize(Pyramid1.getImageFromPyramid(k),temp2,0,1,CV_MINMAX);
		//cvNormalize(Pyramid2.getImageFromPyramid(k),temp3,0,1,CV_MINMAX);

		//toolsKit::cvShowManyImages("warpImage2,image1,image2",3, temp1,temp2,temp3);
		
		start = std::clock();
		/*cout<<"img1:"<<endl;
		//toolsKit::IPL_print( Pyramid1.getImageFromPyramid(k));
		cout<<"img2:"<<endl;
		//toolsKit::IPL_print( Pyramid2.getImageFromPyramid(k));
*/

		SmoothFlowPDE( Pyramid1.getImageFromPyramid(k),Pyramid2.getImageFromPyramid(k),WarpImage2,vx,vy,alpha,gamma,nOuterFPIterations,nInnerFPIterations);//,nCGIterations);	
		diff = ( std::clock() - start ) / (double)CLOCKS_PER_SEC;
		std::cout<<"printf: "<< diff <<'\n';

	}
	
}

void coarse2FineCompute::computePsidashFS_brox(IplImage* iterU,IplImage* iterV,int width,int height,int channels,flowUV* UV){	
	//init masks
	double a[] = {1,1};
	CvMat* matOnes = &cvMat( 1, 2, CV_64FC1, a ); // 64FC1 for double
	CvMat* matOnesT= &cvMat( 2, 1, CV_64FC1, a );
	cvTranspose(matOnes,matOnesT);
	
	double c[] = {0.5,0.5};
	CvMat* matHalf = &cvMat( 1, 2, CV_64FC1, c ); // 64FC1 for double
	CvMat* matHalfT= &cvMat( 2, 1, CV_64FC1, c );
	cvTranspose(matHalf,matHalfT);
	
	double b[] = {1,-1};
	CvMat* matOneNegOne = &cvMat( 1, 2, CV_64FC1, b ); // 64FC1 for double
	CvMat* matOneNegOneT= &cvMat( 2, 1, CV_64FC1, b );;
	cvTranspose(matOneNegOne,matOneNegOneT);
	//init temp params
	IplImage* ux=cvCreateImage(cvSize( width, height ),_imageDepth,channels);
	IplImage* uy=cvCreateImage(cvSize( width, height ),_imageDepth,channels);
	IplImage* vx=cvCreateImage(cvSize( width, height ),_imageDepth,channels);
	IplImage* vy=cvCreateImage(cvSize( width, height ),_imageDepth,channels);
	IplImage* uxd=cvCreateImage(cvSize( width+1, height ),_imageDepth,channels);
	IplImage* vxd=cvCreateImage(cvSize( width+1, height ),_imageDepth,channels);
	IplImage* uyd=cvCreateImage(cvSize( width, height+1 ),_imageDepth,channels);
	IplImage* vyd=cvCreateImage(cvSize( width, height+1 ),_imageDepth,channels);
	IplImage* t  =cvCreateImage(cvSize( width, height+1 ),_imageDepth,channels);
	IplImage* t2  =cvCreateImage(cvSize( width+1, height ),_imageDepth,channels);
	IplImage* uxpd=cvCreateImage(cvSize( width, height+1 ),_imageDepth,channels);
	IplImage* uypd=cvCreateImage(cvSize( width+1, height ),_imageDepth,channels);
	IplImage* vxpd=cvCreateImage(cvSize( width, height+1 ),_imageDepth,channels);
	IplImage* vypd=cvCreateImage(cvSize( width+1, height ),_imageDepth,channels);

	IplImage* temp=cvCreateImage(cvSize( width, height+1 ),_imageDepth,channels);
	IplImage* temp2=cvCreateImage(cvSize( width+1, height ),_imageDepth,channels);;


	//compute psidashFS

	cvFilter2D(iterU,ux,matOneNegOne);// x and y derivatives of u by 2d convolution
	toolsKit::cvMulScalar(ux,-1);	
	cvFilter2D(iterU,uy,matOneNegOneT);// x and y derivatives of u by 2d convolution
	toolsKit::cvMulScalar(uy,-1);
	cvFilter2D(iterV,vx,matOneNegOne);// x and y derivatives of v by 2d convolution
	toolsKit::cvMulScalar(vx,-1);
	cvFilter2D(iterV,vy,matOneNegOneT);	// x and y derivatives of v by 2d convolution
	toolsKit::cvMulScalar(vy,-1);
	/////////////////////////////////////////////////////	
	toolsKit::increaseImageSize(ux,temp2,0);
	cvFilter2D(temp2,uxd,matHalf);
	
	toolsKit::increaseImageSize(vx,temp2,0);
	cvFilter2D(temp2,vxd,matHalf);	
	
	toolsKit::increaseImageSize(uy,temp,1);
	cvFilter2D(temp,uyd,matHalfT);
	
	toolsKit::increaseImageSize(vy,temp,1);
	cvFilter2D(temp,vyd,matHalfT);	
	
	//uxpd============================================================================================
	cvFilter2D(uyd,t,matHalf);// Computes the delta u(i+1/2, j) and delta u(i-1/2, j).
	cvPow(ux,ux,2);//ux^2
	cvPow(t,t,2);//t^2
	toolsKit::increaseImageSize(ux,temp,1);
	toolsKit::IPL_add_bottom(temp,t,uxpd);//uxpd = ux^2 + t^2  last line on uxpd shoud be deleted
	
	//uypd============================================================================================
	cvFilter2D(uxd,t2,matHalfT);//Computes the delta u(i, j+1/2) and delta u(i, j-1/2).
	cvPow(uy,uy,2);//uy^2
	cvPow(t2,t2,2);//t2^2
	toolsKit::increaseImageSize(uy,temp2,0);
	//cvAdd(uy,t2,uypd);//uypd = uy^2 + t2^2	
	toolsKit::IPL_add_right(temp2,t2,uypd);//uypd = uy^2 + t2^2 last line on uxpd shoud be deleted

	//vxpd============================================================================================
	cvFilter2D(vyd,t,matHalf);// Computes the delta v(i+1/2, j) and delta v(i-1/2, j).
	cvPow(vx,vx,2);//vx^2
	cvPow(t,t,2);//t^2
	toolsKit::increaseImageSize(vx,temp,1);
	toolsKit::IPL_add_bottom(temp,t,vxpd);//vxpd = vx^2 + t^2 last line on uxpd shoud be deleted
	
	//vypd============================================================================================
	cvFilter2D(vxd,t2,matHalfT);// Computes the delta v(i+1/2, j) and delta v(i-1/2, j).
	cvPow(vy,vy,2);//vy^2
	cvPow(t2,t2,2);//t2^2
	toolsKit::increaseImageSize(vy,temp2,0);
	toolsKit::IPL_add_right(temp2,t2,vypd);//vypd=vy^2 + t2^2 last line on uxpd shoud be deleted

	//cout<<"uxpd(add bottom)"<<endl;
	//toolsKit::IPL_print(uxpd);
	//cout<<"uypd(add right)"<<endl;
	//toolsKit::IPL_print(uypd);
	//cout<<"vxpd(add bottom)"<<endl;
	//toolsKit::IPL_print(vxpd);
	//cout<<"vypd(add right)"<<endl;
	//toolsKit::IPL_print(vypd);
	
	toolsKit::IPL_add_different_sizes2(uypd,vypd,UV->getPsidashFSAns1());
	toolsKit::IPL_add_different_sizes3(uxpd,vxpd,UV->getPsidashFSAns2());

	//toolsKit::cvShowManyImages("after:uypd,vypd,ans1",3,uypd,vypd,UV->getPsidashFSAns1());			
	//toolsKit::cvShowManyImages("after:vxpd,vxpd,ans2",3,vxpd,vypd,UV->getPsidashFSAns2());

	toolsKit::psiDerivative(UV->getPsidashFSAns1(),_ERROR_CONST);
	toolsKit::psiDerivative(UV->getPsidashFSAns2(),_ERROR_CONST);
	
	//toolsKit::cvZeroRight(UV->getPsidashFSAns1());
	//toolsKit::cvZeroBottomLeft(UV->getPsidashFSAns2());

	cvReleaseImage( &ux ); 
	cvReleaseImage( &uy ); 
	cvReleaseImage( &vx ); 
	cvReleaseImage( &vy ); 
	cvReleaseImage( &uxd ); 
	cvReleaseImage( &vxd ); 
	cvReleaseImage( &uyd ); 
	cvReleaseImage( &vyd ); 
	cvReleaseImage( &t ); 
	cvReleaseImage( &t2 ); 
	cvReleaseImage( &uxpd ); 
	cvReleaseImage( &uypd ); 
	cvReleaseImage( &vxpd ); 
	cvReleaseImage( &vypd ); 
	cvReleaseImage( &temp ); 
	cvReleaseImage( &temp2 ); 

//	return ans;
}

flowUV* coarse2FineCompute::SmoothFlowPDE(  IplImage* Im1, 
											IplImage* Im2, 
											IplImage* warpIm2, 
											IplImage* uinit, 
											IplImage* vinit, 
											double alpha,
											double gamma,
											int nOuterFPIterations, 
											int nInnerFPIterations){
		
		//dimentions
		int height=Im1->height;
		int width=Im1->width;
		int channels=Im1->nChannels;
		//this will hold the optical flow
		flowUV* UV=new flowUV(uinit,vinit);
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
	
		//create the different DX of the pictures
	
		getDXsCVSobel(Im1,Ikx,Iky);	
		getDXsCVSobel(Im2,Ikx2,Iky2);		
		//by brox we need to take the gradient of the gradient:
		getDXsCVSobel(Ikx,Ixx,Ixy);
		getDXsCVSobel(Iky,Iyx,Iyy);
		
		//DXT of original images and their x&y gradiants
	 	//cvAbsDiff(Im1,Im2,Ikt_Org);//IKz
		//cvAbsDiff(Ikx,Ikx2,IXt_axis);
		//cvAbsDiff(Iky,Iky2,IYt_ayis);
		
		cvSub(Im2,Im1,Ikt_Org);
		cvSub(Ikx2,Ikx,IXt_axis);
		cvSub(Iky2,Iky,IYt_ayis);
		//////////////////////////
		/*cvZero(Ikx); cvZero(Iky); cvZero(Ikt_Org); cvZero(Ixx); cvZero(Ixy); cvZero(Iyy); cvZero(IXt_axis); cvZero(IYt_ayis);
		//Ikx, Iky, Ikt_Org, Ixx, Ixy, Iyy, IXt_axis, IYt_ayis
		
		Ikx=toolsKit::IplFromFile("c:\\a\\Ix.txt");
		Iky=toolsKit::IplFromFile("c:\\a\\Iy.txt");
		Ikt_Org=toolsKit::IplFromFile("c:\\a\\Iz.txt");
		Ixx=toolsKit::IplFromFile("c:\\a\\Ixx.txt");
		Ixy=toolsKit::IplFromFile("c:\\a\\Ixy.txt");
		Iyy=toolsKit::IplFromFile("c:\\a\\Iyy.txt");
		IXt_axis=toolsKit::IplFromFile("c:\\a\\Ixz.txt");
		IYt_ayis=toolsKit::IplFromFile("c:\\a\\Iyz.txt");*/
		//////////////////////////
		cout<<"Ikx"<<endl;
		toolsKit::IPL_print(Ikx);
		cout<<"Iky"<<endl;
		toolsKit::IPL_print(Iky);
		
		//cout<<"Ikx2"<<endl;
		//toolsKit::IPL_print(Ikx2);
		cout<<"Iky2"<<endl;
		toolsKit::IPL_print(Iky2);	
	
	
		////cout<<"Iyx"<<endl;
		////toolsKit::IPL_print(Iyx);
		//cout<<"Iyy"<<endl;
		//toolsKit::IPL_print(Iyy);
		
		cout<<"Ixy"<<endl;
		toolsKit::IPL_print(Ixy);
		//cout<<"Ixx"<<endl;
		//toolsKit::IPL_print(Ixx);
		//cout<<"Iyy"<<endl;
		//toolsKit::IPL_print(Iyy);

		//cout<<"Ikt_Org(ikz)"<<endl;
		//toolsKit::IPL_print(Ikt_Org);
		//cout<<"IXt_axis"<<endl;
		//toolsKit::IPL_print(IXt_axis);
		//cout<<"IYt_ayis"<<endl;
		//toolsKit::IPL_print(IYt_ayis);
		
		//outer fixed point iteration
		for(int iter=0;iter<nOuterFPIterations;iter++){
			
			if (iter==0){
				iter=1;
				computePsidashFS_brox(toolsKit::IplFromFile("c:\\a\\u_double.txt"),toolsKit::IplFromFile("c:\\a\\v_double.txt"),width,height,channels,UV);
			}else
				computePsidashFS_brox(UV->getU(),UV->getV(),width,height,channels,UV);
		
			
			

			toolsKit::cvMulScalar(UV->getPsidashFSAns1(),alpha);
			toolsKit::cvMulScalar(UV->getPsidashFSAns2(),alpha);


			cout<<"fs1*alpha"<<endl;
			toolsKit::IPL_print(UV->getPsidashFSAns1());
			cout<<"fs2*alpha"<<endl;
			toolsKit::IPL_print(UV->getPsidashFSAns2());

			vector<float> * dUdV = constructMatrix_brox::constructMatrix_b(Ikx, Iky, Ikt_Org, Ixx, Ixy, Iyy, IXt_axis, IYt_ayis, UV->getPsidashFSAns1(),UV->getPsidashFSAns2(), UV->getU(), UV->getV(),Du,Dv, gamma ,alpha, _ERROR_CONST,nInnerFPIterations);
			
			
			//[duv, err, it, flag] = sor( A, duv, b, omega, inner_iter, tol ) ;
			IplImageIterator<float> DUit(Du);
			IplImageIterator<float> DVit(Dv);
			int i=0;
			cout<<"dUdV size = "<<dUdV->size()<<endl;
			cout<<"Du size is: "<<Du->height<<","<<Du->width<<endl;
			cout<<"Dv size is: "<<Dv->height<<","<<Dv->width<<endl;
			for (vector<float>::iterator it = dUdV->begin(); it!= dUdV->end(); it++, i++)
				if (i<dUdV->size()/2){
						*DUit = *it;
						DUit++;
					}
				else{
						*DVit = *it;
						DVit++;
					}
			cout<<"end;"<<endl;
			delete dUdV;
			cout<<"freed"<<endl;


			
			//cout<<"u"<<endl;
			//toolsKit::IPL_print(UV->getU());
			//cout<<"v"<<endl;
			//toolsKit::IPL_print(UV->getV());


			cvAdd(UV->getU(),Du,UV->getU());
			cvAdd(UV->getV(),Dv,UV->getV());
			IplImage* color_img = cvCreateImage( cvSize(UV->getU()->height,UV->getU()->width), IPL_DEPTH_8U, 3 );
			CvMat mathdr, *tempU = cvGetMat( UV->getU(), &mathdr ), *tempV = cvGetMat(UV->getV(),&mathdr);
				
			MotionToColor( tempU,  tempV,  color_img,  0.1f);
			
			

			toolsKit::cvShowManyImages("plot flow",1,color_img);
			
			cvWaitKey(0);
			cvReleaseImage(&color_img);

		}


	//clean temp vars
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
		
	return UV;

}

