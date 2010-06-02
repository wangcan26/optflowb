#include "constructMatrix_brox.h"



constructMatrix_brox::constructMatrix_brox(void)
{
}

constructMatrix_brox::~constructMatrix_brox(void)
{
}

//theta = 1/(x^2+y^2+epsilon);
void computeTheta(IplImage* theta,IplImage* x,IplImage* y,IplImage* epsilon){
	IplImage* tempx=cvCreateImage(cvSize( x->width, x->height ),IPL_DEPTH_32F,x->nChannels);
	IplImage* tempy=cvCreateImage(cvSize( x->width, x->height ),IPL_DEPTH_32F,x->nChannels);
	tempx=cvCloneImage(x);
	tempy=cvCloneImage(y);	
	double two=2.0;
	//x^2	
	toolsKit::cvMulScalar(tempx,two);
	//y^2
	toolsKit::cvMulScalar(tempy,two);
	//theta=x^2+y^2
	toolsKit::IPL_add(tempx,tempy,theta);
	//theta=theta+epsilon
	
	toolsKit::IPL_add(theta,epsilon,theta);
	//cvShowImage("theta-before",theta);
	//IPL_print(theta);

	toolsKit::IPL_mul_inverse(theta,1);
//	IPL_print(theta);
//	cvShowImage("theta-after",theta);
//	toolsKit::cvShowManyImages("computeTheta",2,tempx,tempy);
	cvReleaseImage( &tempx ); 
	cvReleaseImage( &tempy ); 
}
//psiDerivative( theta0 * ( Ikz + Ikx * du + Iky * dv ) ^ 2 );
void computePsidashBCA(IplImage* psidashBCA,IplImage* theta0,IplImage* Ikz,IplImage* Ikx,IplImage* du,
					   IplImage* Iky,IplImage* dv,double epsilon){
	toolsKit::costumeLineCompute(psidashBCA,Ikz,Ikx,du,Iky,dv);
	cvMul(theta0,psidashBCA,psidashBCA);
	toolsKit::psiDerivative(psidashBCA,epsilon);
}

//psiDerivative( gamma * (  theta1 *  ( Ixz + Ixx * du + Ixy * dv ) ^ 2 + 
	//						theta2 *  ( Iyz + Ixy * du + Iyy * dv ) ^ 2 ) ) ;
void computepsidashGCA(IplImage* psidashGCA,int gamma,IplImage* theta1,IplImage* Ixz,IplImage* Ixx,
					   IplImage* du,IplImage* Ixy,IplImage* dv,IplImage* theta2,
					   IplImage* Iyz,IplImage* Iyy,double epsilon){
	
	IplImage* temp=cvCreateImage(cvSize( psidashGCA->width, psidashGCA->height ),psidashGCA->depth,psidashGCA->nChannels);
	//theta1 *  ( Ixz + Ixx * du + Ixy * dv ) ^ 2
	toolsKit::costumeLineCompute(psidashGCA,Ixz,Ixx,du,Ixy,dv);
	cvMul(theta1,psidashGCA,psidashGCA);
	//theta2 *  ( Iyz + Ixy * du + Iyy * dv ) ^ 2 )
	toolsKit::costumeLineCompute(temp,Iyz,Ixy,du,Iyy,dv);
	cvMul(theta2,temp,temp);
	cvAdd(psidashGCA,temp,psidashGCA);
	toolsKit::cvMulScalar(psidashGCA,gamma);
	toolsKit::psiDerivative(psidashGCA,epsilon);
}


void computeDiagonalPdfSum(IplImage* ans,IplImage* psidashBCA,IplImage* theta0,IplImage* Ikx,double gamma,
						   IplImage* psidashGCA,IplImage* theta1,IplImage* Ixx,
						   IplImage* theta2,IplImage* Ixy,IplImage* pdfsum){
}
void computeDiagonalReg(IplImage* ans,IplImage* psidashBCA,IplImage* theta0,IplImage* Ikx,IplImage* Iky,double gamma,
						IplImage* psidashGCA,IplImage* theta1,IplImage* Ixx,IplImage* Ixy,
						IplImage* theta2,IplImage* Iyy){



}
void constructMatrix_brox::constructMatrix_b(IplImage* Ikx,
											 IplImage* Iky,
											 IplImage* Ikz,
											 IplImage* Ixx,
											 IplImage* Ixy,
											 IplImage* Iyy,
											 IplImage* Ixz,
											 IplImage* Iyz,
											 IplImage* psidash,
											 IplImage* psidashFS1,
											 IplImage* psidashFS2,
											 IplImage* u,
											 IplImage* v,
											 IplImage* du,
											 IplImage* dv,
											 double gamma,
											 double alpha,
											 double _ERROR_CONST ){
	
	//init IPLs											
	IplImage* theta0=cvCreateImage(cvSize(Ikx->width, Ikz->height ),IPL_DEPTH_32F,Ikz->nChannels);
	IplImage* theta1=cvCreateImage(cvSize(Ikx->width, Ikz->height ),IPL_DEPTH_32F,Ikz->nChannels);
	IplImage* theta2=cvCreateImage(cvSize(Ikx->width, Ikz->height ),IPL_DEPTH_32F,Ikz->nChannels);
	IplImage* psidashBCA=cvCreateImage(cvSize(Ikx->width, Ikz->height ),IPL_DEPTH_32F,Ikz->nChannels);
	IplImage* psidashGCA=cvCreateImage(cvSize(Ikx->width, Ikz->height ),IPL_DEPTH_32F,Ikz->nChannels);
	IplImage* epsilon=cvCreateImage(cvSize(Ikx->width, Ikz->height ),Ikz->depth,Ikz->nChannels);
	IplImage* uapp=cvCreateImage(cvSize(Ikx->width, Ikz->height ),Ikz->depth,Ikz->nChannels);
	IplImage* vapp=cvCreateImage(cvSize(Ikx->width, Ikz->height ),Ikz->depth,Ikz->nChannels);
	IplImage* uvapp=cvCreateImage(cvSize(Ikx->width, Ikz->height ),Ikz->depth,Ikz->nChannels);
	IplImage* vuapp=cvCreateImage(cvSize(Ikx->width, Ikz->height ),Ikz->depth,Ikz->nChannels);
	IplImage* pdfSum=cvCreateImage(cvSize(Ikx->width, Ikz->height ),Ikz->depth,Ikz->nChannels);
	//epsilon = 1e-3*ones(size(Ikx))==>zeroing and adding instead
	cvZero(epsilon);
	cvAddS(epsilon,cvScalarAll(_ERROR_CONST),epsilon);

	int width=u->width;
	int height=u->height;

	//theta0 = 1/(Ikx^2+Iky^2+epsilon);
	computeTheta(theta0,Ikx,Iky,epsilon);
	//theta1 = 1/(Ixx^2+Ixy^2+epsilon);
	computeTheta(theta1,Ixx,Ixy,epsilon);
	//theta2 = 1/(Iyy^2+Ixy^2+epsilon);
	computeTheta(theta2,Iyy,Ixy,epsilon);

	// First compute the values of the data  term
    //the brightness constancy assumption


	
	//the brightness constancy assumption
	computePsidashBCA(psidashBCA,theta0,Ikz,Ikx,du,Iky,dv,_ERROR_CONST);
	
    //and the Gradient Constancy Assumption
	computepsidashGCA(psidashGCA,gamma,theta1,Ixz,Ixx,du,Ixy,dv,theta2,Iyz,Iyy,_ERROR_CONST);

	//now compute the  smoothness term
	//Compute new psidashFS(it was computed erlier just mul by alpha here)
	toolsKit::cvMulScalar(psidashFS1,alpha);
	toolsKit::cvMulScalar(psidashFS2,alpha);
	
	//compute pdfSum

	/*uapp= psidashBCA * 
	   (theta0 * ( Ikx .^ 2 ))+ 
	   gamma * psidashGCA * 
	   ( theta1 *  Ixx ^ 2 +
	     theta2 * Ixy ^ 2 )  + 
	   pdfsum */
	computeDiagonalPdfSum(uapp, psidashBCA,theta0,Ikx,gamma,psidashGCA,theta1,Ixx,theta2,Ixy,pdfSum);
	computeDiagonalPdfSum(vapp, psidashBCA,theta0,Iky,gamma,psidashGCA,theta2,Iyy,theta1,Ixy,pdfSum) ;
	computeDiagonalReg   (uvapp,psidashBCA,theta0,Ikx,Iky,gamma,psidashGCA,theta1,Ixx,Ixy,theta2,Iyy);
	computeDiagonalReg   (vuapp,psidashBCA,theta0,Ikx,Iky,gamma,psidashGCA,theta1,Ixx,Ixy,theta2,Iyy);
	/*
	uapp  = psidashBCA * theta0 * ( Ikx ^ 2) +   gamma * psidashGCA * (theta1 *  Ixx ^ 2 +  theta2 * Ixy ^ 2 )  + pdfsum ;
	vapp  = psidashBCA * theta0 * ( Iky ^ 2) +   gamma * psidashGCA * (theta2 *  Iyy ^ 2 +  theta1 * Ixy ^ 2 )  + pdfsum ;
	uvapp = psidashBCA * theta0 * ( Ikx * Iky) + gamma * psidashGCA * (theta1 * Ixx * Ixy + theta2 * Iyy * Ixy ) ;
	vuapp = psidashBCA * theta0 * ( Ikx * Iky) + gamma * psidashGCA * (theta1 * Ixx * Ixy + theta2 * Iyy * Ixy  ) ;
	*/
	toolsKit::cvShowManyImages("constructMatrix_b:psidashBCA,psidashGCA,theta0,theta1,theta2",5,psidashBCA,psidashGCA,theta0,theta1,theta2);
}
