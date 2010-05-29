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
	cvShowImage("theta-before",theta);
	//IPL_print(theta);

	toolsKit::IPL_mul_inverse(theta,1);
//	IPL_print(theta);
//	cvShowImage("theta-after",theta);
//	toolsKit::cvShowManyImages("computeTheta",2,tempx,tempy);
	cvReleaseImage( &tempx ); 
	cvReleaseImage( &tempy ); 
}
//psiDerivative( theta0 * ( Ikz + Ikx * du + Iky * dv ) ^ 2 );
void computePsidashBCA(IplImage* psidashBCA,IplImage* theta0,IplImage* Ikz,IplImage* Ikx,IplImage* du,IplImage* Iky,IplImage* dv){
	//x=Ikx * du

	//y=Iky * dv

	//z=Ikz +x+y
	
	//t=z^2

	//ans=theta0*t

}

//psiDerivative( gamma * (  theta1 *  ( Ixz + Ixx * du + Ixy * dv ) ^ 2 + 
	//						theta2 *  ( Iyz + Ixy * du + Iyy * dv ) ^ 2 ) ) ;
void computepsidashGCA(IplImage* psidashGCA,int gamma,IplImage* theta1,IplImage* Ixz,IplImage* Ixx,IplImage* du,IplImage* Ixy,
					   IplImage* dv,IplImage* theta2,IplImage* Iyz,IplImage* Iyy){

}

void constructMatrix_brox::constructMatrix_b(IplImage* Ikx,IplImage* Iky,IplImage* Ikz,IplImage* Ixx,IplImage* Ixy,IplImage* Iyy,IplImage* Ixz,
											 IplImage* Iyz,IplImage* psidash,IplImage* psidashFS1,IplImage* psidashFS2,IplImage* u,IplImage* v,double gamma,int _ERROR_CONST ){
	

	IplImage* theta0=cvCreateImage(cvSize(Ikx->width, Ikz->height ),IPL_DEPTH_32F,Ikz->nChannels);
	IplImage* theta1=cvCreateImage(cvSize(Ikx->width, Ikz->height ),IPL_DEPTH_32F,Ikz->nChannels);
	IplImage* theta2=cvCreateImage(cvSize(Ikx->width, Ikz->height ),IPL_DEPTH_32F,Ikz->nChannels);
	IplImage* psidashBCA=cvCreateImage(cvSize(Ikx->width, Ikz->height ),IPL_DEPTH_32F,Ikz->nChannels);
	IplImage* psidashGCA=cvCreateImage(cvSize(Ikx->width, Ikz->height ),IPL_DEPTH_32F,Ikz->nChannels);
	IplImage* epsilon=cvCreateImage(cvSize(Ikx->width, Ikz->height ),Ikz->depth,Ikz->nChannels);
	//epsilon = 1e-3*ones(size(Ikx))==>zeroing and adding instead
	cvZero(epsilon);
	cvAddS(epsilon,cvScalarAll(_ERROR_CONST),epsilon);

	//theta0 = 1./(Ikx.^2+Iky.^2+epsilon);
	computeTheta(theta0,Ikx,Iky,epsilon);
	//theta1 = 1./(Ixx.^2+Ixy.^2+epsilon);
	computeTheta(theta1,Ixx,Ixy,epsilon);
	//theta2 = 1./(Iyy.^2+Ixy.^2+epsilon);
	computeTheta(theta2,Iyy,Ixy,epsilon);

	// First compute the values of the data  term
    //the brightness constancy assumption

	//TODO:du or u?? dv or v??
	computePsidashBCA(psidashBCA,theta0,Ikz,Ikx,u,Iky,v);
	//psidashBCA = 
    //and the Gradient Constancy Assumption
    //psidashGCA = psiDerivative( gamma * (  theta1 .*  ( Ixz + Ixx .* du + Ixy .* dv ) .^ 2 + 
	//							theta2 .*   ( Iyz + Ixy .* du + Iyy .* dv ) .^ 2 ) ) ;
	
	//TODO:du or u?? dv or v??
	computepsidashGCA(psidashGCA,gamma,theta1,Ixz,Ixx,u,Ixy,v,theta2,Iyz,Iyy);

}
