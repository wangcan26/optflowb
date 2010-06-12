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


	double two=2.0;
	//x^2
	tempx=cvCloneImage(x);
	toolsKit::cvMulScalar(tempx,two);
	//y^2
	tempy=cvCloneImage(y);	
	toolsKit::cvMulScalar(tempy,two);
	//theta=x^2+y^2
	cvAdd(tempx,tempy,theta);
	//theta=theta+epsilon

	cvAdd(theta,epsilon,theta);
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

/*uapp= psidashBCA * 
(theta0 * ( Ikx ^ 2 ))+ ==>temp
gamma * psidashGCA * ==>temp2
(theta1 *  Ixx ^ 2 + theta2 * Ixy ^ 2 )  + ==>temp3
pdfsum */
void computeDiagonalPdfSum(IplImage* ans,IplImage* psidashBCA,IplImage* theta0,IplImage* Ikx,double gamma,
						   IplImage* psidashGCA,IplImage* theta1,IplImage* Ixx,
						   IplImage* theta2,IplImage* Ixy,IplImage* pdfsum){
							   IplImage* temp=cvCreateImage(cvSize( ans->width, ans->height ),ans->depth,ans->nChannels);
							   IplImage* temp2=cvCreateImage(cvSize( ans->width, ans->height ),ans->depth,ans->nChannels);
							   IplImage* temp3=cvCreateImage(cvSize( ans->width, ans->height ),ans->depth,ans->nChannels);
							   IplImage* temp4=cvCreateImage(cvSize( ans->width, ans->height ),ans->depth,ans->nChannels);
							   //(theta0 * ( Ikx ^ 2 ))==>temp
							   cvPow(Ikx,Ikx,2);
							   cvMul(theta0,Ikx,temp);
							   //psidashBCA * temp==>temp
							   cvMul(psidashBCA,temp,temp);
							   // gamma * psidashGCA==>temp2
							   temp2=cvCloneImage(psidashGCA);
							   toolsKit::cvMulScalar(temp2,gamma);
							   //(theta1 *  Ixx ^ 2 + theta2 * Ixy ^ 2 )==>temp3
							   cvPow(Ixx,temp3,2);
							   cvPow(Ixy,temp4,2);
							   //theta1 * Ixx ^ 2==>temp3
							   cvMul(theta1,temp3,temp3);
							   //theta2 * Ixy ^ 2==>temp4
							   cvMul(theta2,temp4,temp4);
							   //temp3+temp4==>temp3
							   cvAdd(temp3,temp4,temp3);

							   //temp2*temp3
							   cvMul(temp2,temp3,temp2);
							   //temp+temp3
							   cvAdd(temp,temp3,temp);
							   cvAdd(temp,pdfsum,ans);
							   cvReleaseImage(&temp);
							   cvReleaseImage(&temp2);
							   cvReleaseImage(&temp3);
							   cvReleaseImage(&temp4);


}

//uvapp = psidashBCA * 
//					   theta0 * 
//							    ( Ikx * Iky) + 
//												gamma * psidashGCA * 
//																	 (theta1 * Ixx * Ixy + theta2 * Iyy * Ixy ) ;
void computeDiagonalReg(IplImage* ans,IplImage* psidashBCA,IplImage* theta0,IplImage* Ikx,IplImage* Iky,double gamma,
						IplImage* psidashGCA,IplImage* theta1,IplImage* Ixx,IplImage* Ixy,
						IplImage* theta2,IplImage* Iyy,IplImage* IxyB){
							IplImage* temp1=cvCreateImage(cvSize( ans->width, ans->height ),ans->depth,ans->nChannels);
							IplImage* temp2=cvCreateImage(cvSize( ans->width, ans->height ),ans->depth,ans->nChannels);
							IplImage* temp3=cvCreateImage(cvSize( ans->width, ans->height ),ans->depth,ans->nChannels);
							IplImage* temp4=cvCreateImage(cvSize( ans->width, ans->height ),ans->depth,ans->nChannels);
							//psidashBCA * theta0 ==>temp1
							cvMul(psidashBCA,theta0,temp1);
							//( Ikx * Iky)==>temp2
							cvMul(Ikx,Iky,temp2);
							//temp1 * temp2 ==>temp1
							cvMul(temp1,temp2,temp1);
							// gamma * psidashGCA==>temp2
							temp2=cvCloneImage(psidashGCA);
							toolsKit::cvMulScalar(temp2,gamma);
							//temp3<==(theta1 * Ixx * Ixy +
							//temp4<==                     theta2 * Iyy * Ixy )
							cvMul(theta1,Ixx,temp3);
							cvMul(temp3,Ixy,temp3);
							cvMul(theta2,Iyy,temp4);
							cvMul(temp4,IxyB,temp4);
							//temp3<==(theta1 * Ixx * Ixy + theta2 * Iyy * IxyB )
							cvAdd(temp3,temp4,temp3);
							//temp2<==gamma * psidashGCA * (theta1 * Ixx * Ixy + theta2 * Iyy * IxyB ) 
							cvMul(temp2,temp3,temp2);
							//temp1+temp2
							cvAdd(temp1,temp2,ans);
							cvReleaseImage(&temp1);
							cvReleaseImage(&temp2);
							cvReleaseImage(&temp3);
							cvReleaseImage(&temp4);



}


//pdfsum = pdfs( 1 : 2 : 2 * ht, 2 : 2 : end ) + pdfs( 3 : 2 : end, 2 : 2 : end ) +
//		   pdfs( 2 : 2 : end, 1 : 2 : 2 * wt ) + pdfs( 2 : 2 : end, 3 : 2 : end ) ;
void computePdfSum(IplImage* pdfSum,IplImage* psidashFS1,IplImage* psidashFS2){
	IplImage* temp1=cvCreateImage(cvSize( psidashFS1->width, psidashFS1->height ),psidashFS1->depth,psidashFS1->nChannels);
	IplImage* temp2=cvCreateImage(cvSize( psidashFS1->width, psidashFS1->height ),psidashFS1->depth,psidashFS1->nChannels);
	toolsKit::IPL_add_left(psidashFS1,psidashFS1,temp1);
	toolsKit::IPL_add_top(psidashFS2,psidashFS2,temp2);
	cvAdd(temp1,temp2,pdfSum);
	cvReleaseImage(&temp1);
	cvReleaseImage(&temp2);
}

/*pdfaltsumu = psidashFS2(2:2:end ,  1:2:2*wt)* ( u(2:ht+1, 1:wt)  - u(2:ht+1, 2:wt+1) ) + //left ,ans2
psidashFS2(2:2:end ,  3:2:end) * ( u(2:ht+1, 3:end) - u(2:ht+1, 2:wt+1) ) + //right ,ans2
psidashFS1(1:2:2*ht,  2:2:end) * ( u(1:ht, 2:wt+1)  - u(2:ht+1, 2:wt+1) ) + //top ,ans1
psidashFS1(3:2:end ,  2:2:end) * ( u(3:end, 2:wt+1) - u(2:ht+1, 2:wt+1) )   //bottom ,ans1
*/
void computeVectBComponents(IplImage* pdfaltSumXX,IplImage* psidashFS1,IplImage* psidashFS2,IplImage* u,IplImage* v){
	//init
	IplImage* tempLeft1=cvCreateImage(cvSize( psidashFS1->width, psidashFS1->height ),psidashFS1->depth,psidashFS1->nChannels);
	IplImage* tempLeft2=cvCreateImage(cvSize( psidashFS1->width, psidashFS1->height ),psidashFS1->depth,psidashFS1->nChannels);
	IplImage* tempLeft3=cvCreateImage(cvSize( psidashFS1->width, psidashFS1->height ),psidashFS1->depth,psidashFS1->nChannels);
	IplImage* tempLeft4=cvCreateImage(cvSize( psidashFS1->width, psidashFS1->height ),psidashFS1->depth,psidashFS1->nChannels);

	//psidashFS2(2:2:end ,  1:2:2*wt)* ( u(2:ht+1, 1:wt)  - u(2:ht+1, 2:wt+1) ) + //left ,ans2
	toolsKit::IPL_sub_left(u,u,tempLeft1);
	toolsKit::IPL_mul_left(tempLeft1,psidashFS2,tempLeft1);
	//psidashFS2(2:2:end ,  3:2:end) * ( u(2:ht+1, 3:end) - u(2:ht+1, 2:wt+1) ) + //right ,ans2
	toolsKit::IPL_sub_right(u,u,tempLeft2);
	toolsKit::IPL_mul_right(tempLeft2,psidashFS2,tempLeft2);
	//psidashFS1(1:2:2*ht,  2:2:end) * ( u(1:ht, 2:wt+1)  - u(2:ht+1, 2:wt+1) ) + //top ,ans1
	toolsKit::IPL_sub_top(u,u,tempLeft3);
	toolsKit::IPL_mul_top(tempLeft3,psidashFS1,tempLeft3);
	//psidashFS1(3:2:end ,  2:2:end) * ( u(3:end, 2:wt+1) - u(2:ht+1, 2:wt+1) )   //bottom ,ans1
	toolsKit::IPL_sub_bottom(u,u,tempLeft4);
	toolsKit::IPL_mul_bottom(tempLeft4,psidashFS1,tempLeft4);

	//temp1+temp2+temp3+temp4
	cvAdd(tempLeft1,tempLeft2,pdfaltSumXX);
	cvAdd(pdfaltSumXX,tempLeft3,pdfaltSumXX);
	cvAdd(pdfaltSumXX,tempLeft4,pdfaltSumXX);
	//- 1*pdfaltsumXX - needed for use in the next computation
	toolsKit::cvMulScalar(pdfaltSumXX,-1);

}
void constructMatrix_brox::constructMatrix_b(IplImage* Ikx,IplImage* Iky,IplImage* Ikz,IplImage* Ixx,
											 IplImage* Ixy,IplImage* Iyy,IplImage* Ixz,IplImage* Iyz,
											 IplImage* psidash,IplImage* psidashFS1,IplImage* psidashFS2,
											 IplImage* u,IplImage* v,IplImage* du,IplImage* dv,double gamma,
											 double alpha, double _ERROR_CONST ){
			 //init IPLs											
			 IplImage* theta0=	 cvCreateImage(cvSize(Ikx->width, Ikz->height ),Ikz->depth,Ikz->nChannels);
			 IplImage* theta1=	 cvCreateImage(cvSize(Ikx->width, Ikz->height ),Ikz->depth,Ikz->nChannels);
			 IplImage* theta2=	 cvCreateImage(cvSize(Ikx->width, Ikz->height ),Ikz->depth,Ikz->nChannels);
			 IplImage* psidashBCA=cvCreateImage(cvSize(Ikx->width, Ikz->height ),Ikz->depth,Ikz->nChannels);
			 IplImage* psidashGCA=cvCreateImage(cvSize(Ikx->width, Ikz->height ),Ikz->depth,Ikz->nChannels);
			 IplImage* epsilon=	 cvCreateImage(cvSize(Ikx->width, Ikz->height ),Ikz->depth,Ikz->nChannels);
			 IplImage* uapp=	 cvCreateImage(cvSize(Ikx->width, Ikz->height ),Ikz->depth,Ikz->nChannels);
			 IplImage* vapp=	 cvCreateImage(cvSize(Ikx->width, Ikz->height ),Ikz->depth,Ikz->nChannels);
			 IplImage* uvapp=	 cvCreateImage(cvSize(Ikx->width, Ikz->height ),Ikz->depth,Ikz->nChannels);
			 IplImage* vuapp=	 cvCreateImage(cvSize(Ikx->width, Ikz->height ),Ikz->depth,Ikz->nChannels);
			 IplImage* pdfSum=	 cvCreateImage(cvSize(Ikx->width, Ikz->height ),Ikz->depth,Ikz->nChannels);
			 IplImage* constu=	 cvCreateImage(cvSize(Ikx->width, Ikz->height ),Ikz->depth,Ikz->nChannels);
			 IplImage* constv=	 cvCreateImage(cvSize(Ikx->width, Ikz->height ),Ikz->depth,Ikz->nChannels);
			 IplImage* pdfaltSumU=cvCreateImage(cvSize(Ikx->width, Ikz->height ),Ikz->depth,Ikz->nChannels);
			 IplImage* pdfaltSumV=cvCreateImage(cvSize(Ikx->width, Ikz->height ),Ikz->depth,Ikz->nChannels);


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
			 //pdfsum = pdfs( 1 : 2 : 2 * ht, 2 : 2 : end ) + pdfs( 3 : 2 : end, 2 : 2 : end ) +...
			 //		   pdfs( 2 : 2 : end, 1 : 2 : 2 * wt ) + pdfs( 2 : 2 : end, 3 : 2 : end ) ;
			 computePdfSum(pdfSum,psidashFS1,psidashFS2);

			 //uapp  = psidashBCA * theta0 * ( Ikx ^ 2) +   gamma * psidashGCA * (theta1 *  Ixx ^ 2 +  theta2 * Ixy ^ 2 )  + pdfsum ;
			 computeDiagonalPdfSum(uapp, psidashBCA,theta0,Ikx,gamma,psidashGCA,theta1,Ixx,theta2,Ixy,pdfSum);
			 //vapp  = psidashBCA * theta0 * ( Iky ^ 2) +   gamma * psidashGCA * (theta2 *  Iyy ^ 2 +  theta1 * Ixy ^ 2 )  + pdfsum ;
			 computeDiagonalPdfSum(vapp, psidashBCA,theta0,Iky,gamma,psidashGCA,theta2,Iyy,theta1,Ixy,pdfSum);
			 //uvapp = psidashBCA * theta0* (Ikx*Iky)+ gamma*psidashGCA*(theta1*Ixx*Ixy + theta2*Iyy*Ixy ) ;
			 computeDiagonalReg   (uvapp,psidashBCA,theta0,Ikx,Iky,gamma,psidashGCA,theta1,Ixx,Ixy,theta2,Iyy,Ixy);
			 //vuapp =   uvapp
			 vuapp=cvCloneImage(uvapp);

			 //insert to diagonals to matrix A
			 vector<float> diag;
			 IplImageIterator<float> vu(vuapp);
			 IplImageIterator<float> uv(uvapp);
			 while (!vu){
				 diag.push_back(*vu);
				 vu++;
				 }
			 while (!uv){
				 diag.push_back(*uv);
				 uv++;
				 }
			//diag is the vector for the main(0) diag
			 //create the Sparse MAtrix:
			 SparseMat<float> * A = new SparseMat<float>(2*Ikx->height, 2*Ikx->width);
			 A->addDiag(0,diag);
			 cout<<"diag size: "<<diag.size()<<endl;
			 cout<<"A size: "<<A->getM()<<", "<<A->getM()<<" diag size: "<<A->getM()*A->getN()<<endl;
			 //////////////////////build vector B//////////////////////

			 // Computing the constant terms for the first of the Euler Lagrange equations
			 computeVectBComponents(pdfaltSumU,psidashFS1,psidashFS2,u,v);
			 // Computing the constant terms for the second of the Euler Lagrange equations
			 computeVectBComponents(pdfaltSumV,psidashFS1,psidashFS2,u,v);


			 //constu = psidashBCA * theta0 * ( Ikx * Ikz ) + gamma * psidashGCA * (theta1 * Ixx * Ixz + theta2 * Ixy * Iyz ) - 1*pdfaltsumu 
			 computeDiagonalReg   (constu,psidashBCA,theta0,Ikx,Ikz,gamma,psidashGCA,theta1,Ixx,Ixz,theta2,Ixy,Iyz);
			 cvAdd(constu,pdfaltSumU,constu);

			 //constv = psidashBCA * theta0 * ( Iky * Ikz ) + gamma * psidashGCA * (theta1 * Ixy * Ixz + theta2 * Iyy * Iyz ) - 1*pdfaltsumv ;
			 computeDiagonalReg   (constv,psidashBCA,theta0,Iky,Ikz,gamma,psidashGCA,theta1,Ixy,Ixz,theta2,Iyy,Iyz);
			 cvAdd(constv,pdfaltSumV,constv);

			 //insert data to B vector:b = [-constu(:) ; -constv(:) ];




			 toolsKit::cvShowManyImages("constructMatrix_b:uapp,vapp,uvapp,vuapp,pdfaltSumU,pdfaltSumV,constu,constv",8,uapp,vapp,uvapp,vuapp,pdfaltSumU,pdfaltSumV,constu,constv);
}
