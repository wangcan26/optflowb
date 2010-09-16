#include "constructMatrix_brox.h"
#include <ctime>



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
	cvPow(tempx,tempx,two);
	//y^2
	tempy=cvCloneImage(y);	
	cvPow(tempy,tempy,two);
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
//psidashBCA=psiDerivative( theta0 * ( Ikz + Ikx * du + Iky * dv ) ^ 2 );
void computePsidashBCA(IplImage* psidashBCA,IplImage* theta0,IplImage* Ikz,IplImage* Ikx,IplImage* du,
					   IplImage* Iky,IplImage* dv,double epsilon){
						   
						   //init psidashBCA
						   cvZero(psidashBCA);
						   toolsKit::costumeLineCompute(psidashBCA,Ikz,Ikx,du,Iky,dv);
						 
						   cvMul(theta0,psidashBCA,psidashBCA);
						  

						 //  cout<<"psidashBCA-before epsilon add"<<endl;
						 //  toolsKit::IPL_print(psidashBCA);

						   toolsKit::psiDerivative(psidashBCA,epsilon);
}

//psidashGCA=psiDerivative( gamma * (  theta1 *  ( Ixz + Ixx * du + Ixy * dv ) ^ 2 + 
//						    theta2 *  ( Iyz + Ixy * du + Iyy * dv ) ^ 2 ) ) ;
void computepsidashGCA(IplImage* psidashGCA,int gamma,IplImage* theta1,IplImage* Ixz,IplImage* Ixx,
					   IplImage* du,IplImage* Ixy,IplImage* dv,IplImage* theta2,
					   IplImage* Iyz,IplImage* Iyy,double epsilon){
						
						   if (gamma!=0){
							   IplImage* temp=cvCreateImage(cvSize( psidashGCA->width, psidashGCA->height ),psidashGCA->depth,psidashGCA->nChannels);
							   //theta1 *  ( Ixz + Ixx * du + Ixy * dv ) ^ 2
							   toolsKit::costumeLineCompute(psidashGCA,Ixz,Ixx,du,Ixy,dv);
							   cvMul(theta1,psidashGCA,psidashGCA);
							   //theta2 *  ( Iyz + Ixy * du + Iyy * dv ) ^ 2 )
							   toolsKit::costumeLineCompute(temp,Iyz,Ixy,du,Iyy,dv);
							   cvMul(theta2,temp,temp);
							   cvAdd(psidashGCA,temp,psidashGCA);
							   toolsKit::cvMulScalar(psidashGCA,gamma);
						   }
						   else{
							cvZero(psidashGCA);
						   }
						   
						  // cout<<"psidashGCA-before epsilon add"<<endl;
						  // toolsKit::IPL_print(psidashGCA);

						   toolsKit::psiDerivative(psidashGCA,epsilon);
}

/*
if gamma=0:
			uapp= psidashBCA * (theta0 * ( Ikx ^ 2 ))+pdfsum ==>ans
else:
			uapp= psidashBCA * (theta0 * ( Ikx ^ 2 ))+ ==>temp
			gamma * psidashGCA * ==>temp2
								(theta1 *  Ixx ^ 2 + theta2 * Ixy ^ 2 )  + ==>temp3
																			pdfsum */
void computeDiagonalPdfSum(IplImage* ans,IplImage* psidashBCA,IplImage* theta0,IplImage* Ikx,double gamma,
						   IplImage* psidashGCA,IplImage* theta1,IplImage* Ixx,
						   IplImage* theta2,IplImage* Ixy,IplImage* pdfsum){
							  // IplImage* temp=cvCreateImage(cvSize( ans->width, ans->height ),ans->depth,ans->nChannels);
							   IplImage* temp2=cvCreateImage(cvSize( ans->width, ans->height ),ans->depth,ans->nChannels);							  
							   //( Ikx ^ 2)==>temp2
							   temp2=cvCloneImage(Ikx);
							   cvPow(Ikx,temp2,2);
							   //(theta0 * ( Ikx ^ 2 ))==>temp
							   cvMul(theta0,temp2,ans);							   							
							   //psidashBCA * temp==>temp
							   cvMul(psidashBCA,ans,ans);
								
							   if (gamma==0){
								    cvAdd(ans,pdfsum,ans);
									cvReleaseImage(&temp2);
									return;
							   }
							   //for any other gamma:
							   IplImage* temp3=cvCreateImage(cvSize( ans->width, ans->height ),ans->depth,ans->nChannels);
							   IplImage* temp4=cvCreateImage(cvSize( ans->width, ans->height ),ans->depth,ans->nChannels);
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
							   cvAdd(ans,temp3,ans);
							   cvAdd(ans,pdfsum,ans);
							  // cvReleaseImage(&temp);
							   cvReleaseImage(&temp2);
							   cvReleaseImage(&temp3);
							   cvReleaseImage(&temp4);


}
/*if gamma=0
 uvapp = psidashBCA * 
 					   theta0 * 
							    ( Ikx * Iky)												
else
 uvapp = psidashBCA * 
 					   theta0 * 
							    ( Ikx * Iky) + 
												gamma * psidashGCA * 
																	 (theta1 * Ixx * Ixy + theta2 * Iyy * Ixy ) ;
*/

void computeDiagonalReg(IplImage* ans,IplImage* psidashBCA,IplImage* theta0,IplImage* Ikx,IplImage* Iky,double gamma,
						IplImage* psidashGCA,IplImage* theta1,IplImage* Ixx,IplImage* Ixy,
						IplImage* theta2,IplImage* Iyy,IplImage* IxyB)
{
							IplImage* temp1=cvCreateImage(cvSize( ans->width, ans->height ),ans->depth,ans->nChannels);
							IplImage* temp2=cvCreateImage(cvSize( ans->width, ans->height ),ans->depth,ans->nChannels);							
							//psidashBCA * theta0 ==>temp1
							cvMul(psidashBCA,theta0,temp1);
							//( Ikx * Iky)==>temp2
							cvMul(Ikx,Iky,temp2);
							//temp1 * temp2 ==>temp1
							cvMul(temp1,temp2,ans);
							if (gamma==0){
								cvReleaseImage(&temp1);
								cvReleaseImage(&temp2);
								return;
							}
							
							//init more temps
							IplImage* temp3=cvCreateImage(cvSize( ans->width, ans->height ),ans->depth,ans->nChannels);
							IplImage* temp4=cvCreateImage(cvSize( ans->width, ans->height ),ans->depth,ans->nChannels);
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
							//temp1(ans)+temp2
							cvAdd(ans,temp2,ans);
							cvReleaseImage(&temp1);
							cvReleaseImage(&temp2);
							cvReleaseImage(&temp3);
							cvReleaseImage(&temp4);



}


//					fs1									fs1+shiftdown
//pdfsum = pdfs( 1 : 2 : 2 * ht, 2 : 2 : end ) + pdfs( 3 : 2 : end, 2 : 2 : end ) +
//					fs2									fs2+shift left
//		   pdfs( 2 : 2 : end, 1 : 2 : 2 * wt ) + pdfs( 2 : 2 : end, 3 : 2 : end ) ;
void computePdfSum(IplImage* pdfSum,  IplImage* psidashFS1, IplImage* psidashFS2,
				   IplImage* fs1_3222,IplImage* fs1_122ht22,IplImage* fs2_2232,IplImage* fs2_22122wt){	

	IplImage* temp1 =cvCreateImage(cvSize( pdfSum->width, pdfSum->height ),pdfSum->depth,pdfSum->nChannels);	
	IplImage* temp2 =cvCreateImage(cvSize( pdfSum->width, pdfSum->height ),pdfSum->depth,pdfSum->nChannels);
	
	cvZero(fs1_3222);
	cvZero(fs2_2232);
	toolsKit::cvZeroTop(fs1_122ht22);
	toolsKit::cvZeroBottom(fs1_122ht22);	

	cvSetImageROI(fs1_122ht22, cvRect(0,1,fs1_122ht22->width, fs1_122ht22->height));	
	cvCopy(fs1_122ht22, fs1_3222, NULL);
	cvResetImageROI(fs1_122ht22);


	fs1_122ht22->height=fs1_122ht22->height-1;
	//toolsKit::shiftImage(fs1_122ht22,fs1_3222,toolsKit::UP);
	cvAdd(fs1_122ht22,fs1_3222,temp1);

	//cout<<"temp1"<<endl;
	//toolsKit::IPL_print(temp1);

   
	toolsKit::cvZeroLeftRight(fs2_22122wt);

	cvSetImageROI(fs2_22122wt, cvRect(1,0,fs2_22122wt->width, fs2_22122wt->height));
	cvCopy(fs2_22122wt, fs2_2232, NULL);
	cvResetImageROI(fs2_22122wt);
	//toolsKit::shiftImage(fs2_22122wt,fs2_2232,toolsKit::RIGHT);
	fs2_22122wt->width=fs2_22122wt->width-1;
	cvAdd(fs2_22122wt,fs2_2232,temp2);

	//cout<<"temp2"<<endl;
	//toolsKit::IPL_print(temp2);


	//toolsKit::IPL_add_different_sizes(temp1,temp2,pdfSum);
	cvAdd(temp1,temp2,pdfSum);

	cvReleaseImage(&temp1);
	cvReleaseImage(&temp2);
}
				//fs2
//pdfaltsumu = pdfs(2:2:end,  1:2:2*wt) * (u(2:ht+1, 1:wt)  - u(2:ht+1, 2:wt+1) ) + 
//			   pdfs( 2:2:end,  3:2:end) * (u(2:ht+1, 3:end) - u(2:ht+1, 2:wt+1) ) + 
//			   pdfs( 1:2:2*ht, 2:2:end) * (u(1:ht, 2:wt+1)  - u(2:ht+1, 2:wt+1) ) + 
//			   pdfs( 3:2:end,  2:2:end) * (u(3:end, 2:wt+1) - u(2:ht+1, 2:wt+1) ) 


void computeVectBComponents(IplImage* pdfaltSumXX,IplImage* fs1_3222,IplImage* fs1_122ht22,
							IplImage* fs2_2232,IplImage* fs2_22122wt,IplImage* UorV_Org){
	//init
	IplImage* UorV =  cvCreateImage(cvSize( UorV_Org->width+2, UorV_Org->height+2 ),UorV_Org->depth,UorV_Org->nChannels);
	IplImage* tempLeft1=cvCreateImage(cvSize( UorV->width, UorV->height ),UorV->depth,UorV->nChannels);
	IplImage* tempLeft1Reduced=cvCreateImage(cvSize( pdfaltSumXX->width, pdfaltSumXX->height ),pdfaltSumXX->depth,pdfaltSumXX->nChannels);
	IplImage* tempLeft2=cvCreateImage(cvSize( UorV->width, UorV->height ),UorV->depth,UorV->nChannels);
	IplImage* tempLeft2Reduced=cvCreateImage(cvSize( pdfaltSumXX->width, pdfaltSumXX->height ),pdfaltSumXX->depth,pdfaltSumXX->nChannels);
	IplImage* tempLeft3=cvCreateImage(cvSize( UorV->width, UorV->height ),UorV->depth,UorV->nChannels);
	IplImage* tempLeft3Reduced=cvCreateImage(cvSize( pdfaltSumXX->width, pdfaltSumXX->height ),pdfaltSumXX->depth,pdfaltSumXX->nChannels);
	IplImage* tempLeft4=cvCreateImage(cvSize( UorV->width, UorV->height ),UorV->depth,UorV->nChannels);
	IplImage* tempLeft4Reduced=cvCreateImage(cvSize( pdfaltSumXX->width, pdfaltSumXX->height ),pdfaltSumXX->depth,pdfaltSumXX->nChannels);	
	cvZero(UorV);
	cvZero(pdfaltSumXX);
	cvZero(tempLeft1);cvZero(tempLeft1Reduced);cvZero(tempLeft2);cvZero(tempLeft2Reduced);
	cvZero(tempLeft3);cvZero(tempLeft3Reduced);cvZero(tempLeft4);cvZero(tempLeft4Reduced);
	
	
	//pad flow with zeros	
	
	//cvSetImageROI(UorV, cvRect(0,0,UorV->width-2,UorV->height-2));	
	toolsKit::increaseImageSize(UorV_Org,UorV,2);

	//cvCopy(UorV_Org, UorV, NULL);		
	//cvResetImageROI(UorV);
	
	//for testing!
	//UorV=toolsKit::IplFromFile("c:\\a\\upad_test.txt");	
			
	//psidashFS2(2:2:end ,  1:2:2*wt)* ( u(2:ht+1, 1:wt)  - u(2:ht+1, 2:wt+1) ) + //left ,ans2
	toolsKit::IPL_sub_right(UorV,UorV,tempLeft1);
	tempLeft1->width=tempLeft1->width-2;
	tempLeft1->height=tempLeft1->height-1;
	
	cvSetImageROI(tempLeft1, cvRect(0,1,tempLeft1->width, tempLeft1->height));	
	cvCopy(tempLeft1, tempLeft1Reduced, NULL);	
	cvMul(tempLeft1Reduced,fs2_22122wt,tempLeft1Reduced);

	
	//psidashFS2(2:2:end ,  3:2:end) * ( u(2:ht+1, 3:end) - u(2:ht+1, 2:wt+1) ) + //right ,ans2
	 toolsKit::IPL_sub_left(UorV,UorV,tempLeft2);
	 tempLeft2->height=tempLeft2->height-1;

	 cvSetImageROI(tempLeft2, cvRect(2,1,tempLeft2->width, tempLeft2->height));	
	 cvCopy(tempLeft2, tempLeft2Reduced, NULL);		 	
	 cvMul(tempLeft2Reduced,fs2_2232,tempLeft2Reduced);
	
	//psidashFS1(1:2:2*ht,  2:2:end) * ( u(1:ht, 2:wt+1)  - u(2:ht+1, 2:wt+1) ) + //top ,ans1
	toolsKit::IPL_sub_bottom(UorV,UorV,tempLeft3);
	tempLeft3->width=tempLeft3->width-1;
	tempLeft3->height=tempLeft3->height-2;

	cvSetImageROI(tempLeft3, cvRect(1,0,tempLeft3->width, tempLeft3->height));	
	cvCopy(tempLeft3, tempLeft3Reduced, NULL);	
	cvMul(tempLeft3Reduced,fs1_122ht22,tempLeft3Reduced);

	//psidashFS1(3:2:end ,  2:2:end) * ( u(3:end, 2:wt+1) - u(2:ht+1, 2:wt+1) )   //bottom ,ans1
	toolsKit::IPL_sub_top(UorV,UorV,tempLeft4);
	tempLeft4->width=tempLeft4->width-1;
	cvSetImageROI(tempLeft4, cvRect(1,2,tempLeft4->width, tempLeft4->height));	
	cvCopy(tempLeft4, tempLeft4Reduced, NULL);	

	cvMul(tempLeft4Reduced,fs1_3222,tempLeft4Reduced);

	//temp1+temp2+temp3+temp4
	cvAdd(tempLeft1Reduced,tempLeft2Reduced,pdfaltSumXX);
	cvAdd(pdfaltSumXX,tempLeft3Reduced,pdfaltSumXX);
	cvAdd(pdfaltSumXX,tempLeft4Reduced,pdfaltSumXX);
	//- 1*pdfaltsumXX - needed for use in the next computation
	toolsKit::cvMulScalar(pdfaltSumXX,-1);

	//release
	cvReleaseImage(&tempLeft1);
	cvReleaseImage(&tempLeft1Reduced);
	cvReleaseImage(&tempLeft2);
	cvReleaseImage(&tempLeft2Reduced);
	cvReleaseImage(&tempLeft3);
	cvReleaseImage(&tempLeft3Reduced);
	cvReleaseImage(&tempLeft4);	
	cvReleaseImage(&tempLeft4Reduced);
	cvReleaseImage(&UorV);
}


void constructMatrix_brox::computePsidashFS_brox(IplImage* iterU,IplImage* iterV,int width,int height,int channels,flowUV* UV,double _ERROR_CONST){	
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
	IplImage* ux=cvCreateImage(cvSize( width, height ),iterU->depth,channels);
	IplImage* uy=cvCreateImage(cvSize( width, height ),iterU->depth,channels);
	IplImage* vx=cvCreateImage(cvSize( width, height ),iterU->depth,channels);
	IplImage* vy=cvCreateImage(cvSize( width, height ),iterU->depth,channels);
	IplImage* uxd=cvCreateImage(cvSize( width+1, height ),iterU->depth,channels);
	IplImage* vxd=cvCreateImage(cvSize( width+1, height ),iterU->depth,channels);
	IplImage* uyd=cvCreateImage(cvSize( width, height+1 ),iterU->depth,channels);
	IplImage* vyd=cvCreateImage(cvSize( width, height+1 ),iterU->depth,channels);
	IplImage* t  =cvCreateImage(cvSize( width, height+1 ),iterU->depth,channels);
	IplImage* t2  =cvCreateImage(cvSize( width+1, height ),iterU->depth,channels);
	IplImage* uxpd=cvCreateImage(cvSize( width, height+1 ),iterU->depth,channels);
	IplImage* uypd=cvCreateImage(cvSize( width+1, height ),iterU->depth,channels);
	IplImage* vxpd=cvCreateImage(cvSize( width, height+1 ),iterU->depth,channels);
	IplImage* vypd=cvCreateImage(cvSize( width+1, height ),iterU->depth,channels);

	IplImage* temp=cvCreateImage(cvSize( width, height+1 ),iterU->depth,channels);
	IplImage* temp2=cvCreateImage(cvSize( width+1, height ),iterU->depth,channels);;

	cvZero( ux ); 	cvZero( uy ); 	cvZero( vx ); 	cvZero( vy ); 	cvZero( uxd ); 
	cvZero( vxd ); 	cvZero( uyd ); 	cvZero( vyd ); 	cvZero( t ); 	cvZero( t2 ); 
	cvZero( uxpd ); 	cvZero( uypd ); 	cvZero( vxpd ); 	cvZero( vypd ); 
	cvZero( temp ); 	cvZero( temp2 ); 

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
	
	//clean psidash:
	UV->clearPsidash();

	toolsKit::IPL_add_different_sizes2(uypd,vypd,UV->getPsidashFSAns1());
	toolsKit::IPL_add_different_sizes3(uxpd,vxpd,UV->getPsidashFSAns2());	

	toolsKit::psiDerivative(UV->getPsidashFSAns1(),_ERROR_CONST);
	toolsKit::psiDerivative(UV->getPsidashFSAns2(),_ERROR_CONST);	

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


vector<float>*  constructMatrix_brox::constructMatrix_b(IplImage* Ikx,IplImage* Iky,IplImage* Ikz,IplImage* Ixx,
														IplImage* Ixy,IplImage* Iyy,IplImage* Ixz,IplImage* Iyz,											
														flowUV* UV,
														IplImage* du,IplImage* dv,
														SparseMat<float> *  A,
														vector<float> * B,
														vector<float>* dUdV,
														double gamma,double alpha, 
														double _ERROR_CONST,int nInnerFPIterations){
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

			 IplImage* tempU=cvCreateImage(cvSize(Ikx->width, Ikz->height ),Ikz->depth,Ikz->nChannels);
			 IplImage* tempV=cvCreateImage(cvSize(Ikx->width, Ikz->height ),Ikz->depth,Ikz->nChannels);
			
			
			 cvAdd(UV->getU(),du,tempU);
			 cvAdd(UV->getV(),dv,tempV);
			 computePsidashFS_brox(tempU,tempV,UV->getU()->width,UV->getU()->height,UV->getU()->nChannels,UV,_ERROR_CONST);
	
			 toolsKit::cvMulScalar(UV->getPsidashFSAns1(),alpha);
			 toolsKit::cvMulScalar(UV->getPsidashFSAns2(),alpha);

			 cvReleaseImage(&tempU);
			 cvReleaseImage(&tempV);

			//fs1 & fs2 break down
			 IplImage* fs1_3222 =  cvCreateImage(cvSize( Ikx->width, Ikx->height ),Ikx->depth,Ikx->nChannels);
		  	 IplImage* fs1_122ht22=cvCreateImage(cvSize( UV->getPsidashFSAns1()->width, UV->getPsidashFSAns1()->height ),UV->getPsidashFSAns1()->depth,UV->getPsidashFSAns1()->nChannels);
			 fs1_122ht22=cvCloneImage(UV->getPsidashFSAns1());				
			 IplImage* fs2_2232=   cvCreateImage(cvSize( Ikx->width, Ikx->height ),Ikx->depth,Ikx->nChannels);	
			 IplImage* fs2_22122wt=cvCreateImage(cvSize( UV->getPsidashFSAns2()->width, UV->getPsidashFSAns2()->height ),UV->getPsidashFSAns2()->depth,UV->getPsidashFSAns2()->nChannels);
			 fs2_22122wt=cvCloneImage(UV->getPsidashFSAns2());			

			 //epsilon = 1e-3*ones(size(Ikx))==>zeroing and adding instead
			 cvZero(epsilon);
			 cvAddS(epsilon,cvScalarAll(_ERROR_CONST),epsilon);

			 int width=UV->getU()->width;
			 int height=UV->getU()->height;

			 //theta0 = 1/(Ikx^2+Iky^2+epsilon);
			 computeTheta(theta0,Ikx,Iky,epsilon);
			 //theta1 = 1/(Ixx^2+Ixy^2+epsilon);
			 computeTheta(theta1,Ixx,Ixy,epsilon);
			 //theta2 = 1/(Iyy^2+Ixy^2+epsilon);
			 computeTheta(theta2,Iyy,Ixy,epsilon);
					

			 // First compute the values of the data  term
			 //the brightness constancy assumption			 
			 //psiDerivative( theta0 * ( Ikz + Ikx * du + Iky * dv ) ^ 2 );
			 computePsidashBCA(psidashBCA,theta0,Ikz,Ikx,du,Iky,dv,_ERROR_CONST);
			
			 //and the Gradient Constancy Assumption

			 //psidashGCA=psiDerivative( gamma * (  theta1 *  ( Ixz + Ixx * du + Ixy * dv ) ^ 2 + 
			//			                            theta2 *  ( Iyz + Ixy * du + Iyy * dv ) ^ 2 ) ) ;
			 computepsidashGCA(psidashGCA,gamma,theta1,Ixz,Ixx,du,Ixy,dv,theta2,Iyz,Iyy,_ERROR_CONST);
 
			 //now compute the  smoothness term
			 			 
			 //compute pdfSum
			 //pdfsum =pdfs( 1 : 2 : 2 * ht, 2 : 2 : end ) + pdfs( 3 : 2 : end, 2 : 2 : end ) +...
			 //		   pdfs( 2 : 2 : end, 1 : 2 : 2 * wt ) + pdfs( 2 : 2 : end, 3 : 2 : end ) ;
			 
			 computePdfSum(pdfSum, UV->getPsidashFSAns1(),UV->getPsidashFSAns2(),fs1_3222,fs1_122ht22,fs2_2232,fs2_22122wt);										 
			 
			 //prepare data for matrix A

			 //uapp  = psidashBCA * theta0 * ( Ikx ^ 2) +   gamma * psidashGCA * (theta1 *  Ixx ^ 2 +  theta2 * Ixy ^ 2 )  + pdfsum ;
			 computeDiagonalPdfSum(uapp, psidashBCA,theta0,Ikx,gamma,psidashGCA,theta1,Ixx,theta2,Ixy,pdfSum);
			 //vapp  = psidashBCA * theta0 * ( Iky ^ 2) +   gamma * psidashGCA * (theta2 *  Iyy ^ 2 +  theta1 * Ixy ^ 2 )  + pdfsum ;
			 computeDiagonalPdfSum(vapp, psidashBCA,theta0,Iky,gamma,psidashGCA,theta2,Iyy,theta1,Ixy,pdfSum);
			 //uvapp = psidashBCA * theta0* (Ikx*Iky)+ gamma*psidashGCA*(theta1*Ixx*Ixy + theta2*Iyy*Ixy ) ;			
			 computeDiagonalReg   (uvapp,psidashBCA,theta0,Ikx,Iky,gamma,psidashGCA,theta1,Ixx,Ixy,theta2,Iyy,Ixy);
			 //vuapp =   uvapp			 			
			 vuapp=cvCloneImage(uvapp);
			
			 //insert to diagonals to matrix A
			 
			 //uu = spdiags( uapp(:),   0, wt*ht, wt*ht);
			 //vector<float> * uuappCol = toolsKit::IplImageToCoulmnVector(uapp);
			 //SparseMat<float> uu(height*width);
			 //uu.addDiag(0,*uuappCol);
			 //delete uuappCol;
			SparseMat<float> * uu= SparseToolKit::creaseSparse(uapp,0);
			 
				
			 //vv = spdiags( vapp(:),   0, wt*ht, wt*ht);
			 //vector<float> * vvappCol = toolsKit::IplImageToCoulmnVector(vapp);
			 //SparseMat<float> vv(height*width);
			 //vv.addDiag(0,*vvappCol);
			 //delete vvappCol;
			SparseMat<float>* vv = SparseToolKit::creaseSparse(vapp,0);
			

			 //uv = spdiags( uvapp(:), 0, wt*ht, wt*ht);
			 //vector<float> * uvappCol = toolsKit::IplImageToCoulmnVector(uvapp);
			 //SparseMat<float> uv(height*width);
			 //uv.addDiag(0,*uvappCol);
			 //delete uvappCol;
			SparseMat<float>* uv = SparseToolKit::creaseSparse(vuapp,0);
			
			 //vu = spdiags( vuapp(:), 0, wt*ht, wt*ht);
			 //vector<float> * vuappCol = toolsKit::IplImageToCoulmnVector(vuapp);
			 //SparseMat<float> vu(height*width);
			 //vu.addDiag(0,*vuappCol);
			 //delete vuappCol;
			SparseMat<float>* vu = SparseToolKit::creaseSparse(vuapp,0);

			//arguments to u(j) in the linear Euler Lagrange equations.
			//arguments to v(j) in the linear Euler Lagrange equations.

			//negfs1_122ht22 = pdfs( 2 : 2 : end, 1 : 2 : 2 * wt ) 	 		 
			//ul1 = spdiags(-tmp(:), ht, wt*ht, wt*ht);
			IplImage* negfs2_22122wt=cvCreateImage(cvSize(fs2_22122wt->width, fs2_22122wt->height ),fs2_22122wt->depth,fs2_22122wt->nChannels);
			negfs2_22122wt=cvCloneImage(fs2_22122wt);
			toolsKit::cvMulScalar(negfs2_22122wt,-1);			
			//vector<float> * tmp2Col  = toolsKit::IplImageToCoulmnVector(negfs2_22122wt);
			//SparseMat<float> ul1(height*width);
			//ul1.addDiag(height,*tmp2Col);
			SparseMat<float>* ul1 = SparseToolKit::creaseSparse(negfs2_22122wt, height);
		
			
			//vl1 = spdiags(-tmp(:), ht, wt*ht, wt*ht);
			//SparseMat<float> vl1(height*width);
			//vl1.addDiag(height,*tmp2Col);
			//delete tmp2Col;
			SparseMat<float>* vl1 = SparseToolKit::creaseSparse(negfs2_22122wt,height);

			
			cvReleaseImage(&negfs2_22122wt);
			//==========================================

			//negfs2_2232 = pdfs( 2 : 2 : end, 3 : 2 : end ) 
			//ur1 = spdiags(-tmp(:), -ht, wt*ht, wt*ht);
			IplImage* negfs2_2232=cvCreateImage(cvSize(fs2_2232->width, fs2_2232->height ),fs2_2232->depth,fs2_2232->nChannels);		
			negfs2_2232=cvCloneImage(fs2_2232);
			toolsKit::cvMulScalar(negfs2_2232,-1);

			//vector<float> * tmp2Col2  = toolsKit::IplImageToCoulmnVector(negfs2_2232);
			//SparseMat<float> ur1(height*width);
			//ur1.addDiag(-height,*tmp2Col2);
			SparseMat<float>* ur1 = SparseToolKit::creaseSparse(negfs2_2232,-height);
		
			//vr1 = spdiags(-tmp(:), -ht, wt*ht, wt*ht);
			//SparseMat<float> vr1(height*width);
			//vr1.addDiag(-height, *tmp2Col2);
			////no need for tmp2Col			
			//delete tmp2Col2;
			SparseMat<float>* vr1 = SparseToolKit::creaseSparse(negfs2_2232,-height);
			
			cvReleaseImage(&negfs2_2232);
			//==========================================
			
			//negfs1_122ht22 = pdfs( 1 : 2 : 2 * ht, 2 : 2 : end )
			IplImage* negfs1_122ht22=cvCreateImage(cvSize(fs1_122ht22->width, fs1_122ht22->height ),fs1_122ht22->depth,fs1_122ht22->nChannels);
			negfs1_122ht22=cvCloneImage(fs1_122ht22);
			toolsKit::cvMulScalar(negfs1_122ht22,-1);
			
			
			//vector<float> * tmp1Col  = toolsKit::IplImageToCoulmnVector(negfs1_122ht22);
			////ul2 = spdiags(-tmp(:), 1, wt*ht, wt*ht);
			//SparseMat<float> ul2(height*width);
			//ul2.addDiag(1,*tmp1Col);
			SparseMat<float>* ul2 = SparseToolKit::creaseSparse(negfs1_122ht22,1);

			
			//vl2 = spdiags(-tmp(:), 1, wt*ht, wt*ht);
			//SparseMat<float> vl2(height*width);						
			//vl2.addDiag(1, *tmp1Col);			
			SparseMat<float>* vl2 = SparseToolKit::creaseSparse(negfs1_122ht22,1);

			//ur2 = spdiags(-tmp(:), -1, wt*ht, wt*ht);
			//delete tmp1Col;					
			cvReleaseImage(&negfs1_122ht22);
			//==================================
			//tmp = pdfs( 3 : 2 : end, 2 : 2 : end )
			IplImage* negfs1_3222=cvCreateImage(cvSize(fs1_3222->width, fs1_3222->height ),fs1_3222->depth,fs1_3222->nChannels);
			negfs1_3222=cvCloneImage(fs1_3222);
			toolsKit::cvMulScalar(negfs1_3222,-1);


			//vector<float> * tmp1Col2  = toolsKit::IplImageToCoulmnVector(negfs1_3222);//-tmp(:)
			//SparseMat<float> ur2(height*width);
			//ur2.addDiag(-1, *tmp1Col2);
			SparseMat<float> * ur2 = SparseToolKit::creaseSparse(negfs1_3222,-1);

			//vr2 = spdiags(-tmp(:), -1, wt*ht, wt*ht);
			//SparseMat<float> vr2(height*width);
			//vr2.addDiag(-1, *tmp1Col2);
			SparseMat<float>* vr2 = SparseToolKit::creaseSparse(negfs1_3222,-1);

			//no need for negfs1_122ht22, tmp1Col2
			cvReleaseImage(&negfs1_122ht22);

			//A =  [uu+ul1+ul2+ur1+ur2, uv; vu, vv+vl1+vl2+vr1+vr2];

			SparseMat<float> * uuul1ul2ur1ur2 = new SparseMat<float>(uu,ul1,ul2,ur1,ur2);
			
			//cout<<"uuul1ul2ur1ur2:"<<endl<<uuul1ul2ur1ur2<<endl;;
			SparseMat<float> * vvvl1vl2vr1vr2 = new SparseMat<float>(vv,vl1,vl2,vr1,vr2);
			
			//SparseMat<float> * A= new SparseMat<float>(uuul1ul2ur1ur2,uv,vu,vvvl1vl2vr1vr2);				
			A->cube(uuul1ul2ur1ur2,uv,vu,vvvl1vl2vr1vr2);				
			cout<<"finished building mat A"<<endl;
			delete vv;delete uu;
			delete uv;delete vu;
			delete ul1;delete ul2;
			delete vl1;delete vl2;
			delete ur1;delete ur2;
			delete vr1;delete vr2;		
			delete uuul1ul2ur1ur2;
			delete vvvl1vl2vr1vr2;
			//////////////////////build vector B//////////////////////
				
			// Computing the constant terms for the first of the Euler Lagrange equations				
			computeVectBComponents(pdfaltSumU,fs1_3222,fs1_122ht22,fs2_2232,fs2_22122wt,UV->getU());			
			// Computing the constant terms for the second of the Euler Lagrange equations
			computeVectBComponents(pdfaltSumV,fs1_3222,fs1_122ht22,fs2_2232,fs2_22122wt,UV->getV());		
			//constu = psidashBCA * theta0 * ( Ikx * Ikz ) + gamma * psidashGCA * (theta1 * Ixx * Ixz + theta2 * Ixy * Iyz ) - 1*pdfaltsumu 		
			computeDiagonalReg   (constu,psidashBCA,theta0,Ikx,Ikz,gamma,psidashGCA,theta1,Ixx,Ixz,theta2,Ixy,Iyz);					
			cvAdd(constu,pdfaltSumU,constu);		
			//constv = psidashBCA * theta0 * ( Iky * Ikz ) + gamma * psidashGCA * (theta1 * Ixy * Ixz + theta2 * Iyy * Iyz ) - 1*pdfaltsumv ;
			computeDiagonalReg   (constv,psidashBCA,theta0,Iky,Ikz,gamma,psidashGCA,theta1,Ixy,Ixz,theta2,Iyy,Iyz);					
			cvAdd(constv,pdfaltSumV,constv);
			///////////////release all temp iplImages////////////
			cvReleaseImage(&theta0);
			cvReleaseImage(&theta1);
			cvReleaseImage(&theta2);
			cvReleaseImage(&psidashBCA);
			cvReleaseImage(&psidashGCA);
			cvReleaseImage(&epsilon);
			cvReleaseImage(&uapp);
			cvReleaseImage(&vapp);
			cvReleaseImage(&uvapp);
			cvReleaseImage(&vuapp);
			cvReleaseImage(&pdfSum);
			cvReleaseImage(&pdfaltSumU);
			cvReleaseImage(&pdfaltSumV);
			cvReleaseImage(&fs1_3222);
			cvReleaseImage(&fs1_122ht22);
			cvReleaseImage(&fs2_2232);
			cvReleaseImage(&fs2_22122wt);
			/////////////////////////////////////////////////////
			cout<<"starting building B"<<endl;
			//insert data to B vector:b = [-constu(:) ; -constv(:) ];
			IplImage * Mconstu  =  cvCreateImage(cvSize(constu->width,constu->height),constu->depth,constu->nChannels);
			cvCopy(constu,Mconstu);
			toolsKit::cvMulScalar(Mconstu,-1);
			vector<float> * MconstuCol = toolsKit::IplImageToCoulmnVector(Mconstu);
			cvReleaseImage(&Mconstu);
			IplImage * Mconstv = cvCreateImage(cvSize(constv->width,constv->height),constv->depth, constv->nChannels);
			cvCopy(constv,Mconstv);
			toolsKit::cvMulScalar(Mconstv,-1);

			vector<float> * MconstvCol = toolsKit::IplImageToCoulmnVector(Mconstv);
			cvReleaseImage(&Mconstv);
			//vector<float> * B = new vector<float>(MconstuCol->size());

			*B= *MconstuCol;
			//for (vector<float>::iterator it = MconstuCol->begin(); it != MconstuCol->end(); it++)
				//B->push_back(*it);
			for (vector<float>::iterator it = MconstvCol->begin(); it != MconstvCol->end(); it++)
				B->push_back(*it);
			delete MconstuCol;
			delete MconstvCol;
			cvReleaseImage(&constu);
			cvReleaseImage(&constv);
			
			//ofstream thefile("c:\\a\\B_cpp9.txt",ios::out & ios::trunc);thefile<<*B<<endl;thefile.close();
			//ofstream thefile2("c:\\a\\A_cpp9.txt",ios::out & ios::trunc);thefile2<<*A<<endl;thefile2.close();			
			
			cout<<"starting SOR"<<endl;
			float start = std::clock();
			dUdV= SparseToolKit::SOR(A,dUdV,B,1.0,nInnerFPIterations);	
			float diff = ( std::clock() - start ) / (double)CLOCKS_PER_SEC;
			std::cout<<"SOR took: "<< diff <<'\n';
			cout<<"SOR ended!"<<endl;
			
			//ofstream thefile3("c:\\a\\dudv9_cpp_our_sor.txt",ios::out & ios::trunc);thefile3<<*dUdV<<endl;thefile3.close();
			
			cvReleaseImage(&constu);
			cvReleaseImage(&constv);

			//B->clear();
			//delete B;
			//delete A;
			return dUdV;


			
			 
}
