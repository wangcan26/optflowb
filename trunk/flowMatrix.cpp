#include "flowMatrix.h"
#include <ctime>


flowMatrix::flowMatrix(void)
{
}


flowMatrix::~flowMatrix(void)
{
}





static void psiDash(IplImage* src,IplImage* ans,double _ERROR_CONST)
{
	// 1/( 2 sqrt(x+eps) )
	cvAddS(src,cvScalarAll(_ERROR_CONST),ans);
	cvPow(ans,ans,0.5);
	cvAdd(ans,ans,ans);							   							
	cvPow(ans,ans,-1);
}

static void computePdfSum(IplImage* pdfSum,  IplImage* psidashFS1, IplImage* psidashFS2,
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



static void computeVectBComponents(IplImage* pdfaltSumXX,IplImage* fs1_3222,IplImage* fs1_122ht22,
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


void computePsidashFS_brox(IplImage* iterU,IplImage* iterV,int width,int height,int channels,flowUV* UV,double _ERROR_CONST){	
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



 vector<float>*  flowMatrix::constructMatrix(IplImage* Ix, IplImage* Iy, 	IplImage* Iz,						
								  flowUV* UV, IplImage* du, IplImage* dv,vector<float>* dUdV,
								  SparseMat<float> *  A,vector<float> * B,
								  double alpha,
								  double _ERROR_CONST)
{
	
	float start = std::clock();
	std::cout<<"Building Matrix A " ; 
	int width=UV->getU()->width;
	int height=UV->getU()->height;


	IplImage* uapp=	 cvCreateImage(cvSize(width, height ),Ix->depth,Ix->nChannels);
	IplImage* vapp=	 cvCreateImage(cvSize(width, height ),Ix->depth,Ix->nChannels);
	IplImage* uvapp= cvCreateImage(cvSize(width, height ),Ix->depth,Ix->nChannels);
	IplImage* pdfSum = cvCreateImage(cvSize(width, height ),Ix->depth,Ix->nChannels);
	

	IplImage* psiDataTerm= cvCreateImage(cvSize(Ix->width, Ix->height ),Ix->depth,Ix->nChannels);
	IplImage* dataTermNorm= cvCreateImage(cvSize(width, height ),Ix->depth,Ix->nChannels);

	IplImage* temp= cvCreateImage(cvSize(width, height ),Ix->depth,Ix->nChannels);
	IplImage* temp1= cvCreateImage(cvSize(width, height ),Ix->depth,Ix->nChannels);

	IplImage* bu= cvCreateImage(cvSize(width, height ),Ix->depth,Ix->nChannels);
	IplImage* bv= cvCreateImage(cvSize(width, height ),Ix->depth,Ix->nChannels);
	IplImage* pdfaltSumU= cvCreateImage(cvSize(width, height ),Ix->depth,Ix->nChannels);
	IplImage* pdfaltSumV= cvCreateImage(cvSize(width, height ),Ix->depth,Ix->nChannels);

	
	/***********************   SMOOTHNESS TERM ************************/

	 cvAdd(UV->getU(),du,temp);
	 cvAdd(UV->getV(),dv,temp1);
     computePsidashFS_brox(temp,temp1,UV->getU()->width,UV->getU()->height,UV->getU()->nChannels,UV,_ERROR_CONST);
	
     toolsKit::cvMulScalar(UV->getPsidashFSAns1(),alpha);
	 toolsKit::cvMulScalar(UV->getPsidashFSAns2(),alpha);

	 //fs1 & fs2 break down
	 IplImage* fs1_3222 =  cvCreateImage(cvSize( width, height ),Ix->depth,Ix->nChannels);
	 IplImage* fs1_122ht22=cvCreateImage(cvSize( UV->getPsidashFSAns1()->width, UV->getPsidashFSAns1()->height ),UV->getPsidashFSAns1()->depth,UV->getPsidashFSAns1()->nChannels);
	 fs1_122ht22=cvCloneImage(UV->getPsidashFSAns1());				
	 IplImage* fs2_2232=   cvCreateImage(cvSize( width, height ),Ix->depth,Ix->nChannels);	
	 IplImage* fs2_22122wt=cvCreateImage(cvSize( UV->getPsidashFSAns2()->width, UV->getPsidashFSAns2()->height ),UV->getPsidashFSAns2()->depth,UV->getPsidashFSAns2()->nChannels);
	 fs2_22122wt=cvCloneImage(UV->getPsidashFSAns2());	
	 computePdfSum(pdfSum, UV->getPsidashFSAns1(),UV->getPsidashFSAns2(),fs1_3222,fs1_122ht22,fs2_2232,fs2_22122wt);

	/***********************     DATA TERM     ************************/
	
	
	//dataTermNorm = 1/(Ix^2+Iy^2+epsilon);
	cvPow(Ix,temp1,2);	
	cvPow(Iy,dataTermNorm,2);
	cvAdd(dataTermNorm,temp1,dataTermNorm);
	cvAddS(dataTermNorm,cvScalarAll(_ERROR_CONST),dataTermNorm);
	cvPow(dataTermNorm,dataTermNorm,-1);
	

	//psi'( dataTermNorm * ( Iz + Ix * du + Iy * dv ) ^ 2 );
	cvMul(Ix,du,psiDataTerm);
	cvMul(Iy,dv,temp);			
	cvAdd(psiDataTerm,temp,psiDataTerm);
	cvAdd(psiDataTerm,Iz,psiDataTerm);
	cvPow(psiDataTerm,psiDataTerm,2);
	cvMul(psiDataTerm,dataTermNorm,psiDataTerm);
	psiDash(psiDataTerm,psiDataTerm,_ERROR_CONST);
	

	//uapp  = psidashBCA * dataTermNorm * ( Ix ^ 2)   + pdfsum 
	cvPow(Ix,uapp,2);
	cvMul(dataTermNorm,uapp,uapp);							   							
	cvMul(psiDataTerm,uapp,uapp);
	cvAdd(uapp,pdfSum,uapp);
	
	 //vapp  = psi'( dataTermNorm * ( Iy ^ 2)  + pdfsum ;
	cvPow(Iy,vapp,2);
	cvMul(dataTermNorm,vapp,vapp);							   							
	cvMul(psiDataTerm,vapp,vapp);
	cvAdd(vapp,pdfSum,vapp);
	
			 
	//uvapp = psi' ( dataTermNorm * (Ix * Iy) ) ;			
	cvMul(Ix, Iy, uvapp);
	cvMul(dataTermNorm,uvapp,uvapp);							   							
	cvMul(psiDataTerm,uvapp,uvapp);
			 
	//vuapp =   uvapp			 			
		

	/********************** Now construct Matrix A *********************/
	SparseMat<float> * uu= SparseToolKit::creaseSparse(uapp,0);
	SparseMat<float>* vv = SparseToolKit::creaseSparse(vapp,0);
	SparseMat<float>* uv = SparseToolKit::creaseSparse(uvapp,0);
	//SparseMat<float>* vu = SparseToolKit::creaseSparse(uvapp,0);



	
	IplImage* negPsiUV=cvCreateImage(cvSize(fs2_22122wt->width, fs2_22122wt->height ),fs2_22122wt->depth,fs2_22122wt->nChannels);
	negPsiUV=cvCloneImage(fs2_22122wt);
	toolsKit::cvMulScalar(negPsiUV,-1);			
	SparseMat<float>* ul1 = SparseToolKit::creaseSparse(negPsiUV, height);
	
			
	negPsiUV=cvCloneImage(fs2_2232);
	toolsKit::cvMulScalar(negPsiUV,-1);
	SparseMat<float>* ur1 = SparseToolKit::creaseSparse(negPsiUV,-height);
			
	
	negPsiUV=cvCloneImage(fs1_122ht22);
	toolsKit::cvMulScalar(negPsiUV,-1);
	SparseMat<float>* ul2 = SparseToolKit::creaseSparse(negPsiUV,1);
	
	
	negPsiUV=cvCloneImage(fs1_3222);
	toolsKit::cvMulScalar(negPsiUV,-1);
	SparseMat<float>* ur2 = SparseToolKit::creaseSparse(negPsiUV,-1);
		
		
			
	SparseMat<float> * u = new SparseMat<float>(uu,ul1,ul2,ur1,ur2);
	SparseMat<float> * v = new SparseMat<float>(vv,ul1,ul2,ur1,ur2);
	//SparseMat<float> * A= new SparseMat<float>(u,uv,uv,v);				
	
	A->cube(u,uv,uv,v);	
	/********************** Now construct Vector  b *********************/
	std::cout<<"and Vector b " ; 
		// Computing the constant terms for the first of the Euler Lagrange equations				
			computeVectBComponents(pdfaltSumU,fs1_3222,fs1_122ht22,fs2_2232,fs2_22122wt,UV->getU());			
			//toolsKit::IplToFile(pdfaltSumU,"c:\\a\\pdfaultSumU_cpp.txt");
			// Computing the constant terms for the second of the Euler Lagrange equations
			computeVectBComponents(pdfaltSumV,fs1_3222,fs1_122ht22,fs2_2232,fs2_22122wt,UV->getV());
			//toolsKit::IplToFile(pdfaltSumV,"c:\\a\\pdfaultSumV_cpp.txt");		
	

	//bu = psi' * ( dataTermNorm * Ix * Iz) - pdfaltsumu
	cvMul(dataTermNorm,psiDataTerm,bu);							   							
	cvMul(Ix,bu,bu);
	cvMul(Iz,bu,bu);
	cvAdd(bu,pdfaltSumU,bu);
	toolsKit::cvMulScalar(bu,-1);


	//bv = psi' * ( dataTermNorm * Ix * Iz) - pdfaltsumv
	cvMul(dataTermNorm,psiDataTerm,bv);							   							
	cvMul(Iy,bv,bv);
	cvMul(Iz,bv,bv);
	cvAdd(bv,pdfaltSumV,bv);
	toolsKit::cvMulScalar(bv,-1);

	//insert data to B vector:b = [-constu(:) ; -constv(:) ];
	//vector<float> * MconstuCol = toolsKit::IplImageToCoulmnVector(bu);
	vector<float> * MconstvCol = toolsKit::IplImageToCoulmnVector(bv);
		

	*B= *toolsKit::IplImageToCoulmnVector(bu);
	for (vector<float>::iterator it = MconstvCol->begin(); it != MconstvCol->end(); it++)
			B->push_back(*it);
		
			

	
	/********************** Clean the memory *********************/	


	delete vv;
	delete uu;
	delete uv;
	delete ul1;
	delete ul2;
	delete ur1;
	delete ur2;		
	delete u;
	delete v;
	//delete MconstuCol;
	delete MconstvCol;

	cvReleaseImage(&negPsiUV);
	cvReleaseImage(&psiDataTerm);
	cvReleaseImage(&dataTermNorm);
	cvReleaseImage(&temp);
	cvReleaseImage(&temp1);
	cvReleaseImage(&uapp);
	cvReleaseImage(&vapp);
	cvReleaseImage(&uvapp);
	cvReleaseImage(&pdfSum);
	cvReleaseImage(&fs1_3222);
	cvReleaseImage(&fs1_122ht22);
	cvReleaseImage(&fs2_2232);
	cvReleaseImage(&fs2_22122wt);
	cvReleaseImage(&bu);
	cvReleaseImage(&bv);
	cvReleaseImage(&pdfaltSumU);
	cvReleaseImage(&pdfaltSumV);
		
	
	


	//delete B;
	//delete A;
	float diff = ( std::clock() - start ) / (double)CLOCKS_PER_SEC;
	std::cout<<" --- "<< diff <<'\n';

	return dUdV;
			 
}

