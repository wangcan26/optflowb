#include "SparseToolKit.h"
#include <limits>


void SparseToolKit::printSparseMat(CvSparseMat* mat){
	
	CvSparseMatIterator it;
	for(CvSparseNode *node = cvInitSparseMatIterator( mat, &it ); node != 0; node = cvGetNextSparseNode( &it )) {
		int* idx = CV_NODE_IDX(mat,node); 
		//double* count=(double*)CV_NODE_VAL(mat,idx);
		//printf( "(i= %d  ,j= %d): %f\n", idx[0], idx[1], *count ); 
		printf( "(i= %d  ,j= %d): %f\n", idx[0], idx[1], cvGet2D(mat,idx[0],idx[1]).val[0]); 
		}	


	}

CvSparseMat * SparseToolKit::lowerTriangle(CvSparseMat * mat){
	CvSparseMat* ans = cvCreateSparseMat(mat->dims,mat->size,mat->type);
	CvSparseMatIterator it;
	for(CvSparseNode *node = cvInitSparseMatIterator( mat, &it ); node != 0; node = cvGetNextSparseNode( &it )) {
		int* idx = CV_NODE_IDX(mat,node); 
		if (idx[0] != idx[1] //not in main diagonal
		&& idx[0]>idx[1])
			cvSet2D(ans,idx[0],idx[1],cvGet2D(mat,idx[0],idx[1]));
		}
	return ans;
	}


static CvSparseMat * upperTriangle(CvSparseMat * mat){
	CvSparseMat* ans = cvCreateSparseMat(mat->dims,mat->size,mat->type);
	CvSparseMatIterator it;
	for(CvSparseNode *node = cvInitSparseMatIterator( mat, &it ); node != 0; node = cvGetNextSparseNode( &it )) {
		int* idx = CV_NODE_IDX(mat,node); 
		if (idx[0] != idx[1] //not in main diagonal
		&& idx[0]<idx[1])
			cvSet2D(ans,idx[0],idx[1],cvGet2D(mat,idx[0],idx[1]));
		}
	return ans;
	}




	float norm(float val){
		float nan1 = sqrt(-1.0f);
		if(toolsKit::AlmostEqualRelativeOrAbsolute(nan1,val,0.00001,0.00001))
			return 0;
		if(toolsKit::AlmostEqualRelativeOrAbsolute(FLT_MAX,val,0.00001,0.00001))
			return FLT_MAX;
		if(toolsKit::AlmostEqualRelativeOrAbsolute(FLT_MIN,val,0.0001,0.0001))
			return FLT_MIN;
		if(toolsKit::IsNan(val))
			return 0;
		return 0;
	}


	SparseMat<float> * SparseToolKit::creaseSparse(IplImage* im, int diag, string filename){
		vector<float> * colVector = toolsKit::IplImageToCoulmnVector(im);
		SparseMat<float> * ans = new SparseMat<float>(im->height * im->width);
		ans->addDiag(diag,*colVector);
		delete colVector;
		if (filename!=""){
			ofstream thefile(filename.c_str(), ios::out & ios::trunc);
			//if (thefile.good()){
				thefile<<*ans<<endl;
			//	}
			thefile.close();
			}
		return ans;
		}

vector<float> * SparseToolKit::SOR(SparseMat<float>* A, vector<float>* x,vector<float>* B, float w, int numOfIterations){
		int k,i,j;
		float e1,e2,e3,Aii,Aij;
		
		//float value = std::numeric_limits<float>::max();
		
		vector<float> * temp;
		vector<float> * oldX = new vector<float>(*x);
		vector<float> * newX = new vector<float>(x->size());
		for (k=0; k<numOfIterations; k++){
			for (i=0; i<x->size(); i++){
					Aii=(*A)(i,i);
					float Bi =(*B)[i];
					e1 = (*B)[i]/(Aii);
					//e1 = B[i];
					e2=0;
					std::map<int, float> row = A->getRow(i);
					//std::map<int,float>::iterator it = row.begin();
					SparseMat<float>::col_iter it = row.begin();
					//for (j=0; j< i-1; j++){
					for(it; it->first <= i-1 ; it++){
						j=it->first;
						Aij = (*A)(i,j);
						float newXj = ((*newX)[j]);
						e2+=Aij*newXj;
						}
					e2 = e2*w/(Aii);// (W/Aii)* Sigma(0,i-1){Aij*newX[j]}
					e3=0;
					if (it->first == i) it++;//skip the Ith element
					//for (j=i+1; j<x.size(); j++){
					for(it; it != row.end() && it->first < x->size(); it++){
						j = it->first;
						Aij = (*A)(i,j);
						float oldXj = ((*oldX)[j]);
						e3 += Aij * oldXj;
						//if (it == A.getRow(i).end()) break;
						}
					e3 = e3*w/(Aii); // (W/Aii)* Sigma(i+1,n){Aij*oldX[j]}
					float iVal = e1 - e2 - e3;
					(*newX)[i] = iVal; //newX[i] = B[i]/Aii - (W/Aii)Sigma(0,i-1){Aij*newX[j]} - (W/Aii)Sigma(i+1,n){Aij*oldX[j]}
				}
			//temp = oldX;
			*oldX = *newX;
			//newX = temp;
			
			}
		return newX;
	}