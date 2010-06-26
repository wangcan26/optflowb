#include "SparseToolKit.h"


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

vector<float> * SparseToolKit::SOR(SparseMat<float> A, vector<float> x,vector<float> B, float w, int numOfIterations){
		int k,i,j;
		float e1,e2,e3,Aii,Aij;
		vector<float> * temp;
		vector<float> * oldX = new vector<float>(x);
		vector<float> * newX = new vector<float>(x.size());
		for (k=0; k<numOfIterations; k++){
			for (i=0; i<x.size(); i++){
					Aii=A(i,i);
					e1 = B[i]/(Aii!=0?Aii:1);
					e2=0;
					std::map<int, float> row = A.getRow(i);
					//std::map<int,float>::iterator it = row.begin();
					SparseMat<float>::col_iter it = row.begin();
					//for (j=0; j< i-1; j++){
					for(it; it->first < i-1 ; it++){
						j=it->first;
						Aij = A(i,j);
						e2+=Aij*((*newX)[j]);
						}
					e2 = e2*w/(Aii!=0?Aii:1);// (W/Aii)* Sigma(0,i-1){Aij*newX[j]}
					e3=0;
					it++;//skip the Ith element
					//for (j=i+1; j<x.size(); j++){
					for(it; it != row.end() && it->first < x.size(); it++){
						j = it->first;
						Aij = A(i,j);
						e3 += Aij * ((*oldX)[j]);
						//if (it == A.getRow(i).end()) break;
						}
					e3 = e3*w/(Aii!=0?Aii:1); // (W/Aii)* Sigma(i+1,n){Aij*oldX[j]}
					(*newX)[i] = e1 - e2 - e3; //newX[i] = B[i]/Aii - (W/Aii)Sigma(0,i-1){Aij*newX[j]} - (W/Aii)Sigma(i+1,n){Aij*oldX[j]}
				}
			//temp = oldX;
			*oldX = *newX;
			//newX = temp;
			
			}
		return newX;
	}