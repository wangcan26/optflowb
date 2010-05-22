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