#include "SparseMat.h"
#include "Utils.h"




void printVector(vector<float> vec){
	unsigned int i;
	printf("Printing vector... \n ");
	for(i=0;i<vec.size();i++){

		printf(" %f ",vec.at(i));
	}
	printf(" \n ");

}
void residual(SparseMat<float>& A,const vector<float>& B,const vector<float>& X,vector<float>& residualVec){
	unsigned int i;
	multByVec(A,X,residualVec);
	for(i=0;i<residualVec.size();i++){
		float value = B[i]- residualVec[i] ;
		residualVec[i] = value;
	}

}
float residualNorm(SparseMat<float>& A,const vector<float>& B,vector<float>& X){
	float B_norm = norm(B),ans;
	vector<float> *residualVec = new vector<float>(A.getM());
	residual(A,B,X,*residualVec);
	if(B_norm != 0){
		ans = norm(*residualVec)/B_norm;
	}
	else{
		ans = norm(*residualVec);
	}
	residualVec->clear();
	return ans;
}
float norm(const vector<float>& vec){//Euclidian norm		
	float sum_of_powers=0;
	unsigned int i;
	for(i=0;i<vec.size();i++){
		sum_of_powers+= pow(vec[i],2);
	}
	return sqrt(sum_of_powers);
}
void multByVec(SparseMat<float>& mat,const vector<float>& vec,vector<float>& result){
//	int i,j;		
	if((int)vec.size()!= mat.getN()){
		printf("Error multByVec not proper size\n");
		exit(-1);
	}
	/*for(i=0;i<mat.getM();i++){
	float temp =0;
	for(j=0;j<mat.getN();j++){
	float elem = mat(i,j);
	temp+=elem * vec[j];
	}
	result[i] =temp;
	}*/
	for(SparseMat<float>::const_row_iter row = mat.begin();row != mat.end();row++){	
		float temp =0;
		for(SparseMat<float>::const_col_iter col = row->second.begin(); col != row->second.end(); col++){
			temp+= mat(row->first,col->first)* vec[col->first];
		}
		result[row->first] =temp;
	}
}
bool achievedTolerance(SparseMat<float>& A,const vector<float>& B,const vector<float>& X,float tol){
	vector<float> XCopy(X);
	bool ans = false;
	if(residualNorm(A,B,XCopy) < tol){
		ans = true;
	}
	XCopy.clear();
	return ans;
}
//sign = 0 add, sign = 1 sub
void sumVectors(const vector<float>& vec1,const vector<float>& vec2,vector<float>& vecResult,int sign){
	unsigned int i;
	for(i=0;i<vec1.size();i++){
		float temp;
		if(sign == 0){
			temp = vec1[i] + vec2[i] ;
		}
		else{
			temp = vec1[i] - vec2[i] ;
		}
		vecResult[i] = temp;
	}
}