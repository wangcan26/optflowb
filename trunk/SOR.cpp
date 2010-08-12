	#include "../h/SOR.h"
	
//	#include "../h/MultiGrid.h"
	int iterations = 0;
//	int getSORIterations(){
//		return iterations;
//	}
	float tol = 0.00001;
	vector<float>* MY_SOR(SparseMat<float> A, vector<float> X,vector<float> B, float w, int numOfIterations){
		iterations =0;
		if(B.size()!= X.size()){
			printf("Error solve:B and X should be vectors of the same size \n");
			exit(-1);
		}
		if(A.getM()!= A.getN() || A.getM()!= (int)X.size()){
			printf("Error solve:A should be nxn matrix and X should be vector of size n \n");
			exit(-1);
		}
		
		bool converged = achievedTolerance(A, B, X, tol);
		int i,j,k;

		for(k=0;!converged && k<numOfIterations;k++){
			for(i=0;i<A.getM();i++){
				float temp = 0;
				for(j=0;j< i;j++){//j<i
				 	temp += A(i,j)*X[j];
				}//end j-loop
				for(j=i+1;j< A.getM();j++){//j>i
				 	temp += A(i,j)*X[j];
				}//end j-loop
			if(A(i,i) == 0){
				printf("SOR Error input matrix A :diagonal value 0\n");
				exit(-1);
			}
			float xNext = (1-w)*X[i] + w*(B[i]-temp)/A(i,i);
			X[i] = xNext;//X k+1
			}//end i-loop
			// Terminate if residual treshold reached
			
			if(residualNorm(A,B,X) < tol){
				converged = true;
			}
		}//end k loop
		iterations = k;
		//printf("total iterations %d\n",k);
		vector<float> * ans = new vector<float>(X);
		return ans;
	}


	