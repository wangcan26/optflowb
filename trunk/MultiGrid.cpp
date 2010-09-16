#include "MultiGrid.h"


int currentLevels =0;
int iterations = 0;
void multigridIteration(vector<float>* (*solver)(SparseMat<float> *A, vector<float> *X,vector<float> *B, float w, int numOfIterations),
	SparseMat<float>& A, vector<float>& B,float w,vector<float>& X,
	int preIterations,int postIterations,int cycleType);


vector<float> * multigrid(vector<float>* (*solver)(SparseMat<float> *A, vector<float> *X,vector<float> *B, float w, int numOfIterations),
	SparseMat<float>& A, vector<float>& B,float w,vector<float>& X,
	int preIterations,int postIterations,int cycleType){
		if(!(cycleType == V_CYCLE || cycleType == W_CYCLE)){
			printf("Error multigrid cycleType  must be V_CYCLE or W_CYCLE \n");
			return false;
		}		
		currentLevels =0;		
		//do until tolerance condition or MAX_ITERATIONS
		for(iterations = 0;iterations < MAX_ITERATIONS;iterations++){
			multigridIteration(solver,A,B,w,X,preIterations,postIterations,cycleType);
		}	
		return &X;	
}

//cycleType = 1 for v-cycle; =2 for w-cycle
void multigridIteration(vector<float>* (*solver)(SparseMat<float> *A, vector<float> *X,vector<float> *B, float w, int numOfIterations),
	SparseMat<float>& A, vector<float>& B,float w,vector<float>& X,
	int preIterations,int postIterations,int cycleType){	


		int h = X.size();

		if(	(h < 2 ) || ((h%2) != 0) ){
			printf("Error multigrid dimension  must be even number >= 2 \n");
			return ;
		}
		currentLevels++;
		//for(int i=0; i< cycleType	;i++){

		//1.perform relax pre- iterations with the basic solver
		vector<float>* temp =(*solver)(&A,&X,&B,w,preIterations);
		X.assign(temp->begin(),temp->end());
		//2.coarse grid correction
		vector<float> *defect_h= new vector<float>(h);
		//compute(Ax-B)
		residual(A,B,X,	*defect_h);			
		int size2h = h/2;

		//3.compute X2h -approximation of the X on the next coarser grid
		vector<float> *X2h = new vector<float>(size2h);//size X2h
		restrictionX(X,*X2h);

		//4.compute B2h
		vector<float> *A2hX2h= new vector<float>(size2h);//size temp2h
		vector<float> *B2h= new vector<float>(size2h);//size B2h
		coarseGridOperator(A,*X2h,*A2hX2h);

		vector<float> *defect_2h= new vector<float>(size2h);
		restrictionV(*defect_h,*defect_2h);

		sumVectors(*A2hX2h,*defect_2h,*B2h,0);
		//solve A2h*X2h = B2h
		vector<float> X2hCopy(*X2h);
		SparseMat<float> *A2h = new SparseMat<float>(size2h,size2h);
		restrictionM(A,*A2h);
		//Solve error equation

		if((size2h >= 2) && ((size2h%2) == 0)){
			int i;	
			for(i=0; i< cycleType && currentLevels< MAX_LEVELS;i++){
				multigridIteration(solver,*A2h,*B2h,w,X2hCopy,preIterations,postIterations,cycleType);
			}					
		}

		vector<float> *defectSolutionDifference2h = new vector<float>(size2h);
		vector<float> *defectSolutionDifference = new vector<float>(X.size());
		//defectSolutionDifference = X2hCopy - X2h
		sumVectors(X2hCopy,*X2h,*defectSolutionDifference2h,1);
		for(int i=0; i< size2h	;i++){
			(*defectSolutionDifference2h)[i] /= pow((float)DEFECT_FACTOR,currentLevels);
		}
		prolongation(*defectSolutionDifference2h,*defectSolutionDifference);

		//X = X + Difference Of defect
		sumVectors(X,*defectSolutionDifference,X,0);


		vector<float>* temp1 = (*solver)(&A,&X,&B,w,postIterations);
		X.assign(temp1->begin(),temp1->end());

		//	}		
		return;
}
int getMULTIGRIDIterations(){
	return currentLevels;
}
