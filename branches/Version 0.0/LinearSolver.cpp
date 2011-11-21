#include "LinearSolver.h"
#define DEBUG_SOLVER false

LinearSolver::LinearSolver(int iterations,int levels,vector<float>*(*sparseSolver)(SparseMat<float> *A, vector<float> *X,vector<float> *B, float w, int numOfIterations))
	:maxIterations(iterations),maxLevels(levels),solver(sparseSolver)
{
}

LinearSolver::LinearSolver(void):maxIterations(1),maxLevels(2),solver(sparseMatSor){	
}

LinearSolver::~LinearSolver(void)
{
}

/*Fas algorithm functions*/
void restrictionX(const vector<float>& Xh,vector<float>& X2h);
void restrictionV(const vector<float>& Xh,vector<float>& X2h);
void restrictionM(SparseMat<float>& Ah,SparseMat<float>& A2h);
void prolongation(const vector<float>& X2h,vector<float>& Xh);
void coarseGridOperator(SparseMat<float>& Ah,const vector<float>& X2h,vector<float>& B2h);
/*help functions*/
void multByVec(SparseMat<float>& mat,const vector<float>& vec,vector<float>& result);
void residual(SparseMat<float>& A,const vector<float>& B,const vector<float>& X,vector<float>& residualVec);
//sign = 0 add, sign = 1 sub
void sumVectors(const vector<float>& vec1,const vector<float>& vec2,vector<float>& vecResult,int sign);
void printVector(vector<float> &vec);
//Euclidian norm
float norm(const vector<float>& vec);

/*multigrid functions*/
int currentLevels =0;
int iterations = 0;
/*multigrid */
vector<float> * LinearSolver::multigrid(SparseMat<float>& A, vector<float>& B,float w,vector<float>& X,	int preIterations,int postIterations,int cycleType){
	if(!(cycleType == V_CYCLE || cycleType == W_CYCLE)){
		printf("Error multigrid cycleType  must be V_CYCLE or W_CYCLE \n");
		return false;
	}		
	currentLevels =0;				
	for(iterations = 0;iterations < this->maxIterations;iterations++){
		multigridIteration(A,B,w,X,preIterations,postIterations,cycleType);
	}	
	return &X;	
}
/*********************************************************************************************************/
/*multigridIteration**************************************************************************/
/*********************************************************************************************************/
//cycleType = 1 for v-cycle; =2 for w-cycle
void LinearSolver::multigridIteration(SparseMat<float>& A, vector<float>& B,float w,vector<float>& X,int preIterations,int postIterations,int cycleType){	
	int h = X.size();
	if(	(h < 2 ) || ((h%2) != 0) ){
		printf("Error multigrid dimension  must be even number >= 2 \n");
		return ;
	}
	currentLevels++;

	//1.perform relax pre- iterations with the basic solver
	vector<float>* temp =this->solver(&A,&X,&B,w,preIterations);
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
		for(i=0; i< cycleType && currentLevels< (this->maxLevels*cycleType);i++){
			multigridIteration(*A2h,*B2h,w,X2hCopy,preIterations,postIterations,cycleType);
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

	vector<float>* temp1 = this->solver(&A,&X,&B,w,postIterations);
	X.assign(temp1->begin(),temp1->end());	
	if(DEBUG_SOLVER){
		printVector(X);
	}


	//release memory
	temp->clear();
	delete temp;
	temp1->clear();
	delete temp1;
	defect_h->clear();
	delete defect_h;
	X2h->clear();
	delete X2h;
	A2hX2h->clear();
	delete A2hX2h;
	B2h->clear();
	delete B2h;
	defect_2h->clear();
	delete defect_2h;
	X2hCopy.clear();	
	A2h->clean();
	delete A2h;
	defectSolutionDifference2h->clear();
	delete defectSolutionDifference2h;
	defectSolutionDifference->clear();
	delete defectSolutionDifference;
	return;
}
/*********************************************************************************************************/
/*end of multigridIteration**************************************************************************/
/*********************************************************************************************************/
/*SOR function*/
vector<float> * LinearSolver::sparseMatSor(SparseMat<float>* A, vector<float>* X,vector<float>* B, float w, int numOfIterations){	
	int iterations = 0;
	if(B->size()!= X->size()){
		printf("Error solve:B and X should be vectors of the same size \n");
		exit(-1);
	}
	if(A->getM()!= A->getN() || A->getM()!= (int)X->size()){
		printf("Error solve:A should be nxn matrix and X should be vector of size n \n");
		exit(-1);
	}
	int i,j,k;
	vector<float> * oldX = new vector<float>(*X);
	for(k=0; k<numOfIterations;k++){
		//	for(i=0;i<A->getM();i++){
		for (SparseMat<float>::const_row_iter row = A->begin(); row != A->end(); row++){
			float temp = 0;
			//for(j=0;j< i;j++){//j<i		
			i=row->first;
			SparseMat<float>::const_col_iter it = row->second.begin();
			for(it; it->first < i ; it++){
				j=it->first;
				temp += (*A)(i,j)*(*X)[j];
			}//end j-loop
			if (it->first == i) it++;//skip the Ith element
			//for (j=i+1; j<x.size(); j++){
			for(it; it != row->second.end() ; it++){//j>i
				j = it->first;
				//	for(j=i+1;j< A.getM();j++){//j>i
				temp += (*A)(i,j)*(*oldX)[j];
			}//end j-loop
			if((*A)(i,i) == 0){
				printf("SOR Error input matrix A :diagonal value 0\n");
				exit(-1);
			}
			float xNext = (1-w)*(*oldX)[i] + w*((*B)[i]-temp)/(*A)(i,i);
			(*X)[i] = xNext;//X k+1
		}//end i-loop		
	}//end k loop
	iterations = k;
	//printf("total iterations %d\n",k);
	oldX->clear();
	delete oldX;
	vector<float> * ans = new vector<float>(*X);
	return ans;
}// end of SOR




/*********************************************************************************************************/
/*FAS algorithm functions**************************************************************************/
/*********************************************************************************************************/

//build coarse grid X2h of Xh using avg of 4 cells, result is saved in X2h(use for coarse X)
void restrictionX(const vector<float>& Xh,vector<float>& X2h){
	if(Xh.size()!= 2*X2h.size() ){
		printf("Error restrictionX size of Xh is 2* size of X2h\n");
		exit(-1);
	}
	unsigned int i;		
	for(i=0;i<X2h.size();i++){									
		float avg = (Xh[2*i]+ Xh[2*i+1])/2;
		X2h[i] = avg;					
	}		
}
//build coarse grid X2h of Xh using sum of 4 cells, result is saved in X2h(use for coarse residual)
void restrictionV(const vector<float>& Xh,vector<float>& X2h){
	if(Xh.size()!= 2*X2h.size() ){
		printf("Error restrictionX size of Xh is 2* size of X2h\n");
		exit(-1);
	}
	unsigned int i;		
	for(i=0;i<X2h.size();i++){									
		float sum = Xh[2*i]+ Xh[2*i+1];
		X2h[i] = sum;					
	}	

}
void restrictionM(SparseMat<float>& Ah,SparseMat<float>& A2h){
	if(Ah.getM()!= 2*A2h.getM() ){
		printf("Error restrictionX size of Xh is 2* size of X2h\n");
		exit(-1);
	}
	int i,j;	
	for (SparseMat<float>::const_row_iter row = Ah.begin(); row != Ah.end(); row++){
		SparseMat<float>::const_row_iter rowTmp = row;
		rowTmp++;
		bool lastRow = false;
		if(rowTmp == Ah.end()){
			lastRow = true;
		}		
		for( SparseMat<float>::const_col_iter col = row->second.begin(); col != row->second.end(); col++){
			/*SparseMat<float>::const_col_iter colTmp = col;
			colTmp++;*/
			float sum;
			int i=row->first;
			int j = col->first;
			if(lastRow){
				if((j%2) == 0){//EVEN col
					sum = Ah(i,j) + Ah(i,j+1);
					col++;
				}
				else{//ODD col
					sum = Ah(i,j);
				}
			}
			else{//not last row
				if((j%2) == 0){//EVEN col
					sum = Ah(i,j) + Ah(i,j+1)+ Ah(i+1,j)+ Ah(i+1,j+1);
					col++;
				}
				else{//ODD col
					sum = Ah(i,j) + Ah(i+1,j);
				}					
			}
			A2h(i/2,j/2) = sum;
		}		
		row = rowTmp;
	}		
}

void prolongation(const vector<float>& X2h,vector<float>& Xh){
	if(Xh.size()!= 2*X2h.size() ){
		printf("Error prolongation size of Xh is 2* size of X2h\n");
		exit(-1);
	}
	unsigned int i;		
	for(i=0;i<X2h.size();i++){									
		float value = X2h[i];
		Xh[2*i] = value;
		Xh[2*i+1] = value;					
	}	
}

/**1:build fine grid Xh from X2h;
*  2:apply matrix Ah on Xh save the result in Bh
*  3:build coarse grid B2h from Bh
**/
void coarseGridOperator(SparseMat<float>& Ah,const vector<float>& X2h,vector<float>& B2h){	
	vector<float> *Xh= new vector<float>(2*X2h.size());
	vector<float> *Bh= new vector<float>(2*X2h.size());
	prolongation(X2h,*Xh);	
	multByVec(Ah,*Xh,*Bh);	
	restrictionV(*Bh,B2h);	
	Xh->clear();
	delete Xh;
	Bh->clear();
	delete Bh;
}
/*********************************************************************************************************/
/*end of FAS algorithm functions**************************************************************************/
/*********************************************************************************************************/

/*********************************************************************************************************/
/*help functions**************************************************************************/
/*********************************************************************************************************/
void multByVec(SparseMat<float>& mat,const vector<float>& vec,vector<float>& result){		
	if((int)vec.size()!= mat.getN()){
		printf("Error multByVec not proper size\n");
		exit(-1);
	}	
	int i,j;	
	for (SparseMat<float>::const_row_iter row = mat.begin(); row != mat.end(); row++){
		float temp =0;
		i=row->first;
		for( SparseMat<float>::const_col_iter col = row->second.begin(); col != row->second.end(); col++){
			j = col->first;
			temp+=mat(i,j) * vec[j];
		}
		result[i] =temp;
	}
}
void residual(SparseMat<float>& A,const vector<float>& B,const vector<float>& X,vector<float>& residualVec){
	unsigned int i;
	multByVec(A,X,residualVec);
	for(i=0;i<residualVec.size();i++){
		float value = B[i]- residualVec[i] ;
		residualVec[i] = value;
	}

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

void printVector(vector<float> &vec){
	unsigned int i;
	printf("Printing vector... \n ");
	for(i=0;i<vec.size();i++){

		printf(" %f ",vec.at(i));
	}
	printf(" \n ");

}


/*********************************************************************************************************/
/*end of help functions**************************************************************************/
/*********************************************************************************************************/
