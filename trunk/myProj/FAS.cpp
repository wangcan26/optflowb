#include "FAS.h"
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
	for(SparseMat<float>::const_row_iter row = Ah.begin();row != Ah.end();row++){	
		SparseMat<float>::const_row_iter rowTemp = row;
		bool lastLine = false;
		rowTemp++;
		if(rowTemp == Ah.end()){
			lastLine = true;
		}
		for(SparseMat<float>::const_col_iter col = row->second.begin(); col != row->second.end(); col++){
			int i = row->first;
			int j = col->first;
			float sum =0;
			if(!lastLine){
				if(j+1 >= Ah.getN() || j%2 !=0){//last column or odd
					sum = Ah(i,j) + Ah(i+1,j);
				}
				else{
					sum = Ah(i,j) + Ah(i,j+1)+ Ah(i+1,j)+ Ah(i+1,j+1);
				}
			}
			else{//last line
				if(j+1 >= Ah.getN()|| j%2 !=0){//last column or odd
					sum = Ah(i,j);
				}
				else{
					sum = Ah(i,j) + Ah(i,j+1);
				}
			}
			A2h(i/2,j/2) = sum;
			SparseMat<float>::const_col_iter colTemp = col;
			colTemp++;
			if((colTemp != row->second.end())&& (j %2 == 0) && (colTemp->first == j+1)){
				col++;
			}
		}
		row++;
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
void prolongationD(vector<float>& X2h,vector<float>& Xh){
	if(Xh.size()!= 2*X2h.size() ){
		printf("Error prolongation size of Xh is 2* size of X2h\n");
		exit(-1);
	}
	unsigned int i;		
	for(i=0;i<X2h.size();i++){

		X2h[i] /= 1000;

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
	//	printf("start coarse prolong \n");
	//	printVector(X2h);
	//	printVector(*Xh);
	multByVec(Ah,*Xh,*Bh);
	//	printf("start coarse mult \n");
	//	printVector(*Bh);
	//	printVector(*Xh);
	//	printMatrix(Ah);
	restrictionV(*Bh,B2h);	
	//	printf("start coarse restr \n");
	//	printVector(*Bh);
	//	printVector(B2h);
	//	printf("end \n");
}