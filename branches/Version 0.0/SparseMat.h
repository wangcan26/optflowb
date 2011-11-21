#pragma once

#include <map>
#include <vector>
#include <iostream>
#include <cv.h>

using namespace std;

template<class T>
class SparseMat
	{

	public:		
		typedef std::map<int, T> col_t;

		typedef std::map<int, col_t > mat_t;
		typedef typename mat_t::iterator row_iter;
		typedef typename mat_t::const_iterator const_row_iter;
		
		typedef typename col_t::iterator col_iter;
		typedef typename col_t::const_iterator const_col_iter;

		SparseMat(){ m=n=1;}
		template <class T2>
		SparseMat (const SparseMat<T2>& c):mat(c),m(c.m),n(c.n){}
		SparseMat(int i) {m=n=i;}
		SparseMat(int i, int j){m=i; n=j;}
		/*
			==========
			| q1 | q2|
			==========
			| q3 | q4|
			==========
		*/
		template<class T2>
		SparseMat<T> (const SparseMat<T2>* q1, const SparseMat<T2>* q2, const SparseMat<T2>* q3, const SparseMat<T2>* q4):mat(q1->mat){
		//	mat=q1.mat;
			m = q1->m + q3->m;
			n = q1->n + q2->n;
			//add q2
			for (const_row_iter row = q2->mat.begin(); row != q2->mat.end(); row++){
				for( const_col_iter col = row->second.begin(); col != row->second.end(); col++){
						mat[row->first][q1->n+ col->first] = col->second;
					}
				}
			//add q3
			for (const_row_iter row = q3->mat.begin(); row != q3->mat.end(); row++){
				for( const_col_iter col = row->second.begin(); col != row->second.end(); col++){
						mat[q1->m+row->first][col->first] = col->second;
					}
				}
			//add q4
			for (const_row_iter row = q4->mat.begin(); row != q4->mat.end(); row++){
				for( const_col_iter col = row->second.begin(); col != row->second.end(); col++){
						mat[q1->m+row->first][q1->n+col->first] = col->second;
					}
				}
			}
		
		void cube(const SparseMat<T>* q1, const SparseMat<T>* q2, const SparseMat<T>* q3, const SparseMat<T>* q4){
			if (m != q1->m + q3->m || n != q1->n + q2->n)
				cout<<" might be something strange here"<<endl;
			this->mat.clear();
			for (const_row_iter row = q1->mat.begin(); row != q1->mat.end(); row++)
				for (const_col_iter col = row->second.begin(); col != row->second.end(); col++)
					this->mat[row->first][col->first] = col->second;	
			//add q2
			for (const_row_iter row = q2->mat.begin(); row != q2->mat.end(); row++){
				for( const_col_iter col = row->second.begin(); col != row->second.end(); col++){
					mat[row->first][q1->n+ col->first] = col->second;
					}
				}
			//add q3
			for (const_row_iter row = q3->mat.begin(); row != q3->mat.end(); row++){
				for( const_col_iter col = row->second.begin(); col != row->second.end(); col++){
					mat[q1->m+row->first][col->first] = col->second;
					}
				}
			//add q4
			for (const_row_iter row = q4->mat.begin(); row != q4->mat.end(); row++){
				for( const_col_iter col = row->second.begin(); col != row->second.end(); col++){
					mat[q1->m+row->first][q1->n+col->first] = col->second;
					}
				}
			}


		//assuming all matrixes are the same size
		template<class T2>
		SparseMat<T> (SparseMat<T2>* A1 ,SparseMat<T2>* A2 , SparseMat<T2>* A3 ,SparseMat<T2>* A4 ,SparseMat<T2>* A5  ) {
				m=A1->m;
				n=A1->n;
				for (row_iter ii = A1->mat.begin(); ii !=A1->mat.end(); ii++)
					for (col_iter jj = ii->second.begin(); jj!= ii->second.end(); jj++)
						mat[ii->first][jj->first] = jj->second;
				for (row_iter ii = A2->mat.begin(); ii !=A2->mat.end(); ii++)
					for (col_iter jj = ii->second.begin(); jj!= ii->second.end(); jj++)
						mat[ii->first][jj->first] += jj->second;
				for (row_iter ii = A3->mat.begin(); ii !=A3->mat.end(); ii++)
					for (col_iter jj = ii->second.begin(); jj!= ii->second.end(); jj++)
						mat[ii->first][jj->first] += jj->second;
				for (row_iter ii = A4->mat.begin(); ii !=A4->mat.end(); ii++)
					for (col_iter jj = ii->second.begin(); jj!= ii->second.end(); jj++)
						mat[ii->first][jj->first] += jj->second;
				for (row_iter ii = A5->mat.begin(); ii !=A5->mat.end(); ii++)
					for (col_iter jj = ii->second.begin(); jj!= ii->second.end(); jj++)
						mat[ii->first][jj->first] += jj->second;
			}

		IplImage* toIpl(){
			IplImage* ans = cvCreateImage(cvSize(n,m),IPL_DEPTH_32F,1);
			for (row_iter ii = mat.begin(); ii!=mat.end(); ii++)
				for (col_iter jj = ii->second.begin(); jj!= ii->second.end(); jj++)
					cvSetReal2D(ans,ii->first,jj->first,jj->second);
			return ans;
			}

		//Operators:
		inline T& operator()(int i, int j)
			{
			if(i>=m || j>=n) throw;
			return mat[i][j];
			}
		inline T operator()(int i, int j) const
			{
			if(i>=m || j>=n) throw;
			return mat[i][j];
			}



	
		friend ostream& operator <<(ostream& lhs, const SparseMat& rhs)
			{
			const_row_iter ii = rhs.mat.begin();
			const_col_iter jj;
			lhs<<"m="<<rhs.m<<" n="<<rhs.n<<endl;
			for (ii; ii != rhs.mat.end(); ii++){
				for (jj=ii->second.begin(); jj != ii->second.end(); jj++){
					lhs<<"("<<(ii->first)<<","<<(jj->first)<<")"<<"="<<jj->second<<"\t";
					}
				lhs<<endl;
				}

			return lhs;
			}

		friend ofstream& operator <<(ofstream& lhs, SparseMat& rhs)
			{
			for (int j=0; j<rhs.m; j++)
				for (int i=0; i<rhs.n; i++){
					if (rhs.mat[i][j]!=0)
						lhs<<"("<<i+1<<","<<j+1<<")"<<"\t"<<rhs.mat[i][j]<<endl;
					}
			return lhs;
			}

			

		//SparseMat operator*(const T& val)const {

		//	SparseMat<T> ans (this->m,this->n);
		//		
		//	for (const_row_iter ii = mat.begin(); ii != mat.end(); ii++){
		//		for (const_col_iter jj = ii->second.begin(); jj != ii->second.end(); jj++){
		//				
		//				ans(ii->first,jj->first) = jj->second * val;
		//			}
		//		}

		//	return ans;
		//	}


		//SparseMat& operator*=(const T& val){
		//	for (row_iter ii = mat.begin(); ii != mat.end(); ii++){
		//		for (col_iter jj = ii->second.begin(); jj != ii->second.end(); jj++){

		//			jj->second *= val;
		//			}
		//		}

		//	return *this;
		//	}


		//SparseMat operator+(const T& val)const {

		//	SparseMat<T> ans (this->m,this->n);

		//	for (const_row_iter ii = mat.begin(); ii != mat.end(); ii++){
		//		for (const_col_iter jj = ii->second.begin(); jj != ii->second.end(); jj++){

		//			ans(ii->first,jj->first) = jj->second + val;
		//			}
		//		}

		//	return ans;
		//	}


		SparseMat& operator+=(const T& val){
			for (row_iter ii = mat.begin(); ii != mat.end(); ii++){
				for (col_iter jj = ii->second.begin(); jj != ii->second.end(); jj++){

					jj->second += val;
					}
				}

			return *this;
			}



		template <class T2>	SparseMat operator+(const SparseMat<T2>& other)const {

			if (this->m != other.m || this->n != other.n) throw;

			SparseMat<T> ans(this->m, this->n);

			//ans.mat = this->mat;
			for (const_row_iter ii = this->mat.begin(); ii != this->mat.end(); ii++){
				for (const_col_iter jj = ii->second.begin(); jj != ii->second.end(); jj++){
					ans(ii->first,jj->first) = jj->second ;
					}
				}



			for (const_row_iter ii = other.mat.begin(); ii != other.mat.end(); ii++){
				for (const_col_iter jj = ii->second.begin(); jj != ii->second.end(); jj++){
					ans(ii->first,jj->first) += jj->second ;
					}
				}

			return ans;
			}



		template <class T2> SparseMat& operator+=(const SparseMat<T2>& other){
			for (const_row_iter ii = other.mat.begin(); ii != other.mat.end(); ii++){
				for (const_col_iter jj = ii->second.begin(); jj != ii->second.end(); jj++){

					mat[ii->first][jj->first] += jj->second;
					}
				}

			return *this;
			}




		/*clean the matrix form 0 valued elements*/
		void clean(){
			for (row_iter ii=mat.begin(); ii != mat.end(); ii++){
				for (col_iter jj = ii->second.begin(); jj != ii->second.end(); ){
					if (jj->second ==0){
						col_iter del = jj++;
						ii->second.erase(del);
						}
					else
						jj++;
					}
				}
			}
		
		
		
	

		/*
			add the vector elements to the matrix,
			pos - the diagonal to add the elements to:
				  positive ro zero, number of elements added, 
				  negative, error.
		*/
		int addDiag(int pos, vector<T> elements){
			int i=0, j=0, ans =0;
			//assert(abs(pos) < n || abs(pos) <m);
			if (abs(pos) > n || abs(pos) >m)
				return -1;
			if (pos >0){
				j+=pos;
				}
			else if (pos <0){
				i+=pos*-1;
				}
			
			vector<T>::iterator it = elements.begin();
			if (pos>0){
				for (int k=0; k<j; k++)
					it++;

				}
			else if (pos < 0){
				for (int k=0; k<i; k++)
					elements.pop_back();
				}


			for (it; it != elements.end() && i<n && j<m; it++, ans++){
				if (*it!=0){
					mat[i++][j++] = *it;
					
					}
				else{
					i++;j++;
					
					}
				}
			//cout<<"in:"<<in<<" out:"<<out<<endl;
			return ans;
			}


		row_iter begin(){
				return mat.begin();
			}

		row_iter end(){
			return mat.end();
			}

		col_t getRow(int row){
				return mat[row];
			}

		virtual ~SparseMat(void){};
		
		//construct master matrix from uu, uv, vu, vv

		int getN(){return n;}
		int getM(){return m;}
		
	private:
		int m,n;
		mat_t mat;

	};
