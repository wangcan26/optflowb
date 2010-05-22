#pragma once

#include <map>
#include <vector>
#include <iostream>

using namespace std;

template<class T>
class SparseMat
	{

	public:		
		typedef std::map<int, std::map<int , T> > mat_t;
		typedef typename mat_t::iterator row_iter;
		typedef typename mat_t::const_iterator const_row_iter;
		typedef std::map<int, T> col_t;
		typedef typename col_t::iterator col_iter;
		typedef typename col_t::const_iterator const_col_iter;

		SparseMat(){ m=n=1;}
		SparseMat(int i) {m=n=i;}
		SparseMat(int i, int j){m=i; n=j;};
		SparseMat (const SparseMat& c){
			mat(c);
			m=c.m;
			n= c.n;
			}
		inline
			T& operator()(int i, int j)
			{
			if(i>=m || j>=n) throw;
			return mat[i][j];
			}
		inline
			T operator()(int i, int j) const
			{
			if(i>=m || j>=n) throw;
			return mat[i][j];
			}

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
				i+=pos;
				}
			if (pos <0){
				j+=pos*-1;
				}

			for (vector<T>::iterator it = elements.begin(); it != elements.end() && i<n && j<m; it++, ans++){
				mat[i++][j++] = *it;
				}
			return ans;
			}



		virtual ~SparseMat(void){};

		friend ostream& operator <<(ostream& lhs, const SparseMat& rhs)
			{
			const_row_iter ii = rhs.mat.begin();
			const_col_iter jj;
			for (ii; ii != rhs.mat.end(); ii++){
				for (jj=ii->second.begin(); jj != ii->second.end(); jj++){
					lhs<<"["<<ii->first<<","<<jj->first<<"]"<<"="<<jj->second<<endl;
					}
				}			
			return lhs;
			}


	private:
		int m,n;
		mat_t mat;




	};
