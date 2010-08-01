//#pragma once
#include <vector>

#ifndef __VECTOR__TOOLS__KIT
#define __VECTOR__TOOLS__KIT
using namespace std;


//vector operators:
//static vector<float> * vectorSub(vector<float> * a, vector<float>* b);
template<class T>
vector<T>  operator-(vector<T> const& a, vector<T> const& b){
	vector<float> ans;
	vector<float>::const_iterator ait = a.begin();
	vector<float>::const_iterator bit = b.begin();

	while(ait != a.end() && bit != b.end()){
		ans.push_back((*ait) - (*bit));
		ait++; bit++;
		}
	return ans;
	}

//static vector<float>* vectorSub(float val, vector<float>* b);
template<class T,class T2>
vector<T>  operator-(T2 val, vector<T> const& b){
	vector<T> ans;
	for (vector<T>::const_iterator it = b.begin(); it != b.end(); it++)
		ans.push_back(val-*it);
	return ans;
	}

//static vector<float>* vectorSub(vector<float>* a, float val);
template<class T,class T2>
vector<T> operator-(vector<T> const& a, T2 val){
	vector<float> ans;
	for (vector<float>::const_iterator it = a.begin(); it != a.end(); it++)
		ans.push_back(*it-val);
	return ans;

	}

//static vector<float>* vectorMul(vector<float>* a, vector<float>* b);
template <class T>
vector<T> operator*(vector<T> const& a, vector<T> const& b){
	vector<T> ans;
	vector<T>::const_iterator ait = a.begin();
	vector<T>::const_iterator bit = b.begin();

	while(ait != a.end() && bit != b.end()){
		ans.push_back((*ait) * (*bit));
		ait++; bit++;
		}
	return ans;
	}

//static vector<float>* vectorMul(float val, vector<float>* b);
template <class T, class T2>
vector<T> operator*(T2 val, vector<T>const& b){
	vector<T> ans;
	for (vector<T>::const_iterator it = b.begin(); it != b.end(); it++)
		ans.push_back((*it) * val);
	return ans;
	}

//static vector<float>* vectorMul(vector<float>* a, float val);
template<class T, class T2>
vector<T> operator*(vector<T> const& a, T2 val){
	return val*a;
	}

//static vector<float>* vectorAdd(vector<float>* a, vector<float>* b);
template<class T>
vector<T> operator+(vector<T>const& a, vector<T>const& b){
	vector<T> ans;
	vector<T>::const_iterator ait = a.begin();
	vector<T>::const_iterator bit = b.begin();

	while(ait != a.end() && bit != b.end()){
		ans.push_back((*ait) + (*bit));
		ait++; bit++;
		}
	return ans;
	}

//static vector<float>* vectorAdd(vector<float>* a, float val);
template<class T, class T2>
vector<T> operator+(vector<T>const& a, T2 val){
	vector<T> ans;
	for (vector<T>::const_iterator it = a.begin(); it != a.end(); it++)
		ans.push_back(*it+val);
	return ans;
	}



//static vector<float>* vectorAdd(float val, vector<float>* b);
template<class T, class T2>
vector<T> operator+(T2 val, vector<T>const& b){
	return b+val;
	}


template<class T>
vector<T>  operator<<=(IplImage * I, vector<T> const& e){
	vector<T> ans;
	vector<T>* colVect = toolsKit::IplImageToCoulmnVector(I);

	for(vector<T>::const_iterator it = e.begin(); it!= e.end(); it++){
		ans.push_back((*colVect)[(*it)-1]);
		}
	delete colVect;
	return ans;
	}

template<class T>
ostream& operator <<(ostream& lhs, const vector<T>& rhs){
	for (vector<T>::const_iterator it = rhs.begin(); it != rhs.end(); it++)
		lhs<<*it<<endl;
	
	return lhs;
	}


#endif