#pragma once

#include <cv.h>
#include <highgui.h>
#include <stdarg.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdarg.h>
#include "SparseToolKit.h"
#include "IplImageIterator.h"

class toolsKit
	{
	public:
		
		enum operations {ADD=1,SUB=2,MUL=3};
		enum directions {UP=1,DOWN=2,LEFT=3,RIGHT=4};
		toolsKit();
		

		template <class PEL>
		static void IPLsqrt_mul2(IplImageIterator<PEL> it){
			while (!it) {
				float temp=(float)*it;

				if(0!=temp){
					float val=1/(2*sqrt(temp));
					*it=val;
					}

				++it;
				}
			}
		template <class PEL>
		static void IPL_mul_inverse_loop(IplImageIterator<PEL> it){
			double one=1;
			while (!it) {      
				if(0!=((float)*it))
					*it= one/((float)*it); 
				++it;
				}
			}
		static vector<float> * IplImageToCoulmnVector(IplImage* img);
		static void ColumnVectorToIplImage(vector<float>* vCol, IplImage* image);

		static void IPL_mul_inverse(IplImage* img,int opType);
		//adders
		static void IPL_add_left(IplImage* img,IplImage* shiftImg2,IplImage* dest);
		static void IPL_add_right(IplImage* img,IplImage* shiftImg2,IplImage* dest);
		static void IPL_add_top(IplImage* img,IplImage* shiftImg2,IplImage* dest);
		static void IPL_add_bottom(IplImage* img,IplImage* shiftImg2,IplImage* dest);
		//subs
		static void IPL_sub_left(IplImage* img,IplImage* shiftImg2,IplImage* dest);
		static void IPL_sub_right(IplImage* img,IplImage* shiftImg2,IplImage* dest);
		static void IPL_sub_top(IplImage* img,IplImage* shiftImg2,IplImage* dest);
		static void IPL_sub_bottom(IplImage* img,IplImage* shiftImg2,IplImage* dest);
		//multipliers
		static void IPL_mul_left(IplImage* img,IplImage* shiftImg2,IplImage* dest);
		static void IPL_mul_right(IplImage* img,IplImage* shiftImg2,IplImage* dest);
		static void IPL_mul_top(IplImage* img,IplImage* shiftImg2,IplImage* dest);
		static void IPL_mul_bottom(IplImage* img,IplImage* shiftImg2,IplImage* dest);

		static void IPL_mul_different_sizes(IplImage* src,IplImage* srcHorizontal,IplImage* dest);
		static void IPL_mul_different_sizes2(IplImage* src,IplImage* srcVertical,IplImage* dest);
		static void IPL_add_different_sizes(IplImage* imgHorizonal,IplImage* imgVertical,IplImage* dest);
		static void IPL_add_different_sizes2(IplImage* imgVertical,IplImage* imgVertical2,IplImage* dest);
		static void IPL_add_different_sizes3(IplImage* imgHorizonal,IplImage* imgHorizonal2,IplImage* dest);
		static void IPL_add(IplImage* img,IplImage* img2,IplImage* dest);	
		static void IPL_sub(IplImage* img,IplImage* img2,IplImage* dest);	
		static void IPL_print(const IplImage *image);
		static void PrintMat(CvMat *A);
		static void cvMulScalar(IplImage* img,float scalar);
		static void cvZeroBottom(IplImage* img);
		static void cvZeroTop(IplImage* img);
		static void cvZeroLeftRight(IplImage* img);
		static void cvZeroBottomLeft(IplImage* img);
		static void cvZeroRight(IplImage* img);
		static void cvNormalizeEdges(IplImage* img);
		static void seperateDuDv(IplImage* du,IplImage* dv,vector<float> * dUdV);
		static void cvZeroNans(IplImage* img);

		static void costumeLineCompute(IplImage* ans,IplImage* var1,IplImage* var2,IplImage* var3,IplImage* var4,IplImage* var5);
		static IplImage* psiDerivative(IplImage* x,double epsilon);
		static void increaseImageSize(IplImage* src,IplImage* dst,int select);
		static bool AlmostEqualRelativeOrAbsolute(float A, float B,float maxRelativeError, float maxAbsoluteError);
		static bool IsNan(float A);
		
		static IplImage*  toolsKit::transposeImage(IplImage* image) ;
		static IplImage*  toolsKit::transposeImage2(IplImage* image) ;
		static void shiftImage(IplImage* src,IplImage* temp,int select);

		static IplImage * IplFromFile(string filename);
		static void IplToFile(IplImage* img, string filename);


		class vectorTools{
		public:
			static void vectorMin(vector<float>* v, float val);

			static void vectorMax(vector<float>* v, float val);

			static vector<float> * vectorFloor(vector<float> * v);

			static vector<float> * vectorCeil(vector<float> * v);

			static vector<float> * vectorSub(vector<float> * a, vector<float>* b);

			static vector<float>* vectorSub(float val, vector<float>* b);

			static vector<float>* vectorSub(vector<float>* a, float val);

			static vector<float>* vectorMul(vector<float>* a, vector<float>* b);

			static vector<float>* vectorMul(float val, vector<float>* b);

			static vector<float>* vectorMul(vector<float>* a, float val);
		
			static vector<float>* vectorAdd(vector<float>* a, vector<float>* b);
			
			static vector<float>* vectorAdd(int count, ...);

			static vector<float>* vectorAdd(vector<float>* a, float val);

			static vector<float>* vectorAdd(float val, vector<float>* b);

			static vector<float>* elementsFromIpl(IplImage* I, vector<float> * p);

			static void vectorToFile(vector<float>* v, string filename);
			
			static vector<float> * elementsForIpl(vector<float> * A, int B, vector<float> * C);
			};
		
		
		
		virtual ~toolsKit(void);
	private:
		static void IPL_operate_left(IplImage* img,IplImage* img2,IplImage* dest,operations oper);
		static void IPL_operate_right(IplImage* img,IplImage* img2,IplImage* dest,operations oper);
		static void IPL_operate_top(IplImage* img,IplImage* img2,IplImage* dest,operations oper);
		static void IPL_operate_bottom(IplImage* img,IplImage* img2,IplImage* dest,operations oper);


	};
