#pragma once

#include <cv.h>
#include <highgui.h>

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
		toolsKit();
		static void cvShowManyImages(char* title, int nArgs, ...);
		
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


		static void IPL_add(IplImage* img,IplImage* img2,IplImage* dest);	
		static void IPL_print(const IplImage *image);
		static void cvMulScalar(IplImage* img,float scalar);
		static void cvNormalizeEdges(IplImage* img);
		static void cvZeroNans(IplImage* img);
		static void costumeLineCompute(IplImage* ans,IplImage* var1,IplImage* var2,IplImage* var3,IplImage* var4,IplImage* var5);
		static IplImage* psiDerivative(IplImage* x,double epsilon);

		static bool AlmostEqualRelativeOrAbsolute(float A, float B,float maxRelativeError, float maxAbsoluteError);
		static bool IsNan(float A);

		static IplImage * IplFromFile(string filename);
		static void IplToFile(IplImage* img, string filename);
	
		virtual ~toolsKit(void);
	private:
		static void IPL_operate_left(IplImage* img,IplImage* img2,IplImage* dest,operations oper);
		static void IPL_operate_right(IplImage* img,IplImage* img2,IplImage* dest,operations oper);
		static void IPL_operate_top(IplImage* img,IplImage* img2,IplImage* dest,operations oper);
		static void IPL_operate_bottom(IplImage* img,IplImage* img2,IplImage* dest,operations oper);
		

	};
