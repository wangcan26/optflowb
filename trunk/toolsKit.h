#pragma once

#include <cv.h>
#include <highgui.h>

#include <stdio.h>
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
				*it= 1/(2*sqrt((double)*it)); 
				++it;
			}
		}
		template <class PEL>
		static void IPL_mul_inverse_loop(IplImageIterator<PEL> it){
			double one=1;
			while (!it) {      
				if(0!=((double)*it))
					*it= one/((double)*it); 
				++it;
			}
		}
		static void IPL_mul_inverse(IplImage* img,int opType);
		static void IPL_add_left(IplImage* img,IplImage* img2,IplImage* dest);
		static void IPL_add_right(IplImage* img,IplImage* img2,IplImage* dest);
		static void IPL_add_top(IplImage* img,IplImage* img2,IplImage* dest);
		static void IPL_add_bottom(IplImage* img,IplImage* img2,IplImage* dest);
		static void IPL_add(IplImage* img,IplImage* img2,IplImage* dest);	
		static void IPL_print(IplImage *image);
		static void cvMulScalar(IplImage* img,double scalar);
		static void costumeLineCompute(IplImage* ans,IplImage* var1,IplImage* var2,IplImage* var3,IplImage* var4,IplImage* var5);
		static IplImage* psiDerivative(IplImage* x,double epsilon);
		static void opt_flow_lk();
		virtual ~toolsKit(void);
	private:
		static void IPL_operate_left(IplImage* img,IplImage* img2,IplImage* dest,operations oper);
		static void IPL_operate_right(IplImage* img,IplImage* img2,IplImage* dest,operations oper);
		static void IPL_operate_top(IplImage* img,IplImage* img2,IplImage* dest,operations oper);
		static void IPL_operate_bottom(IplImage* img,IplImage* img2,IplImage* dest,operations oper);		

	};
