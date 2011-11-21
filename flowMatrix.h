#pragma once

#include "cv.h"
#include "COOrdSparseMat.h"
#include "CRSSparseMat.h"
#include "flowUV.h"
#include "FArray.h"

class flowMatrix
{
public:
	static float _ERROR_CONST;
	static void constructMatrix(const cv::Mat& Ix, const cv::Mat& Iy, const cv::Mat& Iz, flowUV* UV, const cv::Mat& du, const cv::Mat& dv, CRSSparseMat* A, const FArray& B, const float alpha);

private:
	static cv::Mat* constructSmoothnessTerm(const int rows, const int cols);
	static void bulidA(CRSSparseMat* const A, const int cols, const int rows, const float* const duu, const float* const dvv, const float* const duv, const float* const minorDiag, const float* const outerDiag);
	static void buildB(float* B, const int rows, const int cols, const int type, const cv::Mat& fs1_3222, const cv::Mat& fs1_122ht22, const cv::Mat& fs2_2232, const cv::Mat& fs2_22122wt, const cv::Mat& dataTermNorm, const cv::Mat& psi2DataTerm, const cv::Mat& Ix, const cv::Mat& Iy, const cv::Mat& Iz, flowUV* UV);
	//static void buildBx(const cv::Mat& flow, const cv::Mat& diag, const cv::Mat& minor, const cv::Mat& outer, const cv::Mat& It, cv::Mat& dest, const int rows, const int cols, const int type);
};

