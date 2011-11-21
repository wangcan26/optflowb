#include "flowMatrix.h"
#include <ctime>
#define MISTAKE 1
#define VERBOSE 0

float flowMatrix::_ERROR_CONST;

void flowMatrix::bulidA (CRSSparseMat* const A,
						const int	cols,
						const int	rows,
						const float* const		duu, 
						const float* const 		dvv,
						const float* const		duv,
						const float* const 		minorDiag,
						const float* const 		outerDiag)
{
#if VERBOSE
	std::clock_t start = std::clock();
#endif

	int majorDiagSize = rows * cols;
	int minorDiagSize = (rows * cols) - 1;
	int outerDiagSize = (rows - 1) * cols;
	int nzAestimate = 4 * (majorDiagSize + minorDiagSize + outerDiagSize);

	float* cooVal = new float[nzAestimate];
	float* cooValPtr = cooVal;

	memcpy((cooValPtr),						outerDiag,	outerDiagSize * sizeof(float));
	memcpy((cooValPtr += outerDiagSize),	minorDiag,	minorDiagSize * sizeof(float));
	memcpy((cooValPtr += minorDiagSize),	duu,		majorDiagSize * sizeof(float));
	memcpy((cooValPtr += majorDiagSize),	minorDiag,	minorDiagSize * sizeof(float));
	memcpy((cooValPtr += minorDiagSize),	outerDiag,	outerDiagSize * sizeof(float));
	memcpy((cooValPtr += outerDiagSize),	duv,		majorDiagSize * sizeof(float));

	memcpy((cooValPtr += majorDiagSize),	duv,		majorDiagSize * sizeof(float));
	memcpy((cooValPtr += majorDiagSize),	outerDiag,	outerDiagSize * sizeof(float));
	memcpy((cooValPtr += outerDiagSize),	minorDiag,	minorDiagSize * sizeof(float));
	memcpy((cooValPtr += minorDiagSize),	dvv,		majorDiagSize * sizeof(float));
	memcpy((cooValPtr += majorDiagSize),	minorDiag,	minorDiagSize * sizeof(float));
	memcpy((cooValPtr += minorDiagSize),	outerDiag,	outerDiagSize * sizeof(float));

	int* cooRow = new int[nzAestimate];
	int* cooCol = new int[nzAestimate];
	int* cooRowPtr = cooRow;
	int* cooColPtr= cooCol;

	//TODO:: //BOAZ: maybe the count forvar (and the other redundent forvars) forvar has a relative big overhead.. needs to be benchmarked
	for (int c = 0,				r = cols;	c < outerDiagSize;	*(cooRowPtr++) = r, *(cooColPtr++) = c, ++c, ++r);	//Outer
	for (int c = 0,				r = 1;		c < minorDiagSize;	*(cooRowPtr++) = r, *(cooColPtr++) = c, ++c, ++r);	//Inner
	for (int c = 0,				r = 0;		c < majorDiagSize;	*(cooRowPtr++) = r, *(cooColPtr++) = c, ++c, ++r);	//duu
	for (int c = 1,				r = 0;		r < minorDiagSize;	*(cooRowPtr++) = r, *(cooColPtr++) = c, ++c, ++r);	//Inner
	for (int c = cols,			r = 0;		r < outerDiagSize;	*(cooRowPtr++) = r, *(cooColPtr++) = c, ++c, ++r);	//Outer
	for (int c = majorDiagSize,	r = 0;		r < majorDiagSize;	*(cooRowPtr++) = r, *(cooColPtr++) = c, ++c, ++r);	//duv
	
	for (int c = 0,						r = majorDiagSize;			c < majorDiagSize; *(cooRowPtr++) = r,	*(cooColPtr++) = c, ++c, ++r);							//duv
	for (int c = majorDiagSize,			r = majorDiagSize + cols,	count = 0;	count < outerDiagSize;	*(cooRowPtr++) = r, *(cooColPtr++) = c, ++c, ++r, ++count);	//Outer
	for (int c = majorDiagSize,			r = majorDiagSize + 1,		count = 0;	count < minorDiagSize;	*(cooRowPtr++) = r, *(cooColPtr++) = c, ++c, ++r, ++count);	//Inner
	for (int c = majorDiagSize,			r = majorDiagSize,			count = 0;	count < majorDiagSize;	*(cooRowPtr++) = r, *(cooColPtr++) = c, ++c, ++r, ++count);	//dvv
	for (int c = majorDiagSize + 1,		r = majorDiagSize,			count = 0;	count < minorDiagSize;	*(cooRowPtr++) = r, *(cooColPtr++) = c, ++c, ++r, ++count);	//Inner
	for (int c = majorDiagSize + cols,	r = majorDiagSize,			count = 0;	count < outerDiagSize;	*(cooRowPtr++) = r, *(cooColPtr++) = c, ++c, ++r, ++count);	//Outer

	COOrdSparseMat* coo = new COOrdSparseMat(2 * majorDiagSize, 2 * majorDiagSize, nzAestimate, cooVal, cooRow, cooCol);
#if VERBOSE
	std::clock_t crsTimer = std::clock();
#endif
	A->build(*coo);
#if VERBOSE
	cout << "A Optimized Build  " << " (CRS build took: " << std::clock() - crsTimer << ") ";
#endif
	delete coo;
#if VERBOSE
	cout << std::clock() - start << endl;
#endif
}

void buildM (CRSSparseMat* const A,
						const int	cols,
						const int	rows,
						const float* const		diag,
						const float* const 		minorDiag,
						const float* const 		outerDiag)
{
#if VERBOSE
	std::clock_t start = std::clock();
#endif

	int majorDiagSize = rows * cols;
	int minorDiagSize = (rows * cols) - 1;
	int outerDiagSize = (rows - 1) * cols;
	int nzAestimate = 2 * majorDiagSize + 4 * minorDiagSize + 4 * outerDiagSize;

	float* cooVal = new float[nzAestimate];
	float* cooValPtr = cooVal;

	memcpy((cooValPtr),						outerDiag,	outerDiagSize * sizeof(float));
	memcpy((cooValPtr += outerDiagSize),	minorDiag,	minorDiagSize * sizeof(float));
	memcpy((cooValPtr += minorDiagSize),	diag,		majorDiagSize * sizeof(float));
	memcpy((cooValPtr += majorDiagSize),	minorDiag,	minorDiagSize * sizeof(float));
	memcpy((cooValPtr += minorDiagSize),	outerDiag,	outerDiagSize * sizeof(float));

	memcpy((cooValPtr += outerDiagSize),	outerDiag,	outerDiagSize * sizeof(float));
	memcpy((cooValPtr += outerDiagSize),	minorDiag,	minorDiagSize * sizeof(float));
	memcpy((cooValPtr += minorDiagSize),	diag,		majorDiagSize * sizeof(float));
	memcpy((cooValPtr += majorDiagSize),	minorDiag,	minorDiagSize * sizeof(float));
	memcpy((cooValPtr += minorDiagSize),	outerDiag,	outerDiagSize * sizeof(float));

	int* cooRow = new int[nzAestimate];
	int* cooCol = new int[nzAestimate];
	int* cooRowPtr = cooRow;
	int* cooColPtr= cooCol;

	//TODO:: //BOAZ: maybe the count forvar (and the other redundent forvars) forvar has a relative big overhead.. needs to be benchmarked
	for (int c = 0,				r = cols;	c < outerDiagSize;	*(cooRowPtr++) = r, *(cooColPtr++) = c, ++c, ++r);	//Outer
	for (int c = 0,				r = 1;		c < minorDiagSize;	*(cooRowPtr++) = r, *(cooColPtr++) = c, ++c, ++r);	//Inner
	for (int c = 0,				r = 0;		c < majorDiagSize;	*(cooRowPtr++) = r, *(cooColPtr++) = c, ++c, ++r);	//duu
	for (int c = 1,				r = 0;		r < minorDiagSize;	*(cooRowPtr++) = r, *(cooColPtr++) = c, ++c, ++r);	//Inner
	for (int c = cols,			r = 0;		r < outerDiagSize;	*(cooRowPtr++) = r, *(cooColPtr++) = c, ++c, ++r);	//Outer
	
	for (int c = majorDiagSize,			r = majorDiagSize + cols,	count = 0;	count < outerDiagSize;	*(cooRowPtr++) = r, *(cooColPtr++) = c, ++c, ++r, ++count);	//Outer
	for (int c = majorDiagSize,			r = majorDiagSize + 1,		count = 0;	count < minorDiagSize;	*(cooRowPtr++) = r, *(cooColPtr++) = c, ++c, ++r, ++count);	//Inner
	for (int c = majorDiagSize,			r = majorDiagSize,			count = 0;	count < majorDiagSize;	*(cooRowPtr++) = r, *(cooColPtr++) = c, ++c, ++r, ++count);	//dvv
	for (int c = majorDiagSize + 1,		r = majorDiagSize,			count = 0;	count < minorDiagSize;	*(cooRowPtr++) = r, *(cooColPtr++) = c, ++c, ++r, ++count);	//Inner
	for (int c = majorDiagSize + cols,	r = majorDiagSize,			count = 0;	count < outerDiagSize;	*(cooRowPtr++) = r, *(cooColPtr++) = c, ++c, ++r, ++count);	//Outer

	COOrdSparseMat* coo = new COOrdSparseMat(2 * majorDiagSize, 2 * majorDiagSize, nzAestimate, cooVal, cooRow, cooCol);
#if VERBOSE
	std::clock_t crsTimer = std::clock();
#endif
	A->build(*coo);
#if VERBOSE
	cout << "A Optimized Build  " << " (CRS build took: " << std::clock() - crsTimer << ") ";
#endif
	delete coo;
#if VERBOSE
	cout << std::clock() - start << endl;
#endif
}



//void flowMatrix::buildBx(	const cv::Mat& flow,
//							const cv::Mat& diag,
//							const cv::Mat& minor,
//							const cv::Mat& outer,
//							const cv::Mat& It,
//							cv::Mat& dest,
//							const int rows,
//							const int cols,
//							const int type)
//{
//	cv::Mat minorFlow(rows, cols, type);
//	cv::multiply(flow, minor, minorFlow);
//	cv::Mat outerFlow(rows, cols, type);
//	cv::multiply(flow, outer, outerFlow);
//	cv::Mat diagFlow(rows, cols, type);
//	cv::multiply(flow, diag, diagFlow);
//
//
//	cv::Mat UmD(rows, cols, type);
//	memcpy(UmD.data, minorFlow.data + 1, cols * rows - 1);
//	*(UmD.data + (cols * rows - 1)) = 0;
//	
//	cv::Mat LmD(rows, cols, type);
//	*LmD.data = 0;
//	memcpy(LmD.data + 1, minorFlow.data, cols * rows - 1);
//
//	cv::Mat UoD(rows, cols, type);
//	memcpy(UoD.data, outerFlow.data + rows, cols * (rows - 1));
//	for(int i = cols * (rows - 1); i < cols * rows; *(UoD.data + i) = 0, ++i);
//
//	cv::Mat LoD(rows, cols, type);
//	for(int i = 0; i < rows; *(LoD.data + i) = 0, ++i);
//	memcpy(LoD.data + rows, outerFlow.data, cols * (rows - 1));
//
//	dest = diagFlow + UmD + LmD + UoD + LoD;
//	//dest *= alpah;
//	dest -= It;
//}
//


#define MAIN_DIAG 0
#define MINOR_DIAG 1
#define OUTER_DIAG 2
//Smoothness Term via the Laplacian operator
cv::Mat* flowMatrix::constructSmoothnessTerm(const int rows, const int cols)
{	
	cv::Mat* retVal = new cv::Mat[3];
	retVal[MAIN_DIAG].create(rows, cols, CV_32FC1);
	/*
	*	A matrix of the form
		-2 -3 -3 -3 -2
		-3 -4 -4 -4 -3
		-3 -4 -4 -4 -3
		-3 -4 -4 -4 -3
		-2 -3 -3 -3 -2
	*/
	float* rowOuter = new float[cols];
	*rowOuter = -2;
	float* pOuter = rowOuter + 1;
	for (int i = 1; i < cols - 1; ++i, (*(pOuter++) = -3));
	*pOuter = -2;
	float* rowInner = new float[cols];
	*rowInner = -3;
	float* pInner = rowInner + 1;
	for (int i = 1; i < cols - 1; ++i, (*(pInner++) = -4));
	*pInner = -3;
	memcpy(retVal[MAIN_DIAG].data, rowOuter, cols * sizeof(float));
	for (int i = 1; i < rows - 1; ++i)
		memcpy(retVal[MAIN_DIAG].ptr<float>(i), rowInner, cols * sizeof(float));
	memcpy(retVal[MAIN_DIAG].ptr<float>(rows - 1), rowOuter, cols * sizeof(float));


	float* p = rowOuter;
	for (int i = 0; i < cols - 1; ++i, (*(p++) = 1));
	*p = 0;
	retVal[MINOR_DIAG].create(rows, cols, CV_32FC1);
	/*
	*	A matrix of the form
		1 1 1 1 0 
		1 1 1 1 0 
		1 1 1 1 0 
		1 1 1 1 0 
		1 1 1 1 0 
	*/
	for (int i = 0; i < rows; ++i){
		memcpy(retVal[MINOR_DIAG].ptr<float>(i), rowOuter, cols * sizeof(float));
	}

	retVal[OUTER_DIAG].create(rows, cols, CV_32FC1);
	retVal[OUTER_DIAG].setTo(cv::Scalar(1));

	/*
	*	A matrix of the form
		1 1 1 1 1 
		1 1 1 1 1 
		1 1 1 1 1 
		1 1 1 1 1 
		1 1 1 1 1 
	*/

	return retVal;
}

void flowMatrix::constructMatrix(	const cv::Mat& Ix,
									const cv::Mat& Iy, 
									const cv::Mat& It, 
									flowUV* UV,
									const cv::Mat& du, 
									const cv::Mat& dv, 
									CRSSparseMat* A,
									const FArray& B,
									const float alpha)
{
	
	/*Ax=b
	*  b = [u11 u12 u13 ... u21 u22 .... v11 v12 .....]^t
	*/


	
	/*DATA TERM
	* this is the data sparse matrix as seen in the presentation slide 15, its split into 3 matrices
	* they should be positioned as [duu duv; duv dvv]. to reduce the costs of sparse arthematic operations
	* they are saved as seperate matrices.
	*/
	cv::Mat duu(Ix.rows, Ix.cols, Ix.type());
	cv::multiply(Ix, Ix, duu);

	cv::Mat dvv(Ix.rows, Ix.cols, Ix.type());
	cv::multiply(Iy, Iy, dvv);

	cv::Mat duv(Ix.rows, Ix.cols, Ix.type());
	cv::multiply(Iy, Ix, duv);

	/*SMOOTHNESS TERM
	* this is the smoothness sparse matrix as seen in the presentation slide 15, its split into 3 matrices
	* the main diagonal, the 2 minor diagonals (they are the same thus just use the minor diagonal matrix twice)
	* and the outer diagonal (again, they are the same, so use it twice).
	* togather they form the sparse Lapalcian based smoothness term (well, its only the upper left quadrent.. but again, they are the same)
	*/
	cv::Mat* sTerm = constructSmoothnessTerm(Ix.rows, Ix.cols);

	/* prepare terms to be merged into A via the following calculation: A = DataTerm - alpha * SmoothnessTerm;
	* since the DataTerm sparse matrix consists of quadratic diagonals alone, and the SmoothnessTerm is a band sparse matrix
	* all we need to do is to solve the substrution is to do: 
	* AdiagUpper =  duu - alpha* sTerm[MAIN_DIAG]
	* AdiagLower =  dvv - alpha* sTerm[MAIN_DIAG]
	* and construct the concreat sparse matrix A using the following: AdiagUpper, AdiagLower, duv, -sTerm[MINOR_DIAG] 1th, 
	*														-sTerm[MINOR_DIAG] 2nd, -sTerm[OUTER_DIAG] 1th, -sTerm[OUTER_DIAG] 2nd
	*/
	cv::Mat diag = sTerm[MAIN_DIAG] * alpha;
	//sTerm[MAIN_DIAG] *= alpha;

	duu -= diag;

	dvv -= diag;


	//sTerm[MINOR_DIAG] *= alpha;
	//sTerm[OUTER_DIAG] *= alpha;

	cv::Mat min = sTerm[MINOR_DIAG] * -alpha;
	cv::Mat out = sTerm[OUTER_DIAG] * -alpha;

	bulidA(A, Ix.cols, Ix.rows, (const float* const)duu.data, (const float* const)dvv.data, (const float* const)duv.data, (const float* const)min.data, (const float* const)out.data);
	
	//B
	cv::Mat Itx(Ix.rows, Ix.cols, Ix.type());
	cv::multiply(It, Ix, Itx);

	cv::Mat Ity(Ix.rows, Ix.cols, Ix.type());
	cv::multiply(It, Iy, Ity);


	cv::Mat Bu(Itx.cols, Itx.rows, Itx.type());
	cv::Mat Bv(Ity.cols, Ity.rows, Ity.type());


	/* this is rather a complicated solution for a simple problem but it saves a lot of overhead from working with sparse matrices
	* The SmoothnessTerm M is split into matrices - M = D + UmD + LmD + UoD + Lod (U - upper, L - lower, m - minor, o - outer, D - diagonal)
	* so M*uv =  (D + UmD + LmD + UoD + Lod) * uv. since all those matrices are diagonal (not neccesrily the main diagonal) we can use 
	* cell by cell multipication to solve it.
	*/

	cv::Mat U = UV->getU();
	cv::Mat V = UV->getV();
	//cv::imshow("abc", Iy);
	//cv::waitKey(0);

	//cv::Mat minor(Ix.rows, Ix.cols, Ix.type());
	//cv::multiply(U, sTerm[MINOR_DIAG], minor);
	//cv::Mat outer = U * sTerm[OUTER_DIAG];  uneccery, outer diagonal is all ones

	//buildBx(U, sTerm[MAIN_DIAG], sTerm[MINOR_DIAG], sTerm[OUTER_DIAG], Itx, Bu, Ix.rows, Ix.cols, Ix.type());
	//buildBx(V, sTerm[MAIN_DIAG], sTerm[MINOR_DIAG], sTerm[OUTER_DIAG], Ity, Bv, Ix.rows, Ix.cols, Ix.type());
	CRSSparseMat* M = new CRSSparseMat();

	cv::Mat diag2 = sTerm[MAIN_DIAG] * alpha;
	cv::Mat min2 = sTerm[MINOR_DIAG] * alpha;
	cv::Mat out2 = sTerm[OUTER_DIAG] * alpha;

	buildM(M, U.cols, U.rows, (const float*)diag2.data, (const float*)min2.data, (const float*)out2.data);
	//ofstream f;
	//f.open ("newM.txt");
	//for(int i = 0; i < M->rows(); ++i){
	//	for(int r = M->rowPtr(i); r < M->rowPtr(i + 1); ++r){
	//		f << i << ", " << M->colIdx(r) << ": " << M->val(r) << endl;
	//	}
	//}
	//f.close();
	//M->MulScalar<float>(0.3);
	float* X = new float[B.size()];
	memcpy(X,					 U.data, U.cols * U.rows * sizeof(float));
	memcpy(X + U.cols * U.rows , V.data, V.cols * V.rows * sizeof(float));

	M->MulVector(X,B.ptr);

	for (int i = 0; i < Itx.cols * Itx.rows; ++i){
		B.ptr[i] -= ((float*)Itx.data)[i];
	}
	for (int i = 0; i < Ity.cols * Ity.rows; ++i){
		B.ptr[Itx.cols * Itx.rows + i] -= ((float*)Ity.data)[i];
	}


	//float* bPtr = B.ptr;
	//memcpy(bPtr, Itx/*Bu*/.data, Bu.cols * Bu.rows * sizeof(float));
	//bPtr += Bu.cols * Bu.rows;
	//memcpy(bPtr, Ity/*Bv*/.data, Bv.cols * Bv.rows * sizeof(float));
	//ofstream f2;
	//f2.open ("newB.txt");
	//for(int i = 0; i < B.size; ++i){
	//	f2 << i << ": " << B.ptr[i] << endl;
	//}
	//f2.close();

	//int rows = Ix.rows;

	//cv::Mat UmD(U.size, U.type());
	//memcpy(UmD.data, minor.data + 1, minor.cols * minor.rows - 1);
	//*(Umd.data + (minor.cols * minor.rows - 1)) = 0;
	//
	//cv::Mat LmD(U.size, U.type());
	//*Lmd.data = 0;
	//memcpy(Lmd.data + 1, minor.data, minor.cols * minor.rows - 1);

	//cv::Mat UoD(U.size, U.type());
	//memcpy(UoD.data, minor.data + rows, minor.cols * minor.rows - rows);
	//for(int i = U.cols * U.rows - rows; i < U.cols * U.rows; ++i, *(UoD.data + i) = 0);

	//cv::Mat LoD(U.size, U.type());
	//for(int i = 0; i < rows; ++i, *(LoD.data + i) = 0);
	//memcpy(LoD.data + rows, U.data, U.cols * U.rows - rows);

	//Bu = sTerm[MAIN_DIAG] + UmD + LmD + UoD + LoD;
	////Bu *= alpah;
	//Bu -= Itx;

}