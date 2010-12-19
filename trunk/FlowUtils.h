struct CError : public exception
{
	CError(const char* msg)                 { strcpy(message, msg); }
	CError(const char* fmt, int d)          { sprintf(message, fmt, d); }
	CError(const char* fmt, float f)        { sprintf(message, fmt, f); }
	CError(const char* fmt, const char *s)  { sprintf(message, fmt, s); }
	CError(const char* fmt, const char *s,
		int d)                          { sprintf(message, fmt, s, d); }
	char message[1024];         // longest allowable message
};

class cFlowUtils
{
public:

	enum eErrorType
	{
		errorAE,
		errorEP,
		errorMagAndOri
	};

	cFlowUtils();
	cFlowUtils(IplImage *U,IplImage *V,IplImage *gtU,IplImage *gtV);
	~cFlowUtils();
	void SetFlow(IplImage *U,IplImage *V);
	void SetGroundTruthFlow(IplImage *U,IplImage *V);
	void SetFlow(CvMat *U,CvMat *V);
	void DisplayGroundTruthFlow(float maxMotion);
	void SetGroundTruthFlow(CvMat *U,CvMat *V);
	void DisplayFlow(float maxMotion);
	void DisplayFlowError(eErrorType errorType);
	void DisplayFlowAndError(eErrorType errorType,float maxMotion);
	void CalculateError(eErrorType errorType,float minValue = 0,float maxValue = 0,float minValue2 = 0,float maxValue2 = 0);
	bool LoadFlow(const string &filename);
	bool LoadGroundTruthFlow(const string &filename);
	double GetAverageAE();
	double GetAverageEP();
	double GetAverageMag();
	double GetAverageOri();
	void CalculateMask();
	void Test();
	static const float minValueAngleError,maxValueAngleError,minValueEndpointError,maxValueEndpointError,minValueOriError,maxValueOriError,minValueMagError,maxValueMagError;
	static bool ReadFlowFile(const string &filename,IplImage **U,IplImage **V,int depth);
	static bool WriteFlowFile(const string &filename,IplImage *U,IplImage *V);
	static bool ReadFlowFile(const string &filename,CvMat **U,CvMat **V);
	static bool WriteFlowFile(const string &filename,CvMat *U,CvMat *V);
	static IplImage *GetFlowImage(IplImage *U,IplImage *V,float maxMotion=0.001f); 
	static void cvShowManyImages(char* title, int nArgs, ...);
	static void DrawFlow(IplImage* U,IplImage* V);
	static void DrawFlow(CvMat* U,CvMat* V);
	static void DrawFlow2(IplImage* du,IplImage* u,IplImage* dv,IplImage* v);
private:	

	IplImage *U,*V,*gtU,*gtV;
	char *mask;
	unsigned int validPixels;
	IplImage *aeImage;
	IplImage *eeImage;
	IplImage *orientationImage,*magnitudeImage;
	double averageAE,averageEE,averageMag,averageOri;
};