class FlowError
{
public:
	FlowError();
	void Release();
	void ShowAE();
	void ShowEE();
	
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
	void CalculateError(eErrorType errorType,float maxMotion);
	bool LoadFlow(const string &filename);
	bool LoadGroundTruthFlow(const string &filename);
	double GetAverageAE();
	double GetAverageEP();
	double GetAverageMag();
	double GetAverageOri();
private:
	static void cvShowManyImages(char* title, int nArgs, ...);
	static bool ReadFlowFile(const string &filename,IplImage **U,IplImage **V,int depth);
	static bool WriteFlowFile(const string &filename,IplImage *U,IplImage *V);
	static bool ReadFlowFile(const string &filename,CvMat **U,CvMat **V);
	static bool WriteFlowFile(const string &filename,CvMat *U,CvMat *V);
	static IplImage *GetFlowImage(IplImage *U,IplImage *V,float maxMotion); 
	IplImage *U,*V,*gtU,*gtV;
	IplImage *aeImage;
	IplImage *eeImage;
	IplImage *orientationImage,*magnitudeImage;
	double averageAE,averageEE,averageMag,averageOri;
};