

class GaussPyramid
{
private:
	IplImage** ImPyramid;
	int nLevels;
public:
	GaussPyramid(void);
	~GaussPyramid(void);
	void ConstructPyramid(const IplImage& image,double ratio=0.8,int minWidth=30);
	//void displayTop(const char* filename);
	inline int nlevels() const {return nLevels;};
	inline IplImage* getImageFromPyramid(int index) {return ImPyramid[index];};
};

#endif



