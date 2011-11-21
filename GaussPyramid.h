#ifndef _GaussPyramid_h
#define _GaussPyramid_h

//#include "Image.h"
#include <iostream>
#include <stdio.h>
#include "highgui.h" 
#include "cv.h"
using namespace std;
class GaussPyramid
{
private:
	IplImage** ImPyramid;
	int nLevels;
public:
	GaussPyramid(void);
	~GaussPyramid(void);
	void ConstructPyramid(const IplImage* image,double ratio=0.8,int minWidth=30);
	//void displayTop(const char* filename);
	inline int getNlevels() const {return nLevels;};
	inline void SetNlevels(int levels) {nLevels=levels;};
	inline IplImage* getImageFromPyramid(int index) {return ImPyramid[index];};
};

#endif