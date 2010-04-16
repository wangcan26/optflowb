#pragma once

#include <cv.h>
#include <highgui.h>

#include <stdio.h>
#include <stdarg.h>

class toolsKit
{
public:
	toolsKit(void);
	static void cvShowManyImages(char* title, int nArgs, ...);

	
public:
	virtual ~toolsKit(void);
};
