#pragma once

#include <cv.h>
#include <highgui.h>

#include <stdio.h>
#include <stdarg.h>
#include "SparseToolKit.h"

class toolsKit
{
public:
	toolsKit();
	static void cvShowManyImages(char* title, int nArgs, ...);
	static void opt_flow_lk();
	virtual ~toolsKit(void);

	
};
