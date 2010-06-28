// flowIO.h
#ifndef FLOW_IO
#define FLOW_IO
#include "optical_flow_demo.h"
namespace  middlebury{

// read a flow file into 2-band image
void ReadFlowFile(CFloatImage& img, const char* filename);

// write a 2-band image into flow file 
void WriteFlowFile(CFloatImage img, const char* filename);

};
#endif

