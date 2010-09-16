IplImage *medianFilter(IplImage *src,int n);
void structure_texture_decomposition_rof(IplImage *in1,IplImage *in2,IplImage *texture1,IplImage *texture2,IplImage *strcuture1,IplImage *strcuture2,float theta,int nIters,float alp);
void ShowImage(char *title,IplImage *img);
void ShowImage2(char *title,IplImage *img);
void ImageResize(IplImage *src,IplImage *dst,bool interp=true);