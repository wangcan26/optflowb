

class flowUV
{
public:

	inline flowUV(IplImage* u,IplImage* v)
	{
		innerU=u;
		innerV=v;
		PsidashFSAns1=cvCreateImage(cvSize( u->width, u->height ),u->depth,u->nChannels); 
		PsidashFSAns2=cvCreateImage(cvSize( u->width, u->height ),u->depth,u->nChannels); 
	
	}
	inline flowUV(int width,int height,int depth,int channels){
		innerU=cvCreateImage(cvSize( width, height ),depth,channels); 
		innerV=cvCreateImage(cvSize( width, height ),depth,channels); 
		PsidashFSAns1=cvCreateImage(cvSize( width, height ),depth,channels); 
		PsidashFSAns2=cvCreateImage(cvSize( width, height ),depth,channels); 
	}
	inline virtual ~flowUV(void)
	{
		innerU=NULL;
		innerV=NULL;
	}
	inline IplImage* getU(){
		return innerU;
	}
	inline IplImage* getV(){
		return innerV;
	}

	inline IplImage* getPsidashFSAns1(){
		return PsidashFSAns1;
	}
	inline IplImage* getPsidashFSAns2(){
		return PsidashFSAns2;
	}
private:
	IplImage* innerU;
	IplImage* innerV;
	IplImage* PsidashFSAns1;
	IplImage* PsidashFSAns2;
};