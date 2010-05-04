

class flowUV
{
public:

	inline flowUV(IplImage* u,IplImage* v)
	{
		innerU=u;
		innerV=v;
	}
	inline flowUV(int width,int height,int depth,int channels){
		innerU=cvCreateImage(cvSize( width, height ),depth,channels); 
		innerV=cvCreateImage(cvSize( width, height ),depth,channels); 
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
private:
	IplImage* innerU;
	IplImage* innerV;
};
