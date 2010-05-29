#pragma once
template <class PEL>
class IplImageIterator {
  int i, j, i0;
  PEL* data;
  int step;
  int nl, nc;
  int nch;
public:

  /* constructor */
  IplImageIterator(IplImage* image, int x=0, int y=0, int dx= 0, int dy=0) : i(x), j(y), i0(0) {

    data= reinterpret_cast<PEL*>(image->imageData);
    step= image->widthStep / sizeof(PEL);
    CvRect rect= cvGetImageROI(image);
    nl= rect.height;
    nc= rect.width;
    x+= rect.x;
    y+= rect.y;
    if ((y+dy)>0 && (y+dy)<nl) nl= y+dy;
    if (y<0 || y>=nl) j=0;
    data+= step*j;
    if ((x+dx)>0 && (x+dx)<nc) nc= x+dx;
    nc*= image->nChannels;
    if (x>0 && x<nc) i0= x*image->nChannels;
    i= i0;
    nch= image->nChannels;
  }

 

  /* has next ? */
  bool operator!() const { return j < nl; }

  /* next pixel or next color component */
  IplImageIterator& operator++() {
      i++;
    if (i >= nc) {
            i=i0;
            j++;
            data+= step;
      }
    return *this;
  }

  const IplImageIterator operator++(int) {
      IplImageIterator<PEL> copy(*this);
      ++(*this);
      return copy;
  }

 

  IplImageIterator& operator+=(int s) {
      i+= s;
    if (i >= nc) {
            i=i0;
            j++;
            data+= step;
      }
    return *this;
  }

  /* pixel access */

  PEL& operator*() {
        return data[i];
  }

  const PEL operator*() const {
        return data[i];
  }

 

  const PEL neighbor(int dx, int dy) const {
        return *(data+dy*step+i+dx*nch);
  }

 

  PEL* operator&() const {
        return data+i;
  }

 

  /* current pixel coordinates */

  int column() const {
        return i/nch;
  }

 

  int line() const {
        return j;
  }

};