#include "dog.hpp"
#define GKERNEL_SIZE 9
/**@bug fuck that. use GIMP or photoshop
*/
//static float gkernel[GKERNEL_SIZE]={0.5,0.9,0.12,0.15,0.16,0.15,0.12,0.9,0.5};
//static float gkernel_norm=1/0.98;
#define IX(x,y) ((x)*width+(y))
//32bit rgba image
unsigned char * dog(const unsigned char * image, int width, int height)
{
  unsigned char * output=new unsigned char[width*height*4];
  for(int ii=0;ii<height;ii++){
    for(int jj=0;jj<width;jj++){
  //    int sum[3]={0,0,0};
      for(int kk=jj-GKERNEL_SIZE/2;kk<=jj+GKERNEL_SIZE/2;kk++){

      }
    }
  }
  return output;
}
