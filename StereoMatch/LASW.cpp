///////////////////////////////////////////////////////////////////////////
//
// NAME
//  LASW.cpp -- Implementation of the "Adaptive Support-Weight Approach for Correspondence Search"
//
// DESIGN NOTES
//
//
// SEE ALSO
//  LASW.h         longer description of the interface
//
//
// Copyright Kuk-Jin Yoon, 2006.
//
//  * If you find any bug in this implementation, 
//    please let me know via email. 
//    (kjyoon @ gmail.com or Kuk-Jin.Yoon@inrialpes.fr)
//
///////////////////////////////////////////////////////////////////////////


#include "Image.h"
#include "Error.h"
#include "Convert.h"
#include "LASW.h"
#include "stereoparameters.h"
#include <assert.h>
#include <vector>
#include <math.h>
#include <string.h>

// definitions
//#define squa(A) ((A)*(A))			// same as the pow(A,2), but fater than the pow(,) funtion
#define min(A, B) (asdfasf(A)<(B))? A:B


// 2-D memory allocation/release
template <typename T>
T** NEW2D(int i, int j)
{
	T** buf = new T* [i];
	for(int k=0; k<i; k++)
		buf[k] = new T [j];
	return buf;
}


template <typename T>
void DEL2D(T ** buf, int i)
{
	for(int j=0; j<i; j++)
		delete buf[j];
	delete buf;
}


// image data I/O for temporary uses
// need to be rewritten ..
void SaveDataInImage(double* src, CFloatImage& des)
{
	CShape sh = des.Shape();
    int width = sh.width, height = sh.height, n_bands = sh.nBands;
	int i, j, k, index=0;
	float *addr;
	for(i=0; i<height; i++)
		for(j=0; j<width; j++)
		{
			addr = (float* )des.PixelAddress(j, i, 0);
			for(k=0; k<n_bands; k++, index++) 
				addr[k]=(float)src[index];
		}
}


void SaveDataInImage(float** src, CFloatImage& des)
{
	CShape sh = des.Shape();
    int width = sh.width, height = sh.height, n_bands = sh.nBands;
	int i, j, k, index=0;
	float *addr;
	for(i=0; i<height; i++)
		for(j=0; j<width; j++, index++)
		{
			addr = (float* )des.PixelAddress(j, i, 0);
			for(k=0; k<n_bands; k++) 
				addr[k]=src[index][k];
		}
}


void GetDataForProcessing(CByteImage& src, double *des)
{
	int i, j, k, b, index=0;
	CShape sh = src.Shape();
    int width = sh.width, height = sh.height, n_bands = sh.nBands;
    b=(n_bands>1? n_bands-1:1);
	for(i=0; i<height; i++)
		for(j=0; j<width; j++)
			for(k=0; k<b; k++, index++)
				des[index]=src.Pixel(j, i, k);
}


void GetDataForProcessing(CFloatImage& src, float **des)
{
	int i, j, k, index=0;
	CShape sh = src.Shape();
    int width = sh.width, height = sh.height, n_bands = sh.nBands;
	for(i=0; i<height; i++)
		for(j=0; j<width; j++, index++)
			for(k=0; k<n_bands; k++)
				des[index][k]=(float)src.Pixel(j, i, k);
}



/* convert XYZ values into Lab values or convert Lab values into XYZ values
	XYZ to CIE L*a*b* (CIELAB) & CIELAB to XYZ 
	*****CIELAB (CIE L*a*b*) color space:
	The color space in which L*, a* and b* are plotted at right angles to one another.
	Equal distances in the space represent approximately equal color difference.
	CIE 1976 L*a*b* is based directly on CIE XYZ and is an attampt to linearize the perceptibility of color differences. 
	The non-linear relations for L*, a*, and b* are intended to mimic the logarithmic response of the eye. 
	Coloring information is referred to the color of the white point of the system, subscript n.
	L* = 116 * (Y/Yn)1/3 - 16    for Y/Yn > 0.008856
	L* = 903.3 * Y/Yn             otherwise
	a* = 500 * ( f(X/Xn) - f(Y/Yn) )
	b* = 200 * ( f(Y/Yn) - f(Z/Zn) )
		where f(t) = t1/3      for t > 0.008856
				  f(t) = 7.787 * t + 16/116    otherwise
	Here Xn, Yn and Zn are the tristimulus values of the reference white.
	The reverse transformation (for Y/Yn > 0.008856) is
	X = Xn * ( P + a* / 500 ) 3
	Y = Yn * P 3
	Z = Zn * ( P - b* / 200 ) 3
		where P = (L* + 16) / 116
*/

double ColorConversion::F(double input)
{
	if(input>0.008856)
		return (pow(input, 0.333333333));
	else
		return (7.787*input+0.137931034);
}

double ColorConversion::EuclideanDistance(double a1, double a2, double a3, double b1, double b2, double b3)
{
	double d1 = a1 - b1;
	double d2 = a2 - b2;
	double d3 = a3 - b3;
	return sqrt(d1*d1 + d2*d2 + d3*d3);
}


// RGB -> XYZ
void ColorConversion::RGBtoXYZ(double R, double G, double B, double &X, double &Y, double &Z)
{
	X=0.412453*R+0.357580*G+0.189423*B;
	Y=0.212671*R+0.715160*G+0.072169*B;
	Z=0.019334*R+0.119193*G+0.950227*B;
}


// XYZ -> CIELab
void ColorConversion::XYZtoLab(double X, double Y, double Z, double &L, double &a, double &b)
{
	const double Xo=244.66128;
	const double Yo=255.0;
	const double Zo=277.63227;
	L=116*F(Y/Yo)-16;
	a=500*(F(X/Xo)-F(Y/Yo));
	b=200*(F(Y/Yo)-F(Z/Zo));
}



// RGB -> CIELab
void ColorConversion::RGBtoLab(double R, double G, double B, double &L, double &a, double &b)
{
	double X, Y, Z;
	RGBtoXYZ(R, G, B, X, Y, Z);
	XYZtoLab(X, Y, Z, L, a, b);
}

// converting rgb values into CIELab values in a whole image
void ColorConversion::ColorSpaceConversionFromRGB2Lab(double *pRGB, double *pLab, int width, int height)
{
	int i1, i2, i3, index;
	index=width*height*3;
	for(i1=0, i2=1, i3=2; i3<index; i1+=3, i2+=3, i3+=3)
		RGBtoLab(pRGB[i1], pRGB[i2], pRGB[i3], pLab[i1], pLab[i2], pLab[i3]);
}


// support aggregation using support-weights in both support windows
void CStereoAdaptiveSupportWeight::SupportAggregation(
							float **pSW_ref,		// support-weights of a reference image
							float **pSW_tar,		// support-weights of a target image 
							float **pRawCost,		// per-pixel raw matching cost
							float **pCost,			// aggregated matching cost
							int range,				// disparity range
							int mask,				// radius of a support window 
							int width, int height	// image size 
							)
{
	int i, j, k, l, d, n, ii, index_ref, index1, index2, index_tar;
	//int size = width*height;
	//int mask_pixel_size = (2*mask+1)*(2*mask+1);
	double weight, weight_sum, sum;

	index_ref=0;
	// aggregation 
	for(i=0, ii=0; i<height; i++, ii+=width)
		for(j=0; j<width; j++, index_ref++)
			for(d=0; d<range; d++)
			{					
				index_tar=j-d;
				//border check
				if(index_tar<0) index_tar=width+index_tar;
				else if(index_tar>=width) index_tar=index_tar-width;
				index_tar+=ii;
				for(weight_sum=0, sum=0, n=0, k=-mask; k<=mask; k++)
				{
					index1=i+k;
					//border check
					if(index1<0) index1=height+index1;
					else if(index1>=height) index1=index1-height;
					index1*=width;
					for(l=-mask; l<=mask; l++, n++)
					{
						index2=j+l;
						//border check
						if(index2<0) index2=width+index2;
						else if(index2>=width) index2=index2-width;
						weight_sum+=weight=pSW_ref[index_ref][n]*pSW_tar[index_tar][n];
						sum+=pRawCost[index1+index2][d]*weight;
					}
				}
				pCost[index_ref][d]=(float)(sum/weight_sum);
			}
}




void CStereoAdaptiveSupportWeight::AdaptiveSupportWeightComputation(
								double *pImage,				// input image	
								float **pSW,					// support-weights 
								CYoonStereoParameters param,	// parameters 
								int mask,						// mask size 
								int width, int height			// image size 
								)
{
	int i, j, l, m, ii, x, y, k, index, tar_index, pos_x, pos_y;
	double color_diff;
	double L1, a1, b1, L2, a2, b2;
	ColorConversion cc;

	/* number of pixels in a image */
	int image_size=width*height;
	/* number of pixels in a window */
	int size=(2*mask+1)*(2*mask+1);

	// color space conversion
	double *pImage_CIELab = new double [image_size*3];
	cc.ColorSpaceConversionFromRGB2Lab(pImage, pImage_CIELab, width, height);

	// parameters for weight computation
	double gamma_proximity=param.gamma_proximity;
	double gamma_similarity=param.gamma_similarity;
	
	// weights for proximity
	double *pweight_prox = new double [size];
	for(k=0, y=-mask; y<=mask; y++)
		for(i=y*y, x=-mask; x<=mask; x++, k++)
			pweight_prox[k]=exp(-sqrt((double)(i+x*x))/gamma_proximity);

	// buffer for final support-weights
	for(i=0; i<image_size; i++)
		memset(pSW[i], 0, sizeof(float)*size);

	// computation
	for(i=0; i<height; i++)
	{
		for(ii=i*width, j=0; j<width; j++)
		{
			index=ii+j;
			l=index*3;
			L1=pImage_CIELab[l];
			a1=pImage_CIELab[l+1];
			b1=pImage_CIELab[l+2];
			for(k=0, y=-mask; y<=mask; y++)
			{
				pos_y=i+y;
				//border check
				if(pos_y<0 || pos_y>=height)
					for(x=-mask; x<=mask; k++, x++)
						pSW[index][k]=0;
				else
				{
					pos_y*=width;
					for(x=-mask; x<=mask; k++, x++)
					{
						//border check
						if(pSW[index][k]>0)
							continue;
						pos_x=j+x;
						if(pos_x<0 || pos_x>=width)
							pSW[index][k]=0;

						else
						{
							tar_index=pos_x+pos_y;
							// color difference
							m=tar_index*3;
							L2=pImage_CIELab[m];
							a2=pImage_CIELab[m+1];
							b2=pImage_CIELab[m+2];
							color_diff=cc.EuclideanDistance(L1, a1, b1, L2, a2, b2);
							pSW[index][k]=(float)(pweight_prox[k]*exp(-color_diff/gamma_similarity));
							pSW[tar_index][size-1-k]=pSW[index][k];
						}
					}
				}
			}
		}
	}
	delete pweight_prox;
	delete pImage_CIELab;
}




extern void LASW(CFloatImage& src, CFloatImage& dst, CByteImage& refimage, CByteImage& tarimage, 
				 int xWidth, int yWidth, float gamma_proximity, float gamma_similarity, 
				 int color_space, int diff_iter)
{
	// raw matching cost
    CShape sh = src.Shape();
    //int width1 = sh.width, height1 = sh.height, 
	int n_bands1 = sh.nBands;
	int w = (int)(xWidth/2.0);	// radius of the mask (local window)

    CShape sh_image = refimage.Shape();
    int width = sh_image.width, height = sh_image.height, n_bands = sh_image.nBands;
	
	int size = width*height;
	int bands = n_bands-1;	
	int i, j;
	
	double *pref = new double [size*3];		// reference image (left)
	double *ptar = new double [size*3];;	// target image	(right)
	if(bands==3){
		// color image
		GetDataForProcessing(refimage, pref);
		GetDataForProcessing(tarimage, ptar);
	}

	else
	{
		double *pimage_temp=new double [size];
		GetDataForProcessing(refimage, pimage_temp);
		for(i=j=0; i<size; i++, j+=3)
			pref[j]=pref[j+1]=pref[j+2]=pimage_temp[i];
		GetDataForProcessing(tarimage, pimage_temp);
		for(i=j=0; i<size; i++, j+=3)
			ptar[j]=ptar[j+1]=ptar[j+2]=pimage_temp[i];
		delete pimage_temp;
	}

	CStereoAdaptiveSupportWeight StereoASW;
	CYoonStereoParameters param;
	param.gamma_proximity=gamma_proximity;
	param.gamma_similarity=gamma_similarity;
	param.color_space=color_space;

	int u, k, mm=xWidth*yWidth;
	float **pRawCost = NULL; 
	pRawCost=NEW2D<float>(size, n_bands1);
	float **pCost = NULL; 
	pCost=NEW2D<float>(size, n_bands1);
	float **pSW_ref = NULL; 
	pSW_ref=NEW2D<float>(size, mm);
	float **pSW_tar = NULL; 
	pSW_tar=NEW2D<float>(size, mm);

	GetDataForProcessing(src, pRawCost);

	StereoASW.AdaptiveSupportWeightComputation(pref, pSW_ref, param, w, width, height);
	StereoASW.AdaptiveSupportWeightComputation(ptar, pSW_tar, param, w, width, height);
	for(u=0; u<diff_iter; u++)
	{
		StereoASW.SupportAggregation(pSW_ref, pSW_tar, pRawCost, pCost, n_bands1, w, width, height);
		for(k=0; k<size; k++)
			memcpy(pRawCost[k], pCost[k], sizeof(float)*n_bands1);
	}
	SaveDataInImage(pCost, dst);

	DEL2D(pSW_ref, size);
	DEL2D(pSW_tar, size);
	DEL2D(pRawCost, size);
	DEL2D(pCost, size);
	delete pref;
	delete ptar;	
}
