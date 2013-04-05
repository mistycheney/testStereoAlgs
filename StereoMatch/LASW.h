///////////////////////////////////////////////////////////////////////////
//
// NAME
//  LASW.h -- Implementation of the "Adaptive Support-Weight Approach for Correspondence Search" 
//
//
// SPECIFICATION
//  void LASW(CFloatImage& src, CFloatImage& dst, CByteImage& refimage, CByteImage& tarimage, 
//				 int xWidth, int yWidth, float gamma_proximity, float gamma_similarity, 
//				 int color_space, int diff_iter)
//
//
// PARAMETERS
//  src                 per-pixel raw matching cost
//  dst                 aggregated matching cost
//	ref					reference image
//  tar					target image
//  xWidth, yWidth      horizontal and vertical box widths
//  gamma_proximity		control parameter in assigning support-weights according to the proximity
//  gamma_similarity	control parameter in assigning support-weights according to the similarity
//  color_space			color space in which the color difference is computed
//  diff_iter           number of iteration
//
//
// DESCRIPTION
//  IEEE Transactions on Pattern Analysis and Machine Intelligence (TPAMI), vol. 28, no. 4, pp. 650-656, 2006. 
//
//
// SEE ALSO
//  LASW.cpp       implementation
//
//
// Copyright Kuk-Jin Yoon, 2006.
//
//  * If you find any bug in this implementation, 
//    please let me know via email. 
//    (kjyoon @ gmail.com or Kuk-Jin.Yoon@inrialpes.fr)
//
///////////////////////////////////////////////////////////////////////////


class CYoonStereoParameters
{
public:
	double gamma_proximity;
	double gamma_similarity;
	int match_max;
	int color_space;
};

// class for color space conversion
// RGB -> XYZ -> CIELab
class ColorConversion
{
private:
	double F(double input);
	void RGBtoXYZ(double R, double G, double B, double &X, double &Y, double &Z);
	void XYZtoLab(double X, double Y, double Z, double &L, double &a, double &b);
	void RGBtoLab(double R, double G, double B, double &L, double &a, double &b);

public:
	void ColorSpaceConversionFromRGB2Lab(double *pRGB, double *pLab, int width, int height);
	double EuclideanDistance(double a1, double a2, double a3, double b1, double b2, double b3);
};


// support-weight computation and support aggregation
class CStereoAdaptiveSupportWeight  
{
public:
	// support-weight computation
	void AdaptiveSupportWeightComputation(double *pImage,				// input image	
										float **pSW,					// support-weights 
										CYoonStereoParameters param,	// parameters 
										int mask,						// mask size 
										int width, int height			// image size 
										);


	// support aggregation
	void SupportAggregation(float **pSW_ref,		// support-weights of a reference image
							float **pSW_tar,		// support-weights of a target image 
							float **pRawCost,		// per-pixel raw matching cost
							float **pCost,			// aggregated matching cost
							int range,				// disparity range
							int mask,				// radius of a support window 
							int width, int height	// image size 
							);
};


// main function
void LASW(CFloatImage& src, CFloatImage& dst, CByteImage& refimage, CByteImage& tarimage, 
				 int xWidth, int yWidth, float gamma_proximity, float gamma_similarity, 
				 int color_space, int diff_iter);
