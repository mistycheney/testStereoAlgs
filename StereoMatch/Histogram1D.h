///////////////////////////////////////////////////////////////////////////
//
// NAME
//  Histogram1D.h -- compute a 1-D histogram
//
// SPECIFICATION
//  int  Histogram1D(CImageOf<T> img, int n_bands, CByteImage mask,
//                   std::vector<int>& counts,
//                   float& min_val, float& max_val, float& step_size, int& n_bins);
//  void Histogram1D(CImageOf<T> img, int n_bands, CByteImage mask,
//                   CByteImage& figure,
//                   float& min_val, float& max_val, float& step_size, int& n_bins,
//                   int& height, float& v_scale);
//
// PARAMETERS
//  img                 image whose histogram is being computed
//  n_bands             number of bands to process (0 = all)
//  mask                optional mask image (only look at non-zero pixels)
//  counts              list of pixel counts
//  figure              gray-scale plot of (scaled) counts
//  min_val             minimum value (given or computed)
//  max_val             maximum value (given or computed)
//  step_size           bin step size in value (given or computed)
//  n_bins              number of bins (given or computed)
//  height              height of image (given or computed)
//  v_scale             scaling from counts to height (given or computed)
//
// DESCRIPTION
//  Compute a 1D histogram of an n-banded image (for RGBA images,
//  set n_bands = 3 to only process RGB).
//
//  The minimum and maximum values will be optionally computed (if equal).
//  The step_size will be computed if n_bins is given, or vice-versa
//  (one or both of these must be non-zero).
//
//  For a figure generation (bar chart), the figure width is given by n_bins.
//  the figure height must be specified, or the vertical scaling v_scale
//  must be given (the other will be computed automatically, if non-zero).
//
//  The returned value is the maximum count in any of the bins.
//
// SEE ALSO
//  Histogram1D.cpp     implementation
//  Image.h             image class definition
//
// Copyright © Richard Szeliski, 2001.
// See Copyright.h for more details
//
///////////////////////////////////////////////////////////////////////////

template <class T>
int  Histogram1D(CImageOf<T> img, int n_bands, CByteImage mask,
                 std::vector<int>& counts,
                 float& min_val, float& max_val, float& step_size, int& n_bins);

template <class T>
void Histogram1D(CImageOf<T> img, int n_bands, CByteImage mask,
                 CByteImage& figure,
                 float& min_val, float& max_val, float& step_size, int& n_bins,
                 int& height, float& v_scale);
