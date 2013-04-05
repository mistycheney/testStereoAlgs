///////////////////////////////////////////////////////////////////////////
//
// NAME
//  BoxFilter.h -- separable box filter (moving average convolution)
//
// SPECIFICATION
//  void BoxFilter(CImageOf<T>& src, CImageOf<T>& dst,
//                 int xWidth, int yWidth, bool average)
//
// PARAMETERS
//  src                 source image
//  dst                 destination image
//  xWidth, yWidth      horizontal and vertical box widths
//  average             scale result down by 1/(xWidth * yWidth)
//
// DESCRIPTION
//  Performs a separable box filtering using efficient running sum code.
//
//  Because a temporary row buffer is used, the src and dst images can be the same
//  (in place convolution).
//
//  If dst is not of the right shape, it is reallocated to the right shape.
//
//  The padding type of src (src.borderMode) determines how pixels are
//  filled.
//
// SEE ALSO
//  BoxFilter.cpp       implementation
//  Image.h             image class definition
//
// Copyright © Daniel Scharstein and Richard Szeliski, 2001.
// See Copyright.h for more details
//
///////////////////////////////////////////////////////////////////////////

template <class T>
void BoxFilter(CImageOf<T>& src, CImageOf<T>& dst,
               int xWidth, int yWidth, bool average);
