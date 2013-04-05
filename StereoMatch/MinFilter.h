///////////////////////////////////////////////////////////////////////////
//
// NAME
//  MinFilter.h -- separable min/max filter
//
// SPECIFICATION
//  void MinFilter(CImageOf<T>& src, CImageOf<T>& dst,
//                 int xWidth, int yWidth)
//  void MaxFilter(CImageOf<T>& src, CImageOf<T>& dst,
//                 int xWidth, int yWidth)
//
// PARAMETERS
//  src                 source image
//  dst                 destination image
//  xWidth, yWidth      horizontal and vertical filter widths
//
// DESCRIPTION
//  Performs a separable box filtering using efficient running sum code.
//
// SEE ALSO
//  MinFilter.cpp       implementation
//  Image.h             image class definition
//
// Copyright © Daniel Scharstein and Richard Szeliski, 2001.
// See Copyright.h for more details
//
///////////////////////////////////////////////////////////////////////////

template <class T>
void MinFilter(CImageOf<T>& src, CImageOf<T>& dst,
               int xWidth, int yWidth);

template <class T>
void MaxFilter(CImageOf<T>& src, CImageOf<T>& dst,
               int xWidth, int yWidth);
