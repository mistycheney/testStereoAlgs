///////////////////////////////////////////////////////////////////////////
//
// NAME
//  Warp1D.h -- forward and inverse 1D (horizontal) warping
//
// SPECIFICATION
//  void ForwardWarp(CImageOf<T>& src, CImageOf<T>& dst, CFloatImage& disp,
//                   float d_scale, bool line_interpolate, float eval_disp_gap)
//  void InverseWarp(CImageOf<T>& src, CImageOf<T>& dst, CFloatImage& disp,
//                   float d_scale, float eval_disp_gap)
//
// PARAMETERS
//  src                 source image
//  dst                 destination image
//  disp                floating point disparity (horizontal motion) image
//  d_scale             scale up disparity before warping
//  line_interpolate    render adjacent pixels as lines (unless gap)
//  eval_disp_gap       don't draw across gaps bigger than this
//
// DESCRIPTION
//  Performs a 1D forward or inverse warp, using the given (scaled) disparity.
//
//  The forward warper...
//
// SEE ALSO
//  Warp1D.cpp          implementation
//  Image.h             image class definition
//
// Copyright © Richard Szeliski, 2001.
// See Copyright.h for more details
//
///////////////////////////////////////////////////////////////////////////

template <class T>
void ForwardWarp(CImageOf<T>& src, CImageOf<T>& dst, CFloatImage& disp,
                 float d_scale, bool line_interpolate, float eval_disp_gap);

template <class T>
void InverseWarp(CImageOf<T>& src, CImageOf<T>& dst, CFloatImage& disp,
                 float d_scale, float eval_disp_gap, int order);

float CubicInterpolate(float x0, float v0, float v1, float v2, float v3);
