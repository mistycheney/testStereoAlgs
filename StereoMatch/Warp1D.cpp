///////////////////////////////////////////////////////////////////////////
//
// NAME
//  Warp1D.cpp -- forward and inverse 1D (horizontal) warping
//
// DESIGN NOTES
//
// SEE ALSO
//  Warp1D.h            more detailed description
//
// Copyright © Richard Szeliski, 2001.
// See Copyright.h for more details
//
///////////////////////////////////////////////////////////////////////////

#include "Image.h"
#include "Error.h"
#include "Convert.h"
#include "ImageIO.h"
#include <math.h>
#include <vector>

#define ROUND(x) ((int)( (x) >= 0 ? (x) + .5 : (x) - .5))

template <class T>
static inline void draw_intensity_line(T *src1, T *src2, T *dst,
                                       float x1, float x2, int w, int n_bands,
                                       float round_offset,
                                       T minVal, T maxVal)
{
    if (x2 < x1)
        return;     // backward facing line
    const bool clip = (minVal < maxVal);
#if 0   // trim inward to only draw included pixels
    int i0 = __max(0, (int)ceil(x1));
    int i1 = __min(w-1, (int)floor(x2));
#else   // round the coordinates (overlap drawing of vertices)
    int i0 = __max(0, ROUND(x1));
    int i1 = __min(w-1, ROUND(x2));
#endif
    float iden = 1.0 / (x2 - x1 + (x2 == x1));
    T *dp = &dst[i0 * n_bands];
    for (int i = i0; i <= i1; i++, dp += n_bands)
    {
        float f = (i - x1) * iden;
        for (int b = 0; b < n_bands; b++)
        {
            float v1 = src1[b], v2 = src2[b];
            float v = v1 + f * (v2 - v1);
            if (clip)
                dp[b] = (T) __min(__max(v + round_offset, minVal), maxVal);
            else
                dp[b] = (T)(v + round_offset);
        }
    }
}

template <class T>
extern void ForwardWarp(CImageOf<T>& src, CImageOf<T>& dst, CFloatImage& disp,
                        float d_scale, bool line_interpolate, float disp_gap)
{
    // Warps src into dst using disparities disp.
    //  Each disparity is scaled by d_scale
    // Note that "empty" pixels are left at their original value

    CShape sh = src.Shape();
    int w = sh.width, h = sh.height, n_bands = sh.nBands;
    float round_offset = (typeid(T) == typeid(float)) ? 0.0f : 0.5f;

    if (! sh.SameIgnoringNBands(disp.Shape()))
        throw CError("ForwardWarp: disparity image has wrong size");

    if (sh != dst.Shape())
        dst.ReAllocate(sh);
    
    // Optional clipping (if necessary)
    CFloatImage flt;
    T minVal = dst.MinVal();
    T maxVal = dst.MaxVal();
    if (minVal <= flt.MinVal() && maxVal >= flt.MaxVal())
        minVal = maxVal = 0;

    for (int y = 0; y < h; y++)
    {
        // determine correct warping direction
        int xstart = (d_scale>0 ? 0   : w-1);
        int xend   = (d_scale>0 ? w   : -1  );
        int xincr  = (d_scale>0 ? 1   : -1  );
        
        float *dp = &disp.Pixel(0, y, 0);
        T *ps = &src .Pixel(0, y, 0);
        T *pd = &dst .Pixel(0, y, 0);

        for (int x = xstart; x != xend; x += xincr)
        {
            // determine if a line should be drawn
            int x2 = x + xincr;
            float d_diff = fabs(dp[x] - dp[x2]);
            bool draw_line = line_interpolate && (x2 != xend) &&
                 (d_diff < disp_gap);

            // scaled disparity:
            float d = d_scale * dp[x];

            // line drawing
            if (draw_line)
            {
                float d2 = d_scale * dp[x2];

                if (xincr > 0)
                    draw_intensity_line(&ps[x * n_bands], &ps[x2 * n_bands], pd,
                                        x - d, x2 - d2, w, n_bands, round_offset,
                                        minVal, maxVal);
                else
                    draw_intensity_line(&ps[x2 * n_bands], &ps[x * n_bands], pd,
                                        x2 - d, x - d2, w, n_bands, round_offset,
                                        minVal, maxVal);
                continue;
            }
            
            // splatting
            int xx = x - ROUND(d);
            if (xx >= 0 && xx < w)
                memcpy(&pd[xx * n_bands], &ps[x * n_bands],
                       n_bands*sizeof(T));
        }
    }
}

extern float CubicInterpolate(float x0, float v0, float v1, float v2, float v3)
{
    // See Szeliski & Ito, IEE Proc 133(6) 1986.
    float x1 = 1.0f - x0;
    float s0 = v2 - v0;     // slope matches central difference
    float s1 = v1 - v3;     // slope matches central difference
    float d1 = v2 - v1;
    float phi0  = d1 * (x0 * x0) * (2.0f * x1 + 1.0f);
    float phi1a = s0 * x0 * (x1 * x1);
    float phi1b = s1 * x1 * (x0 * x0);
    float v  = v1 + phi0 + phi1a + phi1b;
    return v;
}

void InverseWarpLine(float src[], float dst[], float disp[], int w, int n_bands,
                     int order, float *fwd_disp, float disp_gap)
{
    // Inverse warp a single line
    for (int x = 0; x < w ; x++, dst += n_bands)
    {
        // Warped (source) location
        float d = disp[x];
        float y = x - d;
        if (y < 0.0 || y > w-1)
            continue;

        // Check for occlusion/gap
        int xx = int(y);
        if (fwd_disp && disp_gap &&
            fabs(d - fwd_disp[xx]) >= disp_gap)
            continue;

        // Resampling
        float *ps0 = &src[xx * n_bands];
        if (order == 0 || xx == y)
            memcpy(dst, ps0, n_bands*sizeof(float));
        else if (order == 1 || xx-1 < 0 || xx+2 > w-1)
        {
            // Linear interpolation
            float f = y - xx;
            float *ps1 = &ps0[n_bands];
            for (int b = 0; b < n_bands; b++)
            {
                float v = ps0[b] + f * (ps1[b] - float(ps0[b]));
                dst[b] = v;
            }
        }
        else if (order == 3)
        {
            // Cubic interpolation
            float f = y - xx;
            float *psp = &ps0[-n_bands];
            float *ps1 = &ps0[n_bands];
            float *ps2 = &ps0[2*n_bands];
            for (int b = 0; b < n_bands; b++)
            {
                float v = CubicInterpolate(f, (float) psp[b], (float) ps0[b],
                                              (float) ps1[b],  (float) ps2[b]);
                dst[b] = v;
            }
        }
        else
            throw CError("InverseWarp: order = %d not implemented", order);
    }
}

template <class T>
extern void InverseWarp(CImageOf<T>& src, CImageOf<T>& dst, CFloatImage& disp,
                        float d_scale, float disp_gap, int order)
{
    // Warps src into dst using disparities disp.
    //  Each disparity is scaled by d_scale
    // Note that "empty" pixels are left at their original value

    CShape sh = src.Shape();
    int w = sh.width, h = sh.height, n_bands = sh.nBands;
    int n = w * n_bands;

    if (! sh.SameIgnoringNBands(disp.Shape()))
        throw CError("InverseWarp: disparity image has wrong size");

    if (sh != dst.Shape())
        dst.ReAllocate(sh);

    // Optional forward warped depth map if checking for visibility
    CFloatImage fwd, fwd_tmp;
    if (disp_gap > 0.0f)
    {
        ScaleAndOffset(disp, fwd_tmp, d_scale, 0.0f);
        fwd.ReAllocate(disp.Shape());
        fwd.FillPixels(-9999.0f);
        ForwardWarp(fwd_tmp, fwd, disp, d_scale, true, disp_gap);
    }

    // Allocate line buffers
    std::vector<float> src_buf, dst_buf, dsp_buf;
    src_buf.resize(n);
    dst_buf.resize(n);
    dsp_buf.resize(n);
    CFloatImage fimg;   // dummy, used for MinVal(), MaxVal()

    for (int y = 0; y < h; y++)
    {
        // Set up (copy) the line buffers
        ScaleAndOffsetLine(&src .Pixel(0, y, 0), &src_buf[0], n,
                           1.0f,       0.0f, fimg.MinVal(), fimg.MaxVal());
        ScaleAndOffsetLine(&disp.Pixel(0, y, 0), &dsp_buf[0], w,
                           d_scale, 0.0f, 0.0f, 0.0f);
        ScaleAndOffsetLine(&dst .Pixel(0, y, 0), &dst_buf[0], n,
                           1.0f,       0.0f, fimg.MinVal(), fimg.MaxVal());

        // Forward warp the depth map
        float *fwd_buf = (disp_gap > 0.0f) ? &fwd.Pixel(0, y, 0) : 0;

        // Process (warp) the line
        InverseWarpLine(&src_buf[0], &dst_buf[0], &dsp_buf[0],
                        w, n_bands, order, fwd_buf, disp_gap);

        // Convert back to native type
        T minVal = dst.MinVal();
        T maxVal = dst.MaxVal();
        float offset = (typeid(T) == typeid(float)) ? 0.0f : 0.5;   // rounding
        if (minVal <= fimg.MinVal() && maxVal >= fimg.MaxVal())
            minVal = maxVal = 0;
        ScaleAndOffsetLine(&dst_buf[0], &dst.Pixel(0, y, 0), n,
                           1.0f, offset, minVal, maxVal);
    }
}

extern void CallWarpers(void)
{
    CFloatImage dimg;
    CByteImage  bimg;
    ForwardWarp(bimg, bimg, dimg, 1.0, true, 1.0);
    InverseWarp(bimg, bimg, dimg, 1.0, 1.0, 3);
    CIntImage   iimg;
    ForwardWarp(iimg, iimg, dimg, 1.0, true, 1.0);
    InverseWarp(iimg, iimg, dimg, 1.0, 1.0, 3);
    CFloatImage fimg;
    ForwardWarp(fimg, fimg, dimg, 1.0, true, 1.0);
    InverseWarp(fimg, fimg, dimg, 1.0, 1.0, 3);
}
