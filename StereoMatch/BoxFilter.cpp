///////////////////////////////////////////////////////////////////////////
//
// NAME
//  BoxFilter.cpp -- separable box filter (moving average convolution)
//
// DESIGN NOTES
//  Two different versions with different computation inner loops and
//  memory access patterns are implemented.  See comments below.
//
// SEE ALSO
//  BoxFilter.h         longer description of the interface
//
// Copyright © Daniel Scharstein and Richard Szeliski, 2001.
// See Copyright.h for more details
//
///////////////////////////////////////////////////////////////////////////

#include "Image.h"
#include "Error.h"
#include "Convert.h"
#include "BoxFilter.h"
#include <assert.h>
#include <vector>

static int borderIndex(int index, int len, EBorderMode m) {
    // TODO:  this code is common with Convolve.cpp:  should be merged...

    // returns correct index depending on border mode
    // cannot be used with eBorderZero since it returns index, not value
    
    if (index < 0) 
    {
        assert(-index < len);   // "double" wrapping not supported
        switch (m) {
        case eBorderReplicate:
            return 0;
        case eBorderReflect:
            return -index;
        case eBorderCyclic:
            return index+len;
        default:
            throw CError("borderIndex: bad borderMode");
        }
    }
    else if (index >= len)
    {
        assert(len+len-index >= 2); // "double" wrapping not supported
        switch (m) {
        case eBorderReplicate:
            return len-1;
        case eBorderReflect:
            return len+len-index-2;
        case eBorderCyclic:
            return index-len;
        default:
            throw CError("borderIndex: bad borderMode");
        }
    }
    else
        return index;
}

///////////////////////////////////////////////////////////////////////////
// box filter functions

template <class T>
static void boxFilterLines(T* src, T* dst, int n_sums,
                           int len, int stride, int w, int pl, int pr, EBorderMode borderMode,
                           bool average)
{
    // Box filter multiple rows or columns of an image
    // w      - window width (e.g., 5)
    // pl     - offset to first cell to left of window  (e.g., -3)
    // pr     - offset to last cell in window (e.g., 2)
    // src    - pointer to original data
    // dst    - pointer to result
    // n_sums - number of running sums to be kept
    // len    - length of each row or column
    // stride - distance between elements in line
    // average- compute an average

   
    int i, k, x, rx;
    
    T scale = (T)((average) ? 1.0 / w : 1);  // factor for computing average

    // compute first value (x==0)
    for (k=0; k<n_sums; k++) 
    {   
        T sum = 0;
        //left half:
        if (borderMode != eBorderZero)
        {
            for (i=pl+1; i<0; i++)
                sum += scale * src[k+stride * borderIndex(i, len, borderMode)];
        }
        // right half:
        for (i=0; i<=pr; i++)
            sum += scale * src[k+stride * i];

        // store average
        dst[k] = sum;
    }

    // compute values in left border area (x=1..-pl+1)
    x = 1;
    rx = x * stride;
    T *dst_current = &dst[rx];
    T *dst_previous = &dst[rx-stride];
    T *src_left = &src[rx+pl*stride];
    T *src_right = &src[rx+pr*stride];

    for (; x < -pl; x++) 
    {
        for (k=0; k<n_sums; k++) 
        {
            dst_current[k] = dst_previous[k] 
                + scale *
                    (-src[k+stride * borderIndex(x+pl, len, borderMode)]    // -src_left[k]
                     + src_right[k]);
        }
        dst_current += stride;
        dst_previous += stride;
        src_left += stride;
        src_right += stride;
    }


    // compute values in center area (x=-pl .. len-pr-1)
    x = -pl;
    rx = x * stride;
    assert (dst_current == &dst[rx]);

    for (; x < len-pr; x++) 
    {
        for (k=0; k<n_sums; k++) 
        {
            dst_current[k] = dst_previous[k] + scale * (-src_left[k] + src_right[k]);
        }
        dst_current += stride;
        dst_previous += stride;
        src_left += stride;
        src_right += stride;
    }

    // compute values in right border area (x=len-pr .. len-1)
    x = len-pr;
    rx = x * stride;
    assert (dst_current == &dst[rx]);

    for (; x < len; x++) 
    {
        for (k=0; k<n_sums; k++) 
        {
            dst_current[k] = dst_previous[k] 
                + scale * 
                    (-src_left[k] 
                     +src[k+stride * borderIndex(x+pr, len, borderMode)]);   // src_right[k];
        }
        dst_current += stride;
        dst_previous += stride;
        src_left += stride;
        src_right += stride;
    }
  
}


template <class T>
extern void BoxFilter(CImageOf<T>& src, CImageOf<T>& dst,
                      int xWidth, int yWidth, bool average)
{
    // Box filter operation
    // second version - operates in memory order to avoid swapping
    if (xWidth != yWidth)
        throw CError("BoxFilter: xWidth != yWidth not implemented yet");
    
    CShape sh = src.Shape();
    int width = sh.width, height = sh.height, n_bands = sh.nBands;

    dst.ReAllocate(sh, false);     // allocate memory for copy of src
    CImageOf<T> tmp;
    tmp.ReAllocate(sh);     // allocate memory for copy of src
    tmp.borderMode = src.borderMode;

    int w = xWidth;  // window size      (e.g.,  5)
    int pr = w / 2;     // pointer to right (e.g.,  2): value to add
    int pl = pr - w;    // pointer to left  (e.g., -3): value to subtract
    
    // aggregate each row, all disparity levels in parallel

    for (int y = 0; y < height; y++)
    {
        T* src_row = &src.Pixel(0, y, 0);
        T* dst_row = &tmp.Pixel(0, y, 0);
        
        boxFilterLines(src_row, dst_row, n_bands, width, n_bands, w, pl, pr, src.borderMode, average);
    }
    
    // aggregate all columns at all disparity levels in parallel
    
    // compute vertical stride
    int stride  = &dst.Pixel(0, 1, 0) - &src.Pixel(0, 0, 0);
    int stride2 = &tmp.Pixel(0, 1, 0) - &tmp.Pixel(0, 0, 0);
    assert(stride == stride2);

    T* src_cols = &tmp.Pixel(0, 0, 0);
    T* dst_cols = &dst.Pixel(0, 0, 0);
        
    boxFilterLines(src_cols, dst_cols, n_bands*width, height, stride, w, pl, pr, dst.borderMode, average);
}

extern void BoxFilter(CByteImage& src, CByteImage& dst,
                      int xWidth, int yWidth, bool average)
{
    // Need to use an integer accumulator for precision
    CIntImage tmp(src.Shape());
    CopyPixels(src, tmp);
    BoxFilter(tmp, tmp, xWidth, yWidth, average);
    dst.ReAllocate(src.Shape(), false);
    CopyPixels(tmp, dst);
}

extern void CallBoxFilters(void)
{
    CByteImage  bimg; BoxFilter(bimg, bimg, 5, 5, false);
    CIntImage   iimg; BoxFilter(iimg, iimg, 5, 5, false);
    CFloatImage fimg; BoxFilter(fimg, fimg, 5, 5, false);
}
