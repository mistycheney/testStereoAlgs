///////////////////////////////////////////////////////////////////////////
//
// NAME
//  MinFilter.cpp -- separable min/max filter
//
// DESCRIPTION
//  Performs a separable box filtering using efficient running sum code.
//
// DESIGN NOTES
//  First version - uses single variable for minimum:
//   potentially inefficent due to large strides in vertical direction
//
// SEE ALSO
//  MinFilter.h         longer description of the interface
//
// Copyright © Daniel Scharstein and Richard Szeliski, 2001.
// See Copyright.h for more details
//
///////////////////////////////////////////////////////////////////////////

#include "Image.h"
#include "Error.h"
#include "Convert.h"
#include "BoxFilter.h"
#include <vector>

// new, faster version: only recompute minimum if old min value leaves current interval

template <class T>
void minFilterLine(T* line, T* buf, int len, int stride,
                          int w, int pl, int pr)
{
    // Min-filter a single row or column of an image
    // line   - original data of length len
    // buf    - buffer of length len+w; buf[0] is offset by -pl
    // len    - length of line
    // stride - distance between elements in line
    // w      - window width (e.g., 5)
    // pl     - offset to first cell to left of window  (e.g., -3)
    // pr     - offset to last cell in window (e.g., 2)
    
    int x, rx;
    
    // copy line into buffer
    for (x = rx = 0; x < len; x++, rx += stride)
        buf[x] = line[rx];
    
    // handle boundaries   
    for (x = pl; x < 0; x++)
        buf[x] = line[0];
    for (x = len; x < len+pr; x++)
        buf[x] = line[stride * (len-1)];
    
    // do the actual minimization
    
    // initialize current minimum
    
    T minval = buf[pl];
    
    for (int x1 = pl+1; x1 < pr; x1++)
        minval = __min(minval, buf[x1]);
    
    for (x = rx = 0; x < len; x++, rx += stride) 
    {
       	T leftval = buf[x+pl];	// value leaving interval
        T rightval = buf[x+pr]; // new value entering interval
        
        if (rightval < minval)
            minval = rightval;			// new value becomes new minimum
        
        if (leftval == minval && leftval < rightval) 
        {
            // leaving value was minimum; need to compute new minimum
            
            minval = rightval;
            
            for (int x1 = x+pl+1; x1 < x+pr; x1++)
                minval = __min(minval, buf[x1]);
        }
        line[rx] = minval;
    }
}

// old slower version
template <class T>
void minFilterLine2(T* line, T* buf, int len, int stride,
                           int w, int pl, int pr)
{
    // Min-filter a single row or column of an image
    // w      - window width (e.g., 5)
    // pl     - offset to first cell to left of window  (e.g., -3)
    // pr     - offset to last cell in window (e.g., 2)
    // line   - original data of length len
    // stride - distance between elements in line
    // buf    - buffer of length len+w; buf[0] is offset by -pl
    
    int x, rx;
    
    // copy line into buffer
    for (x = rx = 0; x < len; x++, rx += stride)
        buf[x] = line[rx];
    
    // handle boundaries   
    for (x = pl; x < 0; x++)
        buf[x] = line[0];
    for (x = len; x < len+pr; x++)
        buf[x] = line[stride * (len-1)];
    
    // do the actual minimization
    
    for (x = rx = 0; x < len; x++, rx += stride) {
        
        T minval = buf[x+pl+1];
        
        for (int x1 = x+pl+2; x1 <= x+pr; x1++)
            minval = __min(minval, buf[x1]);
        
        line[rx] = minval;
    }
}

template <class T>
void MinMaxFilter(CImageOf<T>& src, CImageOf<T>& dst,
                  int xWidth, int yWidth, bool doMax)
{
    // Min filter operation (shiftable windows)
    
    // In each band, replace the value by the minimum value in a square window
    // of size aggr_minfilter.
    // This is identical to using a "shiftable window" to select best disparities.
    
    // first version - uses single variable for minimum
    // potentially inefficent due to large strides in vertical direction
    
    CShape sh = src.Shape();
    int width = sh.width, height = sh.height, bands = sh.nBands;
    
    int w = xWidth;     // window size      (e.g.,  5)
    int pr = w / 2;     // pointer to right (e.g.,  2): value to add
    int pl = pr - w;    // pointer to left  (e.g., -3): value to subtract
    int x, y, d;
    
    // if maximizing, negate all of the values
    int scale  = (! doMax) ? 1 : -1;
    int offset = (! doMax) ? 0 :
                 (typeid(T) == typeid(unsigned char)) ? 255 :
                 (typeid(T) == typeid(int)) ? -1 : 0;

    // Re-allocate dst if necessary
    dst.ReAllocate(sh, false);
    
    // Allocate a buffer for a row plus padding of -pl and pr on left and right
    std::vector<T> buffer;
    buffer.resize(width + w);
    
    for (y = 0; y < height; y++)
    {
        // Copy or negate the source row to the destination
        ScaleAndOffsetLine(&src.Pixel(0, y, 0), &dst.Pixel(0, y, 0),
                           width * bands, scale, offset, (T) 0, (T) 0);

        // Minimize each band separately
        for (d = 0; d < bands; d++)
        {
            T* row = &dst.Pixel(0, y, d); // original row
            T* buf = &buffer[-pl];           // copy is offset by -pl to left

            minFilterLine(row, buf, width, bands, w, pl, pr);
        }
    }
    
    // minimize columns
    w = yWidth;     // window size      (e.g.,  5)
    pr = w / 2;     // pointer to right (e.g.,  2): value to add
    pl = pr - w;    // pointer to left  (e.g., -3): value to subtract
    
    // Allocate a buffer for a column plus padding of -pl and pr on top and bottom
    buffer.resize(height + w);
    
    // compute vertical stride
    int stride = &src.Pixel(0, 1, 0) - &src.Pixel(0, 0, 0);
    
    for (x = 0; x < width; x++)
    {
        // Minimize each band separately
        for (d = 0; d < bands; d++)
        {
            T* col = &dst.Pixel(x, 0, d); // original column
            T* buf = &buffer[-pl];           // copy is offset by -pl to left
            
            minFilterLine(col, buf, height, stride, w, pl, pr);

            // Negate the results if necessary
            if (doMax)
            {
                T* cp = col;
                for (y = 0; y < height; y++, cp += stride)
                    cp[0] = offset - cp[0];
            }
        }
    }
}

template <class T>
extern void MinFilter(CImageOf<T>& src, CImageOf<T>& dst,
                      int xWidth, int yWidth)
{
    MinMaxFilter(src, dst, xWidth, yWidth, false);
}

template void MinFilter<uchar>(CImageOf<uchar>&, CImageOf<uchar>&, int, int);
template void MinFilter<float>(CImageOf<float>&, CImageOf<float>&, int, int);

template <class T>
extern void MaxFilter(CImageOf<T>& src, CImageOf<T>& dst,
                      int xWidth, int yWidth)
{
    MinMaxFilter(src, dst, xWidth, yWidth, true);
}

template void MaxFilter<uchar>(CImageOf<uchar>&, CImageOf<uchar>&, int, int);

extern void CallMinFilters(void)
{
    CByteImage  bimg; MinFilter(bimg, bimg, 5, 5); MaxFilter(bimg, bimg, 5, 5);
    CIntImage   iimg; MinFilter(iimg, iimg, 5, 5); MaxFilter(iimg, iimg, 5, 5);
    CFloatImage fimg; MinFilter(fimg, fimg, 5, 5); MaxFilter(fimg, fimg, 5, 5);
}
