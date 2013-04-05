///////////////////////////////////////////////////////////////////////////
//
// NAME
//  Histogram1D.cpp -- compute a 1-D histogram
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
// SEE ALSO
//  Histogram1D.h       longer description
//  Image.h             image class definition
//
// Copyright © Richard Szeliski, 2001.
// See Copyright.h for more details
//
///////////////////////////////////////////////////////////////////////////

#include "Image.h"
#include "Error.h"
#include <vector>
#include <algorithm>
#include <math.h>
#include <assert.h>
#include "Histogram1D.h"

template <class T>
extern int  Histogram1D(CImageOf<T> img, int n_bands, CByteImage mask,
                        std::vector<int>& counts,
                        float& min_val, float& max_val, float& step_size, int& n_bins)
{
    // Set up the dimensions
    CShape sh = img.Shape(), msh = mask.Shape();
    bool valid_mask = (sh.width == msh.width) && (sh.height == msh.height);
    n_bands = (n_bands > 0) ? n_bands : sh.nBands;

    // Compute the min and max values
    if (min_val >= max_val)
    {
        T min_v = img.MaxVal();
        T max_v = img.MinVal();
        for (int y = 0; y < sh.height; y++)
        {
            T*     p = &img.Pixel(0, y, 0);
            uchar* m = (valid_mask) ? &mask.Pixel(0, y, 0) : 0;
            for (int x = 0; x < sh.width; x++)
            {
                if (m && m[x] == 0)
                    break;
                for (int b = 0; b < n_bands; b++, p++)
                {
                    min_v = __min(min_v, p[0]);
                    max_v = __max(max_v, p[0]);
                }
            }
        }
        max_val = max_v;
        min_val = min_v;
    }

    // Compute the step size and number of bins
    if (n_bins <= 0 && step_size > 0.0)
        n_bins = int(ceil((max_val - min_val) / step_size));
    else if (step_size <= 0 && n_bins > 0)
        step_size = (max_val - min_val) / (float) n_bins;
    else if (n_bins <= 0 && step_size <= 0)
        throw CError("Histogram1D: both step_size and n_bins can't be 0");
    assert(n_bins > 0);

    // Fill in the bins
    counts.resize(n_bins);
    std::fill_n(counts.begin(), n_bins, 0);
    float inv_step_size = 1.0f / step_size;
    float min_v = min_val;
    for (int y = 0; y < sh.height; y++)
    {
        T*     p = &img.Pixel(0, y, 0);
        uchar* m = (valid_mask) ? &mask.Pixel(0, y, 0) : 0;
        for (int x = 0; x < sh.width; x++)
        {
            if (m && m[x] == 0)
                break;
            for (int b = 0; b < n_bands; b++, p++)
            {
                float val = p[0];
                int bin = int((val - min_v) * inv_step_size);
                bin = __min(__max(0, bin), n_bins-1);
                counts[bin] += 1;
            }
        }
    }

    // Return the largest bin size
    int max_count = (std::max_element(counts.begin(), counts.end()))[0];
    return max_count;
}

template <class T>
extern void Histogram1D(CImageOf<T> img, int n_bands, CByteImage mask,
                        CByteImage& figure,
                        float& min_val, float& max_val, float& step_size, int& n_bins,
                        int& height, float& v_scale)
{
    // Compute the histogram counts
    std::vector<int> counts;
    int max_count =
        Histogram1D(img, n_bands, mask, counts,
                    min_val, max_val, step_size, n_bins);

    // Compute the height scaling
    if (height <= 0 && v_scale > 0.0)
        height = int(ceil(max_count * v_scale));
    else if (v_scale <= 0 && height > 0)
        v_scale = height / (float) max_count;
    else if (height <= 0 && step_size <= 0.0)
        throw CError("Histogram1D: both height and v_scale can't be 0");
    assert(height > 0);

    // Re-normalize the counts to be in the range [0,height)
    for (int i = 0; i < n_bins; i++)
        counts[i] = __min(height-1, int(counts[i] * v_scale + 0.5));

    // Draw the figure
    CShape sh(n_bins, height, 1);
    figure.ReAllocate(sh);
    figure.ClearPixels();
    for (int y = 0; y < height; y++)
    {
        uchar* p = &figure.Pixel(0, y, 0);
        int y2 = height-1 - y;
        for (int x = 0; x < n_bins; x++)
            p[x] = (y2 >= counts[x]) ? 255 : 0;
    }
}

extern void Histogram1DInstantiate()
{
    CByteImage mask, figure;
    int n_bins, height;
    float mn, mx, st, v_scale;
    Histogram1D(CByteImage(),  0, mask, figure, mn, mx, st, n_bins, height, v_scale);
    Histogram1D(CIntImage(),   0, mask, figure, mn, mx, st, n_bins, height, v_scale);
    Histogram1D(CFloatImage(), 0, mask, figure, mn, mx, st, n_bins, height, v_scale);
}
