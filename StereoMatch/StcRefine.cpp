///////////////////////////////////////////////////////////////////////////
//
// NAME
//  StcRefine.cpp -- refine the matching disparity to get a sub-pixel match
//
// DESCRIPTION
//  Input:  m_disparity        selected integer disparities in range [0 .. m_disp_n-1]
//  Output: m_float_disparity  refined floating point disparities
//
//  The mapping between the integral disparities k (ranging from 0 to m_disp_n-1)
//  and the floating point
//  disparity d is given by:
//      d = disp_min + k * m_disp_num / m_disp_den
//
// SEE ALSO
//  StereoMatcher.h
//  StereoMatcher.cpp
//
// Copyright © Richard Szeliski, 2001.
// See Copyright.h for more details
//
///////////////////////////////////////////////////////////////////////////

#include "StereoMatcher.h"
#include "Convert.h"
#include <math.h>

void CStereoMatcher::Refine()
{
    // Refine the matching disparity to get a sub-pixel match

    // skip this step if no optimization took place and disparity image
    // was passed through
    float d_offset = disp_min;
    if (opt_fn != eNoOpt)
    {
        ScaleAndOffset(m_disparity, m_float_disparity,
                       m_disp_step, d_offset);
    }

    // optionally, do sub-pixel fit
    if (! refine_subpix || m_disp_n < 3)
        return;

    CShape sh = m_cost.Shape();
    int width = sh.width, height = sh.height;

    // If we already fitted after the aggr stage, use those results
    //  (these results would be meaningless, because costs have been perturbed)
    if (aggr_subpixel || aggr_collapse && disp_step < 1.0f)
    {
        for (int y = 0; y < height; y++)
        {
            float* cost  = &m_cost.Pixel(0, y, 0);
            int*   disp  = &m_disparity.Pixel(0, y, 0);
            float* fdisp = &m_float_disparity.Pixel(0, y, 0);
            float* mind  = &m_sub_pixel_min.Pixel(0, y, 0);

            for (int x = 0; x < width; x++, cost += m_disp_n, mind += m_disp_n)
            {
                int d_min   = disp[x];
                float x0    = mind[d_min];
                float d_new = m_disp_step * (d_min + x0) + d_offset;
                fdisp[x]    = d_new;
            }
        }
        return;
    }

    if (verbose == eVerboseSummary)
        fprintf(stderr, ", subPix");
    if (verbose >= eVerboseProgress)
        fprintf(stderr, "- refining: subpixel parabola fit\n");

    for (int y = 0; y < height; y++)
    {
        float*  cost = &m_cost.Pixel(0, y, 0);
        int*    disp = &m_disparity.Pixel(0, y, 0);
        float* fdisp = &m_float_disparity.Pixel(0, y, 0);

        for (int x = 0; x < width; x++, cost += m_disp_n)
        {
            // Get minimum, but offset by 1 from ends
            int d_min = disp[x] + (disp[x] == 0) - (disp[x] == m_disp_n-1);

            // Compute the equations of the parabolic fit
            float c0 = cost[d_min-1], c1 = cost[d_min], c2 = cost[d_min+1];
            float a = 0.5 * (c0 - 2.0 * c1 + c2);
            float b = 0.5 * (c2 - c0);
            if (a <= 0.0 || a < 0.5 * fabs(b))
                continue;

            // Solve for minimum
            float x0 = - 0.5 * b / a;
            float d_new = m_disp_step * (d_min + x0) + d_offset;
            fdisp[x] = d_new;
        }
    }

}
