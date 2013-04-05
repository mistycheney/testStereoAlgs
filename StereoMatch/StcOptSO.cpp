///////////////////////////////////////////////////////////////////////////
//
// NAME
//  StcOptSO.cpp -- select the best matches using scan-line optimization
//
// DESCRIPTION
//  Input:  m_cost:      DSI (disparity space image)
//  Output: m_disparity  best disparities in range [0 .. n_disp]
//
//  "Scan-line optimization" using dynamic programming: perform the same energy minimization
//  as other global algorithms (e.g., graph cuts), but on each scanline separately.
//  So we get a global minimum on each scan-line, but no inter-scanline consistency
//
//  Unlike other DP approaches, this method is asymmetric and does not compute occlusion
//  information.  Rather, it simply finds the best match in each DSI column, while at the
//  same time minimizing the smoothness cost.
//
// SEE ALSO
//  StereoMatcher.h
//  StereoMatcher.cpp
//
// Copyright © Daniel Scharstein, 2001.
// See Copyright.h for more details
//
///////////////////////////////////////////////////////////////////////////

#include "StereoMatcher.h"
#include <assert.h>

void CStereoMatcher::OptSO()
{

    CShape sh = m_cost.Shape();
    int width = sh.width, height = sh.height, n_disp = sh.nBands;

    CFloatImage sumcostIm(width, 1, n_disp);   
    // lowest cumulative cost for current scanline

    CIntImage transIm(width, 1, n_disp);      
    // transition (i.e., previous d) corresponding to lowest cost

    int x, y, d;

    // process each scanline separately

    for (y = 0; y < height; y++)
    {
        // PART ONE: compute cumulative costs, going left to right

        float*    cost = &m_cost.Pixel(0, y, 0);
        int*     trans = &transIm.Pixel(0, 0, 0);
        float* sumcost = &sumcostIm.Pixel(0, 0, 0);

        // initialize first column

        for (d = 0; d < n_disp; d++)
        {
            trans[d] = -1;  // safety check (should never be on best path)
            sumcost[d] = cost[d];
        }
        cost += n_disp;
        trans += n_disp;
        sumcost += n_disp;

        // keep pointer to m_smooth lag one behind, to reference position x-1
        float* smoothcost = &m_smooth.Pixel(0, y, 1); // band 1 is horizontal cost
       
        // now, do dynamic programming step and build up costs

        // for each disparity column
        for (x = 1; x < width; x++)
        {
            // for each disparity
	    int d;
            for (d = 0; d < n_disp; d++)
            {
                float best_cost = COST_MAX;
                int best_d1 = -1;

                // look at all pixels in previous column
            
                for (int d1 = 0; d1 < n_disp; d1++)
                {
                    // index into previous column
                    int dprev = d1 - n_disp;

                    // candidate cost
                    float c = sumcost[dprev];

                    // add smoothness cost if at different disparity
                    if (d1 != d)
                        c += smoothcost[0];

                    // update min
                    if (c < best_cost) {
                        best_cost = c;
                        best_d1 = d1;
                    }
                }

                sumcost[d] = best_cost + cost[d];
                trans[d] = best_d1;

            } // end for d

            cost += n_disp;
            trans += n_disp;
            sumcost += n_disp;
            smoothcost += 2;

        } // end for x



        // PART TWO: backtrack path of least cost right-to-left and store disparities

        int* disp = &m_disparity.Pixel(0, y, 0);

        // find best sumcost in last column

        x = width-1;
        trans = &transIm.Pixel(x, 0, 0);
        sumcost = &sumcostIm.Pixel(x, 0, 0);

        float best_cost = COST_MAX;
        int best_d = 0;
        
        for (d = 0; d < n_disp; d++) {
            if (sumcost[d] < best_cost) {
                best_cost = sumcost[d];
                best_d = d;
            }
        }

        d = best_d;

#if 0   // old debugging code
        if (y<3)
        fprintf(stderr, "\n\nlowest cost = %g (%g), d = %d, final trans = %d\n",
            best_cost, sumcost[d], d, trans[d]);
#endif

        // backtrack
        while (x >= 0)
        {
            assert(d >= 0 && d < n_disp);
            disp[x] = d;
            d = trans[d];
            x--;
            trans -= n_disp;
        }
        
    } // end for y
}
