///////////////////////////////////////////////////////////////////////////
//
// NAME
//  StcOptDP.cpp -- select the best matches using dynamic programming
//
// DESCRIPTION
//  Input:  m_cost:      DSI (disparity space image)
//  Output: m_disparity  best disparities in range [0 .. n_disp]
//
//  Dynamic programming 1: like Intille & Bobick, no ground control points, but
//  encourage disparity discontinuities along intensity edges
//
//    Each cell in x-d slice has one of three states:
//      0: M - matched
//      1: L - only visible in left image
//      2: R - only visible in right image
//
//    Here is a sample path through a DSI that uses the left image as reference:
//
//    d=3    ' ' ' ' ' ' ' ' M M M ' ' ' ' ' ' 
//    d=2    ' ' ' ' ' M M L ' ' R M M ' ' ' ' 
//    d=1    ' ' ' ' L ' ' ' ' ' ' ' R ' ' ' ' 
//    d=0    M M M L ' ' ' ' ' ' ' ' R M M M M
//
//         x=0 1 2 3 4 5 6 7 8  10  12  14  16
//    
//    The path is traversed from left to right.
//    It can be seen that there are seven possible transitions (current state is upper
//    case, predecessor is lower case):
//
//      0      1     2     3    4    5     6
//
//                    L     M   m    r
//                   /     /    |    |
//    m--M   m--L   l     l     R    R   r--M
//
//    At each cell, we keep three cost values that correspond to the cumulative cost along the 
//    lowest-cost path to the cell, ending in the current state.  These values are stored in
//    the 3-band image sumcost (indexed by d, x rather than x, y),
//    where we use bands 0, 1, 2 for state M, L, R.  The costs are 
//    accumulated by adding a value to the sumcost value of the predecessor cell:
//    M cells get charged the matching cost; L and R cells get charged occlusion costs
//    ocL, ocR.  In addition, for transitions 3 and 6, we charge the horizontal smoothness
//    cost to encourage disparity jumps to align with high intensity gradients.
//  
//
//    NOTE: This implementation makes explicit assumptions about matches "shadowing" other
//    matches, not only in the same DSI column, but also in the diagonal corresponding
//    to the match-frame's line of sight.  In particular, it is assumed that each 
//    d-level in the DSI represents a shift of one pixel to the right.  This is only true
//    if m_disp_step == 1.
//
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

#define N_STATES 3
#define N_TRANS 7

// distribute disparities from certain matches to other pixels
// on each scanline, fill "holes" with closest disparity value on the left
// on left boundary, use closest value to the right
// if revdir==1, go right-to-left instead of left-to-right


// fill occluded pixels in disparity map to get disparity estimates everywhere
// simply fill holes on each scanline from the left (at left edge, fill from right)
// use other direction if revdir==1
static void fill_occluded_pixels(CIntImage dispimg, int occLabel, int revdir)
{
    CShape sh = dispimg.Shape();
    int x, y, w = sh.width, h = sh.height;
    int xstart = (revdir ? w-1 : 0);
    int xend   = (revdir ?  -1 : w);
    int xincr  = (revdir ?  -1 : 1);
    
    for (y = 0; y < h; y++)
    {
        int *disp = &dispimg.Pixel(0, y, 0);
        
        // find first nonoccluded pixel
        for (x = xstart; x != xend; x += xincr)
            if (disp[x] != occLabel)
                break;
            
            int oldd = disp[x];
            
            // fill in all occluded pixels
            for (x = xstart; x != xend; x += xincr)
            {
                int d = disp[x];
                if (disp[x] == occLabel)
                    disp[x] = oldd;
                else
                    oldd = d;
            }
    }
}


void CStereoMatcher::OptDP()
{
    if (m_disp_step != 1.0f)
        fprintf(stderr, 
        "WARNING: dynamic programming will not work properly for m_disp_step != 1\n");

    float ocL = opt_occlusion_cost;
    float ocR = ocL;

    int occLabel = -9999; // marker for occluded pixels 
    // (use 0 if you want to leave occluded pixels black)

    CShape sh = m_cost.Shape();
    int width = sh.width, height = sh.height, n_disp = sh.nBands;

    CFloatImage sumcostIm(n_disp, width, N_STATES);   
    // lowest cumulative cost for current scanline, indexed (d, x, state)

    CIntImage transIm(n_disp, width, N_STATES);      
    // transition corresponding to lowest cost, indexed (d, x, state)
    
    // compute strides (not necessarily == n_disp * N_STATES!)
    int stride = &sumcostIm.Pixel(0, 1, 0) - &sumcostIm.Pixel(0, 0, 0);
    int stride2 = &transIm.Pixel(0, 1, 0) - &transIm.Pixel(0, 0, 0);
    assert(stride == stride2);

    // store current and predecessor state for each transition:
    int cstate[N_TRANS] = {0, 1, 1, 0, 2, 2, 0};
    int pstate[N_TRANS] = {0, 0, 1, 1, 0, 2, 2};
   
    // store "disparity predecessor" for each transition:    
    int dleft = -n_disp, ddiag = -n_disp - 1, dup = 1;
    int pdisp[N_TRANS] = {dleft, dleft, ddiag, ddiag, dup, dup, dleft};
    // store index offset of predecessor for each transition:    
    int left = -stride, diag = -stride - N_STATES, up = N_STATES;
    int pindex[N_TRANS] = {left, left, diag, diag, up, up, left};
    
    // store which transitions don't work for border cases:
    int border0[N_TRANS] = {0, 0, 1, 1, 0, 0, 0}; // for d==0 can't have diag predecessor
    int border1[N_TRANS] = {0, 0, 0, 0, 1, 1, 0}; // for d==max can't have up predecessor
    int x, y, d;

    // process each scanline separately

    for (y = 0; y < height; y++)
    {
        float*    cost = &m_cost.Pixel(0, y, 0);
        int*     trans = &transIm.Pixel(0, 0, 0);
        float* sumcost = &sumcostIm.Pixel(0, 0, 0);

        // initialize first column

        for (d = 0; d < n_disp; d++, trans += N_STATES, sumcost += N_STATES) 
        {
            trans[0] = 0;   // match
            trans[1] = -1;  // safety checks (should never be on best path)
            trans[2] = -1;
            sumcost[0] = cost[d];
            sumcost[1] = COST_MAX; // force first pixel to be non-occluded
            sumcost[2] = COST_MAX;
        }
        cost += n_disp;

        // keep pointer to m_smooth lag one behind, to reference position x-1
        float* smoothcost = &m_smooth.Pixel(0, y, 1); // band 1 is horizontal cost
       
        // now, do dynamic programming step and build up costs

        // for each disparity column
        for (x = 1; x < width; x++)
        {
            trans = &transIm.Pixel(n_disp-1, x, 0);
            sumcost = &sumcostIm.Pixel(n_disp-1, x, 0);

            // for each cell, going down
            for (d = n_disp-1; d >=0; d--, trans -= N_STATES, sumcost -= N_STATES) 

            {
                sumcost[0] = COST_MAX;
                sumcost[1] = COST_MAX;
                sumcost[2] = COST_MAX;
                trans[0] = -1;
                trans[1] = -1;
                trans[2] = -1;


                // try out all transitions
                for (int t=0; t<N_TRANS; t++)
                {
                    // skip transition if predecessor beyond border
                    if ((d == 0 && border0[t]) || 
                        (d == n_disp-1 && border1[t])) {
                        continue;
                    }
                  
                    int current_state = cstate[t];
                    int pred_state = pstate[t];
                    // index offset: "geometric" offset + band offset
                    int pred_index = pindex[t]  + pred_state;

                    // compute candidate cost

                    // increase:
                    float cinc = (current_state == 1 ? ocL :
                                 (current_state == 2 ? ocR :
                                 cost[d]));
                    // add smoothness costs when transitioning from occluded to matched state:
                    if (t==3 || t==6)
                        cinc += smoothcost[0];

                    // cumulative:
                    float c = sumcost[pred_index] + cinc;

                    // update min cost
                    if (c < sumcost[current_state]) 
                    {
                        sumcost[current_state] = c;
                        trans[current_state] = t;
                    }
                }
            
            } // end for d

            cost += n_disp;
            smoothcost += 2;

        } // end for x
        

#if 0   // TODO: remove this debugging code?
        fprintf(stderr, "\n");
        for (d = n_disp-1; d >=0; d--) {
            for (x = width-12; x<width; x++) {
//            for (x = 0; x<10; x++) {
                fprintf(stderr, "%d,%d,%d ", 
                    transIm.Pixel(d, x, 0), transIm.Pixel(d, x, 1), transIm.Pixel(d, x, 2));
            }
            fprintf(stderr, "\n");
        }
     //   exit(1);
#endif

        // now, find best sumcost in last column and backtrack path, storing best disp's

        int* disp = &m_disparity.Pixel(0, y, 0);

        x = width-1;
        trans = &transIm.Pixel(0, x, 0);
        sumcost = &sumcostIm.Pixel(0, x, 0);

        float best_cost = COST_MAX;
        int best_d = 0;
        int st = 0;

        for (d = 0; d < n_disp; d++, sumcost += N_STATES, trans += N_STATES) {
           // for (int s=0; s<N_STATES; s++) {
            // allow only matched rightmost pixels:
            {int s = 0;
                if (sumcost[s] < best_cost) {
                    best_cost = sumcost[s];
                    st = s;
                    best_d = d;
                }
            }
        }

        d = best_d;

        // found end of best path in last column
        trans = &transIm.Pixel(d, x, st);
        sumcost = &sumcostIm.Pixel(d, x, st);

#if 0   // old debugging code
        if (y<3)
        fprintf(stderr, "\n\nlowest cost = %g (%g), d = %d, final trans = %d\n",
            best_cost, sumcost[0], d, trans[0]);
#endif

        int* trans0 = &transIm.Pixel(0, 0, 0);
        

        // backtrack
        while (trans >= trans0)
        {
            int t = trans[0];
                    
#if 0   // old debugging code
            if (y<3 && (x < 5 || x > width-15))
                fprintf(stderr, "cost = %g, trans = %d, x = %d, d = %d\n",
                    sumcost[0], t, x, best_d);
#endif

            int current_state = cstate[t];
            int pred_state = pstate[t];

            // record current disparity
            if (current_state == 0) // matched cell
                disp[x] = d;
            else
                disp[x] = occLabel; // indicator for occluded pixel

            // update pointers
            int pred_index = pindex[t] - current_state + pred_state;
            trans += pred_index;
            sumcost += pred_index; // only for printing

            // update disparity d
            d += pdisp[t];
            if (d < 0) {
                d += n_disp;
                x--;
            }
        }
            
    } // end for y

    if (occLabel != 0)
        fill_occluded_pixels(m_disparity, occLabel, 0);




}
