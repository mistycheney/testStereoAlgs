///////////////////////////////////////////////////////////////////////////
//
// NAME
//  StcAggregate.cpp -- use spatial aggregation to improve the matching scores
//
// DESCRIPTION
//  aggregate disparity space image (DSI) m_cost
//
// SEE ALSO
//  StereoMatcher.h
//  StereoMatcher.cpp
//
// Copyright © Richard Szeliski and Daniel Scharstein, 2001.
// See Copyright.h for more details
//
///////////////////////////////////////////////////////////////////////////

#include "Error.h"
#include "StereoMatcher.h"
#include "Convert.h"
#include "Convolve.h"
#include "ImageIO.h"
#include "BoxFilter.h"
#include "MinFilter.h"
#include "LASW.h"
#include <math.h>

void CStereoMatcher::WriteCosts(CFloatImage& cost_img, char *filename)
{
#ifdef CAN_VIEW_PMFS    // a PMF viewer exists
    // Write this as a multi-band float image
    char tmp_file[1024];
    strcpy(tmp_file, filename);
    char *pct = strrchr(tmp_file, '%');
    strcpy(pct, ".pmf");
    WriteImage(cost_img, tmp_file);
#else
    // Write out the different disparity images
    static float scaleUp = 2; // makes it easier to see DSI
    CShape sh = cost_img.Shape();
    int w = sh.width, h = sh.height, n_disp = sh.nBands;
    CFloatImage cplane(w, h, 1);
    CByteImage  iplane(w, h, 1);
    for (int d = 0; d < n_disp; d++)
    {
        char tmp_file[1024];
        sprintf (tmp_file, filename, d);
        BandSelect(cost_img, cplane, d, 0);
        ScaleAndOffset(cplane, iplane, scaleUp, 0);
        WriteImage(iplane, tmp_file);
    }
#endif
}

void CStereoMatcher::AggrBox()
{
    // Box filter operation: see BoxFilter.h
    BoxFilter(m_cost, m_cost, aggr_window_size, aggr_window_size, true);
}


void CStereoMatcher::AggrMin()
{
    MinFilter(m_cost, m_cost, aggr_minfilter, aggr_minfilter);
}

void CStereoMatcher::AggrAW(int iter)
{
    // Adaptive Support Weight Approach: see LASW.h
    LASW(m_cost,		// initial matching cost
	m_cost,			// aggregated matching cost
	m_reference,		// reference image
	m_matching,		// target image
	aggr_window_size,	// window size - x
	aggr_window_size,	// window size - y
	aggr_gamma_proximity,	// gamma_p
	aggr_gamma_similarity,	// gamma_c
	aggr_color_space,	// color space
	iter			// iteration number (aggregation)
	);
}

void CStereoMatcher::AggrSubPixelFit()
{
    // Replace each aggregated cost with the minimum around that point
    //  (within a half-level).  Note that this operation changes m_cost.
    //  The raw (unaggregated) costs are still available in m_cost0.

    // This code is modeled after Refine(), but doesn't just evaluate
    //  at the finally selected disparity

    // Location of lowest interpolated value
    CShape s = m_cost.Shape();
    m_sub_pixel_min.ReAllocate(s);
    m_sub_pixel_cert.ReAllocate(s);

    for (int y = 0; y < s.height; y++)
    {
        float* cost = &m_cost.Pixel(0, y, 0);
        float* mind = &m_sub_pixel_min .Pixel(0, y, 0);
        float* cert = &m_sub_pixel_cert.Pixel(0, y, 0);

        for (int x = 0, l = 0; x < s.width; x++, l += m_disp_n)
        {
            float c0 = cost[l], c1 = cost[l], c2 = cost[l];

            for (int d = 0; d < m_disp_n; d++)
            {
                // Update the three cost values
                c0 = c1, c1 = c2, c2 = cost[l+__min(m_disp_n-1, d+1)];
                mind[l+d] = 0.0f;
                cert[l+d] = 0.0f;

                // Check for valud values
                if (c0 == m_match_outside ||
                    c1 == m_match_outside ||
                    c2 == m_match_outside)
                    continue;

                // If it's a local minimum, fit a parabola
                if (c1 <= c0 && c1 <= c2)
                {
                    // Compute the equations of the parabolic fit
                    float a = 0.5 * (c0 - 2.0 * c1 + c2);
                    float b = 0.5 * (c2 - c0);
                    if (a <= 0.0 || a < 0.5 * fabs(b))
                        continue;   // more than a 1/2 pixel correction

                    // Solve for minimum
                    float dn  = - 0.5 * (b / a);
                    float cn  = c1 + 0.5 * b * dn;
                    if (cn < 0.0)
                        continue;   // must be a bad parabolic fit (aliased?)
                    cost[l+d] = __max(0.0f, cn);
                    mind[l+d] = dn;
                    cert[l+d] = a;
                }

                // Else find the minimum half-value
                else
                {
                    cost[l+d] = 0.5f * (c1 + __min(c0, c2));
                    mind[l+d] = (c0 < c2) ? -0.5f : 0.5f;
                    cert[l+d] = 0.0f;
                }
            }
        }
    }
}

void CStereoMatcher::AggrCollapse()
{
    // Collapse the DSI to integer disparities
    if (m_disp_step >= 1.0f)
        return;                 // no need to collapse
    int df = int(m_disp_step_inv + 0.5), df2 = df / 2;
    if (df != m_disp_step_inv)
        throw CError("AggrCollapse:  disparity step %g is not a pure fraction",
                     m_disp_step);

    // Recompute the new step size and n_steps
    m_disp_step = m_disp_step_inv = 1.0f;
    m_disp_n = (disp_max - disp_min) + 1;
    CShape s1 = m_cost.Shape(), s2 = s1;
    s2.nBands = m_disp_n;

    // Allocate the new images
    CFloatImage cost(s2);           // collapsed matching costs
    CFloatImage sub_pixel_min(s2);  // collapsed location  of lowest value
    CFloatImage sub_pixel_cert(s2); // collapsed certainty in lowest value

    // Collapse the DSI to integer disparities
    for (int y = 0; y < s2.height; y++)
    {
        for (int x = 0, l = 0; x < s2.width; x++, l += m_disp_n)
        {
            float* cost1 = &m_cost.Pixel(x, y, 0);
            float* mind1 = &m_sub_pixel_min .Pixel(x, y, 0);
            float* cert1 = &m_sub_pixel_cert.Pixel(x, y, 0);
            float* cost2 = &cost.Pixel(x, y, 0);
            float* mind2 = &sub_pixel_min .Pixel(x, y, 0);
            float* cert2 = &sub_pixel_cert.Pixel(x, y, 0);

            // Process each new (integral)disparity
            for (int d1 = 0, d2 = 0; d2 < m_disp_n; d2++)
            {
                int d1_end = __min(disp_n, d2 * df + df - df2);
                int d1_bst = d1;

                // Find the minimum cost
                //  TODO:  consider searching all the way up to including d1_end
                //    (but then remember to re-set d1-- after the loop)
                for (d1 = d1+1; d1 < d1_end; d1++)
                {
                    if (cost1[d1] < cost1[d1_bst])
                        d1_bst = d1;
                }

                // Copy the best value
                cost2[d2] = cost1[d1_bst];
                if (aggr_subpixel)          // already computed
                {
                    mind2[d2] = (mind1[d1_bst] + d1_bst - d2*df)*disp_step;
                    cert2[d2] = cert1[d1_bst];
                }
                else
                {
                    mind2[d2] = (d1_bst - d2*df)*disp_step;
                    cert2[d2] = 0.0;

#if 0   // This code is currently disabled, since it seems to me that if we
        // want sub-pixel fitting of greater accuracy than disp_step, we should
        // turn on aggr_subpixel.
        // TODO:  is this really the best decisions???
                    // Try to fit a parabola to the minimum
                    if (d1_bst > 0 && d1_bst < disp_n-1)
                    {
                        float c0 = cost1[d1_bst-1];
                        float c1 = cost1[d1_bst];
                        float c2 = cost1[d1_bst+1];
                        if (c1 <= c0 && c1 <= c2)
                        {
                            float a = 0.5 * (c0 - 2.0 * c1 + c2);
                            float b = 0.5 * (c2 - c0);
                            if (a > 0.0 && a >= 0.5 * fabs(b))
                                continue;   // under than a 1/2 step correction

                            // Solve for minimum
                            float dn   = - 0.5 * (b / a);
                            mind2[d2] += dn * disp_step;
                            cert2[d2]  = a;
                        }
                    }
#endif
                }
            }
            // TODO:  go through the values in cost2[] and look for duplicate
            //  minima?  If found, can bump the one with a larger fabs(mind2[])
            //  up by a little bit...
        }
    }

    // Clobber the old images
    m_cost           = cost;
    m_sub_pixel_min  = sub_pixel_min;
    m_sub_pixel_cert = sub_pixel_min;
}

///////////////////////////////////////////////////////////////////////////
//
// utility routines for step timing
//

void CStereoMatcher::StartTiming()
{
    m_start_time = clock();    // record start time
}

void CStereoMatcher::PrintTiming()
{
    clock_t end_time = clock();    // record end time
    float m_elapsed_time = (end_time-m_start_time)/(float)CLOCKS_PER_SEC;
    if (verbose >= eVerboseTiming)
        fprintf(stderr, "  * time: %gs\n", m_elapsed_time);
}

///////////////////////////////////////////////////////////////////////////
//
// main dispatch function for aggregation
//

void CStereoMatcher::Aggregate()
{
    // Use spatial aggregation to improve the matching scores
    
    StartTiming();

    // Save the raw matching costs in m_cost0;
    CopyPixels(m_cost, m_cost0);

    // Perform given number of iteration steps
    for (int iter = 0; iter < aggr_iter; iter++) {

        switch (aggr_fn) {
            
        case eBox:        // box filter (square window)
            
           if (verbose == eVerboseSummary && iter < 1)
                fprintf(stderr, ", box=%d", aggr_window_size);
           if (verbose >= eVerboseProgress && iter < 3)
                fprintf(stderr, "- aggregating: box filter, size %d\n", aggr_window_size);
            
            AggrBox();
            
            break;
            
        case eBinomial:   // binomial filter (1 4 6 4 1)

            if (verbose == eVerboseSummary && iter < 1)
                fprintf(stderr, ", binomial filter");
            if (verbose >= eVerboseProgress && iter < 3)
                fprintf(stderr, "- aggregating: binomial 14641 filter\n");

            ConvolveSeparable(m_cost, m_cost, ConvolveKernel_14641, ConvolveKernel_14641,
                 1.0f, 0.0f, 1, 1);
            break;


            // The following three options implement three diffusion algorithms described
            // in "Stereo Matching with Nonlinear Diffusion", Scharstein & Szeliski,
            // IJCV 28(2), 1998.  Each needs to be run for several iterations.
			// 
			// (see stcDiffusion.cpp)

        case eDiffusion:   // simple diffusion: eqn (5)

            if (verbose == eVerboseSummary && iter < 1)
                fprintf(stderr, ", diffusion (lambda=%g)", diff_lambda);
            if (verbose >= eVerboseProgress && iter < 3)
                fprintf(stderr, "- aggregating: diffusion, lambda = %g\n", diff_lambda);
            
            diff_beta = 0.0f;  // make sure beta is 0 to get regular diffusion
            AggrDiffusion(iter);
            break;
            
        case eMembrane:    // membrane model: eqn (7)

            if (verbose == eVerboseSummary && iter < 1)
                fprintf(stderr, ", membrane (lambda=%g, beta=%g)", diff_lambda, diff_beta);
            if (verbose >= eVerboseProgress && iter < 3)
                fprintf(stderr, "- aggregating: membrane model, lambda = %g, beta = %g\n",
                diff_lambda, diff_beta);

            AggrDiffusion(iter);
            break;
            
        case eBayesian:   // Bayesian model (mean field): eqns (26), (27), (30), (29)
            
            if (verbose == eVerboseSummary && iter < 1)
                fprintf(stderr, ", Bayesian (scale=%g, mu=%g, sigmaP=%g, epsP=%g)",
                diff_scale_cost, diff_mu, diff_sigmaP, diff_epsP);
            if (verbose >= eVerboseProgress && iter < 3)
                fprintf(stderr, 
                "- aggregating: Bayesian model, scale=%g, mu=%g, sigmaP=%g, epsP=%g\n",
                diff_scale_cost, diff_mu, diff_sigmaP, diff_epsP);
            
            AggrBayesian(iter);
            break;

        case eASWeight:  
			
            if (verbose == eVerboseSummary && iter < 1)
                fprintf(stderr, ", AdaptiveWeight (box=%d gamma_p=%g gamma_s=%g color_space=%d )", 
				aggr_window_size, aggr_gamma_proximity, aggr_gamma_similarity, aggr_color_space);
            if (verbose >= eVerboseProgress && iter < 3)
                fprintf(stderr, "- aggregating: Adaptive Suport Weight, box=%d gamma_p=%g gamma_s=%g color_space=%d \n", 
				aggr_window_size, aggr_gamma_proximity, aggr_gamma_similarity, aggr_color_space);

			AggrAW(aggr_iter);
    	    iter=aggr_iter;
            break;           
          
		default:
            throw CError("CStereoMatcher::Aggregate(): unknown aggregation function");
            
        }
        
        if (verbose == eVerboseSummary && iter == 1)
            fprintf(stderr, ", %d iterations", aggr_iter);
        if (verbose >= eVerboseProgress && iter == 3)
            fprintf(stderr, "...  (doing %d iterations total)\n", aggr_iter);
        
    }

    PrintTiming();
    if (verbose >= eVerboseTiming && aggr_iter > 1) {
        fprintf(stderr, "    time per iteration: %gs\n", m_elapsed_time / aggr_iter);
    }

    // Perform the optional min-filtering
    if (aggr_minfilter > 1) {
        if (verbose == eVerboseSummary)
            fprintf(stderr, ", minf=%d", aggr_minfilter);
        if (verbose >= eVerboseProgress)
            fprintf(stderr, "- min filter (shiftable windows): filter size %d\n", aggr_minfilter);
        
        StartTiming();
        AggrMin();
        PrintTiming();
    }
    
    // Pad the outside costs back up to bad values
    // TODO:  this should be guarded by a global flag corresponding to 
    //  undefined_cost in StcRawCosts.cpp
    PadCosts();

    // Perform the optional sub-pixel fitting
    if (aggr_subpixel) {
        if (verbose == eVerboseSummary)
            fprintf(stderr, ", aggr_subpixel");
        if (verbose >= eVerboseProgress)
            fprintf(stderr, "- sub-pixel fit\n");

#if 0   // temporary debugging (to see effects of fit)
        if (verbose >= eVerboseDumpFiles)
            WriteImage(m_cost, "reprojected/DSIbf.pmf");
#endif
        StartTiming();
        AggrSubPixelFit();
        PrintTiming();
#if 0   // temporary debugging (to see effects of fit)
        if (verbose >= eVerboseDumpFiles)
            WriteImage(m_sub_pixel_min, "reprojected/DSIspm.pmf");
#endif
    }

    // Perform the optional collapsing to integer disparities
    if (aggr_collapse) {
        if (verbose == eVerboseSummary)
            fprintf(stderr, ", aggr_collapse");
        if (verbose >= eVerboseProgress)
            fprintf(stderr, "- collapse to integral DSI\n");
#if 0   // temporary debugging (to see effects of collapse)
        if (verbose >= eVerboseDumpFiles)
            WriteImage(m_cost, "reprojected/DSIbc.pmf");
#endif
        StartTiming();
        AggrCollapse();
        PrintTiming();
    }

    // Write out the aggregated disparity images
    if (verbose >= eVerboseDumpFiles)
        WriteCosts(m_cost, "reprojected/DSIa_%03d.pgm");
}
