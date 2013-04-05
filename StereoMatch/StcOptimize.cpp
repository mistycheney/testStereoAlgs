///////////////////////////////////////////////////////////////////////////
//
// NAME
//  StcOptimize.cpp -- select the best matches using local or global optimization
//
// DESCRIPTION
//  Input:  m_cost:      DSI (disparity space image)
//  Output: m_disparity  best disparities in range [0 .. n_disp]
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
#include "ImageIO.h"
#include <math.h>


void CStereoMatcher::OptWTA()
{
    // simple minimization (winner takes all)

    CShape sh = m_cost.Shape();
    int width = sh.width, height = sh.height;

    for (int y = 0; y < height; y++)
    {
        float* cost = &m_cost.Pixel(0, y, 0);
        int*   disp = &m_disparity.Pixel(0, y, 0);

        for (int x = 0; x < width; x++, cost += m_disp_n)
        {
            int best_d = 0;
            float best_cost = cost[0];
   
            for (int d = 1; d < m_disp_n; d++)
            {
                float current_cost = cost[d];

                if (current_cost < best_cost) 
                {
                    best_cost = current_cost;
                    best_d = d;
                }
            }
            disp[x] = best_d;
        }
    }
}
 

// Smoothness cost computation used by global optimization algorithms (DP, SO, GC, SA)

static float ComputeNCost(uchar I0[], uchar I1[], int nB, float opt_smoothness,
                          float opt_grad_thresh, float opt_grad_penalty, int& idIa)
{
    float dI2 = 0.0;
    for (int b = 0; b < nB; b++)
    {
        float dI = int(I0[b]) - int(I1[b]);
        dI2 += dI * dI;
    }
    dI2 /= (nB - (nB > 1));     // normalize by color channels (ignore A)
    float dIa = sqrt(dI2);
    idIa = int(dIa + 0.5);      // returned for histogram only!

    //    float s = opt_smoothness / (grad_term * dIa + 1.0f);  // old formula

    // New formula (same as used by Olga Veksler)
    float s = opt_smoothness;
    if (dIa < opt_grad_thresh)
        s *= opt_grad_penalty;
    
    return s;
}

void CStereoMatcher::ComputeSmoothnessCosts()
{
    // Set up the smoothness cost function for global optimization algorithms
    CShape sh = m_cost.Shape();
    sh.nBands = 2;      // vertical and horizontal smoothness costs (N-4)
    m_smooth.ReAllocate(sh, false);
    int H = sh.height;
    int W = sh.width;

    // Compute the histogram of gradients (analysis only)
    static bool compute_gradient_histogram = false;  // little overhead to do this...
    bool compute_hist = compute_gradient_histogram && 
        m_float_disparity.Shape() == m_true_disparity.Shape();
    int g_hist[2][256], i;
    for (i = 0; i < 256; i++)
        g_hist[0][i] = g_hist[1][i] = 1;    // minimum sampling

    // Fill in the values, using the reference image for gradients
    int nB = m_reference.Shape().nBands;
    for (int y = 0; y < H; y++)
    {
        uchar *I0 = &m_reference.Pixel(0, y, 0);
        uchar *I1 = &m_reference.Pixel(0, y+(y<H-1), 0);
        float *s  = &m_smooth.   Pixel(0, y, 0);
        for (int x = 0; x < W; x++, I0 += nB, I1 += nB, s += 2)
        {
            // TODO:  write a better smoothness function...
            int dIv, dIh;
            s[0] = (y < H-1) ? ComputeNCost(&I0[0], &I1[0 ], nB, opt_smoothness,
                                    opt_grad_thresh, opt_grad_penalty, dIv) : 0;
            s[1] = (x < W-1) ? ComputeNCost(&I0[0], &I0[nB], nB, opt_smoothness,
                                    opt_grad_thresh, opt_grad_penalty, dIh) : 0;

            // Increment the histogram
            if (compute_hist && y < H-1)
            {
                float d_diff = m_true_disparity.Pixel(x, y, 0) -
                               m_true_disparity.Pixel(x, y+1, 0);
                g_hist[fabs(d_diff) > eval_disp_gap][dIv]++;
            }
            if (compute_hist && x < W-1)
            {
                float d_diff = m_true_disparity.Pixel(x, y, 0) -
                               m_true_disparity.Pixel(x+1, y, 0);
                g_hist[fabs(d_diff) > eval_disp_gap][dIh]++;
            }
        }
    }

    // Write out smoothness cost maps for debugging
    static bool dump_smoothness = false;     // reset in debugger
    if (dump_smoothness && verbose >= eVerboseDumpFiles)
        WriteImage(m_cost, "reprojected/smoothness.pmf");

    // Write out the disparity histogram and posterior distribution
    if (compute_hist)
    {
        // Compute the log posterior distribution and print it out
        static char* log_file = "reprojected/tmp_gHist.txt";
        fprintf(stderr, " ComputeSmoothnessCosts: writing %s\n", log_file);
        FILE *stream = fopen(log_file, "w");
        fprintf(stream, "  dI\t  D=0\t  D=1\t  p(D)\tlog P(d)\n");
        for (int i = 0; i < 200; i++)
        {
            int tot_count = g_hist[0][i] + g_hist[1][i];
            float prob = g_hist[1][i] / (float) tot_count;
            float log_prob = log(prob);
            fprintf(stream, "%4d\t%7d\t%7d\t%7.4f\t%7.2f\n",
                i, g_hist[0][i], g_hist[1][i], prob, log_prob);
        }
        fclose(stream);
    }
}

//
// main dispatch function
//

void CStereoMatcher::Optimize()
{
    // Select the best matches using local or global optimization
    
    StartTiming();

    // set up the smoothness cost function for the methods that need it
    if (opt_fn == eDynamicProg || 
        opt_fn == eScanlineOpt || 
        opt_fn == eGraphCut ||
        opt_fn == eSimulAnnl ||
        opt_fn == eBPAccel ||
        opt_fn == eBPSync)
    {
        if (verbose == eVerboseSummary)
            fprintf(stderr, ", smooth=%g, grad_thres=%g, penalty=%g",
            opt_smoothness, opt_grad_thresh, opt_grad_penalty);
        if (verbose >= eVerboseProgress)
            fprintf(stderr, 
            "- computing smoothness costs (smoothness = %g, grad_thresh = %g, penalty=%g)\n",
            opt_smoothness, opt_grad_thresh, opt_grad_penalty);
        ComputeSmoothnessCosts();       
    }
    
    switch (opt_fn) {
    case eNoOpt:      // no optimization (pass through input depth maps)
        
        if (verbose == eVerboseSummary)
            fprintf(stderr, ", NO OPT");
        if (verbose >= eVerboseProgress)
            fprintf(stderr, "- optimizing: NONE\n");

        break;

    case eWTA:        // winner-take-all (local minimum)
        
        if (verbose == eVerboseSummary)
            fprintf(stderr, ", WTA");
        if (verbose >= eVerboseProgress)
            fprintf(stderr, "- optimizing: WTA (local minimum)\n");

        OptWTA();

        break;

    case eGraphCut:     // graph-cut global minimization

        if (verbose == eVerboseSummary)
            fprintf(stderr, ", GC");
        if (verbose >= eVerboseProgress)
            fprintf(stderr, "- optimizing: Graph Cut (near global minimum)\n");

        OptWTA();       // get an initial labelling (or just set to 0???)
        OptGraphCut();  // run the optimization

        break;

    case eDynamicProg:  // scanline dynamic programming
        
        if (verbose == eVerboseSummary)
            fprintf(stderr, ", DP (occl_cost=%d)", opt_occlusion_cost);
        if (verbose >= eVerboseProgress)
            fprintf(stderr, "- optimizing: Dynamic Programming (occlusion_cost = %d)\n",
            opt_occlusion_cost);

        OptDP();        // see StcOptDP.cpp

        break;

    case eScanlineOpt:  // scanline optimization
        
        if (verbose == eVerboseSummary)
            fprintf(stderr, ", SO");
        if (verbose >= eVerboseProgress)
            fprintf(stderr, "- optimizing: Scanline Optimization\n");

        OptSO();       // see StcOptSO.cpp

        break;

    case eSimulAnnl:  // simulated annealing

        if (verbose == eVerboseSummary)
            fprintf(stderr, ", SA");
        if (verbose >= eVerboseProgress)
            fprintf(stderr, "- optimizing: Simulated Annealing\n");

        OptWTA();           // initialize to reasonable starting point
                                // (for low-T gradient descent)
        OptSimulAnnl();    // see StcSimulAnn.cpp
        break;
 
    case eSymmetric:    // symmetric winner-take-all
                        //   *NEW RESEARCH CODE, NOT CURRENTLY DISTRIBUTED*
        
        if (verbose == eVerboseSummary)
            fprintf(stderr, ", Symm (margin=%g)", opt_min_margin);
        if (verbose >= eVerboseProgress)
            fprintf(stderr, "- optimizing: symmetric matching (margin=%g)\n", opt_min_margin);

        OptSymmetric();     // see StcOptSym.cpp

        break;

    case eBPAccel:
        OptBP();  // run the optimization
        break;

    case eBPSync:
        OptBPSync();  // run the optimization
        break;

    default:
        throw CError("CStereoMatcher::Optimize(): unknown optimization function");

    }
    PrintTiming();

    static bool always_compute_final_energy = true;
    if (final_energy < 0.0f && always_compute_final_energy &&
        !evaluate_only) {    // data costs not valid if evaluate_only
        if (! m_cost.Shape().SameIgnoringNBands(m_smooth.Shape()))
            ComputeSmoothnessCosts(); 
        float finalEd, finalEn;
        CStereoMatcher::ComputeEnergy(finalEd, finalEn);
        final_energy = finalEd + finalEn;
    }

}

#ifndef OPT_SYMMETRIC
void CStereoMatcher::OptSymmetric()
{
    throw CError("Optimize(eSymmetric) not currently implmented");
}
#endif

