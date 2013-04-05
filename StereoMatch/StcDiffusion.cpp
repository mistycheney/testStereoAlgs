///////////////////////////////////////////////////////////////////////////
//
// NAME
//  StcDiffusion.cpp -- non-linear diffusion
//
// DESCRIPTION
//  Implementation of three diffusion algorithms described
//  in "Stereo Matching with Nonlinear Diffusion", Scharstein & Szeliski,
//  IJCV 28(2), 1998.  Each needs to be run for several iterations.
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
#include "Convert.h"
#include <math.h>
#include <assert.h>

///////////////////////////////////////////////////////////////////////////
// diffusion functions

void CStereoMatcher::AggrDiffusion(int iter)
{
    // single iteration of regular diffusion or membrane model based on parameters
    // diff_lambda and diff_beta (lambda and beta in the following discussion)
    
    // if beta==0, we get regular diffusion
    
    // DSI[x, y, d] :=   (1 - 4*lambda) * DSI[x, y, d]
    //                 + lambda * sum_N4{DSI[x', y', d]}
    
    // if beta > 0, we get the membrane model

    // DSI[x, y, d] :=   (1 - lambda*(beta+4)) * DSI[x, y, d]
    //                 + lambda * beta * DSI_0[x, y, d]
    //                 + lambda * sum_N4{DSI[x', y', d]}
    
    // where DSI_0 is the initial (unaggregate DSI)
    
    // here, space for DSI_0 is allocated the first time this method is called,
    // but only if it is needed (e.g. if beta > 0)

    CShape sh = m_cost.Shape();
    int width = sh.width, height = sh.height;

    if (diff_beta > 0.0 && iter == 0) {
        // allocate and initialize DSI_0 in first iteration

        ScaleAndOffset(m_cost, m_cost0, 1.0, 0.0);
    }

    m_cost2.ReAllocate(sh);     // allocate memory for copy of m_cost

    // swap m_cost and m_cost2:
    CFloatImage m_cost_tmp = m_cost2;
    m_cost2 = m_cost;
    m_cost = m_cost_tmp;

    // aggregate each row, all disparity levels in parallel

    // boundaries: simply use copy of border element for now
    // TODO:  use borderMode

    int dx = &m_cost2.Pixel(1, 0, 0) - &m_cost2.Pixel(0, 0, 0);
    assert(dx == m_disp_n);
    int dy = &m_cost2.Pixel(0, 1, 0) - &m_cost2.Pixel(0, 0, 0);

    int y;
    for (y = 1; y < height-1; y++)
    {
        float* src = &m_cost2.Pixel(0, y, 0);
        float* dst = &m_cost.Pixel(0, y, 0);
        
        // left border
        int x;
        for (x = 0; x < dx; x++)
        {
            dst[x] = (1 - diff_lambda*(diff_beta+4)) * src[x]
                    + diff_lambda * (src[x] + src[x+dx] + src[x-dy] + src[x+dy]);
        }

        // center
        for (x = dx; x < dx*(width-1); x++)
        {
            dst[x] = (1 - diff_lambda*(diff_beta+4)) * src[x]
                    + diff_lambda * (src[x-dx] + src[x+dx] + src[x-dy] + src[x+dy]);
        }

        // right border
        for (x = dx*(width-1); x < dx*width; x++)
        {
            dst[x] = (1 - diff_lambda*(diff_beta+4)) * src[x]
                    + diff_lambda * (src[x-dx] + src[x] + src[x-dy] + src[x+dy]);
        }

    }

    // first and last row:
    for (y = 0; y < height; y += height-1)
    {
        float* src = &m_cost2.Pixel(0, y, 0);
        float* dst = &m_cost.Pixel(0, y, 0);
        
        for (int x = 0; x < dx*width; x++) 
        {
            int left  = (x-dx >= 0 ?       x-dx : x);
            int right = (x+dx < dx*width ? x+dx : x);
            int up    = (y > 0 ?           x-dy : x);
            int down  = (y < height-1 ?    x+dy : x);

            dst[x] = (1 - diff_lambda*(diff_beta+4)) * src[x]
                    + diff_lambda * (src[left] + src[right] + src[up] + src[down]);
        }
    }

    if (diff_beta > 0.0) {
        // add beta term to each element
            
        for (y = 0; y < height; y ++)
        {
            float* src = &m_cost0.Pixel(0, y, 0);
            float* dst = &m_cost.Pixel(0, y, 0);
            
            for (int x = 0; x < dx*width; x++) 
            {
                dst[x] += diff_lambda*diff_beta * src[x];
            }
        }
        
    }
}

void CStereoMatcher::AggrBayesian(int iter)
{
    // Perform a single iteration of Bayesian diffusion using parameters
    //   diff_scale_cost  - scale for m_cost values
    //   diff_mu          - speed of diffusion
    //   diff_sigmaP      - sigma of robust prior
    //   diff_epsP        - epsilon of robust prior
    // (mu, sigmaP, and epsP in discussion below)

    // One iteration consists of 4 steps:
    
    //  1. Probability
    //      P(i,j,d) := 1/Z * exp( - E(i,j,d) )
    //
    //      with normalizing Z(i,j) such that
    //      sum[d](P(i,j,d)) = 1    for all (i,j)
    //
    //  2. "smoothed" Probability
    //      Ps(i,j,d) := sum[d'](P(i,j,d')*w(d,d'))
    //
    //      with w(d,d') = rhoP(d - d')
    //      and rhoP(x) = (1-epsP) * exp(- x^2 / (2 * sigmaP^2)) + epsP
    //
    //  3. "smoothed" Energy
    //      Es(i,j,d) := - ln Ps(i,j,d)
    //
    //  4. updated Energy
    //      E(i,j,d) := E0(i,j,d) + mu * [ Es(i,j,d) +
    //        + Es(i+1,j,d) + Es(i-1,j,d) + Es(i,j+1,d) + Es(i,j-1,d) ]
    //
    
    // we use m_cost to hold the current energies E, m_cost2 to hold the smoothed
    // energies Es, and m_cost0 to hold the original energies E0

    // NOTE: the paper uses a contaminated Gaussian as measurement term.
    // this can be simulated with a truncated quadratic, but overall SCALING
    // of the values is still necessary.  To simulate the values given in the
    // paper for real images, epsM = 0.1 and sigmaM = 5, use parameters
    // match_fn 2 (SSD), match_max 12 (which will get squared to 144), and
    // cost scale factor diff_scale_cost = 0.016 (see also approx-rho.xls):

    CShape sh = m_cost.Shape();

    if (iter == 0) {
        // fprintf(stderr, "initializing m_cost0\n");
        // scale m_cost properly the first time
        ScaleAndOffset(m_cost, m_cost, diff_scale_cost, 0.0);     
        // allocate and initialize E0
        ScaleAndOffset(m_cost, m_cost0, 1.0, 0.0);  

        if (sh != m_cost2.Shape())
            m_cost2.ReAllocate(sh);         // allocate memory for Es
    }

    int width = sh.width, height = sh.height;
    int x, y, d;

    // TODO:
    // need to think of a way of not having to recompute weights every time, and
    // need to get rid of MAXDISP - use a vector for p and an image (?) for w.

#define MAXDISP 100

    assert(m_disp_n < MAXDISP);

    double p[MAXDISP];
    double w[MAXDISP][MAXDISP];

    // compute weights
    if (verbose>eVerboseAllMessages) fprintf(stderr, "weights:\n");
    for (d = 0; d < m_disp_n; d++) {
        float s=0.0;
        int d2;
        for (d2 = 0; d2 < m_disp_n; d2++) {
            float diff = (d-d2);
            s += (w[d][d2] = ((1-diff_epsP)*exp(-diff*diff/(2*diff_sigmaP*diff_sigmaP)) 
                               + diff_epsP));
        }
        for (d2 = 0; d2 < m_disp_n; d2++) {
            w[d][d2] /= s;
            if (verbose>eVerboseAllMessages) fprintf(stderr, "%1.2f ",w[d][d2]);
        }
        if (verbose>eVerboseAllMessages) fprintf(stderr, "\n");
    }

    // run one iteration of Bayesian diffusion

    for (y = 0; y < height; y ++)
    {
        float* pE = &m_cost.Pixel(0, y, 0);
        float* pEs = &m_cost2.Pixel(0, y, 0);
        
        for (x = 0; x < width; x++, pE += m_disp_n, pEs += m_disp_n) 
        {
            // step 1: convert to probabilities and normalize
            double s = 0.;
            for (d = 0; d < m_disp_n; d++) {
                s += (p[d] = exp( - pE[d]));
            }
            if (s==0) {
                for (d = 0; d < m_disp_n; d++)
                    p[d] = 1.0/m_disp_n;
            } else {
                for (d = 0; d < m_disp_n; d++)
                    p[d] /= s;
            }
            
            // step 2: smooth probabilities
            for (d = 0; d < m_disp_n; d++) {
                double ps = 0;
                for (int d2 = 0; d2 < m_disp_n; d2++)
                    ps += (w[d][d2] * p[d2]);

                // step 3: convert back to energies
                pEs[d] = - log(__max(1e-16, ps));
            }
        }
    }


    // step 4: diffuse smoothed energies

    int dx = &m_cost2.Pixel(1, 0, 0) - &m_cost2.Pixel(0, 0, 0);
    assert(dx == m_disp_n);
    int dy = &m_cost2.Pixel(0, 1, 0) - &m_cost2.Pixel(0, 0, 0);

    for (y = 1; y < height-1; y++)
    {
        float* pE  = &m_cost.Pixel(0, y, 0);
        float* pEs = &m_cost2.Pixel(0, y, 0);
        float* pE0 = &m_cost0.Pixel(0, y, 0);
        
        // left border
        for (x = 0; x < dx; x++)
        {
            pE[x] = pE0[x] + diff_mu * (pEs[x] + pEs[x   ] + pEs[x+dx] + pEs[x-dy] + pEs[x+dy]);
        }

        // center
        for (x = dx; x < dx*(width-1); x++)
        {
            pE[x] = pE0[x] + diff_mu * (pEs[x] + pEs[x-dx] + pEs[x+dx] + pEs[x-dy] + pEs[x+dy]);
        }

        // right border
        for (x = dx*(width-1); x < dx*width; x++)
        {
            pE[x] = pE0[x] + diff_mu * (pEs[x] + pEs[x-dx] + pEs[x   ] + pEs[x-dy] + pEs[x+dy]);
        }

    }

    // first and last row:
    for (y = 0; y < height; y += height-1)
    {
        float* pE  = &m_cost.Pixel(0, y, 0);
        float* pEs = &m_cost2.Pixel(0, y, 0);
        float* pE0 = &m_cost0.Pixel(0, y, 0);
        
        for (int x = 0; x < dx*width; x++) 
        {
            int left  = (x-dx >= 0 ?       x-dx : x);
            int right = (x+dx < dx*width ? x+dx : x);
            int up    = (y > 0 ?           x-dy : x);
            int down  = (y < height-1 ?    x+dy : x);

            pE[x] = pE0[x] + diff_mu * (pEs[x]+ pEs[left] + pEs[right] + pEs[up] + pEs[down]);
        }
    }
}
