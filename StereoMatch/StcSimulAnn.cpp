///////////////////////////////////////////////////////////////////////////
//
// NAME
//  StcSimulAnn.cpp -- optimize the MRF using Simulated Annealing
//
// DESCRIPTION
//  Input:  m_cost:      DSI (disparity space image)
//          m_smooth:    smoothness field weights
//          m_rho_s:     table of smoothness penalties vs. disp difference
//  Output: m_disparity  best disparities in range [0 .. n_disp]
//
// SEE ALSO
//  StereoMatcher.h
//  StereoMatcher.cpp
//
// Copyright © Richard Szeliski, 2001.
// See Copyright.h for more details
//
///////////////////////////////////////////////////////////////////////////

#include "Error.h"
#include "StereoMatcher.h"
#include "Convert.h"
#include "ImageIO.h"
#include <algorithm>
#include <math.h>

static float ComputeEnergySingle(CFloatImage& dcost, CFloatImage& ncost,
                                 CIntImage& label,
                                 int x, int y, int mydisp)
{
    // See ComputeEnergy() in StcGraphCut for definition of global energy
    CShape sh = dcost.Shape();
    int W = sh.width;
    int H = sh.height;

	// Compute the data term energy and smoothness energy
    float dSum = dcost.Pixel(x, y, mydisp);
    float nSum = 0.0f;
	if (y > 0   && mydisp != label.Pixel(x, y-1, 0))    // vertical cost
        nSum += ncost.Pixel(x, y-1, 0);
	if (y < H-1 && mydisp != label.Pixel(x, y+1, 0))    // vertical cost
        nSum += ncost.Pixel(x, y  , 0);

    if (x > 0   && mydisp != label.Pixel(x-1, y, 0))    // horizontal cost
        nSum += ncost.Pixel(x-1, y, 1);
    if (x < W-1 && mydisp != label.Pixel(x+1, y, 0))    // horizontal cost
        nSum += ncost.Pixel(x  , y, 1);
    // assert(nSum >= 0.0);
    float tSum = dSum + nSum;
    return tSum;
}

static float urand()
{
    // Uniform random number in [0,1]
    int i = rand();
    const float scale = 1.0f / RAND_MAX;
    float r = scale * i;
    return r;
}

static int SACycle(CFloatImage& dcost, CFloatImage& ncost,
                   CIntImage& label, float kT_inv, EStereoSAVariant sampler,
                   int randomize_pixels, int cycle, int num_cycles,
                   EVerbosityLevel verbose, float& finalE)
{
    // Visit all possible sites
    CShape sh = dcost.Shape();
    int numLabel = sh.nBands;
    int numPixel = sh.width * sh.height;
    if (sh.width > (1 << 14) || sh.height > (1 << 14))
        throw CError("SimulAnnealCycle: image width or height is too large");

    // Make a list of pixels (sites)
    std::vector<int> ranPixelList;
    ranPixelList.resize(numPixel);
    int x, y, l;
    for (y = 0, l = 0; y < sh.height; y++) {
        for (x = 0; x < sh.width; x++, l++)
            ranPixelList[l] = (x << 16) | (y << 0);
    }

    // TODO: should we try to be consistent with local variable names
    //  i.e., either use underscores or mixed case?

    // Make a list of disparities, energies, and probabilities
    int nCand = (sampler == eFullGibbs) ? numLabel : 2;
    std::vector<int> dList;             // possible disparities
    std::vector<float> eList, pList;    // resulting energies and probabilities
    dList.resize(nCand);
    eList.resize(nCand);
    pList.resize(nCand);

    // Optionally randomize the sites
    if (randomize_pixels)
        std::random_shuffle(ranPixelList.begin(), ranPixelList.end());

    // Compute the data energy and smoothness energy
    float oldEd, oldEn, oldE, sumDeltaE = 0.0f;
    CStereoMatcher::ComputeEnergy( dcost, ncost, label, oldEd, oldEn );
    oldE = oldEd + oldEn;
    float min_valid_E = log(FLT_MIN) + 1.0;     // to prevent exp() underflow
    if (verbose >= eVerboseDumpFiles)
        CStereoMatcher::DumpDisparity(label, "reprojected/disp_before_SA.pgm", 8);

    // One loop of simulated annealing
    for (l = 0; l < numPixel; l++)
    {
        int k = ranPixelList[l];
        x = k >> 16;
        y = k & ((1 << 16) - 1);
        int dOld = label.Pixel(x, y, 0);

        // Choose which labels to evaluate
        // Make up a list of disparities
        if (sampler == eFullGibbs)
        {
            for (int d = 0; d < numLabel; d++)
                dList[d] = d;
        }
        else
        {
            dList[0] = dOld;
            int ran1 = rand() % (numLabel - 1);
            int ran2 = (dOld + ran1 + 1) % numLabel;    // not current
            dList[1] = ran2;
        }

        // Compute all of the energies and find minimum
        int d, d_picked = 0;
        float minE = FLT_MAX;
        for (d = 0; d < nCand; d++)
        {
            eList[d] = ComputeEnergySingle(dcost, ncost, label, x, y, dList[d]);
            minE = __min(minE, eList[d]);
        }

        // Pick the new state
        if (sampler == eMetropolis)
        {
            // Metropolis:  always go downhill, sometime go up
            if (eList[1] < eList[0])
                d_picked = 1;
            else
            {
                float EUpHill = kT_inv * (eList[1] - eList[0]);
                float pUpHill = (-EUpHill < min_valid_E) ? 0.0 : exp(-EUpHill);
                float r = urand();
                d_picked = (r <= pUpHill);
            }
        }
        else
        {
            // Gibbs Sampler:  compute cumulative distribution
            float pSum = 0.0;
            for (d = 0; d < nCand; d++)
            {
                float deltaE = kT_inv * (eList[d] - minE);
                float pD = (-deltaE < min_valid_E) ? 0.0 : exp(-deltaE);
                pSum += pD;
                pList[d] = pSum;
            }

            // Pick a disparity based on the cumulative p.d.f.
            float r = urand() * pSum;
            for (d = 0; d < nCand; d++)
            {
                if (r <= pList[d] && pList[d] > 0.0)
                {
                    d_picked = d;
                    break;
                }
            }
        }

        // Update the current site
        int dNew = dList[d_picked];
        label.Pixel(x, y, 0) = dNew;

        // Update the energy change (sanity check)
        int d_orig = (sampler == eFullGibbs) ? dOld : 0;
        float deltaE = (eList[d_picked] - eList[d_orig]);
        if (deltaE > 0.0)
            deltaE = deltaE;    // debugging breakpoint
        sumDeltaE += deltaE;
    }


	// Ccompute the energy after all sites have been visited
    float newEd, newEn, newE;
	CStereoMatcher::ComputeEnergy(dcost, ncost,label,newEd,newEn);
	newE = newEd + newEn;
    finalE = newE;
    if (verbose >= eVerboseProgress)
    {
        if (cycle < 3 || cycle > num_cycles-3) // reduce output
        {
            fprintf(stderr, " simulated annealing: cycle=%d, kT=%g\n", cycle, 1.0 / kT_inv);
            fprintf(stderr, "  old E = %.2f (%.2f + %.2f)\n", oldE, oldEd, oldEn);
            fprintf(stderr, "  new E = %.2f (%.2f + %.2f)\n", newE, newEd, newEn);
        } else {
            if (cycle == 3)
                fprintf(stderr, "...\n");
        }
    }
    if (verbose >= eVerboseDumpFiles)
    {
        char filename[1024];
        sprintf(filename, "reprojected/disp_after_SA_%03d.pgm", cycle);
        CStereoMatcher::DumpDisparity(label, filename, 8);
    }

	int success = (newE < oldE);
    return success;
}

void CStereoMatcher::OptSimulAnnl(void)
{
    // Set up the annealing schedule
    float kT = opt_sa_start_T, kT_delta;
    if (opt_sa_schedule == eSALinear)
    {
        kT_delta = (opt_sa_start_T - opt_sa_end_T) / (opt_max_iter - (opt_max_iter != 1));
    }
    else
    {
        throw CError("OptSimulAnnl: opt_sa_schedule = %d not yet implemented", opt_sa_schedule);
    }

    // Optimize using stochastic gradient descent
    for (int iter = 0; iter < opt_max_iter; iter++)
    {
        // Perform one cycle (visit all pixels / sites)
        int downhill= SACycle(m_cost, m_smooth, m_disparity, 1.0f/kT, opt_sa_var,
                              opt_random, iter, opt_max_iter, verbose, final_energy);
        if (verbose >= eVerboseAllMessages)
            putchar(downhill ? '-' : '+');      // show downhill-uphill progress

        // Update the temperature
        if (opt_sa_schedule == eSALinear)
            kT -= kT_delta;
        else
            kT = 0;     // not yet implemented
        kT = __max(kT, opt_sa_end_T);       // to prevent roundoff error
    }
}
