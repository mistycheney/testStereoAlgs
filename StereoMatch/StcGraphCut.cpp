///////////////////////////////////////////////////////////////////////////
//
// NAME
//  StcGraphCut.cpp -- optimize the MRF using a Graph Cut algorithm
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
extern "C"{
#include "maxflow/maxflow.h"
}
#include <algorithm>
#include <math.h>

// TODO:  this should be a member variable,
//  computed in ComputeEnergy based on the (1 << 30) / (largest data/neighbor cost)
static float GC_scale = (1 << 30) / (256 * 256);  // hope opt_smoothness < 255^2...


void CStereoMatcher::ComputeEnergy(CFloatImage& dcost, CFloatImage& ncost, 
                   CIntImage& label, float& denergy, float& nenergy)
{
    // This version takes all its parameters explicitly
    //  It can therefore be called to evaluate a sub-region of the solution...
    CShape sh = dcost.Shape();
    int W = sh.width;
    int H = sh.height;
    int nD = sh.nBands;
    int nN = ncost.Shape().nBands;
    float dSum = 0.0f, nSum = 0.0f;

	// Compute the data term energy and smoothness energy 
    for (int y = 0; y < H; y++)
    {
        int *curLabel  = &label.Pixel(0, y, 0);
        int *nextLabel = &label.Pixel(0, y+(y<H-1), 0);
        float *dc = &dcost.Pixel(0, y, 0);
        float *nc = &ncost.Pixel(0, y, 0);

        for (int x = 0; x < W; x++, dc += nD, nc += nN)
        {
	        int mydisp = curLabel[x];
			dSum += dc[mydisp];

			if (y < H-1 && mydisp != nextLabel[x])    // vertical cost
				nSum += nc[0];
			if (x < W-1 && mydisp != curLabel[x+1])   // horizontal cost
				nSum += nc[1];
            if (nSum < 0.0)
                throw CError("ComputeEnergy: negative neighborhood costs");
		}
	}
	denergy = dSum;
    nenergy = nSum;

    // Compute a good scaling value, so that total integer cost is still
    //  within the range of a long.
    GC_scale = (1 << 30) / (dSum + nSum);
}

void CStereoMatcher::ComputeEnergy(float& denergy, float& nenergy)
{
    ComputeEnergy(m_cost, m_smooth, m_disparity, denergy, nenergy);
}

/////////////////////////////////////////////////////////////////////////////////////////
// This function is the basic function to implement the alpha and beta swap            //									
// It is used in the disparity labeling and plane labeling							   //
// INPUT: the input of algorithm												       //
//		- dcost: the data cost image whose band number depends on the label num        //
//		- ncost: the neighboring cost image with 4-bands if it is 4-connected neighbor //
//		- label: the updated label after we do the alpha-beta swap					   //
//		- alpha, beta: a pair of label to swap										   //
/////////////////////////////////////////////////////////////////////////////////////////

static inline int nodenum(int x, int y,int width)
{
  return ((y*width) + x + 2);
}


static int  SwapEnergyImprove(CFloatImage& dcost, CFloatImage& ncost,
                              CIntImage& label, int alpha, int beta)
{
#define LIVE(disp) (((disp) == alpha)||((disp) == beta) )

    CShape sh = dcost.Shape();
    int W = sh.width;
    int H = sh.height;
    int nD = sh.nBands;
    int nN = ncost.Shape().nBands;
    struct Graph *g = 0;
    int x, y;
    int source = 0, sink = 1;
    int n_nodes = 0;
    std::vector<int> Cut;
    Cut.resize(W*H+2);

    for (y = 0; y < H; y++)
    {
        int *curLabel  = &label.Pixel(0, y, 0);
        int *nextLabel = &label.Pixel(0, y+(y<(H-1)), 0);

        float *dc = &dcost.Pixel(0, y, 0);
        float *nc = &ncost.Pixel(0, y, 0);
        for (x = 0; x < W; x++, dc += nD, nc += nN)
        {
	        int mydisp = curLabel[x];
	        int mynode = nodenum(x,y,W);

	        if LIVE(mydisp)
	        {
		        if (n_nodes == 0)
                {
			        g = init_graph(source,sink);
                    if (g == 0)
                    {
                        // Maxflow code not downloaded/compiled properly
                        throw CError("CStereoMatcher::OptGraphCut():\
 maxflow code not implemented;\
 please see the README.txt file in the maxflow subdirectory");
                    }
                }
	        n_nodes++;

                // Add D-links
		        add_edge(g, source, mynode, (long int)(dc[alpha] * GC_scale));
		        add_edge(g, mynode, sink, (long int)(dc[beta] * GC_scale));

                // Add 4-connected N-links 
                if (y < H-1 && LIVE(nextLabel[x])) {
			        add_edge(g, mynode, nodenum(x,y+1,W), (long int)(nc[0] * GC_scale)); // North
			        add_edge(g, nodenum(x,y+1,W), mynode, (long int)(nc[0] * GC_scale)); // South
                }
                if (x < W-1 && LIVE(curLabel[x+1])) {
			        add_edge(g, mynode, nodenum(x+1,y,W), (long int)(nc[1] * GC_scale)); // West
			        add_edge(g, nodenum(x+1,y,W), mynode, (long int)(nc[1] * GC_scale)); // East
                }
	        }
        }
	}
	if (n_nodes == 0) return 0;

	long flow = maxflow(g,&Cut[0]);
	
	for (y = 0; y < H; y++)
	{
		int *curLabel = &label.Pixel(0, y, 0);
		for (x = 0; x < W; x++)
		{
			int mydisp = curLabel[x], mynode = nodenum(x,y,W);
			if ((mydisp == alpha) || (mydisp == beta))
			{
				if (Cut[mynode] == 1)		// Link to source is cut
					curLabel[x] = alpha;
				else						// Link to sink is cut
					curLabel[x] = beta;
			}
		}
	}
    return (int) flow;
}

void CStereoMatcher::DumpDisparity(CIntImage& disp, const char* filename, float scale)
{
    CByteImage bdisp(disp.Shape());
    ScaleAndOffset(disp, bdisp, scale, 0.0);
    WriteImage(bdisp, filename);
}

static int CycleAll(CFloatImage& dcost, CFloatImage& ncost,
                    CIntImage& label, int randomize_labels, 
                    EVerbosityLevel verbose, float& finalE)
{
    // Cycle through all possible alpha-beta swaps
    int numLabel = dcost.Shape().nBands;
    static bool randomize_pairings = true;
    int numTotal = (randomize_pairings) ? numLabel*numLabel : numLabel;

    // Make a list of labels
    std::vector<int> ranLabelList;
    ranLabelList.resize(numTotal);
    for (int l = 0; l < numTotal; l++)
        ranLabelList[l] = l;

    // Optionally randomize the labels
    if (randomize_labels)
        std::random_shuffle(ranLabelList.begin(), ranLabelList.end());

    // Compute the data energy and smoothness energy
    float oldEd, oldEn, oldE;
    CStereoMatcher::ComputeEnergy( dcost, ncost, label, oldEd, oldEn );
    oldE = oldEd + oldEn;
    int success = 0;

    // One loop of graph cut algorithm
    for (int label1 = 0; label1 < numTotal; label1++)
    {
        // The second loop is only used if we don't randomize_pairings
        for (int label2 = (randomize_pairings) ? numLabel-1 : label1+1;
             label2 < numLabel; label2++)
        {
			int alpha = ranLabelList[label1]; 
			int beta  = ranLabelList[label2];
            if (randomize_pairings)
            {
                int product = alpha;    // encoded value is product of alpha*beta
                alpha = product % numLabel;
                beta  = product / numLabel;
                if (alpha <= beta)
                    continue;           // wasted loop index, but who cares?
            }

            if (verbose >= eVerboseDumpFiles)
                CStereoMatcher::DumpDisparity(label, "reprojected/disp_before_swap.pgm", 8);
            SwapEnergyImprove(dcost,ncost,label,alpha,beta);
            if (verbose >= eVerboseDumpFiles)
                CStereoMatcher::DumpDisparity(label, "reprojected/disp_after_swap.pgm", 8);

			// Step 2.3.3: compute the energy after we do the swap operation
            float newEd, newEn, newE;
			CStereoMatcher::ComputeEnergy(dcost, ncost,label,newEd,newEn);
			newE = newEd + newEn;
			if (verbose >= eVerboseProgress
                && label1==0 && label2==1 // create less output 
                || verbose >= eVerboseInnerLoops
                )
			{
				fprintf(stderr, " graph cut: alpha=%d, beta=%d\n", alpha, beta);
                fprintf(stderr, "  old E = %.2f (%.2f + %.2f)\n", oldE, oldEd, oldEn);
                fprintf(stderr, "  new E = %.2f (%.2f + %.2f)\n", newE, newEd, newEn);

            }

			if (newE < oldE)
				success = 1;
			oldEd = newEd;
			oldEn = newEn;
			oldE  = newE;
            finalE = newE;
	    }
    }
    return success;
}

void CStereoMatcher::OptGraphCut()
{
    // Optimize using graph cuts
    for (int iter = 0, progress = 1; iter < opt_max_iter && progress; iter++)
    {
        progress = CycleAll(m_cost, m_smooth, m_disparity, 
                        opt_random, verbose, final_energy);
    }
    // ??? Should use m_cost0 (unaggregated cost) for Graph cut and Simulated Annealing??
    // Only aggregate for WTA for start values, but energy should be defined on original
    // costs?
}
