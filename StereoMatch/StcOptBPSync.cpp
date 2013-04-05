// Copyright 2002,2003 Marshall Tappen
//  This program is free software; you can redistribute it and/or modify
//     it under the terms of the GNU General Public License as published by
//     the Free Software Foundation; either version 2 of the License, or
//     (at your option) any later version.

//     This program is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details.

//     You should have received a copy of the GNU General Public License
//     along with this program; if not, write to the Free Software
//     Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
// Adds Belief Propagation to Stereo Matcher Framework 
// Code by Marshall Tappen, MIT AI Lab- 2002
#include "StereoMatcher.h"
#include "ImageIO.h"

#include "node.h"
#include <math.h>
#include <assert.h>

extern double *globalPsiMat;
// char *output_filename;

void CStereoMatcher::OptBPSync()
{
  CShape sh = m_cost.Shape();
  int H = sh.height;
  int W = sh.width;
  const int numRows = H,
    numCols = W,
    numNodes = W *H;
  const int numThingsToMalloc = 9;
  unsigned long locs[numThingsToMalloc];
  Node *nodeArray;
  const int numStates = m_disp_n;
  int maxNumNeighbors = 4;
  double alpha = 0.8;
  printf("H:%d W:%d\n",H,W);
  locs[0] = 0;
  locs[1] = maxNumNeighbors * sizeof(int); /*  //Neighbors */
  locs[2] = locs[1] +numStates* sizeof(double);/*  //localEvidence */
  locs[3] = locs[2] + maxNumNeighbors * sizeof(double *);/*  // psiMats */
  locs[4] = locs[3] + maxNumNeighbors * sizeof(double *); /*  //prevMessages */
  locs[5] = locs[4] + maxNumNeighbors * sizeof(double *);/*  //currMessages */
  locs[6] = locs[5] + maxNumNeighbors * 2 * sizeof(double); /*  // psiMats data */
  locs[7] = locs[6] + maxNumNeighbors * numStates * sizeof(double); /*  //prevMessages data */
  locs[8] = locs[7] + maxNumNeighbors * numStates * sizeof(double);/*  //currMessages data */
  
  nodeArray = new Node[numNodes];
  for (int i = 0; i < numNodes; i++)
  {
    
    unsigned char *memChunk = new unsigned char[locs[8]+10]; /*  // Eliminate this malloc later to be really clever */
    nodeArray[i].neighbors = (int *)(memChunk + locs[0]);
    nodeArray[i].localEvidence = (double *)(memChunk + locs[1]);
    //    nodeArray[i].psiMats = (double **)(memChunk + locs[2]);
    nodeArray[i].pottsProbs = (double**)(memChunk + locs[2]);
    nodeArray[i].prevMessages = (double **)(memChunk + locs[3]);
    nodeArray[i].currMessages = (double **)(memChunk + locs[4]);

    nodeArray[i].numNeighbors = 0;
    nodeArray[i].alpha = alpha;
    nodeArray[i].one_minus_alpha = 1 - alpha;
    nodeArray[i].numStates = numStates;
    nodeArray[i].nodeArray = nodeArray;
    nodeArray[i].maxNumNeighbors = maxNumNeighbors;

    for(int j = 0; j < maxNumNeighbors; j++)
    {
      nodeArray[i].pottsProbs[j] = (double *)(memChunk + locs[5] + (j * 2 * sizeof(double)));
      nodeArray[i].prevMessages[j] =  (double *)(memChunk + locs[6] + j * numStates*sizeof(double));
      nodeArray[i].currMessages[j] =  (double *)(memChunk + locs[7] + j * numStates * sizeof(double));
      for(int k = 0; k < numStates;k++)
	nodeArray[i].currMessages[j][k]=nodeArray[i].prevMessages[j][k]=1;
    }
  }

  // This section of code implements something similar to the potential function used
  // in ECCV-2002 "Stereo Matching with Belief Propagation"
//   globalPsiMat = new double[numStates * numStates];
//   const double ep = 0.05, sigmap = 0.6;
//   for(int i = 0; i < numStates; i++)
//     for(int j = 0; j < numStates; j++)
//     {
//       double di = i;
//       double dj = j;
//       globalPsiMat[i * numStates + j] = (1 - ep) * exp(-1.0 * fabs(di -dj)/sigmap) + ep;
//     }

  printf("Memory Allocated\n");
  for(int m = 0; m < numRows; m++)
  {
    float* smoothcost_vert = &m_smooth.Pixel(0, m, 0); // band 1 is horizontal cost
    float* smoothcost_horz = &m_smooth.Pixel(0, m, 1); // band 1 is horizontal cost
    float* local_cost = &m_cost.Pixel(0, m, 0);
    
    for(int n = 0; n < numCols; n++)
    {
      int on,om;
      int psi_index[2];
      double psi_tensor[2];

      // The constants here only affect the problem is you are using sum-product
      // to compute marginals.

      const double vsmooth_cost = exp(-1*smoothcost_vert[0]/50.0),
	hsmooth_cost = exp(-1*smoothcost_horz[0]/50.0);

      const int cind = m * numCols + n;

      for (int i = 0; i < m_disp_n; i++)
      {
	nodeArray[cind].localEvidence[i] = exp(-1 *local_cost[i]/50.0);
      }

      if (n != numCols-1)
      {
	om = m;
	on = n+1;

	psi_index[0] =m*numCols +n;
	psi_index[1] = om*numCols+on;

	psi_tensor[0] = 1;
	psi_tensor[1] = hsmooth_cost;

	
	addNeighbor(&nodeArray[m*numCols+n],psi_index,psi_tensor,1);
	psi_index[1] =m*numCols +n;
	psi_index[0] = om*numCols+on;

	addNeighbor(&nodeArray[om*numCols+on],psi_index,psi_tensor,1);

      }

      if(m != numRows-1)
      {
	om = m+1;
	on = n;

	psi_index[0] =m*numCols +n;
	psi_index[1] = om*numCols+on;
	
	psi_tensor[0] = 1;
	psi_tensor[1] = vsmooth_cost;

	addNeighbor(&nodeArray[m*numCols+n],psi_index,psi_tensor,1);
	psi_index[1] =m*numCols +n;
	psi_index[0] = om*numCols+on;

	addNeighbor(&nodeArray[om*numCols+on],psi_index,psi_tensor,1);

      }

      smoothcost_horz += 2;
      smoothcost_vert += 2;
      local_cost += m_disp_n;
    }
  }
  

  int numIter;
  if (H>W)
    numIter = H;
  else numIter = W;
  //  numIter = 5;
  //Initialization Completed Time for the good stuff
  for (int iter =0; iter < numIter; iter++)
  {

    printf("Iteration %d / %d\n", iter,numIter-1);

    for (int j = 0; j < numNodes; j++)
      doIteration(&nodeArray[j]);
    for (int j = 0; j < numNodes; j++)
      finishIteration(&nodeArray[j]); /* You can call get Belief any time, 
				     but you make sure you do a finishIteration first */

    double tmpBeliefVec[numStates];
    double *beliefPtr;

    for(int m = 0; m < numRows; m++)
    {
      int *disp = &m_disparity.Pixel(0, m, 0);    
      
      for(int n = 0; n < numCols; n++)
      {
	double best_prob = 0;
	int best_d = 0;
	
	getBelief(&nodeArray[m * numCols + n],tmpBeliefVec);
	beliefPtr =tmpBeliefVec;
	for(int j =0; j <numStates;j++)
	{
	  
	  if (*beliefPtr > best_prob)
	  {
	    best_prob = *beliefPtr;
	    best_d = j;
	  }
	  beliefPtr++;
	}
	
	disp[n] = best_d;
      }
    }
    float locEn, smoothEn;
    ComputeEnergy(m_cost, m_smooth, m_disparity, locEn, smoothEn);
    printf("Energy:%.2f Local %.2f Smooth: %.2f\n",locEn+smoothEn, locEn,smoothEn);
    
  }
  double tmpBeliefVec[numStates];
  double *beliefPtr;
  
  for(int m = 0; m < numRows; m++)
  {
    int *disp = &m_disparity.Pixel(0, m, 0);    
    
    for(int n = 0; n < numCols; n++)
    {
      double best_prob = 0;
      int best_d = 0;
      
      getBelief(&nodeArray[m * numCols + n],tmpBeliefVec);
      beliefPtr =tmpBeliefVec;
      //      beliefPtr = nodeArray[m * numCols + n].localEvidence;
      for(int j =0; j <numStates;j++)
      {
	
	if (*beliefPtr > best_prob)
	{
	  best_prob = *beliefPtr;
	  best_d = j;
	}
	beliefPtr++;
      }
      
      disp[n] = best_d;
    }
  }


  
  for(int i =0; i < numNodes; i++)
    delete nodeArray[i].neighbors;

  delete nodeArray;
  delete globalPsiMat;
}
