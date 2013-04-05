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
//This uses new BP code

#include "StereoMatcher.h"
#include "ImageIO.h"

#include "bpregions.h"
#include <math.h>
#include <assert.h>
#include <sys/timeb.h>

int OneNodeCluster::numStates = 2;
int TwoNodeCluster::numStates = 2;

void CStereoMatcher::OptBP()
{
  CShape sh = m_cost.Shape();

  int H = sh.height;
  int W = sh.width;

  const int numRows = H,
    numCols = W,
    numNodes = W *H;

  OneNodeCluster *nodeArray;
  const int numStates = m_disp_n;
    FLOATTYPE alpha = 0.95;

  unsigned int totalMem = 0;

  nodeArray = new OneNodeCluster[numNodes];
  FLOATTYPE *oneNodeMsgArray = new FLOATTYPE[numNodes * 8 * numStates];
  FLOATTYPE *oneNodeLocalEv = new FLOATTYPE[numStates * numNodes],
    *currLocalEvPtr = oneNodeLocalEv;

  totalMem += sizeof(OneNodeCluster) * numNodes;
  totalMem += sizeof(FLOATTYPE) * numNodes * 8 * numStates;
  totalMem += sizeof(FLOATTYPE) * numNodes * numStates;
  
  OneNodeCluster::numStates = numStates;
  TwoNodeCluster::numStates = numStates;
  printf("NumStates: %d\n",nodeArray[0].numStates);

  for (int i = 0; i < numNodes * 8 * numStates; i++)
  {
    oneNodeMsgArray[i] = 1.0f/numStates;
  }

  initOneNodeMsgMem(nodeArray,oneNodeMsgArray,numNodes, numStates);

  TwoNodeCluster dummyNode;
  FLOATTYPE dummyMessages[numStates * numStates*4];

  for(int i = 0; i < numStates * numStates*4; i++)
    dummyMessages[i] = 1.0f;

  
  
  printf("Msg Space Allocated and Initialized\n");

  CShape sh2 = m_cost.Shape();
  printf("%d %d\n",m_disp_n, sh2.nBands);
    for(int m = 0; m < numRows; m++)
    {
      float* smoothcost_vert = &m_smooth.Pixel(0, m, 0); // band 1 is horizontal cost
      float* smoothcost_horz = &m_smooth.Pixel(0, m, 1); // band 1 is horizontal cost
      float* local_cost = &m_cost.Pixel(0, m, 0);
    
      for(int n = 0; n < numCols; n++)
      {
	const double div_factor = 50;
        const double vsmooth_cost = exp(-1*smoothcost_vert[0]/div_factor),
	  hsmooth_cost = exp(-1*smoothcost_horz[0]/div_factor);
	

        const int cind = m * numCols + n;

        nodeArray[cind].localEv = currLocalEvPtr;
        currLocalEvPtr += numStates;
        for (int i = 0; i < m_disp_n; i++)
        {
	  nodeArray[cind].localEv[i] = exp(-1 *local_cost[i]/div_factor);
        }

        nodeArray[cind].psiData_pottsSameProb[OneNodeCluster::UP]=1;
        nodeArray[cind].psiData_pottsSameProb[OneNodeCluster::DOWN]=1;
        nodeArray[cind].psiData_pottsSameProb[OneNodeCluster::LEFT]=1;
        nodeArray[cind].psiData_pottsSameProb[OneNodeCluster::RIGHT]=1;

        if (n == 0)
	  nodeArray[cind].psiData_pottsDiffProb[OneNodeCluster::LEFT] = 1;

        if (n < numCols-1)
        {
	  nodeArray[cind+1].psiData_pottsDiffProb[OneNodeCluster::LEFT] = hsmooth_cost;
	  nodeArray[cind].psiData_pottsDiffProb[OneNodeCluster::RIGHT] = hsmooth_cost;
	}
        else
        {
	  nodeArray[cind].psiData_pottsDiffProb[OneNodeCluster::RIGHT]=1;
        }

        if (m == 0)
	  nodeArray[cind].psiData_pottsDiffProb[OneNodeCluster::UP] = 1;

        if(m != numRows -1)
	{	
	  nodeArray[cind+numCols].psiData_pottsDiffProb[OneNodeCluster::UP] = vsmooth_cost;
	  nodeArray[cind].psiData_pottsDiffProb[OneNodeCluster::DOWN] = vsmooth_cost;
	}
        else
        {
	  nodeArray[cind].psiData_pottsDiffProb[OneNodeCluster::DOWN]=1;
        }


        smoothcost_horz += 2;
        smoothcost_vert += 2;
        local_cost += m_disp_n;
      }
    }
  
    int numIter;
    //Initialization Completed Time for the good stuff
    
    numIter = 50;
    printf("Structures Initialized\n");
    
    int iter = 0;
    
    int beliefIm[H * W];
    for (int i = 0; i < H * W; i++)
      beliefIm[i] = 0;

    while (iter < numIter)//( !converged)
    {

      printf("Iter %03d\n",iter);
      struct timeb start, endtime;
      ftime(&start);
      for (int m = 0; m < numRows; m++)
      {
	passOneNodeMsgsLeft(&nodeArray[m *numCols], numCols, 0,0,dummyNode,alpha);
	passOneNodeMsgsRight(&nodeArray[m *numCols], numCols, 0,0,dummyNode,alpha);
      }
      
      for(int n = 0; n < numCols; n++)
      {
	passOneNodeMsgsDown(nodeArray,(TwoNodeCluster **)0,dummyNode,n,numRows,numCols,alpha);

	passOneNodeMsgsUp(nodeArray,(TwoNodeCluster **)0,dummyNode,n,numRows,numCols,alpha);
      }
      ftime(&endtime);
      float etime = (float)(endtime.time - start.time) + 0.001 * (float) (endtime.millitm - start.millitm);
      printf ("%f Seconds ",etime);

      float tmpBeliefVec[numStates];
      float *beliefPtr;
    
      for(int m = 0; m < numRows; m++)
      {
	int *disp = &m_disparity.Pixel(0, m, 0);    
      
	for(int n = 0; n < numCols; n++)
	{
	  double best_prob = 0;
	  int best_d = 0;
	  
	  nodeArray[m * numCols + n].getBelief(tmpBeliefVec);
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
      iter++;
    }

    float tmpBeliefVec[numStates];
    float *beliefPtr;
    
    for(int m = 0; m < numRows; m++)
    {
      int *disp = &m_disparity.Pixel(0, m, 0);    
      
      for(int n = 0; n < numCols; n++)
      {
	double best_prob = 0;
	int best_d = 0;
	
	nodeArray[m * numCols + n].getBelief(tmpBeliefVec);
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
    
  
    
    printf("\nTotalMem:%ud\n",totalMem/1024/1024);
    delete nodeArray;
    delete oneNodeMsgArray;
    delete oneNodeLocalEv;

}

