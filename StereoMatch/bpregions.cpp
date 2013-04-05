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

#define FLOATTYPE double
#include "bpregions.h"
#include <stdio.h>
int disabledClusters;
int numIterRun;
// SOme of the GBP code has been disabled here
char output_filename[10];
void getPsiMat(OneNodeCluster &cluster, FLOATTYPE *destMatrix, int direction)
{
  FLOATTYPE *currPtr = destMatrix;
  const int numStates = cluster.numStates;
  for (int i = 0; i < numStates * numStates; i++)
  {
    *currPtr = cluster.psiData_pottsDiffProb[direction];
    currPtr++;
  }
  currPtr = destMatrix;
  for(int i = 0; i < numStates; i++)
  {
    *currPtr = cluster.psiData_pottsSameProb[direction];
    currPtr+= numStates+1;
  }

}

void initOneNodeMsgMem(OneNodeCluster *nodeArray, FLOATTYPE *memChunk, 
		       const int numNodes, const int msgChunkSize)
{
  FLOATTYPE *currPtr = memChunk;
  OneNodeCluster *currNode = nodeArray;
  for(int i = 0; i < numNodes; i++)
  {

    currNode->receivedMsgs[0] = currPtr; currPtr+=msgChunkSize;
    currNode->receivedMsgs[1] = currPtr; currPtr+=msgChunkSize;
    currNode->receivedMsgs[2] = currPtr; currPtr+=msgChunkSize;
    currNode->receivedMsgs[3] = currPtr; currPtr+=msgChunkSize;
    currNode->nextRoundReceivedMsgs[0] = currPtr; currPtr +=msgChunkSize;
    currNode->nextRoundReceivedMsgs[1] = currPtr; currPtr +=msgChunkSize;
    currNode->nextRoundReceivedMsgs[2] = currPtr; currPtr +=msgChunkSize;
    currNode->nextRoundReceivedMsgs[3] = currPtr; currPtr +=msgChunkSize;
  
    for(int j = 0; j < 4;j++)
    {
      currNode->receivedMsgsb[j] = currNode->receivedMsgs[j];
      currNode->nextRoundReceivedMsgsb[j] = currNode->nextRoundReceivedMsgs[j];
    }
    currNode++;
  }
}

void initTwoNodeMsgMem(TwoNodeCluster *nodeArray, FLOATTYPE *memChunk, 
		       const int numNodes, const int msgChunkSize)
{
  FLOATTYPE *currPtr = memChunk;
  TwoNodeCluster *currNode = nodeArray;
  for(int i = 0; i < numNodes; i++)
  {

    currNode->receivedMsgs[0] = currPtr; currPtr+=msgChunkSize;
    currNode->receivedMsgs[1] = currPtr; currPtr+=msgChunkSize;
    currNode->nextRoundReceivedMsgs[0] = currPtr; currPtr +=msgChunkSize;
    currNode->nextRoundReceivedMsgs[1] = currPtr; currPtr +=msgChunkSize;
    for(int j = 0; j < 2;j++)
    {
      currNode->receivedMsgsb[j] = currNode->receivedMsgs[j];
      currNode->nextRoundReceivedMsgsb[j] = currNode->nextRoundReceivedMsgs[j];
    }

    currNode++;
  }
}

void OneNodeCluster::ComputeMsgRight(TwoNodeCluster &cluster,
				     FLOATTYPE *msgDest, FLOATTYPE *lastMsg, 
				     const FLOATTYPE alpha)
{
  FLOATTYPE *nodeLeftMsg = receivedMsgs[LEFT],
    *nodeDownMsg = receivedMsgs[DOWN],
    *nodeUpMsg =    receivedMsgs[UP];


  FLOATTYPE psiMat[numStates * numStates];

  getPsiMat(*this,psiMat,OneNodeCluster::RIGHT);
  
  FLOATTYPE *cmessage = msgDest,
    *cPrevMsg = lastMsg,
    total = 0;
  
  for(int rightNodeInd = 0; rightNodeInd < numStates; rightNodeInd++)
  {

    *cmessage = 0;
    for(int leftNodeInd = 0; leftNodeInd < numStates; leftNodeInd++)
    {

	FLOATTYPE tmp = nodeLeftMsg[leftNodeInd] * 
	nodeUpMsg[leftNodeInd] * 
	nodeDownMsg[leftNodeInd] * 
	psiMat[leftNodeInd * numStates + rightNodeInd] * 
	localEv[leftNodeInd];

	if (tmp > *cmessage)
	  *cmessage = tmp;
    }
    total += *cmessage;
    cmessage++;
  }

  for(int i = 0; i < numStates; i++)
  {
    msgDest[i] /= total;
    msgDest[i] = alpha * msgDest[i] + (1-alpha) * *cPrevMsg;
    cPrevMsg++;
  }
}

// This means, "Compute the message to send left."

//Be Careful about orientation of Messages for two-cluster message!
void OneNodeCluster::ComputeMsgLeft(TwoNodeCluster &cluster,
				    FLOATTYPE *msgDest, FLOATTYPE *lastMsg, const FLOATTYPE alpha)
{
  FLOATTYPE *nodeRightMsg = receivedMsgs[RIGHT],
    *nodeDownMsg = receivedMsgs[DOWN],
    *nodeUpMsg =    receivedMsgs[UP];


  FLOATTYPE psiMat[numStates * numStates];

  getPsiMat(*this,psiMat,OneNodeCluster::LEFT);
  
  FLOATTYPE *cmessage = msgDest,
    *cPrevMsg = lastMsg,
    total = 0;
  
  for(int leftNodeInd = 0; leftNodeInd < numStates; leftNodeInd++)
  {

    *cmessage = 0;
    for(int rightNodeInd = 0; rightNodeInd < numStates; rightNodeInd++)
    {
      FLOATTYPE tmp = 	nodeRightMsg[rightNodeInd] * 
	nodeUpMsg[rightNodeInd] * 
	nodeDownMsg[rightNodeInd] * 
	psiMat[leftNodeInd * numStates + rightNodeInd] * 
	localEv[rightNodeInd];

	if (tmp > *cmessage)
	  *cmessage = tmp;
    }
    total += *cmessage;
    cmessage++;
  }

  for(int i = 0; i < numStates; i++)
  {
    msgDest[i] /= total;
    msgDest[i] = alpha * msgDest[i] + (1-alpha) * *cPrevMsg;
    cPrevMsg++;
  }
}

void OneNodeCluster::ComputeMsgUp(TwoNodeCluster &cluster,
				  FLOATTYPE *msgDest, FLOATTYPE *lastMsg, const FLOATTYPE alpha)
{
  FLOATTYPE *nodeRightMsg = receivedMsgs[RIGHT],
    *nodeDownMsg = receivedMsgs[DOWN],
    *nodeLeftMsg =    receivedMsgs[LEFT];
 

  FLOATTYPE psiMat[numStates * numStates];

  getPsiMat(*this,psiMat,OneNodeCluster::UP);
  
  FLOATTYPE *cmessage = msgDest,
    *cPrevMsg = lastMsg,
    total = 0;
  
  for(int upNodeInd = 0; upNodeInd < numStates; upNodeInd++)
  {

    *cmessage = 0;
    for(int downNodeInd = 0; downNodeInd < numStates; downNodeInd++)
    {
      FLOATTYPE tmp = nodeRightMsg[downNodeInd] * 
	nodeLeftMsg[downNodeInd] * 
	nodeDownMsg[downNodeInd] * 
	psiMat[upNodeInd * numStates + downNodeInd] * 
	localEv[downNodeInd];
	if (tmp > *cmessage)
	  *cmessage = tmp;

    }

    total += *cmessage;
    cmessage++;
  }

  for(int i = 0; i < numStates; i++)
  {
    msgDest[i] /= total;
    msgDest[i] = alpha * msgDest[i] + (1-alpha) * *cPrevMsg;
    cPrevMsg++;
    if (msgDest[i] != msgDest[i]) // Check for NaN
    { printf("Break Here2 %f\n",total); }
  }
}

void OneNodeCluster::ComputeMsgDown(TwoNodeCluster &cluster,
				    FLOATTYPE *msgDest, FLOATTYPE *lastMsg, const FLOATTYPE alpha)
{
  FLOATTYPE *nodeRightMsg = receivedMsgs[RIGHT],
    *nodeUpMsg = receivedMsgs[UP],
    *nodeLeftMsg =    receivedMsgs[LEFT];
 

  FLOATTYPE psiMat[numStates * numStates];

  getPsiMat(*this,psiMat,OneNodeCluster::DOWN);
  
  FLOATTYPE *cmessage = msgDest,
    *cPrevMsg = lastMsg,
    total = 0;
  
  for(int downNodeInd = 0; downNodeInd < numStates; downNodeInd++)
  {

    *cmessage = 0;
    for(int upNodeInd = 0; upNodeInd < numStates; upNodeInd++)
    {
      FLOATTYPE tmp = 	nodeRightMsg[upNodeInd] * 
	nodeLeftMsg[upNodeInd] * 
	nodeUpMsg[upNodeInd] * 
	psiMat[upNodeInd * numStates + downNodeInd] * 
	localEv[upNodeInd];

	if (tmp > *cmessage)
	  *cmessage = tmp;

    }
    total += *cmessage;
    cmessage++;
  }

  for(int i = 0; i < numStates; i++)
  {
    msgDest[i] /= total;

    msgDest[i] = alpha * msgDest[i] + (1-alpha) * *cPrevMsg;
    cPrevMsg++;
    if (msgDest[i] != msgDest[i])
    { printf("Break Here2 %f\n",total);}
  }

}







void OneNodeCluster::deliverMsgs()
{

	for(int i = 0; i < numStates; i++)
	{
	  receivedMsgs[UP][i] = 0.5 * receivedMsgs[UP][i] + 0.5 * nextRoundReceivedMsgs[UP][i];
	}
	for(int i = 0; i < numStates; i++)
	{
	  receivedMsgs[DOWN][i] = 0.5 * receivedMsgs[DOWN][i] + 0.5 * nextRoundReceivedMsgs[DOWN][i];
	}
	for(int i = 0; i < numStates; i++)
	{
	  receivedMsgs[LEFT][i] = 0.5 * receivedMsgs[LEFT][i] + 0.5 * nextRoundReceivedMsgs[LEFT][i];
	}
	for(int i = 0; i < numStates; i++)
	{
	  receivedMsgs[RIGHT][i] = 0.5 * receivedMsgs[RIGHT][i] + 0.5 * nextRoundReceivedMsgs[RIGHT][i];
	}

}

void OneNodeCluster::getBelief(float *beliefVec)
{
	FLOATTYPE sum=0;
	for(int i = 0; i < numStates; i++)
	{
		beliefVec[i] = receivedMsgs[UP][i] * receivedMsgs[DOWN][i] *
				receivedMsgs[LEFT][i] * receivedMsgs[RIGHT][i] *
				localEv[i];
		sum += beliefVec[i];
	}

	for(int i = 0; i < numStates; i++)
	{
		beliefVec[i] /= sum;
	}
}

int checkBeliefConvergence(OneNodeCluster *nodeArray, const int numRows, const int numCols, 
		       int *destImage)
{
  const int numStates = OneNodeCluster::numStates;
  float currPtr[numStates];
  int cnt =0;
  int *destPtr = destImage;

  for (int m = 0; m < numRows; m++)
  {
    for(int n = 0; n < numCols; n++)
    {
      nodeArray[m * numCols + n].getBelief(currPtr);
      float best_prob = currPtr[0];
      int best_state = 0;
      for(int i = 1; i < numStates; i++)
      {
	if (currPtr[i] > best_prob)
	{
	  best_prob = currPtr[i];
	  best_state = i;
	}
      }
      if(best_state != *destPtr)
	cnt++;
      *destPtr = best_state;
      destPtr++;
    }
  }
  return cnt;
}


void passOneNodeMsgsLeft(OneNodeCluster *nodeArray, const int numCols,
			 TwoNodeCluster *topHorzRow, TwoNodeCluster *vertRow, 
			 TwoNodeCluster &dummyNode, const float alpha)
{
  const int numStates = nodeArray[0].numStates;
  //nodeArray should have 2 rows!
  int cind = numCols-1;
  for(int n = numCols-1; n > 0; n--)
  {
    nodeArray[cind].ComputeMsgLeft(dummyNode,//topHorzRow[cind-1],
				     nodeArray[cind-1].nextRoundReceivedMsgs[OneNodeCluster::RIGHT],
				     nodeArray[cind-1].receivedMsgs[OneNodeCluster::RIGHT],
				     alpha);

      FLOATTYPE *tmp = nodeArray[cind-1].nextRoundReceivedMsgs[OneNodeCluster::RIGHT];
      nodeArray[cind-1].nextRoundReceivedMsgs[OneNodeCluster::RIGHT] = nodeArray[cind-1].receivedMsgs[OneNodeCluster::RIGHT];
      nodeArray[cind-1].receivedMsgs[OneNodeCluster::RIGHT] = tmp;


      cind--;
  }
  for(int i = 0; i < numStates; i++)
    nodeArray[0].nextRoundReceivedMsgs[OneNodeCluster::LEFT][i] = 1;

}


void passOneNodeMsgsRight(OneNodeCluster *nodeArray, const int numCols,
			  TwoNodeCluster *topHorzRow, TwoNodeCluster *vertRow, 
			  TwoNodeCluster &dummyNode, const float alpha)
{
  const int numStates = OneNodeCluster::numStates;
  //nodeArray should have 2 rows!
  int cind = 0;
  for(int n = 0; n < numCols-1; n++)
  {

    nodeArray[cind].ComputeMsgRight(dummyNode,
				    nodeArray[cind+1].nextRoundReceivedMsgs[OneNodeCluster::LEFT],
				    nodeArray[cind+1].receivedMsgs[OneNodeCluster::LEFT],
				    alpha);
    FLOATTYPE *tmp = nodeArray[cind+1].nextRoundReceivedMsgs[OneNodeCluster::LEFT];
    nodeArray[cind+1].nextRoundReceivedMsgs[OneNodeCluster::LEFT] = nodeArray[cind+1].receivedMsgs[OneNodeCluster::LEFT];
    nodeArray[cind+1].receivedMsgs[OneNodeCluster::LEFT] = tmp;

    cind++;

  }
  
  for(int i = 0; i < numStates; i++)
    nodeArray[cind].nextRoundReceivedMsgs[OneNodeCluster::RIGHT][i] = 1;

}

void passOneNodeMsgsDown(OneNodeCluster *nodeArray,
			 TwoNodeCluster **vertRows, TwoNodeCluster &dummyNode,
			 const int col, const int numRows, const int numCols, const FLOATTYPE alpha)
{
  OneNodeCluster *currOneNodePtr = &nodeArray[col],
    *nextRowPtr = &nodeArray[col+numCols];
  for(int m = 0; m < numRows-1; m++)
  {
    currOneNodePtr->ComputeMsgDown(dummyNode,
				   nextRowPtr->nextRoundReceivedMsgs[OneNodeCluster::UP],
				   nextRowPtr->receivedMsgs[OneNodeCluster::UP],
				   alpha);
    FLOATTYPE *tmp = nextRowPtr->nextRoundReceivedMsgs[OneNodeCluster::UP];
    nextRowPtr->nextRoundReceivedMsgs[OneNodeCluster::UP] = nextRowPtr->receivedMsgs[OneNodeCluster::UP];
    nextRowPtr->receivedMsgs[OneNodeCluster::UP] = tmp;
    
    currOneNodePtr += numCols;
    nextRowPtr += numCols;
    
  }

  
 
}

void passOneNodeMsgsUp(OneNodeCluster *nodeArray,
			 TwoNodeCluster **vertRows, TwoNodeCluster &dummyNode,
			 const int col, const int numRows, const int numCols, const FLOATTYPE alpha)
{
  OneNodeCluster *currOneNodePtr = &nodeArray[numCols * (numRows-1)+col],
    *nextRowPtr = currOneNodePtr-numCols;

  for(int m = numRows-1; m > 0; m--)
  {
    currOneNodePtr->ComputeMsgUp(dummyNode,
				   nextRowPtr->nextRoundReceivedMsgs[OneNodeCluster::DOWN],
				   nextRowPtr->receivedMsgs[OneNodeCluster::DOWN],
				   alpha);
    FLOATTYPE *tmp = nextRowPtr->nextRoundReceivedMsgs[OneNodeCluster::DOWN];
    nextRowPtr->nextRoundReceivedMsgs[OneNodeCluster::DOWN] = nextRowPtr->receivedMsgs[OneNodeCluster::DOWN];
    nextRowPtr->receivedMsgs[OneNodeCluster::DOWN] = tmp;
    
    currOneNodePtr -= numCols;
    nextRowPtr -= numCols;
    
  }
 
}


TwoNodeCluster::TwoNodeCluster()
{
  staySame[0] = 0;
  staySame[1] = 0;
  disabled=0;
}

OneNodeCluster::OneNodeCluster()
{
  staySame = 0;
}

