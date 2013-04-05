// Copyright 2001,2002,2003 Marshall Tappen
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

#include "node.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

double *globalPsiMat;

double multiplyAndSum(double *src1, double *src2, int numEl)
{
  int i;
  double res,*p1,*p2;
  res = 0;

  p1 = src1;
  p2 = src2;
  for (i = 0; i < numEl; i++)
  {
    res += *p1 * *p2;;
    ++p1;
    ++p2;
  }
  
  return res;
}

double multiplyAndSum_Potts(double *src1, const double pottsDiff, const double pottsSame, const int sameIndex, int numEl)
{
  int i;
  double res,*p1;
  res = 0;

  p1 = src1;
  for (i = 0; i < numEl; i++)
  {
    if (i != sameIndex)
    {
      res += *p1 * pottsDiff;
    }
    else
    {
      res += *p1 * pottsSame;
    }
    ++p1;
  }
  
  return res;
}


double multiplyAndMax(double *src1, double *src2, int numEl)
{
  double tmp,res,*p1,*p2;
  int i;

  res = -1;

  p1 = src1;
  p2 = src2;
  for ( i = 0; i < numEl; i++)
  {
    tmp= *p1 * *p2;;
    if(tmp > res)
      res = tmp;
    ++p1;
    ++p2;
  }
  
  return res;
}

void multiply(double *src1, double *src2, double *dst, int numEl)
{
  double *p1,*p2,*pd;
  int i;

  p1 = src1;
  p2 = src2;
  pd = dst;
  for (i = 0; i < numEl; i++)
  {
    *pd = *p1 * *p2;
    ++p1;
    ++p2;
    ++pd;
  }
  

}



int addNeighbor(Node *node,int *psiIndices, double *psiMatrix, int firstNode)
/*    //first node means whether the node is the first node in psiIndices or not */
{
//    double *oldPsiMat = psiMatrix;
//    double *newPsiMat = (double *)0;
  int error = 0;
  int numNeighbors;

/*    for(i = 0; i < 4;i++) */
/*      printf("%lf ",psiMatrix[i]); */
/*    printf("\n"); */
/*    numStates = node->numStates; */
  numNeighbors = node->numNeighbors;

  node->neighbors[numNeighbors]=psiIndices[firstNode];
  node->myIndex = psiIndices[!firstNode];

/*    if(!firstNode) */
/*    { */
/*      newPsiMat = (double *)malloc(sizeof(double)*node->numStates * node->numStates); */
/*      for (i = 0; i < numStates; i++) */
/*      { */
/*        for (j =0;j < numStates; j++) */
/*        { */
/*  	newPsiMat[i * numStates + j] = psiMatrix[j * numStates +i]; */
/*        } */
/*      }   */
/*      psiMatrix = newPsiMat; */
/*    } */
  
  //  memcpy(node->psiMats[numNeighbors], psiMatrix,sizeof(double)*numStates*numStates);
  memcpy(node->pottsProbs[numNeighbors], psiMatrix,sizeof(double)*2);
  node->numNeighbors=node->numNeighbors+1;
  numNeighbors = node->numNeighbors;
//    if(newPsiMat)
//    {
//      free(newPsiMat);
//      psiMatrix = oldPsiMat;
//    }
/*    for(i = 0; i < 4;i++) */
/*      printf("%lf ",node->psiMats[numNeighbors-1]); */
/*    printf("\n"); */
/*    printf("\n"); */

  return error;
}


int getMessage(Node *node,int destNode, double *messageVec)
{
  const int numStates = node->numStates;
  int neighborIndex;
  int error = 0;
  int *neighbors = node->neighbors;
  int i,j, numNeighbors = node->numNeighbors;
  double **prevMessages = node->prevMessages;
  double sum = 0;
  double *msgVecPtr = messageVec,
    *psiColPtr = 0;
  double beliefVec[numStates]; 
  double psiMat[node->numStates * node->numStates];

  neighborIndex = numNeighbors +1;
  for (i = 0; i < numNeighbors; i++)
  {
    if (neighbors[i] == destNode)
    {
      neighborIndex = i;
      break;
    }
  }
  

  if (neighborIndex > numNeighbors)
  {
    error =1;
    return error;
  }
  
  const double pottsProb = node->pottsProbs[neighborIndex][1];
  for (int m = 0; m < numStates * numStates; m++)
  {
    psiMat[m] = pottsProb;
  }
  for (int m = 0; m < numStates * numStates; m+=numStates+1)
  {

    psiMat[m] =  node->pottsProbs[neighborIndex][0];
  }

//    for (int m = 0; m < numStates; m++)
//    {
//      for (int n = 0; n < numStates; n++)
//      {
//        //      printf("%f ",psiMat[m *numStates + n]);

//      }
//      //    printf("\n");
//    }
/*    // Matrices are in column-major order */
  psiColPtr = psiMat;
//    psiColPtr = globalPsiMat;
  for (i = 0; i < numStates; i++)
  {
    memcpy(beliefVec,node->localEvidence,sizeof(double)*numStates);
    for (j = 0; j < numNeighbors; j++)
    {
      if (j != neighborIndex)
	multiply(beliefVec,prevMessages[j],beliefVec,numStates);
    }
    

//      *msgVecPtr = multiplyAndSum(psiColPtr, beliefVec, numStates);
    *msgVecPtr = multiplyAndMax(psiColPtr, beliefVec, numStates);
    psiColPtr += numStates;
    //    *msgVecPtr = multiplyAndSum_Potts(beliefVec, pottsProb,node->pottsProbs[neighborIndex][0],i, numStates);
    ++msgVecPtr;
  }


  for(i = 0; i < numStates; i++)
    sum += messageVec[i];

  for(i = 0; i < numStates; i++)
    messageVec[i] /= sum;


  return error;
}

void doIteration(Node *node)
{
  int i,j,  numNeighbors = node->numNeighbors,
    *neighbors = node->neighbors;
  double **currMessages = node->currMessages,alpha = node->alpha,
    one_minus_alpha = node->one_minus_alpha, **prevMessages = node->prevMessages;
  Node *nodeArray = node->nodeArray;
  const int numStates = node->numStates;

  double tmpMessage[numStates];

  for (i = 0; i < numNeighbors; i++)
  {
    getMessage(&nodeArray[neighbors[i]],node->myIndex, tmpMessage);
    for(j = 0; j < numStates; j++)
    {
      currMessages[i][j] = alpha* tmpMessage[j] + one_minus_alpha * prevMessages[i][j];
    }
  }
}


void finishIteration(Node *node)
{
  double **tmpPtr;
 

  tmpPtr = node->prevMessages;
  node->prevMessages = node->currMessages;
  node->currMessages = tmpPtr;
}


void getBelief(Node *node,double *beliefVec)
{
  
  const int numStates = node->numStates;
  int i,  numNeighbors = node->numNeighbors;
  double **prevMessages = node->prevMessages;
  double sum=0;

  memcpy(beliefVec,node->localEvidence,sizeof(double)*numStates);

  for (i = 0; i < numNeighbors; i++)
  {
    multiply(beliefVec,prevMessages[i],beliefVec,numStates);
  }
  for (i = 0; i < numStates; i++)
  {
    sum += beliefVec[i];
  }
  for (i = 0; i < numStates; i++)
  {
    beliefVec[i] = beliefVec[i]/sum;
  }
}


double *getLocalEvidencePtr(Node *node)
{
  return node->localEvidence;
}


double getMMSE(Node *node,double *states)
{

  const int numStates = node->numStates;
  double beliefVec[numStates];

  getBelief(node,beliefVec);
  return multiplyAndSum(states, beliefVec,node->numStates);
}

void freeNode(Node *node)
{
  int i;
  free(node->neighbors);
  for (i = 0; i < node->maxNumNeighbors; i++)
  {
    free(node->psiMats[i]);
    free(node->prevMessages[i]);
    free(node->currMessages[i]);
  }
  free(node->belief);
  free(node->psiMats);
  free(node->prevMessages);
  free(node->currMessages);
  free(node->localEvidence);
}

/*
void Node::writeMessages()
{
  for (int i = 0; i < numNeighbors; i++)
  {
    printf("    Message from %d: ",neighbors[i]);
    for(int j = 0; j < numStates; j++)
    {
      printf("%.4f ",prevMessages[i][j]);
    }
    printf("\n");
  }
}
void Node::writeNode()
{
  printf("NumNeigbors: %d\n Local Evidence: ",numNeighbors);
  for (int i = 0; i < numStates; i++)
    printf("%.3f ",localEvidence[i]);
  printf("\n");
  
  for (int i = 0; i < numNeighbors; i++)
  {
    printf("Neighbor: %d ",neighbors[i]);
    for(int j = 0; j < numStates; j++)
    {
      if (j>0) printf("             ");
      for(int k = 0; k < numStates; k++)
      {
	printf("%.3f ",psiMats[i][k* numStates + j]);
      }
      printf("\n");
    }

    for(int j = 0; j < numStates * numStates; j++)
      printf("%.3f ",psiMats[i][j]);
    printf("\n");
  }
}


*/
