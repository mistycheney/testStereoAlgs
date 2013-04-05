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

#ifndef _NODE_H
#define _NODE_H
struct Node;

struct Node {
  int maxNumNeighbors,
    numNeighbors,
    numStates,
    myIndex,
    *neighbors;

  double alpha,
    one_minus_alpha,
    *localEvidence,
    *belief,
    **psiMats, //ELiminated for Potts Model
    **prevMessages,
    **currMessages,
    **pottsProbs;
  
  
    Node *nodeArray;
  
};

void initNode(Node *node, int newMaxNumNeighbors, int newNumStates, Node *nodeArrayPtr, double alpha);
int getMessage(Node *node, int destNode, double *messageVec); /*Returns message from this node to destNode*/
int addNeighbor(Node *node, int *psiIndices, double *psiMatrix, int firstNode);
/* if psiIndices = [3 4] and you are working on node 3, firstNode =1*/
/* if you are working on node 4 firstNode =0;*/
void doIteration(Node *node);
void finishIteration(Node *node);
void getBelief(Node *node, double *beliefVec); /*You must call finishIteration before calling getBelief*/
double getMMSE(Node *node, double *states);
double *getLocalEvidencePtr(Node *node);
void writeNode(Node *node);
void writeMessages(Node *node);
Node *createMRFStruct(int numNodes, int numTensors, double *localEvidence[], 
		      double *psi_indices[], double *psi_mats[],
		      double alpha, int newNumStates, int maxNumNeighbors);


#endif
