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

class TwoNodeCluster;
class OneNodeCluster
{
public:
  OneNodeCluster();
  const static int UP=0, DOWN = 1, LEFT = 2, RIGHT = 3;
  static int numStates;int staySame,disabled;
  
  FLOATTYPE *localEv;
  FLOATTYPE   *receivedMsgs[4],
              *nextRoundReceivedMsgs[4];

  FLOATTYPE   *receivedMsgsb[4],
              *nextRoundReceivedMsgsb[4];

  FLOATTYPE psiData_pottsSameProb[4], 
            psiData_pottsDiffProb[4];  

  void ComputeMsgRight(TwoNodeCluster &cluster,
		       FLOATTYPE *msgDest, FLOATTYPE *lastMsg, 
		       const FLOATTYPE alpha);

  void ComputeMsgUp(TwoNodeCluster &cluster,
		    FLOATTYPE *msgDest, FLOATTYPE *lastMsg, const FLOATTYPE alpha);

  void ComputeMsgLeft(TwoNodeCluster &cluster,
		      FLOATTYPE *msgDest, FLOATTYPE *lastMsg, const FLOATTYPE alpha);

  void ComputeMsgDown(TwoNodeCluster &cluster,
		      FLOATTYPE *msgDest, FLOATTYPE *lastMsg, const FLOATTYPE alpha);

  void getBelief(float *beliefVec);
  
  int msgChange(FLOATTYPE thresh);

  void deliverMsgs();
    
};

class TwoNodeCluster
{
public:
  const static int UP=0,DOWN=1, LEFT=0, RIGHT=1;
  static int numStates;
  int staySame[2],thisStaySame,disabled;
  FLOATTYPE *receivedMsgs[2],
    *nextRoundReceivedMsgs[2];

  FLOATTYPE *receivedMsgsb[2],
    *nextRoundReceivedMsgsb[2];
  
  TwoNodeCluster();
  void ComputeMsgUp(TwoNodeCluster &leftCluster, TwoNodeCluster &rightCluster,
		    OneNodeCluster &leftNode, OneNodeCluster &rightNode,
		    FLOATTYPE *msgDest, FLOATTYPE *lastMsg, const FLOATTYPE alpha);

  void ComputeMsgDown(TwoNodeCluster &leftCluster, TwoNodeCluster &rightCluster,
		      OneNodeCluster &leftNode, OneNodeCluster &rightNode,
		      FLOATTYPE *msgDest, FLOATTYPE *lastMsg, const FLOATTYPE alpha);
  
  void ComputeMsgRight(TwoNodeCluster &upCluster, TwoNodeCluster &downCluster,
		       OneNodeCluster &upNode, OneNodeCluster &downNode,
		       FLOATTYPE *msgDest, FLOATTYPE *lastMsg, const FLOATTYPE alpha);
  
  void ComputeMsgLeft(TwoNodeCluster &upCluster, TwoNodeCluster &downCluster,
		      OneNodeCluster &upNode, OneNodeCluster &downNode,
		      FLOATTYPE *msgDest, FLOATTYPE *lastMsg, const FLOATTYPE alpha);

  void deliverMsgs();
  int msgChange(FLOATTYPE thresh);

  int canDisableHorz(const int r, const int c, 
		     const int numRows, const int numCols,
		     OneNodeCluster *nodeArray,const TwoNodeCluster *upRow,
		     const TwoNodeCluster *downRow);
};

void initOneNodeMsgMem(OneNodeCluster *nodeArray, FLOATTYPE *memChunk, const int numNodes, 
		       const int msgChunkSize);
void initTwoNodeMsgMem(TwoNodeCluster *nodeArray, FLOATTYPE *memChunk, const int numNodes, 
		       const int msgChunkSize);

void computeOneNodeMessages(OneNodeCluster *nodeArray, const int numCols,
			    TwoNodeCluster *topHorzRow, TwoNodeCluster *vertRow, 
			    TwoNodeCluster &dummyNode, const float alpha);

void computeTwoNodeMessages(OneNodeCluster *nodeArray, const int numCols,
			    TwoNodeCluster *topHorzRow, TwoNodeCluster *vertRow, 
			    TwoNodeCluster *botHorzRow, TwoNodeCluster &dummyNode, 
			    const float alpha);

void getPsiMat(OneNodeCluster &cluster, FLOATTYPE *destMatrix, int direction);

void createBeliefImage(OneNodeCluster *nodeArray, const int numRows, 
		       const int numCols, FLOATTYPE *destImage);
void passOneNodeMsgsDown(OneNodeCluster *nodeArray,
			 TwoNodeCluster **vertRows, TwoNodeCluster &dummyNode,
			 const int col, const int numRows, const int numCols, const FLOATTYPE alpha);
void passOneNodeMsgsRight(OneNodeCluster *nodeArray, const int numCols,
			  TwoNodeCluster *topHorzRow, TwoNodeCluster *vertRow, 
			  TwoNodeCluster &dummyNode, const float alpha);
void passOneNodeMsgsLeft(OneNodeCluster *nodeArray, const int numCols,
			 TwoNodeCluster *topHorzRow, TwoNodeCluster *vertRow, 
			 TwoNodeCluster &dummyNode, const float alpha);

int checkBeliefConvergence(OneNodeCluster *nodeArray, const int numRows, const int numCols, 
			   int *destImage);



void passOneNodeMsgsUp(OneNodeCluster *nodeArray,
			 TwoNodeCluster **vertRows, TwoNodeCluster &dummyNode,
		       const int col, const int numRows, const int numCols, const FLOATTYPE alpha);

void computeOneNodeMessagesTopRow(OneNodeCluster *nodeArray, TwoNodeCluster *vertRow,
				  const int numCols, TwoNodeCluster &dummyNode, 
				  const float alpha);
void computeOneNodeMessagesLastRow(OneNodeCluster *nodeArray, TwoNodeCluster *vertRow,
				  const int numCols, TwoNodeCluster &dummyNode, 
				   const float alpha);
void computeTwoNodeMessagesTopRow(OneNodeCluster *nodeArray, const int numCols,
				  TwoNodeCluster *vertRow, TwoNodeCluster *botHorzRow, 
				  TwoNodeCluster &dummyNode, const float alpha);

