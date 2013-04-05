CPPFLAGS = -O2 -Wall
#CPPFLAGS = -O2 -Wall -mcpu=pentium3 -msse -funroll-loops -ffast-math -mfpmath=sse 
#CXX=g++-2.95
#CPPFLAGS = -O2 -Wall -mcpu=athlon-xp -msse -funroll-loops -ffast-math -mfpmath=sse -m3dnow

#CPPFLAGS = -g
SRC = BoxFilter.cpp Convert.cpp Convolve.cpp Histogram1D.cpp Image.cpp ImageIO.cpp main.cpp MinFilter.cpp ParameterIO.cpp RefCntMem.cpp StcAggregate.cpp StcDiffusion.cpp StcEvaluate.cpp StcGraphCut.cpp StcOptDP.cpp StcOptimize.cpp StcOptSO.cpp StcPreProcess.cpp StcRawCosts.cpp StcRefine.cpp StcSimulAnn.cpp StereoIO.cpp StereoMatcher.cpp StereoParameters.cpp Warp1D.cpp StcOptBP.cpp StcOptBPSync.cpp node.cpp bpregions.cpp LASW.cpp

HDR = BoxFilter.h Convert.h Convolve.h Copyright.h Error.h Histogram1D.h Image.h ImageIO.h MinFilter.h ParameterIO.h RefCntMem.h StereoIO.h StereoMatcher.h StereoParameters.h Verbose.h Warp1D.h  LASW.h

OBJ = $(SRC:.cpp=.o) maxflow/maxflow.o
 
all: StereoMatch

StereoMatch: $(OBJ)
	g++ -o StereoMatch $(OBJ) -lm -O2
#	g++-2.95 -o StereoMatch $(OBJ) -lm -O2 

clean:
	rm -f $(OBJ)
