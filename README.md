testStereoAlgs
==============

Scripts for testing various stereo matching algorithms.

Intro
-------------
This repo contains the scripts that build a common testing interface to several stereo matching algorithms. Currently these include:

1.  Middlebury Markov Random Field (mrf) (http://vision.middlebury.edu/MRF/code/)
2.  Middlebury Taxonomy Paper Test Suite (sm) (http://vision.middlebury.edu/stereo/code/)
3.  ELAS - Efficient Large Scale Stereo Matching (elas) (http://www.cvlibs.net/software/libelas.html)

Usage
--------------

First, set variables in the run.sh. `PROGSPATH` is the directory containing the three libraries. 

Then, in the terminal, run `source run.sh algorithm image_id parameters_file`

For example: `source run.sh elas 5 paramsELAS.txt` runs ELAS on image pair 5, using the parameters set in paramsELAS.txt.