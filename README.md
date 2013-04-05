testStereoAlgs
==============

Scripts for testing various stereo matching algorithms.

Intro
-------------
This repo contains the scripts that build a common testing interface to several stereo matching algorithms. Currently these include:

1.  Middlebury Markov Random Field (mrf) (http://vision.middlebury.edu/MRF/code/)
2.  Middlebury Taxonomy Paper Test Suite, with all patches and add-ons integrated (sm) (http://vision.middlebury.edu/stereo/code/)
3.  ELAS - Efficient Large Scale Stereo Matching (elas) (http://www.cvlibs.net/software/libelas.html)

Usage
--------------

1.  Set variables in run.sh.

    `IMGPATH` contains *horizontally rectified* image pairs. 
    Vertically rectified pairs (such as our staircase scenes) must be rotated 90 degrees before passing to stereo matching programs. You can use the command-line tools of [ImageMagick](http://www.imagemagick.org/Usage/) to do these kinds of simple transformation. For example: `convert $IMG -rotate 90 $IMG_ROT`. Some algorithms are very slow for high-resolution images, in that case you can shrink the images using, for example, `convert $IMG -resize 25% $IMG_SMALL`.
    
    `OUTPUTPATH` contains the output disparity maps. Each output is suffixed by the algorithm name used to generate it.
    
    `PROGSPATH` is the directory containing the three libraries.

2.  Edit parameters file.

    The parameters for ELAS and MRF can be set by editing respective txt file.
    Middlebury StereoMatcher relies on a directory of control scripts. The parameters are in stairs/param_in.txt and exp*_*.txt. The txt files in the corresponding output path seem to contain the complete parameter list.
    
3.  Then, in the terminal, run `source run.sh algorithm image_id parameters_file_or_folder`

    For example: `source run.sh elas 5 paramsELAS.txt` runs ELAS on image pair 5, using the parameters set in paramsELAS.txt. Or, 
    `source run.sh sm 5 StereoMatcherScript` runs Middlebury StereoMatcher on image pair 5, using the control scripts in StereoMatcherScript.
