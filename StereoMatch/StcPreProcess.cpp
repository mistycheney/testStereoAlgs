///////////////////////////////////////////////////////////////////////////
//
// NAME
//  StcPreProcess.cpp -- pre-process the images to clean them up or normalize them
//
// DESCRIPTION
//  Currently, we only support iterated binomial blur, to clean up the images
//  a little.  This should help sub-pixel fitting work better, by making
//  image shifts closer to a Taylor series expansion, but will result in worse
//  performance near discontinuity regions and in finely textured regions.
// 
//  Other potential pre-processing operations (currently not implemented),
//  might include:
//    bias and gain normalization
//    histogram equalization (global or local)
//    rank statistics pre-processing
//
// SEE ALSO
//  StereoMatcher.h         longer description of this class
//
// Copyright © Richard Szeliski, 2001.
// See Copyright.h for more details
//
///////////////////////////////////////////////////////////////////////////

#include "StereoMatcher.h"
#include "Convolve.h"
#include "ImageIO.h"
#include "Error.h"

void CStereoMatcher::PreProcess()
{
    if (preproc_addnoise_stddev > 0.0f)
    {
        // TODO: implement adding noise to input images
        throw CError("PreProcess: additive noise not yet implemented");
    }

    if (preproc_blur_iter <= 0)
        return;

    if (verbose == eVerboseSummary)
        fprintf(stderr, "pre-process, ");
    if (verbose >= eVerboseProgress)
        fprintf(stderr, "- pre-process: binomial 121 filter iterated %d times\n",
                preproc_blur_iter);

    StartTiming();
    for (int iter = 0; iter < preproc_blur_iter; iter++)
    {
        ConvolveSeparable(m_reference, m_reference, ConvolveKernel_121, ConvolveKernel_14641,
             1.0f, 0.0f, 1, 1);
        ConvolveSeparable(m_matching, m_matching, ConvolveKernel_121, ConvolveKernel_14641,
             1.0f, 0.0f, 1, 1);
    }
    PrintTiming();

    // Write out the pre-processed images
    static bool dump_pre_processed = false;     // reset in debugger
    if (dump_pre_processed && verbose >= eVerboseDumpFiles)
    {
        WriteImage(m_reference, "reprojected/tmp_m_reference.ppm");
        WriteImage(m_matching,  "reprojected/tmp_m_matching.ppm");
    }
}
