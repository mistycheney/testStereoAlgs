///////////////////////////////////////////////////////////////////////////
//
// NAME
//  StereoMatcher.cpp -- stereo correspondence algorithms their evaluation
//
// DESCRIPTION
//  The CStereoMatcher class implements a number of popular stereo
//  correspondence algorithms (currently for 2-frame rectified image pairs).
//
// The definition of disparity is as follows:
//  Disparity is the (floating point) displacement
//  of a pixel between reference frame and match frame multi-frame stereo pair.
//  The images are always specified left-to-right.
//
//  When we store the floating point disparities in a gray_level image, we
//  use the formulas
//      g_d = (d - disp_min) * disp_scale
//        d =  disp_min  + g_d / disp_scale
//  to do the conversions.

// SEE ALSO
//  StereoMatcher.h         longer description of this class
//  Stc*.cpp                implementations of components (member functions)
//
// Copyright © Richard Szeliski and Daniel Scharstein, 2001.
// See Copyright.h for more details
//
///////////////////////////////////////////////////////////////////////////

#include "Error.h"
#include "StereoMatcher.h"
#include "Convert.h"
#include <math.h>
#include <time.h>

#define ROUND(x) ((int)( (x) >= 0 ? (x) + .5 : (x) - .5))

void CStereoMatcher::ComputeCorrespondence()
{

    // Check that frame numbers are valid
    if (frame_ref < 0 || frame_ref >= (int)m_frame.size())
        throw CError("ComputeCorrespondence: invalid reference frame number %d",
        frame_ref);
    if (frame_match < 0 || frame_match >= (int)m_frame.size())
        throw CError("ComputeCorrespondence: invalid matching frame number %d",
        frame_match);

    // Check that the input images are valid
    CByteImage& ref = m_frame[frame_ref  ].input_image;  // reference image
    CByteImage& mat = m_frame[frame_match].input_image;  // matching (target) image
    CShape sh0 = ref.Shape();
    CShape sh1 = mat.Shape();
    if (sh0.width * sh0.height == 0)
        throw CError("ComputeCorrespondence: invalid reference image");
    if (sh1.width * sh1.height == 0)
        throw CError("ComputeCorrespondence: invalid matching image");
    if (sh0 != sh1)
        throw CError("ComputeCorrespondence: reference and matching not the same size");
    m_frame_diff = frame_match - frame_ref;         // frame difference
    m_frame_diff_sign = (m_frame_diff > 0) ? 1 : -1;

    // compute m_disp_num, m_disp_den, and disp_n

    const float min_precision = 1e-3f;
    if (disp_step <= 0.0 || 
        (disp_step < 1.0 && fabs(1.0/disp_step - ROUND(1.0/disp_step)) > min_precision) ||
        (disp_step > 1.0 && fabs(disp_step - ROUND(disp_step)) > min_precision))
    {
        throw CError("ComputeCorrespondence: disp_step must integer N or 1.0/N");
    }
    m_disp_num = disp_step < 1.0f ? 1 : ROUND(disp_step);
    m_disp_den = disp_step < 1.0f ? ROUND(1.0/disp_step) : 1;

    // recompute disp_step and its inverse
    disp_step       = m_disp_num / (float) m_disp_den;
    m_disp_step_inv = m_disp_den / (float) m_disp_num;
    disp_n = int(m_disp_step_inv * (disp_max - disp_min)) + 1;

    // store in internal variables (may be collapsed later)
    m_disp_step     = disp_step;
    m_disp_n        = disp_n;

    if (verbose >= eVerboseProgress)
    {
        if (evaluate_only)
           fprintf(stderr, "\nEvaluating existing depth map\n");
        else
        {
            fprintf(stderr, "\nComputing correspondence between frames %d and %d\n",
                frame_ref, frame_match);

            // TODO: this should be a different verbosity level (e.g., 
            // eVerboseDispParams)
            fprintf(stderr, 
            "disp_n = %d, disp_min = %d, m_disp_num = %d, m_disp_den = %d, disp_scale = %g\n",
            disp_n, disp_min, m_disp_num, m_disp_den, disp_scale);
        }
    }

    // Copy the selected reference and matching frames
    //  (since they may end up being pre-processed)
    m_reference.ReAllocate(sh0);
    CopyPixels(ref, m_reference);
    m_matching.ReAllocate(sh1);
    CopyPixels(mat, m_matching);

    // Copy the depth map and reallocate if necessary
    CShape shD = sh0;   // shape of disparity image (1 band)
    shD.nBands = 1;
    CByteImage& depth_image = m_frame[frame_ref].depth_image;
    bool already_initialized = (depth_image.Shape() == shD);
    depth_image.ReAllocate(shD);        // keep data if same shape (for evaluate_only)
    if (! already_initialized)
        depth_image.ClearPixels();      // set to 0 if we just reallocated
    m_float_disparity.ReAllocate(shD);  // sub-pixel disparities

    ScaleAndOffset(depth_image, m_float_disparity, 1.0f/disp_scale, disp_min);

    // Copy the ground truth image if it exists
    CByteImage& truth_image = m_frame[frame_ref].truth_image;
    m_true_disparity.ReAllocate(shD);   // ground truth
    if (truth_image.Shape() == shD)
    {
        ScaleAndOffset(truth_image, m_true_disparity, 1.0f/disp_scale, disp_min);
    }
    else
    {
        m_true_disparity.ClearPixels();     // set to 0, can't think of anything better
    }

    if (evaluate_only) {

        // I put my old code back in.  Without it, evaluate_only didn't work.
        // See test_evaluate_only.txt

        float d_offset = disp_min;
        
        // formulas for converting between DSI disparities (k) and float disparities (d):
        //      d = disp_min + k * m_disp_num / m_disp_den
        //      k = (d - disp_min) / (m_disp_num/m_disp_den)

        // here we need the second formula:
        ScaleAndOffset(m_float_disparity, m_disparity,
            m_disp_step_inv, - d_offset * disp_step);

        // Don't run the rest of the matcher
        return;

        // Note: if you want to compute energies, don't use evaluate_only.
        // Rather, use opt_fn = eNoOpt
    }
    
    // Allocate the integer disparity and cost images
    m_disparity.ReAllocate(shD);        // winning disparities
    CShape shC = shD;
    shC.nBands = m_disp_n;                // number of disparity levels
    if (m_disp_n < 2)
        throw CError("ComputeCorrespondence: too few disparity levels (%d)", m_disp_n);
    m_cost.ReAllocate(shC);             // raw matching costs (# bands = # disparities)

    clock_t time0 = clock();    // record start time (can't use global timer...)

    // Invoke all of the stereo matching processing steps
    PreProcess();   // see StcPreProcess.cpp
    RawCosts();     // see StcRawCosts.cpp
    Aggregate();    // see StcAggregate.cpp
    Optimize();     // see StcOptimize.cpp
    Refine();       // see StcRefine.cpp
        
    clock_t time1 = clock();    // record end time
    total_time =  (float)(time1-time0)/(float)CLOCKS_PER_SEC;
    if (verbose >= eVerboseTiming)
      fprintf(stderr, "* total time: %gs\n", total_time);

    // Convert the final disparities back into the scaled depth map (round)
    ScaleAndOffset(m_float_disparity, depth_image, disp_scale, -disp_min * disp_scale + 0.5);

    // Restore the value of m_reference for evaluation purposes (undo PreProcessing)
    CopyPixels(ref, m_reference);
}
