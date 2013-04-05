///////////////////////////////////////////////////////////////////////////
//
// NAME
//  StereoParameters.cpp -- parameter values controlling the stereo matcher
//
// DESCRIPTION
//  The CStereoParameters class encapsulates all of the parameters
//  (e.g., window size, ...) necessary to control the CStereoMatcher class.
//
//  The parameters are broken down into several simple structs for easier
//  legibility and documentation.
//
//  See the comments after each parameter for a description of its use.
//
// SEE ALSO
//  StereoParameters.h      longer description
//
// Copyright © Richard Szeliski and Daniel Scharstein, 2001.
// See Copyright.h for more details
//
///////////////////////////////////////////////////////////////////////////

#include "StereoParameters.h"


CStereoParameters::CStereoParameters()
{
    ReInitializeSeqParams();
    ReInitializeAlgParams();
    ResetOutputParams();
}

void CStereoParameters::ReInitializeSeqParams()
{
    // (Re-)Initialize Parameters specific to each image sequence
    // Note: these parameters are not affected by a "reset" command.  They should
    // be specified explicitly in the sequence's parameter file param_in.txt

    frame_ref = 0;              // reference frame
    frame_match = 1;            // matching image frame
    disp_min = 0;               // smallest disparity
    disp_max = 7;               // largest disparity
    disp_step = 1.0f;           // disparity step size (must be int or 1/int)
    disp_n = 0;                 // number of disparity levels (computed from above 3)
    disp_scale = 1.0f;          // disparity scaling factor for gray depthmap

    eval_ignore_border = 0;     // number of border pixels to ignore in ground truth image
    eval_disp_gap = 2.0f;       // don't interpolate across disparity jumps bigger than this

    // put these parameters here so they won't get changed by "reset"
    verbose = eVerboseProgress; // verbosity level (see Verbose.h)
    evaluate_only = 0;          // read specified depth map and evaluate only
}

void CStereoParameters::ReInitializeAlgParams()
{
    // (Re-)Initialize Parameters specific to each algorithm
    // Note: this method is called upon a "reset" command

    // Parameters controlling the pre-processing of input images:
    preproc_addnoise_stddev = 0.0f;   // additive noise standard deviation
    preproc_blur_iter = 0;            // number of pre-processing blur iterations

    // Parameters controlling the per-pixel matching:
    match_fn = eAD;             // matching function
    match_interp = eCubic;      // interpolation function
    match_max = 1000;           // maximum difference for truncated SAD/SSD
    match_interval = 0;         // search for best 1/2 disparity match
    match_interpolated = 0;     // interpolate both lines when matching

    // Parameters controlling the spatial aggregation:
    aggr_fn = eBox;             // aggregation function
    aggr_window_size = 7;       // size of window
    aggr_iter = 1;              // number of aggregation iterations
    aggr_minfilter = 0;         // spatial min-filter after aggregation (shiftable window)
    aggr_subpixel = 0;          // do local fits around minima to get lower value
    aggr_collapse = 0;          // collapse DSI back to integer disp sampling
    diff_lambda = 0.15f;        // lambda parameter for diffusion algorithms
    diff_beta = 0.5f;           // beta parameter for membrane model diffusion
    diff_scale_cost = 0.01f;    // scale for m_cost values (necessary for Bayesian diffusion)
    diff_mu = 0.5f;             // mu parameter for Bayesian diffusion
    diff_sigmaP = 0.4f;         // sigma for robust prior (0.1 / synthetic, 0.4 / real in IJCV)
    diff_epsP = 0.01f;          // epsilon for robust prior (0.01 in IJCV paper)
    // adaptive support weight parameters (see LASW.h and LASW.cpp)
    aggr_gamma_proximity = 5.5f;  // sigma value in the distance term
    aggr_gamma_similarity = 5.0f; // gamma value in the color term
    aggr_color_space = eCIELab;       // color space

    // Parameters controlling the optimization algorithms:
    opt_fn = eWTA;              // optimization function
    opt_smoothness = 1.0f;      // smoothness penalty multiplier
    opt_grad_thresh = 5.0f;     // threshold for magnitude of intensity gradient
    opt_grad_penalty = 1.0f;    // smoothness penalty factor if gradient is too small
    opt_occlusion_cost = 20;    // cost for occluded pixels in DP algorithm
    opt_max_iter = 100;         // maximum number of optimization iterations
    opt_random = 1;             // randomize optimization (disparity and/or pixel sites)
    opt_sa_var = eFullGibbs;    // simulated annealing variant (update rule)
    opt_sa_start_T = 10.0f;     // starting temperature
    opt_sa_end_T = 0.01f;       // ending temperature
    opt_sa_schedule = eSALinear;// annealing schedule
    opt_min_margin = 0.7f;      // neccessary margin of "clear" min in symmetric matcher
    opt_sym_passes = 1;         // number of outer iterations with further aggregation 

    // Parameters controlling the final sub-pixel refinement step:
    refine_subpix = 0;          // fit sub-pixel value to local correlation

    // Parameters controlling the quality / fitness evaluation
    eval_bad_thresh = 1.0f;     // acceptable disparity error (for counting "bad points")
    eval_error_scale = 2.0f;    // scale disparity errors and write out (along with bad pixels)
    eval_lin_interp = 1;        // use linear interpolation in forward warping
    eval_predict_type = ePredictForward; // type of prediction error to compute
    eval_textureless_width = 3;          // width of box filter for summing squared horiz. gradients
    eval_textureless_thresh = 4.0f;      // threshold applied to summed (h-grad)^2
    eval_discont_width = 9;     // width of discontinuity region (box filter)
    eval_predict_diff = 0;      // write out difference images instead of resampled
    eval_empty_color = 0x00ffc0ff;  // color of empty pixels: light magenta
    eval_partial_shuffle = 0.f; // use interval analysis for prediction error
    eval_match_quality = 0;     // evaluate match quality (final cost and certainty)
    eval_certain_matches_only = 0;  // only compute statistics for matches with status "certain"
}

void CStereoParameters::ResetOutputParams()
{
    // Reset output parameters that measure performance of algorithm

    rms_error_all = -1.0f;          // RMS disparity error (all pixels)
    rms_error_nonocc = -1.0f;       // RMS disparity error (non-occluded pixels only)
    rms_error_occ = -1.0f;          // RMS disparity error (occluded pixels only)
    rms_error_textured = -1.0f;     // RMS disparity error (textured pixels only)
    rms_error_textureless = -1.0f;  // RMS disparity error (textureless pixels only)
    rms_error_discont = -1.0f;      // RMS disparity error (near depth discontinuities)
    bad_pixels_all = -1.0f;         // fraction of bad points (all pixels)
    bad_pixels_nonocc = -1.0f;      // fraction of bad points (non-occluded pixels only)
    bad_pixels_occ = -1.0f;         // fraction of bad points (occluded pixels only)
    bad_pixels_textured = -1.0f;    // fraction of bad points (textured pixels only)
    bad_pixels_textureless = -1.0f; // fraction of bad points (textureless pixels only)
    bad_pixels_discont = -1.0f;     // fraction of bad points (near depth discontinuities)

    // prediction error for 4 reference positions
    // initialize to -1 because they are not used in all image sequences
    predict_err_near = -1.0f;      // prediction error near   (e.g., frame 0)
    predict_err_middle = -1.0f;    // prediction error middle (e.g., frame 4)
    predict_err_match = -1.0f;     // prediction error match  (e.g., frame 6)
    predict_err_far = -1.0f;       // prediction error far    (e.g., frame 8)

    final_energy = -1.0f;          // final energy of solution
    total_time = -1.0f;            // total time taken by algorithm
};

void CStereoParameters::PIOInitialize(CParameterIO& o)
{
    // Initialize parameter input/output object

    // Parameters controlling the pre-processing of input images:
    o.PPPF(preproc_addnoise_stddev);  // additive noise standard deviation
    o.PPPD(preproc_blur_iter);        // number of pre-processing blur iterations

    // Parameters controlling the per-pixel matching:
    o.PPPD(frame_ref);          // reference frame
    o.PPPD(frame_match);        // matching image frame
    o.PPPD(disp_min);           // smallest disparity
    o.PPPD(disp_max);           // largest disparity
    o.PPPF(disp_step);          // disparity step size (must be int or 1/int)
    o.PPPD(disp_n);             // disparity step size (must be int or 1/int)
    o.PPPF(disp_scale);         // disparity scaling factor for gray depthmap
    o.PPPD(match_fn);           // matching function
    o.PPPD(match_interp);       // interpolation function
    o.PPPD(match_max);          // maximum difference for truncated SAD/SSD
    o.PPPD(match_interval);     // search for best 1/2 disparity match
    o.PPPD(match_interpolated); // interpolate both lines when matching

    // Parameters controlling the spatial aggregation:
    o.PPPD(aggr_fn);            // aggregation function
    o.PPPD(aggr_window_size);   // size of window
    o.PPPD(aggr_iter);          // number of aggregation iterations
    o.PPPD(aggr_minfilter);     // spatial min-filter after aggregation (shiftable window)
    o.PPPD(aggr_subpixel);      // do local fits around minima to get lower value
    o.PPPD(aggr_collapse);      // collapse DSI back to integer disp sampling
    o.PPPF(diff_lambda);        // lambda parameter for diffusion algorithms
    o.PPPF(diff_beta);          // beta parameter for membrane model diffusion
    o.PPPF(diff_scale_cost);    // scale for m_cost values (necessary for Bayesian diffusion)
    o.PPPF(diff_mu);            // mu parameter for Bayesian diffusion
    o.PPPF(diff_sigmaP);        // sigma of robust prior for Bayesian diffusion
    o.PPPF(diff_epsP);          // epsilon of robust prior for Bayesian diffusion

    o.PPPD(aggr_color_space);     // color space
    o.PPPF(aggr_gamma_proximity); // gamma value for proximity
    o.PPPF(aggr_gamma_similarity);    // gamma value for similarity

    // Parameters controlling the optimization algorithms:
    o.PPPD(opt_fn);             // optimization function
    o.PPPF(opt_smoothness);     // smoothness penalty multiplier
    o.PPPF(opt_grad_thresh);    // threshold for magnitude of intensity gradient
    o.PPPF(opt_grad_penalty);   // smoothness penalty factor if gradient is too small
    o.PPPD(opt_occlusion_cost); // cost for occluded pixels in DP algorithm
    o.PPPD(opt_max_iter);       // maximum number of optimization iterations
    o.PPPD(opt_random);         // randomize optimization (disparity and/or pixel sites)
    o.PPPD(opt_sa_var);         // simulated annealing variant (update rule)
    o.PPPF(opt_sa_start_T);     // starting temperature
    o.PPPF(opt_sa_end_T);       // ending temperature
    o.PPPD(opt_sa_schedule);    // annealing schedule
    o.PPPF(opt_min_margin);     // neccessary margin of "clear" min in symmetric matcher
    o.PPPD(opt_sym_passes);     // number of outer iterations with further aggregation 

    // Parameters controlling the final sub-pixel refinement step:
    o.PPPD(refine_subpix);      // fit sub-pixel value to local correlation

    // Parameters controlling the quality / fitness evaluation
    o.PPPD(eval_ignore_border); // number of border pixels to ignore in ground truth image
    o.PPPF(eval_bad_thresh);    // acceptable disparity error (for counting "bad points")
    o.PPPF(eval_error_scale);   // scale disparity errors and write out (along with bad pixels)
    o.PPPD(eval_lin_interp);    // use linear interpolation in forward warping
    o.PPPF(eval_disp_gap);      // don't interpolate across disparity jumps bigger than this
    o.PPPD(eval_predict_type);  // type of prediction error to compute
    o.PPPD(eval_textureless_width);     // width of box filter for summing squared horiz. gradients
    o.PPPF(eval_textureless_thresh);    // threshold applied to summed (h-grad)^2
    o.PPPD(eval_discont_width); // width of discontinuity region (box filter)
    o.PPPD(eval_predict_diff);  // write out difference images instead of resampled
    o.PPPH(eval_empty_color);   // color of empty pixels
    o.PPPF(eval_partial_shuffle);   // use interval analysis for prediction error
    o.PPPD(eval_match_quality); // evaluate match quality (final cost and certainty)
    o.PPPD(eval_certain_matches_only);  // only compute statistics for "certain" matches

    // Output parameters: measure performance of algorithm
    o.PPPF(rms_error_all);          // RMS disparity error (all pixels)
    o.PPPF(rms_error_nonocc);       // RMS disparity error (non-occluded pixels only)
    o.PPPF(rms_error_occ);          // RMS disparity error (occluded pixels only)
    o.PPPF(rms_error_textured);     // RMS disparity error (textured pixels only)
    o.PPPF(rms_error_textureless);  // RMS disparity error (textureless pixels only)
    o.PPPF(rms_error_discont);      // RMS disparity error (near depth discontinuities)
    o.PPPF(bad_pixels_all);         // fraction of bad points (all pixels)
    o.PPPF(bad_pixels_nonocc);      // fraction of bad points (non-occluded pixels only)
    o.PPPF(bad_pixels_occ);         // fraction of bad points (occluded pixels only)
    o.PPPF(bad_pixels_textured);    // fraction of bad points (textured pixels only)
    o.PPPF(bad_pixels_textureless); // fraction of bad points (textureless pixels only)
    o.PPPF(bad_pixels_discont);     // fraction of bad points (near depth discontinuities)
    o.PPPF(fraction_matched);       // fraction of pixels with match status "certain"
    o.PPPF(predict_err_near);     // prediction error near   (e.g., frame 0)
    o.PPPF(predict_err_middle);   // prediction error middle (e.g., frame 4)
    o.PPPF(predict_err_match);    // prediction error match  (e.g., frame 6)
    o.PPPF(predict_err_far);      // prediction error far    (e.g., frame 8)

    o.PPPF(final_energy);         // final energy of solution
    o.PPPF(total_time);           // total time taken by algorithm

    // Miscellaneous parameters:
    o.PPPD(verbose);            // verbosity level
    o.PPPD(evaluate_only);      // read specified depth map and evaluate only
};
