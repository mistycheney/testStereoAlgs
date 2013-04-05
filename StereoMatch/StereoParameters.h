///////////////////////////////////////////////////////////////////////////
//
// NAME
//  StereoParameters.h -- parameter values controlling the stereo matcher
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
//  StereoParameters.cpp    implementation
//
// Copyright © Richard Szeliski and Daniel Scharstein, 2001.
// See Copyright.h for more details
//
///////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <vector>
#include "ParameterIO.h"
#include "Verbose.h"

//  Note:  we could also fold all of the parameters into a single struct,
//      if that turns out to be easier for people to follow...

enum EStereoMatchFn
{
    // Per-pixel matching function:

    eAD        = 1,    // absolute differences (for SAD)
    eSD        = 2,    // squared differences (for SSD)
    // we DON'T do normalized cross-correlation, since only good at feature points...
};

enum EStereoInterpFn
{
    // Resampling / interpolation function for sub-pixel shifts:

    eLinear     = 1,    // linear interpolation
    eCubic      = 3     // cubic interpolation
};

struct CStereoPreProcParameters
{
    // Parameters controlling the pre-processing of input images:

    float preproc_addnoise_stddev;  // additive noise standard deviation
    int preproc_blur_iter;          // number of pre-processing blur iterations
    
    // TODO: do we want any other kinds of pre-processing (hist.eq., rank stats.)?
};

struct CStereoMatchParameters
{
    // Parameters controlling the per-pixel matching:

    int frame_ref;              // reference image frame
    int frame_match;            // matching image frame

    int disp_min;               // smallest disparity
    int disp_max;               // largest disparity
    float disp_step;            // disparity step size (must be int or 1/int)
    int disp_n;                 // number of disparity levels (computed from above 3)
    float disp_scale;           // disparity scaling factor for gray depthmap

    EStereoMatchFn match_fn;      // matching function
    EStereoInterpFn match_interp; // interpolation function

    int match_max;              // maximum difference for truncated SAD/SSD
    int match_interval;         // use Birchfield/Tomasi matching measure
    int match_interpolated;     // interpolate both lines when matching
};

enum EStereoAggrFn
{
    // Spatial aggregation function:

    eBox        = 1,    // box filter (square window)
    eBinomial   = 2,    // binomial filter (1 4 6 4 1)

    // nonlinear diffusion [Scharstein & Szeliski, IJCV 98]
    // these are included here because they should be iterated (aggr_iter > 1), and followed
    // by winner-take-all optimization

    eDiffusion  = 3,    // simple diffusion
    eMembrane   = 4,    // membrane model
    eBayesian   = 5,    // Bayesian model
    eASWeight   = 6,    // locally adaptive support weight (see LASW.h and LASW.cpp)  
    // how about double exponential (IIR), ...?
};

enum EStereoAggrCS
{
   // color space
   eRGB = 1,
   eCIELab  = 2,
};

struct CStereoAggregateParameters
{
    // Parameters controlling the spatial aggregation:

    EStereoAggrFn aggr_fn;  // aggregation function

    int aggr_window_size;   // size of window
    int aggr_iter;          // number of aggregation iterations

    int aggr_minfilter;     // spatial min-filter after aggregation (shiftable window)
                            // if 0, don't min-filter; otherwise use the value of the
                            // variable min_filter as the filter window size

    int aggr_subpixel;      // do local fits around minima to get lower value
    int aggr_collapse;      // collapse DSI back to integer disp sampling

    float diff_lambda;      // lambda parameter for diffusion algorithms
    float diff_beta;        // beta parameter for membrane model diffusion

    float diff_scale_cost;  // scale for m_cost values (necessary for Bayesian diffusion)
    float diff_mu;          // mu parameter for Bayesian diffusion
    float diff_sigmaP;      // sigma of robust prior for Bayesian diffusion
    float diff_epsP;        // epsilon of robust prior for Bayesian diffusion

    // adaptive support weight parameters (see LASW.h and LASW.cpp)
    EStereoAggrCS aggr_color_space;   // color space
    float aggr_gamma_proximity;       // gamma value in the proximity term
    float aggr_gamma_similarity;  // gamma value in the similarity term

};

enum EStereoOptimizeFn
{
    // Optimization technique used:

    eNoOpt       = 0,    // no optimization (pass through input depth maps)
    eWTA         = 1,    // winner-take-all (local minimum)
    eDynamicProg = 2,    // scanline dynamic programming
    eScanlineOpt = 3,    // iterative scanline DP with inter-scanline energies
    eGraphCut    = 4,    // graph-cut global minimization
    eSimulAnnl   = 5,    // simulated annealing
	eSymmetric   = 6,	 // symmetric winner-take-all *NEW*
    eBPSync      = 7,
    eBPAccel     = 8,
};

enum EStereoSAVariant
{
    // Update rule used
    eMetropolis = 1,    // Metropolis algorithm: accept downhill, sometime uphill
    eFlipGibbs  = 2,    // accept move with Gibbs distribution probability
    eFullGibbs  = 3,    // compute all possible disparity changes at a pixel
};

enum EStereoSASchedule
{
    // Annealing schedule used
    eSALinear   = 1,    // linear schedule (start..end)
    eSALog      = 2,    // logarithmic schedule
};

struct CStereoOptimizeParameters
{
    // Parameters controlling the optimization algorithms:

    EStereoOptimizeFn opt_fn;       // optimization function

    float opt_smoothness;           // smoothness penalty multiplier (lambda)
    float opt_grad_thresh;          // threshold for magnitude of intensity gradient
    float opt_grad_penalty;         // smoothness penalty factor if gradient is too small
    int opt_occlusion_cost;         // cost for occluded pixels in DP algorithm
    int opt_max_iter;               // maximum number of optimization iterations
    int opt_random;                 // randomize optimization (disparity and/or pixel sites)

    EStereoSAVariant opt_sa_var;        // simulated annealing variant (update rule)
    float opt_sa_start_T;               // starting temperature
    float opt_sa_end_T;                 // ending temperature
    EStereoSASchedule opt_sa_schedule;  // annealing schedule

    float opt_min_margin;           // neccessary margin of "clear" min in symmetric matcher
    int opt_sym_passes;             // number of outer iterations with further aggregation 
};

struct CStereoRefineParameters
{
    // Parameters controlling the final sub-pixel refinement step:

    int refine_subpix;              // fit sub-pixel value to local correlation
};

enum EStereoPredictionType
{
    // What type of prediction error to compute
    ePredictNone    = 0,    // skip prediction error computation
    ePredictForward = 1,    // compute forward prediction error
    ePredictInverse = 2,    // compute inverse prediction error
};

struct CStereoEvaluateParameters
{
    // Parameters controlling the quality / fitness evaluation

    int eval_ignore_border;         // number of border pixels to ignore in ground truth
                                    // image (set to 18 for Ohta data set)
    float eval_bad_thresh;          // acceptable disparity error (for counting "bad points")
    float eval_error_scale;         // scale disparity errors and write out (along with bad pixels)

    int eval_lin_interp;            // use linear interpolation in forward warping
    float eval_disp_gap;            // don't interpolate across disparity jumps bigger than this
    EStereoPredictionType eval_predict_type; // type of prediction error to compute
    int eval_textureless_width;              // width of box filter for summing squared horiz. gradients
    float eval_textureless_thresh;           // threshold applied to summed (h-grad)^2
    int eval_discont_width;         // width of discontinuity region (box filter)
    int eval_empty_color;           // color of empty pixels
    int eval_predict_diff;          // write out difference images instead of resampled
    float eval_partial_shuffle;     // use interval analysis for prediction error

    int eval_match_quality;         // evaluate match quality (final cost and certainty)
    int eval_certain_matches_only;  // only compute statistics for matches with status "certain"
};

struct CStereoOutputParameters
{
    // Output parameters: measure performance of algorithm

    float rms_error_all;            // RMS disparity error (all pixels)
    float rms_error_nonocc;         // RMS disparity error (non-occluded pixels only)
    float rms_error_occ;            // RMS disparity error (occluded pixels only)
    float rms_error_textured;       // RMS disparity error (textured pixels only)
    float rms_error_textureless;    // RMS disparity error (textureless pixels only)
    float rms_error_discont;        // RMS disparity error (near depth discontinuities)
    float bad_pixels_all;           // fraction of bad points (all pixels)
    float bad_pixels_nonocc;        // fraction of bad points (non-occluded pixels only)
    float bad_pixels_occ;           // fraction of bad points (occluded pixels only)
    float bad_pixels_textured;      // fraction of bad points (textured pixels only)
    float bad_pixels_textureless;   // fraction of bad points (textureless pixels only)
    float bad_pixels_discont;       // fraction of bad points (near depth discontinuities)

    float fraction_matched;         // fraction of pixels with match status "certain"

    // prediction error for 4 reference positions: -0.5, 0.5, 1.0 (match frame position), 1.5
    // for example, if frame_ref=2 and frame_match=6, we predict frames 0, 4, 6, 8
    float predict_err_near;         // prediction error near   (e.g., frame 0)
    float predict_err_middle;       // prediction error middle (e.g., frame 4)
    float predict_err_match;        // prediction error match  (e.g., frame 6)
    float predict_err_far;          // prediction error far    (e.g., frame 8)

    float final_energy;             // final energy of solution
    float total_time;               // total time taken by algorithm

};

struct CStereoParameters : public
    CStereoPreProcParameters,
    CStereoMatchParameters,
    CStereoAggregateParameters,
    CStereoOptimizeParameters,
    CStereoRefineParameters,
    CStereoEvaluateParameters,
    CStereoOutputParameters
{
    // miscellaneous parameters:

    EVerbosityLevel verbose;    // verbosity level (see Verbose.h)

    int evaluate_only;          // read specified depth map and evaluate only

    CStereoParameters();
        // Constructor
    void ReInitializeSeqParams();
        // (Re-)Initialize Parameters specific to each image sequence
    void ReInitializeAlgParams();
        // (Re-)Initialize Parameters specific to each algorithm
    void ResetOutputParams();
        // Reset output parameters that measure performance of algorithm

    void PIOInitialize(CParameterIO& prw);
        // Initialize parameter input/output object
};

