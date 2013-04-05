///////////////////////////////////////////////////////////////////////////
//
// NAME
//  StcEvaluate.cpp -- evaluate the quality of the matching score
//
// DESCRIPTION
//  Evaluate the quality of the computed disparity map m_float_disparities using
//  the ground truth m_true_disparity.  Compute statistics and print them to stderr.
//  Also compute occlusion map from the ground truth disparities and include
//  statistics for non-occluded pixels.
// 
//  Relationship between float disparities d and "actual" (scaled) disparities s_d
//  for occlusion computation and warping:
//      s_d  =  m_frame_diff_sign * d
//  The coordinate of a pixel x_m in the matching frame corresponding
//  to a reference pixel x_r in the reference frame is given by
//      x_m = x_r - s_d
//  (this is because input images are ordered left to right, so that pixel
//  motion is leftward).
//
// SEE ALSO
//  StereoMatcher.h
//  StereoMatcher.cpp
//
// Copyright © Richard Szeliski and Daniel Scharstein, 2001.
// See Copyright.h for more details
//
// *** fixed bug with occlusion masks
//
///////////////////////////////////////////////////////////////////////////

#include "Error.h"
#include "StereoMatcher.h"
#include "Convert.h"
#include "BoxFilter.h"
#include "MinFilter.h"
#include "Warp1D.h"
#include "Histogram1D.h"
#include "ImageIO.h"
#include <math.h>

CByteImage CStereoMatcher::ComputeOcclusion(int frame)
{
    // Compute an occlusion map
    float fractional_shift = (frame - frame_ref) / fabs((float)m_frame_diff);

    // Simply forward warp a white image, then inverse warp it.
    //  TODO:  another approach would be to forward warp the disparity image,
    //      then to splat it again and compare the destination.
    CShape sh = m_true_disparity.Shape();

    // Forward-warped depth map (with marked empty pixels)
    CFloatImage fwd_depth(sh);
    static float invalid_depth = -9999.0f;
    fwd_depth.FillPixels(invalid_depth);
    ForwardWarp(m_true_disparity, fwd_depth, m_true_disparity,
                fractional_shift, true, eval_disp_gap);
    
    // for "empty" pixels in fwd_depth, need to make pixels black in "white" image.
    // otherwise, these have disparity 0 and don't move, creating white areas in
    // occluded areas
	CByteImage white(sh);
    white.FillPixels(255);
    int x, y, w=sh.width, h=sh.height;
    for (y = 0; y < h; y++) {
        uchar* wp = &white.Pixel(0, y, 0);
        float* fp = &fwd_depth.Pixel(0, y, 0);
        for (x = 0; x < w; x++)
            if (fp[x] == invalid_depth) // empty pixel
                wp[x] = 0;
    }

	// inverse warp marked white image to get occlusion mask
	CByteImage occlusion(sh);
	occlusion.ClearPixels();	// pixels that nothing maps to will stay 0
    ForwardWarp(white, occlusion, fwd_depth,
                -fractional_shift, true, eval_disp_gap);
    return occlusion;
}

void CStereoMatcher::ComputeOcclusions()
{
    // Call the general-purpose code
    m_occlusion = ComputeOcclusion(frame_match);

    // make borders black to indicate that statistics are not collected there
    CShape sh = m_true_disparity.Shape();
    int x, y, w=sh.width, h=sh.height;
    int b = eval_ignore_border;
    for (y = 0; y < h; y++) {
        uchar* occ = &m_occlusion.Pixel(0, y, 0);
        for (x = 0; x < w; x++)
            if (x < b || x >= w-b || y < b || y >= h-b)
                occ[x] = 0; // CHANGE: make borders black! (was 255);
    }
    // occlusion now contains a visibility map (white = visible)

    // write out warped depth map and occlusion map for debugging
    if (verbose >= eVerboseDumpFiles) { 
        CByteImage fwd(sh);
        // ScaleAndOffset(fwd_depth, fwd, disp_scale, -disp_min * disp_scale);
        // fprintf(stderr, "Writing fwd_depth map fwddepth.pgm\n");
        // WriteImage(fwd, "reprojected/fwddepth.pgm");
        fprintf(stderr, "Writing occlusion map occl.pgm\n");
        WriteImage(m_occlusion, "reprojected/occl.pgm");
    }
    ScaleAndOffset(m_occlusion, m_occlusion, -1, 255);
}

void CStereoMatcher::ComputeTextureless()
{
    // Compute textureless regions:
    //  1. evaluate squared intensity/color gradient
    //  2. apply a box filter of eval_textureless_width
    //  3. threshold using eval_textureless_thresh
    CShape sh = m_reference.Shape();
    int w = sh.width, h = sh.height;
    int nB = sh.nBands, nC = nB - (nB > 1);
    sh.nBands = 1;  // single channel for m_textureless
    m_textureless.ReAllocate(sh);
    CFloatImage sq_grad(sh);    // filtered squared horizontal gradient

    // Compute the squared horizontal gradient
    int x, y, max_sum2 = 0;     // max_sum2 for debugging only
    for (y = 0; y < h; y++)
    {
        uchar* p = &m_reference.Pixel(0, y, 0);
        float* s = &sq_grad.Pixel(0, y, 0);
        for (x = 0; x < w-1; x++, p += nB)
        {
            float sum2 = 0;
            for (int b = 0; b < nC; b++) {
                float diff = int(p[b]) - int(p[b+nB]);// this is NOT centered
                sum2 += diff * diff;
            }
            sum2 /= (float)nC;
            s[x+1] = sum2;
            if (x == 0)
                s[x] = sum2;
            s[x] = __max(sum2, s[x]);  // makes it centered
            max_sum2 = __max((int)sum2, max_sum2);
        }
    }

    // Aggregate with a box filter
    if (eval_textureless_width > 0)
        BoxFilter(sq_grad, sq_grad, eval_textureless_width, eval_textureless_width, true);

    // Threshold to get final map
    float squared_thresh = eval_textureless_thresh * eval_textureless_thresh;
    for (y = 0; y < h; y++)
    {
        float* s = &sq_grad.Pixel(0, y, 0);
        uchar* p = &m_textureless.Pixel(0, y, 0);
        for (x = 0; x < w; x++)
            p[x] = (s[x] < squared_thresh) ? 255 : 0;
    }

    // Write out textureless map for debugging
    if (verbose >= eVerboseDumpFiles) { 
        CByteImage tl;
        // Write out combined textureless and occlusion map
        // black where occluded, grey where textured
        // -> only white areas are counted in textureless statistics!
        ScaleAndOffset(m_textureless, tl, 0.5, 128); // make textured pixels grey
        for (y = 0; y < h; y++)
        {
            uchar* t = &tl.Pixel(0, y, 0);
            uchar* o = &m_occlusion.Pixel(0, y, 0);
            for (x = 0; x < w; x++)
                if (o[x])
                    t[x] = 0;        // make occluded pixels black
        }
        fprintf(stderr, "Writing occ_and_textl.pgm\n");
        WriteImage(tl, "reprojected/occ_and_textl.pgm");

#if 0   // probably won't need the following anymore (old code):
        ScaleAndOffset(m_textureless, tl, -1, 255);    // black on white (for printing)
        fprintf(stderr, "Writing textureless.pgm\n");
        WriteImage(tl, "reprojected/textureless.pgm");
        ScaleAndOffset(sq_grad, tl, 1, 0);
        fprintf(stderr, "Writing sqHorizDiff.pgm\n");
        WriteImage(tl, "reprojected/sqHorizDiff.pgm");
#endif
    }
}

// older version - newer version below gives link error in second
// BoxFilter call under linux for some reason
void CStereoMatcher::ComputeDisparityDiscont()
{
    // Compute disparity discontinuities
    //  1. threshold horiz. and vert. depth discontinuities with eval_disp_gap
    //  2. apply a box filter of eval_discont_width
    //  3. re-threshold above 0
    CShape sh = m_true_disparity.Shape();
    int w = sh.width, h = sh.height;
    sh.nBands = 1;  // single channel for m_textureless
    m_depth_discont.ReAllocate(sh, false);
    CIntImage d_disc(sh);    // filtered squared horizontal gradient

    // Compute the squared horizontal gradient
    int x, y;
    int bor = eval_ignore_border+1;
    for (y = 0; y < h-1; y++)
    {
        float* p0 = &m_true_disparity.Pixel(0, y, 0);
        float* p1 = &m_true_disparity.Pixel(0, y+1, 0);
        int*   d0 = &d_disc.Pixel(0, y, 0);
        int*   d1 = &d_disc.Pixel(0, y+1, 0);

        // Clear to 0
        if (y == 0)
            memset(d0, 0, w*sizeof(int));
        memset(d1, 0, w*sizeof(int));

        // ignore borders
        if (y < bor || y >= h-bor)
            continue;

        // Find the discontinuities
        for (x = bor; x < w-bor-1; x++)
        {
            float h_diff = fabs(p0[x] - p0[x+1]);
            if (h_diff >= eval_disp_gap)
                d0[x] = d0[x+1] = 255;
            float v_diff = fabs(p0[x] - p1[x]);
            if (v_diff >= eval_disp_gap)
                d0[x] = d1[x] = 255;
        }
    }

    // Aggregate with a box filter
    if (eval_discont_width > 0)
        BoxFilter(d_disc, d_disc, eval_discont_width, eval_discont_width, false);

    // Threshold to get final map
    for (y = 0; y < h; y++)
    {
        int*   d = &d_disc.Pixel(0, y, 0);
        uchar* p = &m_depth_discont.Pixel(0, y, 0);
        for (x = 0; x < w; x++)
            p[x] = (d[x] != 0) ? 255 : 0;
    }

    // Write out depth discontinuity map for debugging
    if (verbose >= eVerboseDumpFiles) { 
        CByteImage dd;
        // Write out combined discontinuity and occlusion map
        // black where occluded, grey where no discontinuities
        // -> only white areas are counted in discontinuity statistics!
        ScaleAndOffset(m_depth_discont, dd, 0.5, 128);
        for (y = 0; y < h; y++)
        {
            uchar* d = &dd.Pixel(0, y, 0);
            uchar* o = &m_occlusion.Pixel(0, y, 0);
            for (x = 0; x < w; x++)
                if (o[x])
                    d[x] = 0;        // make occluded pixels black
        }
        fprintf(stderr, "Writing occ_and_discont.pgm\n");
        WriteImage(dd, "reprojected/occ_and_discont.pgm");

#if 0   // probably won't need the following anymore (old code):
        fprintf(stderr, "Writing depth_discont.pgm\n");
        ScaleAndOffset(m_depth_discont, dd, -1, 255);   // black on white (for printing)
        WriteImage(dd, "reprojected/depth_discont.pgm");
#endif
    }
}

// this is the newer version but it doesn't want to link under linux...
/*
void CStereoMatcher::ComputeDisparityDiscont_new()
{
    // Compute disparity discontinuities
    //  1. threshold horiz. and vert. depth discontinuities with eval_disp_gap
    //  2. apply a box filter of eval_discont_width
    //  3. re-threshold above 0
    CShape sh = m_true_disparity.Shape();
    int w = sh.width, h = sh.height;
    sh.nBands = 1;  // single channel for m_textureless
    m_depth_discont.ReAllocate(sh);
    m_depth_discont.ClearPixels();

    int x, y;
    // ignore borders
    int bor = eval_ignore_border+1;
    for (y = bor; y < h-bor-1; y++)
    {
        float* p0 = &m_true_disparity.Pixel(0, y, 0);
        float* p1 = &m_true_disparity.Pixel(0, y+1, 0);
        uchar* d0 = &m_depth_discont.Pixel(0, y, 0);
        uchar* d1 = &m_depth_discont.Pixel(0, y+1, 0);
        uchar* occ = &m_occlusion.Pixel(0, y, 0);

        // Find the discontinuities
        for (x = bor; x < w-bor-1; x++)
        {
            float h_diff = fabs(p0[x] - p0[x+1]);
            if (h_diff >= eval_disp_gap)
                d0[x] = d0[x+1] = 255;
            float v_diff = fabs(p0[x] - p1[x]);
            if (v_diff >= eval_disp_gap)
                d0[x] = d1[x] = 255;

            // also mark all occluded areas 
			// TODO: this may make sense to always turn on!
            static bool mark_occlusions = false;
            if (mark_occlusions && occ[x])
                d0[x] = 255;

        }
    }

    // Aggregate with a box filter
    if (eval_discont_width > 0) {
        BoxFilter(m_depth_discont, m_depth_discont, eval_discont_width, eval_discont_width, false);

    }
    // Threshold to get final map
    for (y = 0; y < h; y++)
    {
        uchar* p = &m_depth_discont.Pixel(0, y, 0);
        for (x = 0; x < w; x++)
            p[x] = (p[x] != 0) ? 255 : 0;
    }

    // Write out depth discontinuity map for debugging
    if (verbose >= eVerboseDumpFiles) { 
        CByteImage dd;
        // Write out combined discontinuity and occlusion map
        // black where occluded, grey where no discontinuities
        // -> only white areas are counted in discontinuity statistics!
        ScaleAndOffset(m_depth_discont, dd, 0.5, 128);
        for (y = 0; y < h; y++)
        {
            uchar* d = &dd.Pixel(0, y, 0);
            uchar* o = &m_occlusion.Pixel(0, y, 0);
            for (x = 0; x < w; x++)
                if (o[x])
                    d[x] = 0;        // make occluded pixels black
        }
        fprintf(stderr, "Writing occ_and_discont.pgm\n");
        WriteImage(dd, "reprojected/occ_and_discont.pgm");

#if 0	// write out textured nonoccluded nondiscont. maps for matching cost paper
		// make textured pixels black, everything else white
        ScaleAndOffset(m_textureless, dd, 1, 0);
        for (y = 0; y < h; y++)
        {
            uchar* d = &dd.Pixel(0, y, 0);
            uchar* o = &m_occlusion.Pixel(0, y, 0);
	        uchar* disc = &m_depth_discont.Pixel(0, y, 0);
            for (x = 0; x < w; x++)
                if (disc[x] || o[x])
                    d[x] = 255;
        }
        fprintf(stderr, "Writing textured.pgm\n");
        WriteImage(dd, "reprojected/textured.pgm");
#endif

#if 0   // probably won't need the following anymore (old code):
        fprintf(stderr, "Writing depth_discont.pgm\n");
        ScaleAndOffset(m_depth_discont, dd, -1, 255);   // black on white (for printing)
        WriteImage(dd, "reprojected/depth_discont.pgm");
#endif
    }
}
*/


void CStereoMatcher::ComputeDisparityErrors()
{
    // Compute disparity errors
    // evaluation of depth map:
    // compare m_float_disparity to m_true_disparity

    CShape sh = m_float_disparity.Shape();

    // check if ground truth is present
    if (m_true_disparity.Shape() != sh)
        throw CError("Evaluate: invalid ground truth\n");

    int w = sh.width, h = sh.height;

    // Allocate error images
    bool error_images = (eval_error_scale > 0.0f);
    if (error_images)
    {
        m_disparity_error.ReAllocate(sh);   // scaled error in disparities
        m_bad_pixels.ReAllocate(sh);        // pixels flagged with bad disparities
        m_disparity_error.FillPixels(128);  // neutral gray (outside boundary)
        m_bad_pixels.FillPixels(255);        // white (outside boundary)

        // Compute the disparity histogram
        static bool compute_histogram_image = false;  // enable in debug mode
        if (compute_histogram_image)
        {
            float d_min = disp_min, d_max = disp_min + 256 / disp_scale;
            float d_step = 0, vscale = 0;
            int width = 256, height = 256;
            Histogram1D(m_float_disparity, 0, CByteImage(), m_disparity_histogram,
                        d_min, d_max, d_step, width, height, vscale);
        }
    }

    // count all pixels, even when only evaluating "certain" matches
    int total_cnt_all = 0;

    // statistics for all (certain) pixels
    float sum_sq_diff = 0.0;
    int bad_pix_cnt = 0;
    int total_cnt = 0;

    // statistics for non-occluded pixels
    float nocc_sum_sq_diff = 0.0;
    int nocc_bad_pix_cnt = 0;
    int nocc_total_cnt = 0;

    // statistics for occluded pixels
    float occ_sum_sq_diff = 0.0;
    int occ_bad_pix_cnt = 0;
    int occ_total_cnt = 0;

    // statistics for textured pixels
    float texd_sum_sq_diff = 0.0;
    int texd_bad_pix_cnt = 0;
    int texd_total_cnt = 0;

    // statistics for textureless pixels
    float texless_sum_sq_diff = 0.0;
    int texless_bad_pix_cnt = 0;
    int texless_total_cnt = 0;
    
    // statistics for depth discontinuity pixels
    float disc_sum_sq_diff = 0.0;
    int disc_bad_pix_cnt = 0;
    int disc_total_cnt = 0;

    // can't evaluate certain matches if no status image available
    if (m_status.Shape().width == 0)
        eval_certain_matches_only = 0;

    for (int y = eval_ignore_border; y < h - eval_ignore_border; y++)
    {
        float* disp = &m_float_disparity.Pixel(0, y, 0);
        float* trud = &m_true_disparity.Pixel(0, y, 0);
        uchar* ocnt = &m_occlusion.Pixel(0, y, 0);
        uchar* texl = &m_textureless.Pixel(0, y, 0);
        uchar* disc = &m_depth_discont.Pixel(0, y, 0);
        uchar* erri = (error_images) ? &m_disparity_error.Pixel(0, y, 0) : 0;
        uchar* badp = (error_images) ? &m_bad_pixels.Pixel(0, y, 0) : 0;
        uchar* stat = eval_certain_matches_only ? &m_status.Pixel(0, y, 0) : 0;

        for (int x = eval_ignore_border; x < w - eval_ignore_border; x++)
        {
            float diff = disp[x] - trud[x];
            if (error_images)
            {
                int v = 128 + int(diff * eval_error_scale * disp_scale + 0.5f);
                erri[x] = __max(0, __min(v, 255));
            }

            total_cnt_all++;

            // collect statistics for certain matches only?
            if (eval_certain_matches_only && stat[x] != eCertainMatch) {
                if (error_images) {
                    erri[x] = 128;  // mark as "no error"
                    badp[x] = 255;
                }
                continue;
            }

            // statistics for all pixels
            sum_sq_diff += diff*diff;
            int bad = (fabs(diff) > eval_bad_thresh);
            bad_pix_cnt += bad;

            // If not textured, don't draw it (optional for debugging & matching cost paper)
            static bool only_textured = false;
            bool skip_bad = (only_textured && (ocnt[x] > 1 || texl[x] != 0));

            // If near discontinuity or occluded, NEITHER COUNT NOR DRAW IT
			// (optional for matching cost paper)
            static bool no_discontinuities = false;
            skip_bad = skip_bad || (no_discontinuities && disc[x] != 0 );

            if (error_images && ! skip_bad)
                badp[x] = (bad) ? 0 : 255;

            total_cnt++;

            // statistics for occluded pixels
            if (ocnt[x] > 1) {
                occ_sum_sq_diff += diff*diff;
                occ_bad_pix_cnt += bad;
                occ_total_cnt++;

                if (error_images)
                    badp[x] = __min(255, badp[x]+200); // "grey out" occluded pixels in bad pixel map
            } else {
            // statistics for nonoccluded pixels

                nocc_sum_sq_diff += diff*diff;
                nocc_bad_pix_cnt += bad;
                nocc_total_cnt++;

				// NOTE: collect all other statistics in non-occluded
				// regions only!
				
				// statistics for textured/textureless pixels
				if (texl[x] != 0) {
					texless_sum_sq_diff += diff*diff;
					texless_bad_pix_cnt += bad;
					texless_total_cnt++;
				} else if (!(no_discontinuities && disc[x] != 0)) {
					texd_sum_sq_diff += diff*diff;
					texd_bad_pix_cnt += bad;
					texd_total_cnt++;
				}
				
				// statistics for depth discontinuity pixels
				if (disc[x] != 0) {
					disc_sum_sq_diff += diff*diff;
					disc_bad_pix_cnt += bad;
					disc_total_cnt++;
				}
			}
        }
    }

    total_cnt += (total_cnt == 0);
    total_cnt_all += (total_cnt_all == 0);

    fraction_matched = (float)total_cnt / (float)total_cnt_all;

    rms_error_all = sqrt(sum_sq_diff / (float)total_cnt);
    bad_pixels_all = (float)bad_pix_cnt / (float)total_cnt;
    
    nocc_total_cnt += (nocc_total_cnt == 0);
    rms_error_nonocc = sqrt(nocc_sum_sq_diff / (float)nocc_total_cnt);
    bad_pixels_nonocc = (float)nocc_bad_pix_cnt / (float)nocc_total_cnt;

    occ_total_cnt += (occ_total_cnt == 0);
    rms_error_occ = sqrt(occ_sum_sq_diff / (float)occ_total_cnt);
    bad_pixels_occ = (float)occ_bad_pix_cnt / (float)occ_total_cnt;
   
    texd_total_cnt += (texd_total_cnt == 0);
    rms_error_textured = sqrt(texd_sum_sq_diff / (float)texd_total_cnt);
    bad_pixels_textured = (float)texd_bad_pix_cnt / (float)texd_total_cnt;
    
    texless_total_cnt += (texless_total_cnt == 0);
    rms_error_textureless = sqrt(texless_sum_sq_diff / (float)texless_total_cnt);
    bad_pixels_textureless = (float)texless_bad_pix_cnt / (float)texless_total_cnt;

    disc_total_cnt += (disc_total_cnt == 0);
    rms_error_discont = sqrt(disc_sum_sq_diff / (float)disc_total_cnt);
    bad_pixels_discont = (float)disc_bad_pix_cnt / (float)disc_total_cnt;

    if (verbose >= eVerboseSummary) {
        if (verbose >= eVerboseProgress) {
            fprintf(stderr, "Results");
            if (eval_ignore_border > 0)
                fprintf(stderr, " (ignoring borders of %d pixels)", eval_ignore_border );
        }
        if (eval_certain_matches_only)
            fprintf(stderr, ":\n  Certain matches only (%5.2f%% of all pixels)", 100.0f*fraction_matched);
        fprintf(stderr, ":\n  ALL   NON OCCL   OCCL   TEXTRD TEXTRLS D_DISCNT\n");
        fprintf(stderr, "%7.2f %7.2f %7.2f %7.2f %7.2f %7.2f  RMS disparity error\n",
                rms_error_all, rms_error_nonocc, rms_error_occ,
                rms_error_textured, rms_error_textureless, rms_error_discont);
        fprintf(stderr, "%7.2f%%%7.2f%%%7.2f%%%7.2f%%%7.2f%%%7.2f%% bad pixels (disp error > %g)\n",
                100.0f * bad_pixels_all, 100.0f * bad_pixels_nonocc, 100.0f * bad_pixels_occ,
                100.0f * bad_pixels_textured, 100.0f * bad_pixels_textureless, 100.0f * bad_pixels_discont,
                eval_bad_thresh);
        if (total_time >= 0.0f)
            fprintf(stderr, "  time = %gs", total_time);
        if (final_energy >= 0.0f)
            fprintf(stderr, "  energy = %.0f", final_energy);

        fraction_matched *= 0.1f; // bring into range of 0..10% so it can be plotted
                                  // with same scale as bad_pixel fractions

    }
}

static void PartialShuffle(CByteImage img_org,
                           CByteImage& img_min, CByteImage& img_max,
                           float shuffle_amt)
{
    // Compute 3x3 min and max images, then compute intervals that include
    //  original (center) pixel partially shifted towards min/max
    MinFilter(img_org, img_min, 3, 3);
    MaxFilter(img_org, img_max, 3, 3);

    // Compute the [min,max] intervals (linear interps of img and img_{min,max})
    CShape sh = img_org.Shape();
    int w = sh.width, h = sh.height, nB = sh.nBands;
    int n = w * nB;
    for (int y = 0; y < h; y++)
    {
        uchar *io = &img_org.Pixel(0, y, 0);
        uchar *in = &img_min.Pixel(0, y, 0);
        uchar *ix = &img_max.Pixel(0, y, 0);
        for (int x = 0; x < n; x++)
        {
            in[x] = int(io[x] + shuffle_amt * (in[x] - (int) io[x]));           // round down
            ix[x] = int(io[x] + shuffle_amt * (ix[x] - (int) io[x]) + 0.99f);   // round up
        }
    }

}
void CStereoMatcher::ComputePredictionError(CByteImage& predicted, CByteImage& original,
                                            float& rms_error, float& fraction_visible)
{
    // Compare two color images (synthesized and original)
    CShape sh = predicted.Shape();
    int w = sh.width, h = sh.height, nB = sh.nBands;
    int nC = nB - (nB > 1);     // actual number of color channels
    float sum2 = 0.0;
    int visible = 0, total = w * h;

    // Optional partial-shuffle transform to compensate for mis-alignments
    CByteImage pred_min, pred_max, orig_min, orig_max;
    bool shuffle = (eval_partial_shuffle > 0.0f);
    if (shuffle)
    {
        PartialShuffle(predicted, pred_min, pred_max, eval_partial_shuffle);
        PartialShuffle(original , orig_min, orig_max, eval_partial_shuffle);
        static bool dump_shuffled = false;
        if (dump_shuffled && verbose >= eVerboseDumpFiles) { 
            fprintf(stderr, "Writing partially shuffled images tmp_*_{min,max}.ppm\n");
            WriteImage(pred_min, "reprojected/tmp_pred_min.ppm");
            WriteImage(pred_max, "reprojected/tmp_pred_max.ppm");
            WriteImage(orig_min, "reprojected/tmp_orig_min.ppm");
            WriteImage(orig_max, "reprojected/tmp_orig_max.ppm");
        }
    }

    // Compute the differences and statistics
    for (int y = 0; y < h; y++)
    {
        uchar *p  = &predicted.Pixel(0, y, 0);
        uchar *o  = &original.Pixel(0, y, 0);
        for (int x = 0; x < w; x++, p += nB, o += nB)
        {
            // Early out:  predicted pixel is not visible (alpha != 255)
            //  note that the A channel in eval_empty_color MUST be 0!
            if (nB > 1 && p[nC] != 255)
                continue;
            visible++;

            // Process each band, accumulating the statistics
            for (int b = 0; b < nC; b++)
            {
                float diff = int(p[b]) - int(o[b]);

                // Partial shuffle transform:  analyze ranges
                if (shuffle)
                {
                    int pn = pred_min.Pixel(x, y, b);
                    int px = pred_max.Pixel(x, y, b);
                    int on = orig_min.Pixel(x, y, b);
                    int ox = orig_max.Pixel(x, y, b);
                    int xn = __max(pn, on); // max of mins
                    int nx = __min(px, ox); // min of maxs
                    if (xn <= nx)
                        diff = 0;   // overlapping ranges -> no error
                    else
                        diff = (pn > ox) ?          // check sign
                               pn - ox : on - px;   // gap between intervals
                }

                sum2 += diff * diff;

                // Optionally replace predicted image with scaled difference (error)
                if (eval_predict_diff && (nB == 1 || b < nB-1))
                {
                    p[b] = __max(0, __min(255, 128 + int(diff * eval_predict_diff)));
                }
            }
        }
    }

    // Compute the final errors
    rms_error = sqrt(sum2 / nC / (visible + (visible == 0)));
    fraction_visible = visible / (float) total;
}

void CStereoMatcher::ComputePredictionErrors()
{
    // Compute frame prediction (re-projection) errors
    static int inverse_warp_order = 3;      // TODO: put in CStereoEvaluateParameters?

    // Cycle through all of the images
    for (int f = 0; f < (int) m_frame.size(); f++)
    {
        CStereoFrame& frame = m_frame[f];
        CByteImage& original = frame.input_image;       // input image (gray or color)
        CByteImage& resampled= frame.resampled_image;   // resampled image (for reprojection error)

        // Fill with the empty color
        CShape sh = original.Shape();
        resampled.ReAllocate(sh, false);
        if (sh.nBands == 1)
            resampled.FillPixels(eval_empty_color);     // gray-level
        else
        {
            // Color image: make integer alias, and then fill (hack!)
            sh.nBands = 1;
            CIntImage aliased;
            aliased.ReAllocate(sh, (int *) &resampled.Pixel(0, 0, 0),
                               false, 0);
            aliased.FillPixels(eval_empty_color);
        }

        // Forward or inverse warp, then compute the errors
        float& predict_err = frame.predict_err;         // prediction error (visible pixels)
        float& predict_visible = frame.predict_visible; // fraction of pixels visible
        float fractional_shift = (f - frame_ref) / fabs((float)m_frame_diff);    // amount of shift
        if (eval_predict_type == ePredictForward)
        {
            ForwardWarp(m_reference, resampled, m_float_disparity,
                        fractional_shift, eval_lin_interp != 0, eval_disp_gap);
            ComputePredictionError(resampled, original, predict_err, predict_visible);
        }
        else
        {
            InverseWarp(original, resampled, m_float_disparity,
                        fractional_shift, eval_disp_gap, inverse_warp_order);
            ComputePredictionError(resampled, m_reference, predict_err, predict_visible);
            }

        // Print the error statistics
        if (verbose >= eVerbosePredictionError) {
            fprintf(stderr, " prediction error for frame %d: RMS error = %.2f, visible = %.2f%%\n",
                f, predict_err, predict_visible * 100);
        }

        // store prediction error in parameter variables
        // store errors for view positions (0=ref, 1=match)
        //   near:   -0.5  (or -1.0)
        //   middle   0.5 
        //   match:   1.0 
        //   far:     1.5  (or 2.0)
        // note: near assumes that only one one of view positions -1, -0.5 is present
        //        far assumes that only one one of view positions 1.5, 2.0 is present
        if (2*f ==  4*frame_ref - 2*frame_match ||
            2*f ==  3*frame_ref - 1*frame_match)  predict_err_near   = predict_err;
        // (2*f ==  2*frame_ref + 0*frame_match)  // reference frame: error == 0
        if (2*f ==  1*frame_ref + 1*frame_match)  predict_err_middle = predict_err;
        if (2*f ==  0*frame_ref + 2*frame_match)  predict_err_match  = predict_err;
        if (2*f == -1*frame_ref + 3*frame_match ||
            2*f == -2*frame_ref + 4*frame_match)  predict_err_far    = predict_err;
    }
}

void CStereoMatcher::ComputeMatchQuality()
{
    // Compute final matching cost and certainty
    CShape sh = m_cost.Shape(), sh1 = sh;
    int width = sh.width, height = sh.height;
    sh1.nBands = 1;

    // Allocate the final results
    m_final_cost.ReAllocate(sh1);   // best (lowest) cost at winning disparity
    m_certainty.ReAllocate(sh1);    // certainty (inv. variance) at winning disparity

    // Compute parabolic fit and store results
    //  (this code patterned after CStereoMatcher::Refine())
    float d_offset = disp_min;
    int nBands = (m_reference.Shape().nBands == 1) ? 1 : 3;
    for (int y = 0; y < height; y++)
    {
        float* cost  = &m_cost.Pixel(0, y, 0);
        int*   disp  = &m_disparity.Pixel(0, y, 0);
        float* fdisp = &m_float_disparity.Pixel(0, y, 0);
        float* fcost = &m_final_cost.Pixel(0, y, 0);
        float* fcert = &m_certainty.Pixel(0, y, 0);
        float* scert = &m_sub_pixel_cert.Pixel(0, y, 0);

        for (int x = 0; x < width; x++, cost += m_disp_n, scert += m_disp_n)
        {
            // Recompute disp[x] in case it's not there (evaluate_only)
            float d_new = fdisp[x];
            float d_sub = (d_new - d_offset) * m_disp_step_inv;
            disp[x] = int(d_sub + 0.5);
            float x0    = d_sub - disp[x];
            if (eval_match_quality == 2)
                x0 = 0.0;   // don't interpolate to sub-pixel

            // If we already computed sub-pixel cost and certainty, use that
            if (aggr_subpixel)
            {
                int d_min = disp[x];
                fcost[x] = cost[d_min];
                fcert[x] = scert[d_min];
            }
            else
            {

            // Get minimum, but offset by 1 from ends
            int d_min = disp[x] + (disp[x] == 0) - (disp[x] == m_disp_n-1);

            // Compute the equations of the parabolic fit
            float c0 = cost[d_min-1], c1 = cost[d_min], c2 = cost[d_min+1];
            float a = 0.5 * (c0 - 2.0 * c1 + c2);
            float b = 0.5 * (c2 - c0);

            // Degenerate parabola
            if (a <= 0.0 || a < 0.5 * fabs(b))
            {
                fcost[x] = c1;
                fcert[x] = 0.0f;
            }

            // Good parabola
            else
            {
                float ffit = a * x0 * x0 + b * x0 + c1;
                fcost[x] = ffit;
                fcert[x] = a;   // make this depend on overall variance??
            }
            }

            // Normalize by bands and optionally take sqrt
            float favg = fcost[x] / nBands;
            float fnew = (match_fn == eSD) ? sqrt(favg) : favg;
            fcost[x] = fnew;
        }
    }

    // Optionally write out the images
    if (verbose >= eVerboseDumpFiles) { 
        fprintf(stderr, "Writing final_cost.pgm and certainty.pgm\n");
        CByteImage fc(m_final_cost.Shape());

        static float cost_scale = 16.0f;
        ScaleAndOffset(m_final_cost, fc, cost_scale, 0);
        WriteImage(fc, "reprojected/final_cost.pgm");
        
        static float cert_scale = 0.5f;
        ScaleAndOffset(m_certainty, fc, cert_scale, 0);
        WriteImage(fc, "reprojected/certainty.pgm");

        // Compute the cost function histograms and write them out
        float c_min = 0.0f, c_max = 32.0f, c_step = 0.0f, vscale = 0.0f;
        int width = 256, height = 256;
        Histogram1D(m_final_cost, 0, CByteImage(), fc,
                    c_min, c_max, c_step, width, height, vscale);
        WriteImage(fc, "reprojected/final_cost_hist_all.pgm");

        // If we don't reset vscale = 0 below, uses _all's vertical scaling
        Histogram1D(m_final_cost, 0, m_occlusion, fc,
                    c_min, c_max, c_step, width, height, vscale);
        WriteImage(fc, "reprojected/final_cost_hist_occluded.pgm");
        Histogram1D(m_final_cost, 0, m_textureless, fc,
                    c_min, c_max, c_step, width, height, vscale);
        WriteImage(fc, "reprojected/final_cost_hist_textureless.pgm");
    }
}

void CStereoMatcher::ComputeStatusErrors()
{
    // Compute disparity errors grouped by status
    // evaluation of depth map:
    // compare m_float_disparity to m_true_disparity

    CShape sh = m_float_disparity.Shape();

    // check if ground truth is present
    if (m_true_disparity.Shape() != sh)
        throw CError("Evaluate: invalid ground truth\n");
 
    // check if status map is present
    if (m_status.Shape() != sh)
        throw CError("Evaluate: no status map available\n");

    int w = sh.width, h = sh.height;

    const int num_status = 1+eOccludedMatch;   // how many different status categories
                                               // TODO: make cleaner (how?)
    int total_cnt_all = 0;

    // statistics for each status
    float sum_sq_diff[num_status];
    int bad_pix_cnt[num_status];
    float bad_pix_percent[num_status];
    int total_cnt[num_status];
    
    int k;
    for (k = 0; k < num_status; k++) {
        sum_sq_diff[k] = 0.0f;
        bad_pix_cnt[k] = 0;
        total_cnt[k] = 0;
    }

    // statistics for occluded pixels (see if correctly identified)
    int occ_total_cnt = 0;
    int occ_false_positive = 0;     // wrongly labeled unoccluded pixel as occluded
    int occ_false_negative = 0;     // wrongly labeled occluded pixel as nonoccluded

    for (int y = eval_ignore_border; y < h - eval_ignore_border; y++)
    {
        float* disp = &m_float_disparity.Pixel(0, y, 0);
        float* trud = &m_true_disparity.Pixel(0, y, 0);
        uchar* ocnt = &m_occlusion.Pixel(0, y, 0);
        uchar* stat = &m_status.Pixel(0, y, 0);

        for (int x = eval_ignore_border; x < w - eval_ignore_border; x++)
        {
            float diff = disp[x] - trud[x];
            int bad = (fabs(diff) > eval_bad_thresh);
            total_cnt_all++;

            int status = stat[x];
            sum_sq_diff[status] += diff*diff;
            bad_pix_cnt[status] += bad;
            total_cnt[status]++;

            if (ocnt[x] > 1) {
                occ_total_cnt++;
                if (status != eOccludedMatch)
                    occ_false_negative++;
            } else {
                if (status == eOccludedMatch)
                    occ_false_positive++;
            }
        }
    }

    for (k = 0; k < num_status; k++) {
        total_cnt[k] += (total_cnt[k] == 0);
        sum_sq_diff[k] = sqrt(sum_sq_diff[k] / (float)total_cnt[k]);
        bad_pix_percent[k] = (float)bad_pix_cnt[k] / (float)total_cnt[k] * 100.0f;
    }
    occ_total_cnt += (occ_total_cnt == 0);
    float occ_fn = (float)occ_false_negative / (float)occ_total_cnt * 100.0f;
    float occ_fp = (float)occ_false_positive / (float)occ_total_cnt * 100.0f;

    if (verbose >= eVerboseSummary) {
        fprintf(stderr, "\n  UNK(%2d) CER(%2d) AMB(%2d) OCC(%2d) FNEG    FPOS\n",
            100*total_cnt[0]/total_cnt_all, 100*total_cnt[1]/total_cnt_all,
            100*total_cnt[2]/total_cnt_all, 100*total_cnt[3]/total_cnt_all);
        fprintf(stderr, "%7.2f %7.2f %7.2f %7.2f                  RMS disparity error\n",
                sum_sq_diff[0], sum_sq_diff[1], sum_sq_diff[2], sum_sq_diff[3]);
        occ_false_negative = 17;
        fprintf(stderr, "%7.2f%%%7.2f%%%7.2f%%%7.2f%%%7.2f%%%7.2f%% bad pixels (disp error > %g)\n",
                bad_pix_percent[0], bad_pix_percent[1], bad_pix_percent[2], bad_pix_percent[3],
                occ_fn, occ_fp, eval_bad_thresh);
    }
}


void CStereoMatcher::Evaluate()
{
    // Evaluate the quality of the matching score

    // Compute simple occlusion map using forward mapping
    ComputeOcclusions();

    // Compute textureless regions
    ComputeTextureless();

    // Compute disparity discontinuities
    ComputeDisparityDiscont();

    // Compute some simple statistics for now
    ComputeDisparityErrors();

    // Compute frame prediction (re-projection) errors
    if (eval_predict_type != ePredictNone)
        ComputePredictionErrors();

    // Compute best matching costs (final errors) and certainty
    if (eval_match_quality && ! evaluate_only)
        ComputeMatchQuality();

    // Compute disparity errors grouped by status
    if (m_status.Shape().width > 0)
        ComputeStatusErrors();
}
