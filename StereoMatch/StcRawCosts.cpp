///////////////////////////////////////////////////////////////////////////
//
// NAME
//  StcRawCosts.cpp -- compute raw per-pixel matching score between a pair of frames
//
// DESCRIPTION
//  The CStereoMatcher class implements a number of popular stereo
//  correspondence algorithms (currently for 2-frame rectified image pairs).
//
// The definition of disparity is as follows:
//  Disparity is the (floating point) inter-frame (instantaneous) displacement
//  of a pixel between successive frames in a multi-frame stereo pair.
//
//  In building the disparity-space image (DSI) m_cost, we enumerate a set of
//  m_disp_n disparitites, k = 0..m_disp_n-1.
//
//  The mapping between the integral disparities k and the floating point
//  disparity d is given by:
//      d = disp_min + k * m_disp_num / m_disp_den
//
//  Because we may be matching frames that are not adjacent in the sequence
//  we also define the frame difference f_d and scaled disparity s_d as:
//      f_d = frame_match - frame_ref
//      s_d = f_d * d
//  The coordinate of a pixel x_m in the matching frame corresponding
//  to a reference pixel x_r in the reference frame is given by
//      x_m = x_r - s_d
//  (this is because input images are ordered left to right, so that pixel
//  motion is leftward).
//
//  When we store the floating point disparities in a gray_level image, we
//  use the formulas
//      g_d = (d - disp_min) * disp_scale
//        d =  disp_min  + g_d / disp_scale
//  to do the conversions.
//
// IMPLEMENTATION
//  For better efficiency, we first interpolate the matching line
//  up by a factor of m_disp_den.
//
// SEE ALSO
//  StereoMatcher.h         longer description of this class
//
// Copyright © Richard Szeliski, 2001.
// See Copyright.h for more details
//
///////////////////////////////////////////////////////////////////////////

#include <assert.h>
#include "Error.h"
#include "StereoMatcher.h"
#include "Convert.h"
#include "ImageIO.h"
#include "Warp1D.h"
#include <time.h>

static void InterpolateLine(int buf[], int s, int w, int nB,
                            EStereoInterpFn match_interp)     // interpolation function
{
    // Interpolate the missing values
    float si = 1.0f / s;
    for (int x = 0; x < w-1; x++)
    {
        for (int b = 0; b < nB; b++)
        {
            int *v = &buf[s*x*nB+b];
            float I0 = v[0];
            float I1 = v[s*nB];
            if (match_interp == eCubic) // cubic interpolation
            {
                float Im = (x > 0) ? v[-s*nB] :
                           (I0 - (I1 - I0));  // extend linearly
                float Ip = (x+1 < w-1) ? v[2*s*nB] :
                           (I1 + (I1 - I0));  // extend linearly
                float sf = si;
                for (int is = 1; is < s; is++, sf += si)
                {
                    v  += nB;
                    float Ii = CubicInterpolate(sf, Im, I0, I1, Ip);
                    v[0] = int(Ii);
                }
            }
            else  // linear interpolation
            {
                float d = (I1 - I0) / (float) s;
                for (int is = 1; is < s; is++)
                {
                    v  += nB;
                    I0 += d;
                    v[0] = int(I0);
                }
            }
        }
    }
}

static void BirchfieldTomasiMinMax(const int* buffer, int* min_buf, int* max_buf,
                                   const int w, const int b)
{
    // Compute for every (interpolated) pixel, the minimum and maximum
    //  values in the two half-intervals before and after it
    //  (see [Birchfield & Tomasi, PAMI 20(40), April 1998, p. 401]).

    // Process each band separaterly
    for (int k = 0; k < b; k++)
    {
        int Ir = buffer[k], b1 = buffer[k];
        for (int x = 0, l = k; x < w; x++, l += b)
        {
            int Il = Ir, b0 = b1;   // shift down previously computed values
            if (x < w-1)
                b1 = buffer[l+b];
            Ir = (b0 + b1 + 1)/2;   // interpolated half-value
            min_buf[l] = __min(Il, __min(b0, Ir));
            max_buf[l] = __max(Il, __max(b0, Ir));
        }
    }
}

static bool undefined_cost = true;     // set this to true to pad with outside_cost

static void MatchLine(int w, int b, int interpolated,
                      int rmn[], int rmx[],     // min/max of ref (ref if rmx == 0)
                      int mmn[], int mmx[],     // min/max of mtc (mtc if mmx == 0)
                      float cost[],
                      int m_disp_n, int disp, int disp_den,
                      EStereoMatchFn match_fn,  // matching function
                      int match_max,            // maximum difference for truncated SAD/SSD
                      float match_outside)        // special value for outside match
{
    // Set up the starting addresses, pointers, and cutoff value
    int n = (w-1)*disp_den + 1;             // number of reference pixels
    int s = (interpolated) ? 1 : disp_den;     // skip in reference pixels
    std::vector<float> cost1;
    cost1.resize(n);
    int cutoff = (match_fn == eSD) ? match_max * match_max : abs(match_max);
    // TODO:  cutoff is not adjusted for the number of bands...
    const float bad_cost = -1;

    // Match valid pixels
    float  left_cost = bad_cost;
    float right_cost = bad_cost;
    int x, y;
    for (x = 0; x < n; x += s)
    {
        // Compute ref and match pointers
        cost1[x] = bad_cost;
        int x_r = x, x_m = x + disp;
        if (x_m < 0 || x_m >= n)
            continue;
        int* rn = &rmn[x_r*b];    // pointer to ref or min pixel(s)
        int* rx = &rmx[x_r*b];    // pointer to ref    max pixel(s)
        int* mn = &mmn[x_m*b];    // pointer to mtc or min pixel(s)
        int* mx = &mmx[x_m*b];    // pointer to mtc    max pixel(s)
        int  diff_sum = 0;        // accumulated error

        // This code could be special-cased for b==1 for more efficiency...
        for (int ib = 0; ib < b; ib++)
        {
            int diff1 = mn[ib] - rn[ib];    // straightforward difference
            if (rmx && mmx)
            {
                // Compare intervals (see partial shuffle code in StcEvaluate.cpp)
                int xn = __max(rn[ib], mn[ib]);     // max of mins
                int nx = __min(rx[ib], mx[ib]);     // min of maxs
                if (xn <= nx)
                    diff1 = 0;          // overlapping ranges -> no error
                else
                    diff1 = (mn[ib] > rx[ib]) ?     // check sign
                            mn[ib]-rx[ib] : 
                            rn[ib]-mx[ib];          // gap between intervals
            }
            int diff2 = (match_fn == eSD) ?    // squared or absolute difference
                            diff1 * diff1 : abs(diff1);
            diff_sum += diff2;
        }
        int diff3 = __min(diff_sum, cutoff);    // truncated difference
        if (left_cost == bad_cost)
            left_cost = diff3;  // first cost computed
        right_cost = diff3;     // last  cost computed
        cost1[x] = diff3;        // store in temporary array
    }

    // Fill in the left and right edges
    if (undefined_cost)
        left_cost = right_cost = match_outside;
    for (x = 0  ; x <  n && cost1[x] == bad_cost; x += s)
        cost1[x] = left_cost;
    for (x = n-1; x >= 0 && cost1[x] == bad_cost; x -= s)
        cost1[x] = right_cost;

    // Box filter if interpolated costs
    int dh = disp_den / 2;
    float box_scale = 1.0 / (2*dh + 1);
    for (x = 0, y = 0; y < w*m_disp_n; x += disp_den, y += m_disp_n)
    {
        if (interpolated && disp_den > 1)
        {
            float sum = 0;
            for (int k = -dh; k <= dh; k++)
            {
                int l = __max(0, __min(n-1, x+k));  // TODO: make more efficient
                sum += cost1[l];
            }
            cost[y] = int(box_scale * sum + 0.5);
        }
        else
            cost[y] = cost1[x];
    }
}

static int gcd(int a, int b)
{
    if (b>a)
        return gcd(b, a);
    if (b==0)
        return a;
    return gcd(b, a % b);
}

void CStereoMatcher::RawCosts()
{
    StartTiming();

    // Compute raw per-pixel matching score between a pair of frames
    CShape sh = m_reference.Shape();
    int w = sh.width, h = sh.height, b = sh.nBands;

    if (verbose >= eVerboseProgress)
        fprintf(stderr, "- computing costs: "); 
    if (verbose >= eVerboseSummary) {
        fprintf(stderr, match_fn == eAD ? "AD" : (match_fn == eSD ? "SD" : "???"));
        if (m_disp_step != 1.0f)
            fprintf(stderr, ", step=%g", m_disp_step);
        if (match_max < 1000)
            fprintf(stderr, ", trunc=%d", match_max);
        if (match_interval)
            fprintf(stderr, ", interval");
        if (match_interpolated)
            fprintf(stderr, ", interpolated");
    }
    if (verbose >= eVerboseProgress)
        fprintf(stderr, "\n");

    // Allocate a buffer for interpolated values
    //  Note that we don't have to interpolate the ref image if we
    //  aren't using match_interpolated, but it's simpler to code this way.
    match_interval = (match_interval ? 1 : 0);  // force to [0,1]
    int n_interp = m_disp_den * (w-1) + 1;
    std::vector<int> buffer0, buffer1, min_bf0, max_bf0, min_bf1, max_bf1;
    buffer0.resize(n_interp * b);
    buffer1.resize(n_interp * b);
    min_bf0.resize(n_interp * b);
    max_bf0.resize(n_interp * b);
    min_bf1.resize(n_interp * b);
    max_bf1.resize(n_interp * b);

    // Special value for border matches
    int worst_match = b * ((match_fn == eSD) ? 255 * 255 : 255);
    int cutoff = (match_fn == eSD) ? match_max * match_max : abs(match_max);
	m_match_outside = __min(worst_match, cutoff);	// trim to cutoff

    // Process all of the lines
    for (int y = 0; y < h; y++)
    {
        uchar* ref = &m_reference.Pixel(0, y, 0);
        uchar* mtc = &m_matching.Pixel(0, y, 0);
        int*  buf0 = &buffer0[0];
        int*  buf1 = &buffer1[0];
        int*  min0 = &min_bf0[0];
        int*  max0 = &max_bf0[0];
        int*  min1 = &min_bf1[0];
        int*  max1 = &max_bf1[0];

        // Fill the line buffers
        int x, l, m;
        for (x = 0, l = 0, m = 0; x < w; x++, m += m_disp_den*b)
        {
            for (int k = 0; k < b; k++, l++)
            {
                buf0[m+k] = ref[l];
                buf1[m+k] = mtc[l];
            }
        }

        // Interpolate the matching signal
        if (m_disp_den > 1)
        {
            InterpolateLine(buf1, m_disp_den, w, b, match_interp);
            InterpolateLine(buf0, m_disp_den, w, b, match_interp);
        }
        if (match_interval) {
            BirchfieldTomasiMinMax(buf1, min1, max1, n_interp, b);
            if (match_interpolated)
                BirchfieldTomasiMinMax(buf0, min0, max0, n_interp, b);
        }

        // Compute the costs, one disparity at a time
        for (int k = 0; k < m_disp_n; k++)
        {
            float* cost = &m_cost.Pixel(0, y, k);
            int disp = -m_frame_diff_sign * (m_disp_den * disp_min + k * m_disp_num);
            MatchLine(w, b, match_interpolated,
                      (match_interval) ? (match_interpolated) ? min0 : buf0 : buf0,
                      (match_interval) ? (match_interpolated) ? max0 : buf0 : 0,
                      (match_interval) ? min1 : buf1,
                      (match_interval) ? max1 : 0,
                      cost, m_disp_n, disp, m_disp_den,
                      match_fn, match_max, m_match_outside);
        }
    }
    PrintTiming();

    // Write out the different disparity images
    if (verbose >= eVerboseDumpFiles)
        WriteCosts(m_cost, "reprojected/RAW_DSI_%03d.pgm");
}

static void PadLine(int w, int b, float cost[],
                    int m_disp_n, int disp, int disp_den,
                    float match_outside)        // special value for outside match
{
    // Set up the starting addresses, pointers, and cutoff value
    int n = (w-1)*disp_den + 1;             // number of reference pixels
    int s = disp_den;                       // skip in reference pixels

	//  Hack: add -(s-1) to disp to make the left boundary 1 pixel wider!
	//  This is to account for possible interpolated pixels having mixed match_outside values
	//  TODO:  find a more principled solution (that also works for true occlusions)
	disp -= (s-1);

    // Fill invalid pixels
    for (int x = disp, y = 0; x < n+disp; x += s, y += b)
    {
        // Check if outside the bounds
        if (x < 0 || x >= n)
            cost[y] = match_outside;
    }
}

void CStereoMatcher::PadCosts()
{
    // Pad the matching scores with the previously computed m_match_outside value
    CShape sh = m_cost.Shape();
    int w = sh.width, h = sh.height, b = sh.nBands;

    // Process all of the lines
    for (int y = 0; y < h; y++)
    {
        // Compute the costs, one disparity at a time
        for (int k = 0; k < m_disp_n; k++)
        {
            float* cost = &m_cost.Pixel(0, y, k);
            int disp = -m_frame_diff_sign * (m_disp_den * disp_min + k * m_disp_num);
            PadLine(w, b, cost, m_disp_n, disp, m_disp_den, m_match_outside);
        }
    }
}
