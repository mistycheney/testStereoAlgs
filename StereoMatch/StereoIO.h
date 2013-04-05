///////////////////////////////////////////////////////////////////////////
//
// NAME
//  StereoIO.h -- stereo matcher file and parameter input/output
//
// DESCRIPTION
//  The CStereoIO class reads and write image and parameter files
//  needed to drive the CStereoMatcher class.
//
//  The parameters can be read from either a parameter file,
//  the command line, or both.
//
//  The input/output image filenames are given in a separate data
//  file (specified as parameters), and the output depth image name
//  can be overriden on the command line.
//
//  Scripting is also supported with the script command.
//  You can also change directories with the cd command.
//
//  Each "command" line either contains a list of parameter/value pairs,
//  which are read and stored for future use, or a command / argument.
//  The current commands that are supported are
//      script script_filename      // execute the commands in the file
//      cd directory                // change directories
//      reset                       // reset CStereoParameters to defaults
//      exit                        // exit this script file
//  If depth_map is specified, the stereo matching algorithm is executed.
//  If output_params is specified, an evaluation is performed and
//    the results are written out.
//
//  The parameters in the CStereoIOParameters class below are "use once"
//  parameters, i.e., they are used in the current command invocation,
//  and then erased (not re-used).  However, other parameters, such as
//  those in the CStereoParameters class, are remembered between commands.
//
//  See the descriptions below for more details on other parameters.
//
// SEE ALSO
//  StereoIO.cpp            implementation of this class
//  StereoParameters.h      description of parameters controlling the matcher
//  StereoMatcher.h         description of stereo matcher class
//  ParameterIO.h           description of the parameter input/output functionality
//
// Copyright © Richard Szeliski and Daniel Scharstein, 2001.
// See Copyright.h for more details
//
///////////////////////////////////////////////////////////////////////////

#include "StereoMatcher.h"

#ifndef FILENAME_MAX
#define FILENAME_MAX 1024
#endif

struct CStereoFrameDescriptor
{
    // An input image and its associated meta-data (depth map, ground truth, ...):

    char input_file[FILENAME_MAX];      // input image (gray or color)
    char depth_file[FILENAME_MAX];      // depth image (input/output)
    char truth_file[FILENAME_MAX];      // ground truth depth image
    char resampled_file[FILENAME_MAX];  // resampled image (for reprojection error)

    CStereoFrameDescriptor();
        // Constructor
    void PIOInitialize(CParameterIO& pio);
        // Initialize parameter input/output object
    void ReadFrame(CStereoFrame& fr, int verbose);
        // Read the constituent images
    void WriteFrame(CStereoFrame& fr, int verbose, bool write_depth);
        // Write the depth and resampled images
};

struct CStereoIOParameters
{
    // Additional parameters for running stereo from command-line or script:

    char input_params[FILENAME_MAX];    // input parameter file
    char output_params[FILENAME_MAX];   // output parameter and result file
    char input_data[FILENAME_MAX];      // original image and depth map filenames
    char output_data[FILENAME_MAX];     // final image and depth map filenames
    char depth_map[FILENAME_MAX];       // depth map image (override input_data)
    char cost_map[FILENAME_MAX];        // cost map (DSI), for debugging

    CStereoIOParameters();
        // Constructor
    void PIOInitialize(CParameterIO& pio);
        // Initialize parameter input/output object
};

class CStereoIO : public
    CStereoMatcher,
    CStereoIOParameters
{
public:
    CStereoIO();
        // Constructor
    int InterpretCommandLine(int argc, const char *argv[]);
        // Interpret the command (or script file) line

protected:
    void InterpretScriptFile(const char* script_file);
        // Interpret the commands in the script file
    void RunMatcher();
        // Execute the command
    void ReadParameters();
        // Read in the algorithm parameters
    void WriteParameters();
        // Write out the algorithm parameters and accuracy measures
    void ReadData();
        // Read in the images and (optional) ground truth
    void WriteData();
        // Write out the resulting disparity maps

private:
    CParameterIO m_pio;   // parameter input/output object
    CParameterIO m_dio;   // data file input/output object
    CStereoFrameDescriptor m_fd0;    // reference input/output frame descriptor
    std::vector<CStereoFrameDescriptor> m_fd;   // list of frame descriptors
};
