///////////////////////////////////////////////////////////////////////////
//
// NAME
//  StereoIO.cpp -- stereo matcher file and parameter input/output
//
// DESCRIPTION
//  The CStereoIO class reads and write image and parameter files
//  needed to drive the CStereoMatcher class.
//
// SEE ALSO
//  StereoIO.h              longer description this class
//
// Copyright Richard Szeliski and Daniel Scharstein, 2001.
// See Copyright.h for more details
//
// *** fixed bug with "cd" command that causes trouble under Unix
//
// 7/20/07 - replaced "\t\n" with "\t\r\n" in the command line parser
//           so that windows script files work under Unix.  Thanks
//           to Mariano Tepper for providing this fix!
//
// 1/21/2010 - fixed bug involving string length in remove_stray_chars.
//	       Thanks to Martin Gernhard for pointing this out!
//
///////////////////////////////////////////////////////////////////////////

#include "Error.h"
#include "StereoIO.h"
#include "ImageIO.h"
#include "Convert.h"

#ifdef WIN32
#include <direct.h>
#else
#include <unistd.h>
#endif

// struct CStereoFrameDescriptor:
//  An input image and its associated meta-data (depth map, ground truth, ...)

CStereoFrameDescriptor::CStereoFrameDescriptor()
{
    memset(this, 0, sizeof(CStereoFrameDescriptor));
}

void CStereoFrameDescriptor::PIOInitialize(CParameterIO& o)
{
    // Initialize parameter input/output object

    // An input image and its associated meta-data (depth map, ground truth, ...):
    o.PPPS(input_file);         // input image (gray or color)
    o.PPPS(depth_file);         // depth image (input/output)
    o.PPPS(truth_file);         // ground truth depth image
    o.PPPS(resampled_file);     // resampled image (for reprojection error)
}

static void ReadIfThere(CByteImage& img, const char* filename,
                        bool ok_if_not_there, int verbose)
{
    if (filename[0] == 0)
        return;
        
    if (verbose >= eVerboseFileIO)
        fprintf(stderr, "reading image %s\n", filename);

    if (! ok_if_not_there)
        ReadImage(img, filename);
    else try
    {
        ReadImage(img, filename);
    }
    catch (CError& err)
    {
      //    err.message;    // Don't do anything, since image may not exist
    }
    // Most of the stereo matching code assumes nbands==1 or nbands==4
    // make sure this is true
    int nb = img.Shape().nBands;
    if (! (nb == 1 || nb == 4))
        throw CError("ReadImage(%s): number of bands (%d) not supported by stereo matcher",
        filename, nb);
    
}

void CStereoFrameDescriptor::ReadFrame(CStereoFrame& fr, int verbose)
{
    // Read the constituent images
    ReadIfThere(fr.input_image, input_file, false, verbose);
    ReadIfThere(fr.depth_image, depth_file, true, verbose); 
    ReadIfThere(fr.truth_image, truth_file, false, verbose);
    // Don't read the resampled image (will be recomputed, if needed)
 // ReadIfThere(fr.resampled_image, resampled_file, true, verbose);
}

static void WriteIfThere(CByteImage& img, const char* filename, int verbose)
{
    if (filename[0] == 0 || img.Shape().height == 0)
        return;

    if (verbose >= eVerboseFileIO)
        fprintf(stderr, "writing image %s\n", filename);

    WriteImage(img, filename);
}

void CStereoFrameDescriptor::WriteFrame(CStereoFrame& fr, int verbose, bool write_depth)
{
    // Write the depth and resampled images
    if (write_depth)
        WriteIfThere(fr.depth_image, depth_file, verbose);
    WriteIfThere(fr.resampled_image, resampled_file, verbose);
}


// struct CStereoIOParameters
//  Additional parameters for running stereo from command-line or script:

CStereoIOParameters::CStereoIOParameters()
{
    memset(this, 0, sizeof(CStereoIOParameters));
}

void CStereoIOParameters::PIOInitialize(CParameterIO& o)
{
    // Initialize parameter input/output object

    // Additional parameters for running stereo from command-line or script:
    o.PPPS(input_params);       // input parameter file
    o.PPPS(output_params);      // output parameter and result file
    o.PPPS(input_data);         // original image and depth map filenames
    o.PPPS(output_data);        // final image and depth map filenames
    o.PPPS(depth_map);          // depth map image (override input_data)
    o.PPPS(cost_map);           // cost map (DSI), for debugging
}

CStereoIO::CStereoIO()
{
    // Constructor
    CStereoParameters::PIOInitialize(m_pio);
    CStereoIOParameters::PIOInitialize(m_pio);

    m_fd0.PIOInitialize(m_dio);
}

// *** new function to strip carriage return characters from filenames
// *** and directory names that cause trouble under Unix
// *** (TODO: needs to be cleaned up)
char * remove_stray_chars(const char * input) {
    char * my_argument;
    int length = strlen(input);
    my_argument = (char *)malloc((length+1)*sizeof(char));

    int i=0;

    while(input[i]==' ')
	i++;

    int t=i;
    while(input[i]!='\r' && input[i]!=' ' && input[i]!='\0'){
	my_argument[i-t]=input[i];
	i++;
    }
    my_argument[i-t]='\0';

    return my_argument;
}

int CStereoIO::InterpretCommandLine(int argc, const char *argv[])
{
    // Interpret the command (or script file) line
    const char *command  = argv[1];
    const char *argument = argv[2];

    // Zero out the parameters in this class before initial read
    memset(&this->input_params, 0, sizeof(CStereoIOParameters));

    // Read in the arguments (may get clobbered by input_params, but fix that later)
    m_pio.ReadFromCommandLine(argc-1, argv+1);
    if (command == 0 || command[0] == 0)
        throw CError("No command is given on command line.\n\
  Please see the README-StereoMatch.txt file for proper usage.\n");

    // Interpret the script file
    if (strcmp("script", command) == 0)
    {
      
        char * my_argument=remove_stray_chars(argument);
        InterpretScriptFile(my_argument);
        return 0;                           // no other commands executed
    }

    // Change the directory
    if (strcmp("cd", command) == 0)
    {
      //int retval = chdir(argument);
        if (verbose >= eVerboseSummary)
        {
	  char buf[1024];
	        getcwd(buf, 1024);
	        //int length=strlen(buf);
	        //fprintf(stderr, "buf length: %d\n", length);
	       	fprintf(stderr, "working dir: %s\n", buf);
	  
        }
        char * my_argument=remove_stray_chars(argument);
	int retval = chdir(my_argument);


    if (retval)
	  throw CError("Could not cd to %s", my_argument);
    //else 
    //  fprintf(stderr, "cd succeed\n");
    return 0;                           // no other commands executed
    }

    // Reset all CStereoParameters that affect the algorithm to their default values
    if (strcmp("reset", command) == 0)
    {
        ReInitializeAlgParams();
    }

    // Exit this script
    if (strcmp("exit", command) == 0)
    {
        return 1;
    }

    // Read in the parameter file if specified
    ReadParameters();

    // Re-read the command-line arguments (to override parameter file)
    m_pio.ReadFromCommandLine(argc-1, argv+1, false);

    // Reset the output parameters
    ResetOutputParams();

    // Read in the data (input images) if specified
    ReadData();

    // Copy the depth_map name into the frame descriptor
    if (frame_ref >= 0 && frame_ref < (int)m_fd.size() && depth_map[0])
        strcpy(m_fd[frame_ref].depth_file, depth_map);

    // Early out if no matching or evaluation
    if (depth_map[0] == 0 && output_params[0] == 0)
        return 0;

    // At this point, a "depth_map" and/or an "output_params" command is present
    // in the current command line
    if (frame_ref   >= (int)m_fd.size() ||
        frame_match >= (int)m_fd.size())
        throw CError("Data has not yet been read in");

    // Read in the depth_map if only evaluating
    if (evaluate_only) {
        ReadIfThere(m_frame[frame_ref].depth_image, depth_map, false, verbose); 
    }

    // If output_params is not specified, derive name from depth_map
    if (output_params[0] == 0)
    {
        strcpy(output_params, depth_map);
        strcpy(&output_params[strlen(output_params)-3], "txt");
    }

    // Execute the command
    RunMatcher();

    // Write out the data (results)
    WriteData();

    // Write out the (updated) parameter file if specified
    WriteParameters();

    if (verbose >= eVerboseSummary)
        fprintf(stderr, "  writing %s\n", output_params[0] ? output_params : depth_map);

    return 0;
}

void CStereoIO::RunMatcher()
{
    // Run the matcher and evaluator

    // TODO: ALWAYS need to call ComputeCorrespondence because that's where
    // m_float_disparity is initialized from depth_image??
    if (depth_map[0])
        ComputeCorrespondence();

    if (output_params[0])
        Evaluate();
}

class CCommandLineParser
{
public:
    CCommandLineParser(char *command_line);
    int argc;
    const char *argv[256];
};

CCommandLineParser::CCommandLineParser(char *command_line)
{
    // Set up the pointers to the parsed line (modifies the command_line string)
    argc = 1;
    argv[0] = 0;    // no valid argv[0] if parsing from string
    int spn = strspn(command_line, " \t\r\n");  // skip whitespace
    command_line += spn;
    while (command_line[0] && argc < 256)
    {
        if (command_line[0] == '#') // comment line or trailing comment
            return;

        argv[argc] = command_line;  // push the next command line argument
        argc += 1;
        char *last = strpbrk(command_line, " \t\r\n");   // find next whitespace
        if (last)
        {
            last[0] = 0;    // null-terminate the argument (destructive)
            command_line = last+1;
            int spn = strspn(command_line, " \t\r\n");  // skip whitespace
            command_line += spn;
        }
        else
            break;
    }
}

void CStereoIO::InterpretScriptFile(const char* script_file)
{
    // Interpret the commands in the script file
    FILE *stream = fopen(script_file, "r");
    if (stream == 0)
        throw CError("InterpretScriptFile: could not open %s", script_file);

    // Echo where we are
    if (verbose >= eVerboseScriptFile)
        fprintf(stderr, "Interpreting script file %s\n", script_file);

    // Iterate over the lines
    const int MAX_LINE = 4096;  // maximum command line length
    char command_line[MAX_LINE];
    while (fgets(command_line, MAX_LINE, stream))
    {
        // Echo if verbose
        if (verbose >= eVerboseScriptFile)
            fprintf(stderr, command_line);

        CCommandLineParser clp(command_line);
        if (clp.argc > 1)
        {
            // Re-use current object, so remember previous state/parameters
            int retval = InterpretCommandLine(clp.argc, clp.argv);
            if (retval)
                break;
        }
    }
    fclose(stream);
}

void CStereoIO::ReadParameters()
{
    // Read in the algorithm parameters
    if (input_params[0])
    {
        m_pio.ReadFromFile(input_params);
    }
}

void CStereoIO::WriteParameters()
{
    // Write out the algorithm parameters and accuracy measures
    if (output_params[0])
    {
        m_pio.WriteToFile(output_params);
    }
}

void CStereoIO::ReadData()
{
    // Read in the images and (optional) ground truth
    if (input_data[0] == 0)
        return;

    FILE *stream = fopen(input_data, "r");
    if (stream == 0)
        throw CError("ReadData: could not open %s", input_data);

    // Iterate over the frames
    const int MAX_LINE = 4096;  // maximum descriptor line length
    char line[MAX_LINE];
    m_fd.clear();
    m_frame.clear();
    while (fgets(line, MAX_LINE, stream))
    {
        CCommandLineParser clp(line);
        if (clp.argc > 1)
        {
            // Read the descriptor
            memset(&m_fd0, 0, sizeof(m_fd0));
            m_dio.ReadFromCommandLine(clp.argc-1, clp.argv+1);

            // Override the depth_map?
            if ((int)m_fd.size() == frame_ref && depth_map[0])
                strcpy(m_fd0.depth_file, depth_map);
            m_fd.push_back(m_fd0);

            // Read the files
            CStereoFrame fr;
            m_fd0.ReadFrame(fr, verbose);
            m_frame.push_back(fr);
        }
    }
    fclose(stream);
}

void CStereoIO::WriteData()
{
    // Write out the resulting disparity maps (unless evaluate_only)
    //  reprojected images (if specified) and
    //  the descriptor file (if specified)
    FILE *stream = 0;
    if (output_data[0])
    {
        stream = fopen(output_data, "w");
        if (stream == 0)
            throw CError("WriteData: could not open %s", output_data);
    }

    // Iterate over the frames
    for (int i = 0, n = m_fd.size(); i < n; i++)
    {
        m_fd0 = m_fd[i];
        CStereoFrame& fr = m_frame[i];
        m_fd0.WriteFrame(fr, verbose, ! evaluate_only);
        if (stream)
            m_dio.WriteToStream(stream, true);
    }
    if (stream)
        fclose(stream);

    // Write out the error images
    if (eval_error_scale > 0.0f)
    {
        // Create the new filenames (suffix _e and _b)
        char filename[FILENAME_MAX];
        char *dot = strrchr(depth_map, '.');
        int n_stem = dot - depth_map;
        strncpy(filename, depth_map, n_stem);
        strcpy(&filename[n_stem], "_e");
        strcat(filename, dot);

        WriteImage(m_disparity_error, filename);
        filename[n_stem+1] = 'b';
        WriteImage(m_bad_pixels, filename);
        if (m_disparity_histogram.Shape().width > 0)
        {
            filename[n_stem+1] = 'h';
            WriteImage(m_disparity_histogram, filename);
        }
    }

    // Write out combined disparity and status color image
    if (m_status_disp.Shape().width > 0)
    {
        // Create the new filename (.ppm instead of .pgm)
        char filename[FILENAME_MAX];
        strcpy(filename, depth_map);
        filename[strlen(filename)-2] = 'p';
        WriteImage(m_status_disp, filename);
    }

}
