///////////////////////////////////////////////////////////////////////////
//
// NAME
//  Verbose.h -- defines verbosity levels for the stereo matcher
//
// DESCRIPTION
//  Defines the verbosity levels used to print varying amounts of information.
//  The higher the level, the more information is printed
//  Messages are printed to stderr.
//
// Copyright © Richard Szeliski and Daniel Scharstein, 2001.
// See Copyright.h for more details
//
///////////////////////////////////////////////////////////////////////////

enum EVerbosityLevel
{
    eVerboseQuiet           = 0,   // no output except for errors
    eVerboseWarning         = 1,   // include warning messages
    eVerboseSummary         = 2,   // print brief summary of algorithm and result 
    eVerboseProgress        = 3,   // report progress through steps of algorithm
    eVerboseFileIO          = 4,   // report reading and writing of files
    eVerboseTiming          = 5,   // print timing information
    eVerbosePredictionError = 6,   // print prediction error for all frames
    eVerboseScriptFile      = 10,  // echo commands from script file
    eVerboseInnerLoops      = 15,  // show inner iterations
    eVerboseDumpFiles       = 25,  // dump intermediate results as image files
    eVerboseAllMessages     = 99   // include all information
};
