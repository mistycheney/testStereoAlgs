///////////////////////////////////////////////////////////////////////////
//
// NAME
// ParameterIO.h : read or write parameter value pairs
//
// DESCRIPTION
//  Read or write program arguments and parameter values.
//
//  The arguments to the program are a list of parameter value pairs.
//  These can be given on the command line, or read from a file
//  or both (the command line overrides the file).
//
//  The program needs to know how large a value to write
//  on output.  Therefore, this class is currently restricted
//  to writing only the following types:
//      %s  strings
//      %f  floats
//      %d  ints
//
// SEE ALSO
//  ParameterIO.cpp       implementation of this class
//
// Copyright © Richard Szeliski, 2001.
// See Copyright.h for more details
//
///////////////////////////////////////////////////////////////////////////

#include <vector>

struct CParameterPair {
    char* name;         // name string
    char* format;       // format string
    void* value_loc;    // location of value
};

class CParameterIO {
    // Read/write parameter name/value pairs to/from file or command line
public:
    void PushParamPair(char *name, char *format, void *value_loc);
    void ReadFromStream(FILE *stream);
    void ReadFromFile(const char *filename);
    void ReadFromCommandLine(int argc, const char *argv[],
                             bool warn_bad_name = true);
    void WriteToStream(FILE *stream, bool single_line = false);
    void WriteToFile(const char *filename);
protected:
    void ReadParamPair(const char *name, const char *value,
                       bool warn_bad_name = true);
    std::vector<CParameterPair> m_pp_list;
};

#define PPPS(x)  PushParamPair(#x, "%s", &x)
#define PPPF(x)  PushParamPair(#x, "%f", &x)
#define PPPD(x)  PushParamPair(#x, "%d", &x)
#define PPPH(v)  PushParamPair(#v, "0x%08x", &v)
