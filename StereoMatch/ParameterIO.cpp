///////////////////////////////////////////////////////////////////////////
//
// NAME
// ParameterIO.cpp : read or write parameter value pairs
//
// DESCRIPTION
//  Read or write program arguments and parameter values.
//
// SEE ALSO
//  ParameterIO.h       class definition
//
// Copyright © Richard Szeliski, 2001.  See Copyright.h for more details
//
///////////////////////////////////////////////////////////////////////////

#include <vector>
#include <stdlib.h>
#include "Error.h"
#include "ParameterIO.h"

void CParameterIO::PushParamPair(char *name, char *format, void *value_loc)
{
    CParameterPair pp;
    pp.name = name;
    pp.format = format;
    pp.value_loc = value_loc;
    m_pp_list.push_back(pp);
}

void CParameterIO::ReadFromStream(FILE *stream)
{
    char name[1024], value[1024];
    while (fscanf(stream, "%s %s", name, value) == 2)
        ReadParamPair(name, value);
}

void CParameterIO::ReadFromFile(const char *filename)
{
    FILE *stream = fopen(filename, "r");
    if (stream == 0)
        throw CError("CParameterIO::ReadFromFile: could not open %s", filename);
    ReadFromStream(stream);
    fclose(stream);
}

void CParameterIO::ReadFromCommandLine(int argc, const char *argv[],
                                       bool warn_bad_name)
{
    for (int i = 0; i+1 < argc; i += 2)
        ReadParamPair(argv[i], argv[i+1], warn_bad_name);
}

void CParameterIO::WriteToStream(FILE *stream, bool single_line)
{
    for (int i = 0, n = m_pp_list.size(); i < n; i++) {
        CParameterPair& pp = m_pp_list[i];
        char* s = (char *) pp.value_loc;
        char fmt_type = pp.format[strlen(pp.format)-1];
        bool is_string =  (fmt_type == 's');
        if (is_string && s[0] == 0)
            continue;
        fprintf(stream, "%s ", pp.name);
        if (is_string) 
            fprintf(stream, pp.format, s);
        else  if (fmt_type == 'f')
            fprintf(stream, pp.format, *(float *) pp.value_loc);
        else  if (fmt_type == 'd' || fmt_type == 'x')
            fprintf(stream, pp.format, *(int *) pp.value_loc);
        else
            throw CError("CParamaterIO: illegal format '%s'\n", pp.format);
        fprintf(stream, (single_line) ? " " : "\n");
    }
    if (single_line)
        fprintf(stream, "\n");  // final newline
}

void CParameterIO::WriteToFile(const char *filename)
{
    bool append_mode = 0;

    // if first character of filename is '+', append
    if (filename[0] == '+') {
        filename++; // skip over the '+'
        append_mode = 1;
    }

    FILE *stream = fopen(filename, append_mode? "a" : "w");
    if (stream == 0)
        throw CError("CParameterIO::WriteToFile: could not open %s", filename);

    WriteToStream(stream);

    if (append_mode)
        fprintf(stream, "\n"); // trailing newline to separate different runs

    fclose(stream);
}

void CParameterIO::ReadParamPair(const char *name, const char *value,
                                 bool warn_bad_name)
{
    for (int i = 0; i < (int) m_pp_list.size(); i++) {
        CParameterPair& pp = m_pp_list[i];
        if (strcmp(pp.name, name) == 0) {
            sscanf(value, pp.format, pp.value_loc);
            return;
        }
    }

    // these names are ok:
    if (strcmp(name, "script") == 0)
        return;
    if (strcmp(name, "cd") == 0)
        return;

    // If we get here, we have an unknown name value
    if (warn_bad_name)
        fprintf(stderr, "Warning: CParameterIO::ReadParamPair: unknown parameter %s\n", name);
}
