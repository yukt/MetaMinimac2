#ifndef METAM_MYVARIABLES_H
#define METAM_MYVARIABLES_H

#include "StringBasics.h"

class UserVariables
{
public:
    String inputFiles;
    String outfile;
    String FileDelimiter;
    bool debug;

    UserVariables()
    {
        inputFiles = "";
        outfile = "MetaMinimac.Output";
        debug=false;
        FileDelimiter=":";
    };
};

#endif //METAM_MYVARIABLES_H
