#ifndef METAM_MYVARIABLES_H
#define METAM_MYVARIABLES_H

#include "StringBasics.h"

using namespace std;

class UserVariables
{
public:
    String inputFiles;
    String outfile;
    String FileDelimiter;
    bool debug;
    int PrintBuffer;
    bool infoDetails;
    String formatStringForVCF;
    bool GT, DS, HDS, GP, SD;

    string CommandLine;

    UserVariables()
    {
        inputFiles = "";
        outfile = "MetaMinimac.Output";
        debug=false;
        FileDelimiter=":";
        PrintBuffer = 100000000;
        infoDetails = true;
        formatStringForVCF = "GT:DS:HDS:GP";
        GT = true;
        DS = true;
        HDS = true;
        GP = true;
    };

    void CreateCommandLine(int argc, char ** argv)
    {
        int len = 0;

        for (int i=0; i<argc; i++)
            len += strlen(argv[i]) + 1;



        char MyCommandLine[len];
        strcpy(MyCommandLine,argv[0]);

        for (int i=1; i<argc; i++)
        {
            strcat(MyCommandLine, " ");
            strcat(MyCommandLine, argv[i]);
        }
        CommandLine=MyCommandLine;
    }
};

#endif //METAM_MYVARIABLES_H
