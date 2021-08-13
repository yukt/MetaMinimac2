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
    String formatString;
    bool debug;
    int PrintBuffer, VcfBuffer;
    bool infoDetails;
    String formatStringForVCF;
    bool GT, DS, HDS, GP, SD;
    bool gzip, nobgzip;
    bool log;

    // check phasing consistency
    bool hapcheck;

    string CommandLine;

    UserVariables()
    {
        inputFiles = "";
        outfile = "MetaMinimac.Output";
        formatString = "GT,DS,HDS";
        debug=false;
        FileDelimiter=":";
        PrintBuffer = 100000000;
        infoDetails = true;
        formatStringForVCF = "";
        GT=false;
        DS=false;
        GP=false;
        HDS=false;
        SD=false;
        debug=false;
        gzip = true;
        nobgzip = false;
        VcfBuffer = 1000;
        log = false;
        hapcheck = false;
    };

    void Status()
    {
        cout << " Command Line Options: " << endl;
        printf( "      --input [%s],\n", inputFiles.c_str());
        printf( "      --output [%s],\n", outfile.c_str());
        printf( "      --format [%s],\n", formatString.c_str());
        printf( "      --phasingCheck %s,\n", hapcheck?"[ON]":"[OFF]");
        printf( "      --skipInfo %s,\n", infoDetails?"[OFF]":"[ON]");
        printf( "      --nobgzip  %s,\n", nobgzip?"[ON]":"[OFF]");
        printf( "      --weight   %s,\n", debug?"[ON]":"[OFF]");
        printf( "      --log      %s.", log?"[ON]":"[OFF]");
        printf("\n\n");
    }

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

    bool CheckValidity()
    {
        string formatPiece,formatTemp=formatString.c_str();
        char *end_str1;

        for(char * pch = strtok_r ((char*)formatTemp.c_str(),",", &end_str1);
            pch!=NULL;
            pch = strtok_r (NULL, ",", &end_str1))
        {

            formatPiece=(string)pch;

            if(formatPiece.compare("GT")==0)
            {
                GT=true;
            }
            else if(formatPiece.compare("DS")==0)
            {
                DS=true;
            }
            else if(formatPiece.compare("GP")==0)
            {
                GP=true;
            }
            else if(formatPiece.compare("HDS")==0)
            {
                HDS=true;
            }
            else if(formatPiece.compare("SD")==0)
            {
                SD=true;
            }
            else
            {
                cout << " ERROR !!! \n Cannot identify handle for -f [--format] parameter : "<<formatPiece<<endl;
                cout << " Available handles GT, DS, HDS and GP (for genotype, dosage, haplotype dosage and posterior probability). \n\n";
                cout<<" Program Exiting ..."<<endl<<endl;
                return false;
            }
        }

        bool colonIndex=false;
        if(GT)
        {
            formatStringForVCF+="GT";
            colonIndex=true;
        }
        if(DS)
        {
            formatStringForVCF+= (colonIndex?":DS":"DS");
            colonIndex=true;
        }
        if(HDS)
        {
            formatStringForVCF+= (colonIndex?":HDS":"HDS");
            colonIndex=true;
        }
        if(GP)
        {
            formatStringForVCF+= (colonIndex?":GP":"GP");
            colonIndex=true;
        }
        if(SD)
        {
            formatStringForVCF+= (colonIndex?":SD":"SD");
            colonIndex=true;
        }


        if(nobgzip)
            gzip=false;


        if (inputFiles == "")
        {
            cout<< " Missing -i [--input], a required parameter.\n\n";
            cout<< " Try -h [--help] for usage ...\n\n";
            cout<< " Program Exiting ...\n\n";
            return false;
        }

        if(PrintBuffer<=100)
        {
            cout << " ERROR !!! \n Invalid input for -b [--buffer] = "<<PrintBuffer<<"\n";;
            cout << " Buffer for writing output files should be at least 1,000 characters long !!! \n\n";
            cout<< " Try -h [--help] for usage ...\n\n";
            cout<<  " Program Exiting ..."<<endl<<endl;
            return false;
        }

        return true;
    };
};

#endif //METAM_MYVARIABLES_H
