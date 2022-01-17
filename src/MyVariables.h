#ifndef METAM_MYVARIABLES_H
#define METAM_MYVARIABLES_H

#include <string>
#include <cstring>
#include <vector>
#include <iostream>
#include <savvy/file.hpp>

using namespace std;

class UserVariables
{
public:
    std::string inputPrefixes;
    std::vector<std::string> empInputFiles;
    std::vector<std::string> doseInputFiles;
    std::string outPrefix;
    std::string FileDelimiter;
    std::string formatString;
    bool debug;
    int PrintBuffer, VcfBuffer;
    bool infoDetails;
    std::string formatStringForVCF;
    bool GT, DS, HDS, GP, SD;
    std::string outFileExtension;
    bool log;

    // check phasing consistency
    bool hapcheck;

    string CommandLine;

    UserVariables()
    {
        inputPrefixes = "";
        outPrefix = "MetaMinimac.Output";
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
        outFileExtension = "vcf.gz";
        VcfBuffer = 1000;
        log = false;
        hapcheck = true;
    };

    void Status()
    {
        cout << " Command Line Options: " << endl;
        printf( "      --input [%s],\n", inputPrefixes.c_str());
        printf( "      --output [%s],\n", outPrefix.c_str());
        printf( "      --outputFormat  %s,\n", outFileExtension.c_str());
        printf( "      --format [%s],\n", formatString.c_str());
        printf( "      --skipPhasingCheck %s,\n", hapcheck?"[OFF]":"[ON]");
        printf( "      --skipInfo %s,\n", infoDetails?"[OFF]":"[ON]");
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
        std::strcpy(MyCommandLine,argv[0]);

        for (int i=1; i<argc; i++)
        {
            std::strcat(MyCommandLine, " ");
            std::strcat(MyCommandLine, argv[i]);
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

        if (inputPrefixes == "" && empInputFiles.empty())
        {
            cout<< " Missing -i [--input], a required parameter.\n\n";
            cout<< " Try -h [--help] for usage ...\n\n";
            cout<< " Program Exiting ...\n\n";
            return false;
        }

        if (!empInputFiles.empty() && empInputFiles.size() != doseInputFiles.size())
        {
            cout<< " Missing input files.\n\n";
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
    }

    int outCompressionLevel() const
    {
      if (outFileExtension == "vcf.gz" || outFileExtension == "bcf" || outFileExtension == "sav")
          return 6;
      return 0;
    }

    savvy::file::format outFileFormat() const
    {
      if (outFileExtension == "bcf" || outFileExtension == "ubcf")
          return savvy::file::format::bcf;
      else if (outFileExtension == "sav" || outFileExtension == "usav")
          return savvy::file::format::sav;
      return savvy::file::format::vcf;
    }
};

#endif //METAM_MYVARIABLES_H
