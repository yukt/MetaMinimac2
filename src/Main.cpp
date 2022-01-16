#include <getopt.h>
#include "MetaMinimac.h"
#include <unistd.h>

void MetaMinimacVersion();
void helpFile();

int main(int argc, char ** argv)
{
    MetaMinimac myAnalysis;
    bool help = false;

    int c;
    static struct option loptions[] =
            {
                    {"input",required_argument,NULL,'i'},
                    {"output",required_argument,NULL,'o'},
                    {"vcfBuffer",required_argument,NULL,'v'},
                    {"format",required_argument,NULL,'f'},
                    {"skipInfo",no_argument,NULL,'s'},
                    {"outputFormat",required_argument,NULL,'O'},
                    {"log",no_argument,NULL,'l'},
                    {"weight",no_argument,NULL,'w'},
                    {"skipPhasingCheck",no_argument,NULL,'p'},
                    {"help",no_argument,NULL,'h'},
                    {NULL,0,NULL,0}
            };

    while ((c = getopt_long(argc, argv, "i:o:O:v:f:slwh",loptions,NULL)) >= 0)
    {
        switch (c) {
            case 'i': myAnalysis.myUserVariables.inputPrefixes = optarg; break;
            case 'o': myAnalysis.myUserVariables.outPrefix = optarg; break;
            case 'w': myAnalysis.myUserVariables.debug=true; break;
            case 'f': myAnalysis.myUserVariables.formatString = optarg; break;
            case 's': myAnalysis.myUserVariables.infoDetails = false; break;
            case 'O': myAnalysis.myUserVariables.outFileExtension = optarg; break;
            case 'v': myAnalysis.myUserVariables.VcfBuffer=atoi(optarg); break;
            case 'p': myAnalysis.myUserVariables.hapcheck=false; break;
            case 'h': help=true; break;
            case 'l': myAnalysis.myUserVariables.log=true; break;
            case '?': helpFile(); return 1;
            default:  printf("[ERROR:] Unknown argument: %s\n", optarg);
        }
    }

    if (help)
    {
        helpFile();
        return(-1);
    }

    int remaining_arg_count = argc - optind;
    if (remaining_arg_count % 2)
    {
        std::cout << "Error: both dosage and empirical dosage files must be passed as input for each reference panel" << std::endl;
        return(-1);
    }

    if (remaining_arg_count)
        myAnalysis.myUserVariables.inputPrefixes.clear();

    if (std::strncmp(myAnalysis.myUserVariables.outFileExtension.c_str(), "vcf", 3) != 0 &&
        myAnalysis.myUserVariables.outFileExtension != "bcf" && myAnalysis.myUserVariables.outFileExtension != "ubcf" &&
        myAnalysis.myUserVariables.outFileExtension != "sav" && myAnalysis.myUserVariables.outFileExtension != "usav")
    {
        std::cout << "Error: unsupported --outputFormat" << std::endl;
        return(-1);
    }

    for (int i = 0; i < remaining_arg_count / 2; ++i)
        myAnalysis.myUserVariables.empInputFiles.emplace_back(argv[optind + i]);

    for (int i = remaining_arg_count / 2; i < remaining_arg_count; ++i)
        myAnalysis.myUserVariables.doseInputFiles.emplace_back(argv[optind + i]);

    int start_time = time(0);
    myAnalysis.myUserVariables.CreateCommandLine(argc,argv);

    FILE *LogFile=NULL;
    if(myAnalysis.myUserVariables.log)
        LogFile=freopen((myAnalysis.myUserVariables.outPrefix +".logfile").c_str(),"w",stdout);
    dup2(fileno(stdout), fileno(stderr));

    MetaMinimacVersion();
    myAnalysis.myUserVariables.Status();

    std::string compStatus;
    std::string MySuccessStatus="Error";

    MySuccessStatus = myAnalysis.Analyze();

    if(MySuccessStatus!="Success")
    {
        compStatus=MySuccessStatus;
        return(-1);
    }

    int time_tot = time(0) - start_time;

    cout<<"\n ------------------------------------------------------------------------------"<<endl;
    cout<<"                                END OF PROGRAM                                 "<<endl;
    cout<<" ------------------------------------------------------------------------------"<<endl;

    printf("\n Meta-Imputation Completed in %d hours, %d mins, %d seconds.\n",
           time_tot / 3600, (time_tot % 3600) / 60, time_tot % 60);
    myAnalysis.summary();

    cout<<"\n Thank You for using MetaMinimac2 !!! "<<endl<<endl;

    compStatus="Success";

    return 0;
}


void MetaMinimacVersion()
{
    printf("\n\n ------------------------------------------------------------------------------ \n");
    printf("             MetaMinimac2 -- An Efficient Tool for Meta-Imputation    \n");
    printf(" ------------------------------------------------------------------------------\n");
    printf(" (c) 2019 - Ketian Yu, Sayantan Das, Goncalo Abecasis \n");
    cout<< " Version : " << VERSION<< ";\n Built   : " << DATE << " by " << USER << endl;
    printf("\n URL = http://genome.sph.umich.edu/wiki/MetaMinimac2");
    printf("\n GIT = https://github.com/yukt/MetaMinimac2.git\n");
    cout << endl;
}

void helpFile()
{
    MetaMinimacVersion();
    printf( "\n About   : Combine GWAS data imputed against different panels  \n");
    printf( " Usage   : MetaMinimac2 [options] empirical_dose_files... dose_files...  \n");
    printf( "\n");
    printf( " Options :\n");
    //printf( "   -i, --input  <prefix1:prefix2 ...>  Colon-separated prefixes of input data to meta-impute.\n");
    printf( "   -o, --output <prefix>               Output prefix [MetaMinimac.Output] \n");
    printf( "   -O, --outputFormat <fmt>           Output file format (vcf, vcf.gz, bcf, ubcf, sav, or usav) [vcf.gz] \n");
    printf( "   -f, --format <string>               Comma-separated FORMAT tags [GT,DS,HDS]\n");
//    printf( "   -v, --vcfBuffer <int>               Maximum number of samples processed at a time [200] \n");
    printf( "   -p, --skipPhasingCheck              If OFF (by default and highly recommended), program will check phasing consistency before analysis.\n");
    printf( "   -s, --skipInfo                      If ON, the INFO fields are removed from the output file.\n");
    printf( "   -w, --weight                        If ON, weights will be saved in $prefix.metaWeights(.gz)\n");
    printf( "   -l, --log                           If ON, log will be written to $prefix.logfile. \n");
    printf( "   -h, --help                          If ON, detailed help on options and usage. \n");
    cout<<endl<<endl;
    return;
}
