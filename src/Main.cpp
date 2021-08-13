#include <getopt.h>
#include "Parameters.h"
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
                    {"nobgzip",no_argument,NULL,'n'},
                    {"log",no_argument,NULL,'l'},
                    {"weight",no_argument,NULL,'w'},
                    {"hapCheck",no_argument,NULL,'c'},
                    {"help",no_argument,NULL,'h'},
                    {NULL,0,NULL,0}
            };

    while ((c = getopt_long(argc, argv, "i:o:v:f:snlwh",loptions,NULL)) >= 0)
    {
        switch (c) {
            case 'i': myAnalysis.myUserVariables.inputFiles = optarg; break;
            case 'o': myAnalysis.myUserVariables.outfile = optarg; break;
            case 'w': myAnalysis.myUserVariables.debug=true; break;
            case 'f': myAnalysis.myUserVariables.formatString = optarg; break;
            case 's': myAnalysis.myUserVariables.infoDetails = false; break;
            case 'n': myAnalysis.myUserVariables.nobgzip=true; break;
            case 'v': myAnalysis.myUserVariables.VcfBuffer=atoi(optarg); break;
            case 'c': myAnalysis.myUserVariables.hapcheck=true; break;
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

    int start_time = time(0);
    myAnalysis.myUserVariables.CreateCommandLine(argc,argv);

    FILE *LogFile=NULL;
    if(myAnalysis.myUserVariables.log)
        LogFile=freopen(myAnalysis.myUserVariables.outfile +".logfile","w",stdout);
    dup2(fileno(stdout), fileno(stderr));

    MetaMinimacVersion();
    myAnalysis.myUserVariables.Status();

    String compStatus;
    String MySuccessStatus="Error";

    MySuccessStatus = myAnalysis.Analyze();

    if(MySuccessStatus!="Success")
    {
        compStatus=MySuccessStatus;
        PhoneHome::completionStatus(compStatus.c_str());
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
    PhoneHome::completionStatus(compStatus.c_str());

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
    printf( " Usage   : MetaMinimac2 [options] \n");
    printf( "\n");
    printf( " Options :\n");
    printf( "   -i, --input  <prefix1:prefix2 ...>  Colon-separated prefixes of input data to meta-impute.\n");
    printf( "   -o, --output <prefix>               Output prefix [MetaMinimac.Output] \n");
    printf( "   -f, --format <string>               Comma-separated FORMAT tags [GT,DS,HDS]\n");
//    printf( "   -v, --vcfBuffer <int>               Maximum number of samples processed at a time [200] \n");
    printf( "   -s, --skipInfo                      If ON, the INFO fields are removed from the output file.\n");
    printf( "   -n, --nobgzip                       If ON, output files will NOT be bgzipped.\n");
    printf( "   -w, --weight                        If ON, weights will be saved in $prefix.metaWeights(.gz)\n");
    printf( "   -l, --log                           If ON, log will be written to $prefix.logfile. \n");
    printf( "   -h, --help                          If ON, detailed help on options and usage. \n");
    cout<<endl<<endl;
    return;
}
