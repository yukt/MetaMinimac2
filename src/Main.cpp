#include <getopt.h>
#include "Parameters.h"
#include "MetaMinimac.h"

void MetaMinimacVersion();
void helpFile();

int main(int argc, char ** argv)
{
    MetaMinimac myAnalysis;
    ParameterList inputParameters;
    bool log = false, help = false;

    int c;
    static struct option loptions[] =
            {
                    {"input",required_argument,NULL,'i'},
                    {"output",required_argument,NULL,'o'},
                    {"log",no_argument,NULL,'l'},
                    {"debug",no_argument,NULL,'d'},
                    {"help",no_argument,NULL,'h'},
                    {NULL,0,NULL,0}
            };

    while ((c = getopt_long(argc, argv, "i:o:hld",loptions,NULL)) >= 0)
    {
        switch (c) {
            case 'i': myAnalysis.myUserVariables.inputFiles = optarg; break;
            case 'o': myAnalysis.myUserVariables.outfile = optarg; break;
            case 'd': myAnalysis.myUserVariables.debug=true; break;
            case 'h': help=true; break;
            case 'l': log=true; break;
            case '?': helpFile(); return 1;
            default:  printf("[ERROR:] Unknown argument: %s\n", optarg);
        }
    }

    int start_time = time(0);
    myAnalysis.myUserVariables.CreateCommandLine(argc,argv);

    FILE *LogFile=NULL;
    if(log)
        LogFile=freopen(myAnalysis.myUserVariables.outfile +".logfile","w",stdout);
    dup2(fileno(stdout), fileno(stderr));

    MetaMinimacVersion();
    if (help)
    {
        helpFile();
        return(-1);
    }

    myAnalysis.Analyze();

    int time_tot = time(0) - start_time;

    cout<<"\n ------------------------------------------------------------------------------"<<endl;
    cout<<"                                END OF PROGRAM                                 "<<endl;
    cout<<" ------------------------------------------------------------------------------"<<endl;


    cout << "\n Program Successfully Implemented... \n ";


    printf("\n Total Run completed in %d hours, %d mins, %d seconds.\n",
           time_tot / 3600, (time_tot % 3600) / 60, time_tot % 60);

    cout<<"\n Thank You for using MetaMinimac !!! "<<endl<<endl;

    return 0;
}


void MetaMinimacVersion()
{
    printf("\n\n -------------------------------------------------- \n");
    printf("                   MetaMinimac     \n");
    printf(" --------------------------------------------------\n");
    printf(" (c) 2019 - Ketian Yu, Sayantan Das, Goncalo Abecasis \n");
    cout<< " Version : " << VERSION<< ";\n Built   : " << DATE << " by " << USER << endl;
    cout << endl;
}

void helpFile()
{
    printf( "\n About   : Combine GWAS data imputed against different panels  \n");
    printf( " Usage   : metaMinimac [options] \n");
    printf( "\n");
    printf( " Options :\n");
    printf( "   -i, --input  <prefix1 prefix2 ...>  Prefixes of input data to meta-impute\n");
    printf( "   -o, --output <prefix>               Output prefix [MetaMinimac.Output] \n");
    printf( "   -f, --format <string>               Comma separated FORMAT tags [GT,DS,HDS]\n");
    printf( "   -s, --skipInfo                      Skip INFO in output [FALSE] \n");
    printf( "   -n, --nobgzip                       Output unzipped file [FALSE]\n");
    printf( "   -b, --buffer                        Print Buffer [1e8] \n");
    printf( "   -d, --debug                         Debug mode [FALSE] \n");
    printf( "   -l, --log                           Save logfile [TRUE] \n");
//    printf("\n URL = http://genome.sph.umich.edu/wiki/MetaMinimac\n");
//    printf(" GIT = https://github.com/Santy-8128/MetaMinimac\n");
//    printf("\n Visit website for more details ...\n");

    cout<<endl<<endl;
    return;
}
