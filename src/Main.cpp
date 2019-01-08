#include <getopt.h>
#include "Parameters.h"
#include "MetaMinimac.h"

int main(int argc, char ** argv)
{
    MetaMinimac myAnalysis;
    ParameterList inputParameters;
    bool log = true;

    int c;
    static struct option loptions[] =
            {
                    {"input",required_argument,NULL,'i'},
                    {"output",required_argument,NULL,'o'},
                    {"log",no_argument,NULL,'l'},
                    {"debug",no_argument,NULL,'d'},
                    {NULL,0,NULL,0}
            };

    while ((c = getopt_long(argc, argv, "i:o:f:c:snhld",loptions,NULL)) >= 0)
    {
        switch (c) {
            case 'i': myAnalysis.myUserVariables.inputFiles = optarg; break;
            case 'o': myAnalysis.myUserVariables.outfile = optarg; break;
            case 'd': myAnalysis.myUserVariables.debug=true; break;
            case 'l': log=true; break;
            default:  printf("[ERROR:] Unknown argument: %s\n", optarg);
        }
    }

    FILE *LogFile=NULL;
    if(log)
        LogFile=freopen(myAnalysis.myUserVariables.outfile +".logfile","w",stdout);
    dup2(fileno(stdout), fileno(stderr));

    myAnalysis.Analyze();

    return 0;
}

