#ifndef METAM_MARKOVPARAMETERS_H
#define METAM_MARKOVPARAMETERS_H

#include "MathVector.h"

class MarkovParameters
{
public:
    int        NoMarker;
    double     Recom,Error;

    MarkovParameters()
    {
        Recom = 1e-5;
        Error = 1e-5;
    };
};

#endif //METAM_MARKOVPARAMETERS_H
