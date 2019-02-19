#ifndef METAM_MARKOVPARAMETERS_H
#define METAM_MARKOVPARAMETERS_H

#include "MathVector.h"

class MarkovParameters
{
public:
    double     Recom,Error, backgroundError;

    MarkovParameters()
    {
        Recom = 1e-5;
        Error = 0.01;
        backgroundError = 1e-5;
    };
};

#endif //METAM_MARKOVPARAMETERS_H
