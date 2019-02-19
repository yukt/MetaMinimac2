#ifndef METAM_MARKOVMODEL_H
#define METAM_MARKOVMODEL_H

//#include "MarkovParameters.h"
#include "MetaMinimac.h"

using namespace std;

class MarkovModel
{
public:
    int NoVariants;
    int NoInPrefix;

    vector< double > InitProb;
    vector< vector<double> > LeftProb;
    vector< double > PrevRightProb;
    vector< double > CurrentRightProb;

    double Recom, backgroundError;
    double JumpFix, JumpThreshold;

    MarkovModel()
    {
        Recom = 1e-5;
        backgroundError = 1e-5;
        JumpThreshold = 1e-10;
        JumpFix = 1e10;

    };

    void initialize(int Sample, MetaMinimac *const ThisStudy);
    void walkLeft(int Sample, MetaMinimac *const ThisStudy);
    void walkRight(int Sample, MetaMinimac *const ThisStudy, int SampleInBatch);
};

class LogOddsModel
{
private:

    int NoMarkers;
    int NoStudies;
    int SampleID;
    vector<vector<double> > LooDosageVal;
    vector<double> ChipGTVal;


public:
    void initialize(int SampleId, MetaMinimac *const ThisStudy);
    double  operator()(vector<double> x);
};






#endif //METAM_MARKOVMODEL_H
