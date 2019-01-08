#ifndef METAM_MARKOVMODEL_H
#define METAM_MARKOVMODEL_H

#include "MarkovParameters.h"
#include "MetaMinimac.h"

using namespace std;

class MarkovModel : public MarkovParameters
{
public:
    vector<double> InitialProb;
    vector< vector<double> > TypedEmissionProb;

    vector< vector<double> > LeftProb;
    vector< vector<double> > RightProb;
    vector< vector<double> > FinalProb;

//    void initialize(int Sample, MetaMinimac *const ThisStudy);
};

//class LogOddsModel
//{
//private:
//
//    int NoMarkers;
//    int NoStudies;
//    int SampleID;
//    vector<vector<double> > LooDosageVal;
//    vector<double> ChipGTVal;
//
//
//public:
//    void initialize(int SampleId, MetaMinimac *const ThisStudy);
//    double  operator()(vector<double> x);
//};
//
//
//void logitTransform(vector<double> &From,
//                    vector<double> &To)
//{
//
//    double sum=1.0;
//    int NoDimensions = (int)To.size();
//    for(int i=0; i < (NoDimensions-1); i++) sum+=exp(From[i]);
//    for(int i=0; i < (NoDimensions-1); i++)  To[i]=exp(From[i])/sum;
//    To[NoDimensions-1]=1.0/sum;
//
//
//    double checkSum=0.0;
//    for(int i=0;i<To.size();i++)
//        checkSum+=To[i];
//    if(checkSum>1.0001)
//        abort();
//}




#endif //METAM_MARKOVMODEL_H
