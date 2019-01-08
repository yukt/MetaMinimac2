#include "MarkovModel.h"


//void MarkovModel::initialize(int Sample, MetaMinimac *const ThisStudy)
//{
//    LeftProb.resize(ThisStudy->NoVariants);
//    RightProb.resize(ThisStudy->NoVariants);
//    FinalProb.resize(ThisStudy->NoVariants);
//
//    for(int i=0; i<ThisStudy->NoVariants; i++)
//    {
//        LeftProb[i].resize(ThisStudy->NoInPrefix, 0.0);
//        RightProb[i].resize(ThisStudy->NoInPrefix, 0.0);
//        FinalProb[i].resize(ThisStudy->NoInPrefix, 0.0);
//    }
//}

//double LogOddsModel::operator()(vector<double> x)
//{
//
//
//    vector<double> tempVar(NoStudies);
//    double sum=0.0;
//
//
//    logitTransform(x,tempVar);
//
//    for(int ThisMarker=0;ThisMarker<NoMarkers;ThisMarker++)
//    {
//
//        double temp=0.0;
//
//        for(int j=0;j<NoStudies;j++)
//        {
//            temp+=((tempVar[j])*(LooDosageVal[j][ThisMarker]));
//        }
//
//        temp-=(ChipGTVal)[ThisMarker];
//        temp=temp*temp;
//        sum+=temp;
//
//    }
//    return sum;
//}
//
//void LogOddsModel::initialize(int SampleId, MetaMinimac *const ThisStudy)
//{
//    NoStudies=ThisStudy->NoInPrefix;
//    NoMarkers=400;
//    SampleID=SampleId;
//
//    LooDosageVal.resize(NoStudies);
//    for(int i=0; i<NoStudies; i++)
//    {
//        LooDosageVal[i].resize(NoMarkers);
//        for(int j=0; j<NoMarkers; j++) {
//            LooDosageVal[i][j] = ThisStudy->InputData[i].LooDosage[SampleID][j];
//        }
//    }
//
//    ChipGTVal.resize(NoMarkers);
//    for(int j=0; j<NoMarkers; j++)
//    {
//        ChipGTVal[j]=ThisStudy->InputData[0].TypedGT[SampleID][j];
//    }
//
//}
//
