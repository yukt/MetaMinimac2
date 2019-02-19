#include "MarkovModel.h"
#include "simplex.h"

using BT::Simplex;

void logitTransform(vector<double> &From,
                    vector<double> &To)
{

    double sum=1.0;
    int NoDimensions = (int)To.size();
    for(int i=0; i < (NoDimensions-1); i++) sum+=exp(From[i]);
    for(int i=0; i < (NoDimensions-1); i++)  To[i]=exp(From[i])/sum;
    To[NoDimensions-1]=1.0/sum;


    double checkSum=0.0;
    for(int i=0;i<To.size();i++)
        checkSum+=To[i];
    if(checkSum>1.0001)
        abort();
}


void MarkovModel::initialize(int Sample, MetaMinimac *const ThisStudy)
{
    NoVariants = ThisStudy->NoVariants;
    NoInPrefix = ThisStudy->NoInPrefix;

    LeftProb.clear();
//    RightProb.clear();
    LeftProb.resize(NoVariants);
//    RightProb.resize(NoVariants);

    for(int i=0; i<NoVariants; i++)
    {
        LeftProb[i].resize(NoInPrefix, 0.0);
//        RightProb[i].resize(NoInPrefix, 0.0);
    }

    LogOddsModel ThisSampleAnalysis;
    ThisSampleAnalysis.initialize(Sample, ThisStudy);
    vector<double> init(NoInPrefix-1, 0.0);
    vector<double> MiniMizer = Simplex(ThisSampleAnalysis, init);
    InitProb.resize(NoInPrefix, 0.0);
    logitTransform(MiniMizer, InitProb);

    for(int j=0; j<NoInPrefix; j++)
    {
        InitProb[j]+=backgroundError;
        LeftProb[0][j]=InitProb[j];
    }
}

void MarkovModel::walkLeft(int Sample, MetaMinimac *const ThisStudy)
{
    double CurrentGT = 0, CurrentLooDosage = 0;
    int NoCommonTypedProcessed = 0;
    double r = Recom*1.0/NoInPrefix, complement = 1-Recom;
    if(ThisStudy->CommonTyped[0])
    {
        CurrentGT = ThisStudy->InputData[0].TypedGT[Sample][0];
        for(int j=0; j<NoInPrefix; j++)
        {
            CurrentLooDosage = ThisStudy->InputData[j].LooDosage[NoCommonTypedProcessed][Sample];
            LeftProb[0][j]*=(CurrentGT==1)?CurrentLooDosage:(1-CurrentLooDosage);
        }
        NoCommonTypedProcessed++;
    }

    for(int i=1; i<NoVariants; i++)
    {
        double sum = 0;

        for (int j=0; j<NoInPrefix; j++)
        {
            for (int k=0; k<NoInPrefix; k++)
            {
                LeftProb[i][j] += LeftProb[i-1][k] * r;
            }
            LeftProb[i][j] += LeftProb[i-1][j] * complement;
            sum += LeftProb[i][j];
        }
        while (sum < JumpThreshold)
        {
            sum = 0;
            for (int j=0; j<NoInPrefix; j++)
            {
                LeftProb[i][j] *= JumpFix;
                sum += LeftProb[i][j];
            }
        }

        if (ThisStudy->CommonTyped[i])
        {
            sum = 0;
            CurrentGT = ThisStudy->InputData[0].TypedGT[Sample][NoCommonTypedProcessed];
            for(int j=0; j<NoInPrefix; j++)
            {
                CurrentLooDosage = ThisStudy->InputData[j].LooDosage[Sample][NoCommonTypedProcessed];
                LeftProb[i][j]*=(CurrentGT==1)?(CurrentLooDosage+backgroundError):(1-CurrentLooDosage+backgroundError);
                sum += LeftProb[i][j];
            }
            while (sum < JumpThreshold)
            {
                sum = 0;
                for (int j=0; j<NoInPrefix; j++)
                {
                    LeftProb[i][j] *= JumpFix;
                    sum += LeftProb[i][j];
                }
            }
            NoCommonTypedProcessed++;
        }

    }

    assert(NoCommonTypedProcessed==ThisStudy->NoCommonTypedVariants);
}


void MarkovModel::walkRight(int Sample, MetaMinimac *const ThisStudy, int SampleInBatch)
{
    vector< vector<double> > & FinalProb = ThisStudy->Posterior[SampleInBatch];
    FinalProb.resize(NoVariants);
    for(int i=0; i<NoVariants; i++)
    {
        FinalProb[i].resize(NoInPrefix, 0.0);
    }


    double CurrentGT = 0, CurrentLooDosage = 0;
    int NoCommonTypedProcessed = 0, NoCommonTypedVariants = ThisStudy->NoCommonTypedVariants;
    double r = Recom*1.0/NoInPrefix, complement = 1-Recom;

    PrevRightProb.clear();
    CurrentRightProb.clear();
    PrevRightProb.resize(NoInPrefix, 1.0);
    CurrentRightProb.resize(NoInPrefix, 0.0);
    for(int j=0; j<NoInPrefix; j++)
    {
        FinalProb[NoVariants-1][j] = LeftProb[NoVariants-1][j];
    }


    for(int i=NoVariants-1; i>0; i--)
    {
        if (ThisStudy->CommonTyped[i])
        {
            NoCommonTypedProcessed++;
            CurrentGT = ThisStudy->InputData[0].TypedGT[Sample][NoCommonTypedVariants-NoCommonTypedProcessed];
            for(int j=0; j<NoInPrefix; j++)
            {
                CurrentLooDosage = ThisStudy->InputData[j].LooDosage[Sample][NoCommonTypedVariants-NoCommonTypedProcessed];
                PrevRightProb[j]*=(CurrentGT==1)?(CurrentLooDosage+backgroundError):(1-CurrentLooDosage+backgroundError);
            }
        }

        double sum = 0.0;

        for(int j=0; j<NoInPrefix; j++)
        {
            for (int k=0; k<NoInPrefix; k++)
            {
                CurrentRightProb[j] += PrevRightProb[k] * r;
            }
            CurrentRightProb[j] += PrevRightProb[j] * complement;
            sum += CurrentRightProb[j];
        }

        if (sum < JumpThreshold)
        {
            for(int j=0; j<NoInPrefix; j++)
            {
                PrevRightProb[j] = JumpFix * CurrentRightProb[j];
                FinalProb[i-1][j] = PrevRightProb[j] * LeftProb[i-1][j];
                CurrentRightProb[j] = 0.0;
            }
        }
        else
        {
            for(int j=0; j<NoInPrefix; j++)
            {
                PrevRightProb[j] = CurrentRightProb[j];
                FinalProb[i-1][j] = PrevRightProb[j] * LeftProb[i-1][j];
                CurrentRightProb[j] = 0.0;
            }
        }
    }

    assert(NoCommonTypedProcessed==ThisStudy->NoCommonTypedVariants);
}


double LogOddsModel::operator()(vector<double> x)
{
    vector<double> tempVar(NoStudies);
    double sum=0.0;


    logitTransform(x,tempVar);

    for(int ThisMarker=0;ThisMarker<NoMarkers;ThisMarker++)
    {

        double temp=0.0;

        for(int j=0;j<NoStudies;j++)
        {
            temp+=((tempVar[j])*(LooDosageVal[j][ThisMarker]));
        }

        temp-=(ChipGTVal)[ThisMarker];
        temp=temp*temp;
        sum+=temp;

    }
    return sum;
}

void LogOddsModel::initialize(int SampleId, MetaMinimac *const ThisStudy)
{
    NoStudies=ThisStudy->NoInPrefix;
    NoMarkers=400;
    SampleID=SampleId;

    LooDosageVal.resize(NoStudies);
    for(int i=0; i<NoStudies; i++)
    {
        LooDosageVal[i].resize(NoMarkers);
        for(int j=0; j<NoMarkers; j++) {
            LooDosageVal[i][j] = ThisStudy->InputData[i].LooDosage[SampleID][j];
        }
    }

    ChipGTVal.resize(NoMarkers);
    for(int j=0; j<NoMarkers; j++)
    {
        ChipGTVal[j]=ThisStudy->InputData[0].TypedGT[SampleID][j];
    }

}

