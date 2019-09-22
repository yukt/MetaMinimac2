#ifndef METAM_METAMINIMAC_H
#define METAM_METAMINIMAC_H

#include "MyVariables.h"
#include "HaplotypeSet.h"

using namespace std;

void logitTransform(vector<double> &From,vector<double> &To);

class MetaMinimac
{
public:
    UserVariables myUserVariables;
    vector<String> InPrefixList;
    int NoInPrefix;
    string finChromosome;

    vector<variant> VariantList;
    vector<variant> CommonTypedVariantList;
    vector<string> CommonGenotypeVariantNameList;
    int NoHaplotypes, NoSamples;
    int NoVariants, NoCommonTypedVariants;

    // Variables for input dosage file stream and records
    vector<VcfFileReader*> InputDosageStream;
    vector<VcfRecord*> CurrentRecordFromStudy;
    vector<int> StudiesHasVariant;
    int CurrentFirstVariantBp;
    int NoStudiesHasVariant;

    vector<HaplotypeSet> InputData;

    // Process Part of Samples each time
    int StartSamId, EndSamId;
    double Recom, backgroundError;
    double JumpFix, JumpThreshold;
    vector<int> NoOffsetThisBlock;
    vector<vector<vector<double>>> Weights;
    vector<vector<vector<double>>> LeftProb;
    vector<vector<double>> PrevLeftProb, PrevRightProb;
    int NoVariantsProcessed, NoCommonVariantsProcessed;

    // Output files
    IFILE vcfdosepartial, vcfweightpartial;
    IFILE vcfsnppartial, vcfrsqpartial;
    IFILE metaWeight;
    char *VcfPrintStringPointer;
    char *WeightPrintStringPointer;
    char *RsqPrintStringPointer;
    char *SnpPrintStringPointer;
    int VcfPrintStringPointerLength, WeightPrintStringPointerLength, RsqPrintStringPointerLength, SnpPrintStringPointerLength;
    int batchNo;
    vector<IFILE> vcfrsqpartialList;

    variant* CurrentVariant;
    variant ThisVariant;
    string NextTypedName;
    int PrevBp, CurrBp, ThisBp;
    vector<vector<double>> *PrevWeights;
    vector<vector<double>> *CurrWeights;
    vector<float> CurrentMetaImputedDosage;
    vector<vector<float>> CurrentPosterior;

    float CurrentHapDosageSum, CurrentHapDosageSumSq;
    vector<float> HapDosageSum, HapDosageSumSq;

    // Buffer
    int NoRecords;
    int BufferBp, BufferNoVariants;
    vector<variant> BufferVariantList;


    MetaMinimac()
    {
        Recom = 1e-3;
        backgroundError = 1e-5;
        JumpThreshold = 1e-10;
        JumpFix = 1e10;
    };



    String Analyze();

    bool ParseInputVCFFiles();
    bool CheckSampleNameCompatibility();
    void OpenStreamInputDosageFiles(bool siteOnly);
    void CloseStreamInputDosageFiles();
    bool OpenStreamOutputDosageFiles();
    string GetDosageFileFullName(String prefix);
    bool doesExistFile(String filename);

    bool LoadEmpVariantInfo();
    void FindCommonGenotypedVariants();
    void FindCurrentMinimumPosition();
    int IsVariantEqual(VcfRecord &Rec1, VcfRecord &Rec2);
    void ReadCurrentVariantInfo();
    void UpdateCurrentRecords();

    void LoadLooDosage();

    String PerformFinalAnalysis();
    void CalculateWeights();
    void InitiateWeights();
    void CalculateLeftProbs();
    void CalculatePosterior();
    void InitiateLeftProb(int SampleInBatch);
    void InitiateRightProb(int SampleInBatch);
    void UpdateOneStepLeft(int SampleInBatch);
    void UpdateOneStepRight(int SampleInBatch);
    void MetaImputeAndOutput();
    void UpdateWeights();
    void OutputPartialVcf();
    void OutputAllVcf();
    void OutputMetaWeights();

//    void GetMetaEstimate(int Sample, int SampleInBatch);
//    void GetMetaEstimate(int SampleInBatch);
    void OpenTempOutputFiles();
    void FlushPartialVcf();
    void FlushAllVcf();

    void InitiateProbs();
    void InitLeftProb(int SampleInBatch);
    void WalkLeft();
    void WalkOneStepLeft();
    void WalkOneStepLeft(int SampleInBatch);
    void ProcessTypedLeftProb();
    void ProcessTypedLeftProb(int SampleInBatch);

    void AppendtoMainVcf();
    void AppendtoMainWeightsFile();

    void MetaImputeCurrentBuffer();
    void MetaImputeCurrentBuffer2();
    void MetaImputeCurrentBuffer3();
    void ClearCurrentBuffer();
    void ReadCurrentDosageData();
    void CreateMetaImputedData(int VariantId);
    void CreateMetaImputedData();
    void CalculateStats();
    void WalkOneStepRight();
    void WalkOneStepRight(int SampleInBatch);
    void UpdateLeftProb();
    void UpdateRightProb();
    void UpdateRightProb(int HapInBatch);
    void MetaImpute(int Sample);
    void PrintMetaImputedData();
    void PrintMetaWeight();
    void PrintVariantInfo();
    void PrintVariantPartialInfo();
    void PrintWeightVariantInfo();
    void PrintMetaImputedRsq();

    string CreateInfo();
    string CreateInfo(int i);
    string CreatePartialInfo();
    string CreateRsqInfo();
    void PrintDiploidDosage(float &x, float &y);
    void PrintHaploidDosage(float &x);
    void PrintWeightForHaplotype(int haploId);

};
#endif //METAM_METAMINIMAC_H
