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
    vector<bool> CommonTyped;
    int NoHaplotypes, NoSamples;
    int NoVariants, NoRead, NoCommonTypedVariants, CommonTypedVariantListCounter;

    // Variables for input dosage file stream and records
    vector<VcfFileReader*> InputDosageStream;
    vector<VcfFileReader*> InputEmpDosageStream;
    vector<VcfRecord*> CurrentRecordFromStudy;
    vector<VcfRecord*> CurrentEmpRecordFromStudy;
    vector<int> StudiesHasVariant;
    int CurrentFirstVariantBp;
    int NoStudiesHasVariant;

    vector<HaplotypeSet> InputData;

    // Process Part of Samples each time
    int StartSamId, EndSamId;
    double Recom, backgroundError;
    double JumpFix, JumpThreshold;
    vector<vector<vector<double>>> Posterior;
    vector<vector<vector<double>>> LeftProb;
    vector<vector<double>> PrevLeftProb, CurrentLeftProb;
    int NoVariantsProcessed, NoCommonVariantsProcessed;

    // Output files
    IFILE vcfdosepartial, vcfweightpartial;
    IFILE metaWeight;
    char *VcfPrintStringPointer;
    char *WeightPrintStringPointer;
    int VcfPrintStringPointerLength, WeightPrintStringPointerLength;
    int batchNo;

    variant* CurrentVariant;
    vector<float> CurrentMetaImputedDosage;
    int NoVariantsImputed;

    double CurrentHapDosageSum, CurrentHapDosageSumSq;
    vector<double> HapDosageSum, HapDosageSumSq;


    MetaMinimac()
    {
        Recom = 1e-5;
        backgroundError = 1e-5;
        JumpThreshold = 1e-10;
        JumpFix = 1e10;
    };



    int Analyze();

    bool ParseInputVCFFiles();
    bool CheckSampleNameCompatibility();
    void OpenStreamInputDosageFiles(bool siteOnly);
    void CloseStreamInputDosageFiles();
    bool OpenStreamOutputDosageFiles();
    string GetDosageFileFullName(String prefix);
    bool doesExistFile(String filename);

    bool LoadVariantInfo();
    void FindCurrentMinimumPosition();
    int IsVariantEqual(VcfRecord &Rec1, VcfRecord &Rec2);
    void ReadCurrentVariantInfo();
    void UpdateCurrentRecords();

    void LoadLooDosage();

    int PerformFinalAnalysis();
//    void GetMetaEstimate(int Sample, int SampleInBatch);
//    void GetMetaEstimate(int SampleInBatch);
    void FlushPartialVcf(int batchNo);

    void InitLeftProb();
    void InitLeftProb(int SampleInBatch);
    void WalkLeft();
    void ProcessTypedLeftProb();
    void ProcessLeftProb(int SampleInBatch);

    void AppendtoMainVcf();
    void AppendtoMainWeightsFile();

    void ReadCurrentDosageData();
    void CreateMetaImputedData();
    void MetaImpute(int Sample);
    void PrintMetaImputedData();
    void PrintMetaWeight();

    string CreateInfo(int i);
    void PrintDiploidDosage(float &x, float &y);
    void PrintHaploidDosage(float &x);
    void PrintWeightForHaplotype(int haploId);

};
#endif //METAM_METAMINIMAC_H
