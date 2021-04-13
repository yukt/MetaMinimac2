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

    // Variables for input weight file stream and records
    VcfFileReader* InputWeightStream;
    VcfRecord* CurrentRecordFromWeight;
    vector<vector<float>> WeightsFromCurrentRecord, WeightsFromPreviousRecord;

    // Process Part of Samples each time
    int StartSamId, EndSamId;
    double lambda;
    vector<double> TransitionProb;
    double Recom, backgroundError;
    double JumpFix, JumpThreshold;
    vector<vector<vector<double>>> Weights;
    vector<vector<double>> PrevRightProb;
    int NoCommonVariantsProcessed;

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
    int PrevBp, CurrBp;
    vector<vector<double>> *PrevWeights;
    vector<vector<double>> *CurrWeights;
    vector<float> CurrentMetaImputedDosage;

    float CurrentHapDosageSum, CurrentHapDosageSumSq;

    // Buffer
    int NoRecords;
    int BufferBp, BufferNoVariants;
    vector<variant> BufferVariantList;


    MetaMinimac()
    {
        lambda = 2e-7;
        Recom = 1e-3;
        backgroundError = 1e-5;
        JumpThreshold = 1e-10;
        JumpFix = 1e10;
    };



    String Analyze();

    bool ParseInputVCFFiles();
    bool CheckSampleNameCompatibility();
    void OpenStreamInputDosageFiles(bool siteOnly);
    void OpenStreamInputWeightFiles();
    void CloseStreamInputDosageFiles();
    void CloseStreamInputWeightFiles();
    bool OpenStreamOutputDosageFiles();
    bool OpenStreamOutputWeightFiles();
    string GetDosageFileFullName(String prefix);
    bool doesExistFile(String filename);

    bool LoadEmpVariantInfo();
    void FindCommonGenotypedVariants();
    void FindCurrentMinimumPosition();
    int IsVariantEqual(VcfRecord &Rec1, VcfRecord &Rec2);
    void UpdateCurrentRecords();

    void LoadLooDosage();


    bool PerformWeightEstimation();
    void CalculateWeights();
    void InitiateWeights();
    void CalculateLeftProbs();
    void CalculatePosterior();
    void InitiateLeftProb(int SampleInBatch);
    void InitiateRightProb(int SampleInBatch);
    void UpdateOneStepLeft(int SampleInBatch);
    void UpdateOneStepRight(int SampleInBatch);
    void OutputWeights();
    void OutputPartialWeights();
    void OutputAllWeights();

    String PerformFinalAnalysis();
    bool InitiateWeightsFromRecord();
    void ReadCurrentWeights();
    void StoreWeightsFromCurrentRecord(int i, char* str);
    void CopyCurrentWeightsToPreviousWeights();
    void UpdateWeights();
    void AppendtoMainWeightsFile();

    void MetaImputeCurrentBuffer();
    void ClearCurrentBuffer();
    void ReadCurrentDosageData();
    void CreateMetaImputedData(int VariantId);
    void MetaImpute(int Sample);
    void PrintMetaImputedData();
    void PrintMetaWeight();
    void PrintVariantInfo();
    void PrintWeightVariantInfo();

    string CreateInfo();
    void PrintDiploidDosage(float &x, float &y);
    void PrintHaploidDosage(float &x);
    void PrintWeightForHaplotype(int haploId);
    void summary()
    {
        for (int i=0; i<NoInPrefix; i++)
        {
            cout << " Study " << i+1 << " #Markers = " << InputData[i].noMarkers << endl;
        }

        cout << " Total #Markers   = " << NoVariants << endl;
    }

};
#endif //METAM_METAMINIMAC_H
