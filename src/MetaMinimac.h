#ifndef METAM_METAMINIMAC_H
#define METAM_METAMINIMAC_H

#include "MyVariables.h"
#include "HaplotypeSet.h"

using namespace std;

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
    vector<vector<vector<double>>> Posterior;

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

    bool LoadLooDosage();

    int PerformFinalAnalysis();
    void GetMetaEstimate(int Sample, int SampleInBatch);
    void FlushPartialVcf(int batchNo);

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
