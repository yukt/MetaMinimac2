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

    vector<variant> VariantList;
    vector<variant> CommonTypedVariantList;
    vector<string> CommonGenotypeVariantNameList;
    int NoVariants, NoCommonTypedVariants;

    // Variables for input dosage file stream and records
    vector<VcfFileReader*> InputDosageStream;
    vector<VcfFileReader*> InputEmpDosageStream;
    vector<VcfRecord*> CurrentRecordFromStudy;
    vector<VcfRecord*> CurrentEmpRecordFromStudy;
    vector<int> StudiesHasVariant;
    int CurrentFirstVariantBp;
    int NoStudiesHasVariant;

    vector<HaplotypeSet> InputData;

    vector<vector<double>> Weights;
    vector<double> CurrentInitialProb;


    int Analyze();

    bool ParseInputVCFFiles();
    bool CheckSampleNameCompatibility();
    void OpenStreamInputDosageFiles(bool siteOnly);
    void CloseStreamInputDosageFiles();
    string GetDosageFileFullName(String prefix);
    bool doesExistFile(String filename);

    bool LoadVariantInfo();
    void FindCurrentMinimumPosition();
    int IsVariantEqual(VcfRecord &Rec1, VcfRecord &Rec2);
    void ReadCurrentVariantInfo();
    void UpdateCurrentRecords();

    bool LoadLooDosage();
    void OpenStreamInputEmpDosageFiles();
    void CloseStreamInputEmpDosageFiles();
    string GetEmpDosageFileFullName(String prefix);

    int PerformFinalAnalysis();
    void GetMetaEstimate(int Sample);

};
#endif //METAM_METAMINIMAC_H
