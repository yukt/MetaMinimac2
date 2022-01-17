#ifndef METAM_METAMINIMAC_H
#define METAM_METAMINIMAC_H

#include "MyVariables.h"
#include "HaplotypeSet.h"

#include <savvy/reader.hpp>
#include <savvy/writer.hpp>

#include <memory>

using namespace std;

void logitTransform(vector<double> &From,vector<double> &To);

class MetaMinimac
{
public:
    UserVariables myUserVariables;
    int NoInPrefix;
    string finChromosome;

    vector<variant> CommonTypedVariantList;
    vector<string> CommonGenotypeVariantNameList;
    int NoHaplotypes, NoSamples;
    int NoVariants, NoCommonTypedVariants;

    // Variables for input dosage file stream and records
    vector<savvy::reader*> InputDosageStream;
    vector<savvy::variant*> CurrentRecordFromStudy;
    vector<int> StudiesHasVariant;
    int CurrentFirstVariantBp;
    int NoStudiesHasVariant;

    vector<HaplotypeSet> InputData;

    // Variables for input weight file stream and records
    savvy::reader* InputWeightStream;
    savvy::variant* CurrentRecordFromWeight;
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
    std::unique_ptr<savvy::writer> metaDoseOut;
    int batchNo;
    const int BinScalar = 1000;

    variant* CurrentVariant;
    int PrevBp, CurrBp;
    vector<vector<double>> *PrevWeights;
    vector<vector<double>> *CurrWeights;
    vector<float> CurrentMetaImputedDosage;
    savvy::compressed_vector<float> sparse_dosages_;
    savvy::compressed_vector<std::int8_t> sparse_gt_;
    std::vector<float> dense_float_vec_;

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



    std::string Analyze();

    bool ParseInputVCFFiles();
    bool CheckSampleNameCompatibility();
    void OpenStreamInputDosageFiles();
    void OpenStreamInputEmpiricalDoseFiles();
    void OpenStreamInputWeightFiles();
    void CloseStreamInputDosageFiles();
    void CloseStreamInputWeightFiles();
    bool OpenStreamOutputDosageFiles();

    bool LoadEmpVariantInfo();
    void FindCommonGenotypedVariants();
    bool CheckPhasingConsistency();
    void FindCurrentMinimumPosition();
    int IsVariantEqual(savvy::variant &Rec1, savvy::variant &Rec2);
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
    void OutputAllWeights(const std::string& out_file_path);

    std::string PerformFinalAnalysis();
    bool InitiateWeightsFromRecord();
    void ReadCurrentWeights();
    void CopyCurrentWeightsToPreviousWeights();
    void UpdateWeights();
    bool AppendtoMainWeightsFile();

    void MetaImputeCurrentBuffer();
    void ClearCurrentBuffer();
    void ReadCurrentDosageData();
    void CreateMetaImputedData(int VariantId);
    void MetaImpute(int Sample);
    void SetMetaImputedData(savvy::variant& out_var);
    void SetMetaWeightVecs(std::vector<std::vector<float>>& wt_vecs);
    void SetWeightVariantInfo(savvy::variant& dest);

    void CreateInfo(savvy::variant& rec);
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
