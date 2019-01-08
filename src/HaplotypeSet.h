#ifndef METAM_HAPLOTYPESET_H
#define METAM_HAPLOTYPESET_H

#include "VcfFileReader.h"
#include "VcfHeader.h"

using namespace std;

class variant
{
public:

    string name;
    int bp;
    string chr;
    string refAlleleString,altAlleleString;
    bool typed;
    vector<int> StudiesHasVariant;
    int NoStudiesHasVariant;

    variant()
    {
        typed=false;
    };
    variant(string &id,string &CHR,int &BP)
    {
        name=id;
        bp=BP;
        chr=CHR;
    };
};

class HaplotypeSet
{


public:

    // File Name Variables
    String InfilePrefix;
    string DoseFileName;
    string EmpDoseFileName;

    // Summary Variables
    int         numHaplotypes,numSamples;
    int         numActualHaps;
    int         numMarkers;
    vector<string> individualName;
    vector<int> SampleNoHaplotypes;
    vector<int> CummulativeSampleNoHaplotypes;
    vector<variant> VariantList;
    vector<variant> TypedVariantList;
    int noTypedMarkers;
    string finChromosome;


    // Dosage Data
    vector<float> CurrentHapDosage;
    vector<vector<double> > LooDosage;
    vector<vector<double> > TypedGT;


    // FUNCTIONS
    void        LoadEmpVariantList                      ();
    void        LoadVariantList                         (string inFile);
    bool        CheckSampleConsistency                  (int tempNoSamples, vector<string> &tempindividualName, vector<int> tempSampleNoHaplotypes, string File1, string File2);
    void        ReadBasedOnSortCommonGenotypeList       (vector<string> &SortedCommonGenoList);
    void        SortCommonGenotypeList                  (std::unordered_set<string> &CommonGenotypeVariantNameList, vector<string> &SortedCommonGenoList, vector<variant> &CommonTypedVariantList);
    bool        CheckSuffixFile                         (string prefix, const char* suffix, string &FinalName);
    void        LoadHapDoseVariant                      (VcfRecordGenotype &ThisGenotype);


    bool        GetSampleInformation                    (string filename);
    bool        GetSampleInformationfromHDS                    (string filename);
    void        LoadLooVariant                          (VcfRecordGenotype &ThisGenotype,int loonumReadRecords);
    bool        LoadSampleNames                         (string prefix);
    bool        doesExistFile                           (string filename);
};



#endif //METAM_HAPLOTYPESET_H
