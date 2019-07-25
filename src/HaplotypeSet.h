#ifndef METAM_HAPLOTYPESET_H
#define METAM_HAPLOTYPESET_H

#include "VcfFileReader.h"
#include "VcfHeader.h"
#include "assert.h"

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
    int         numSamples;
    int         numActualHaps;
    vector<string> individualName;
    vector<int> SampleNoHaplotypes;
    vector<int> CummulativeSampleNoHaplotypes;
    vector<variant> VariantList;
    vector<variant> TypedVariantList;
    int noMarkers;
    int noTypedMarkers;
    string finChromosome;


    // Dosage Data
    vector<float> CurrentHapDosage;
    vector<vector<float> > LooDosage;
    vector<vector<float> > TypedGT;


    // FUNCTIONS
    bool        CheckSampleConsistency                  (int tempNoSamples, vector<string> &tempindividualName, vector<int> tempSampleNoHaplotypes, string File1, string File2);
    void        ReadBasedOnSortCommonGenotypeList       (vector<string> &SortedCommonGenoList, int StartSamId, int EndSamId);
    bool        CheckSuffixFile                         (string prefix, const char* suffix, string &FinalName);
    void        LoadHapDoseVariant                      (VcfRecordGenotype &ThisGenotype, int StartSamId, int EndSamId);


    bool        GetSampleInformation                    (string filename);
    bool        GetSampleInformationfromHDS                    (string filename);
    void        LoadLooVariant                          (VcfRecordGenotype &ThisGenotype,int loonumReadRecords, int StartSamId, int EndSamId);
    bool        LoadSampleNames                         (string prefix);
    bool        doesExistFile                           (string filename);
};



#endif //METAM_HAPLOTYPESET_H
