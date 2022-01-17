#ifndef METAM_HAPLOTYPESET_H
#define METAM_HAPLOTYPESET_H

#include <map>
#include <cassert>

#include <savvy/site_info.hpp>

using namespace std;

class variant
{
public:

    string name;
    int bp;
    string chr;
    string refAlleleString,altAlleleString;
    vector<int> StudiesHasVariant;
    int NoStudiesHasVariant;

    variant() {};
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

    // Buffer Data
    int BufferNoVariants;
    vector<vector<float> > BufferHapDosage;
    map<int,int> VariantId2Buffer;



    // FUNCTIONS
    bool        CheckSampleConsistency                  (int tempNoSamples, vector<string> &tempindividualName, vector<int> tempSampleNoHaplotypes, string File1, string File2);
    void        ReadBasedOnSortCommonGenotypeList       (vector<string> &SortedCommonGenoList, int StartSamId, int EndSamId);
    bool        CheckSuffixFile                         (string prefix, const char* suffix, string &FinalName);

    bool        GetSampleInformation                    (string filename);
    bool        GetSampleInformationfromHDS             (string filename);
    void        LoadEmpVariantList                      ();
    void        ClearEmpVariantList                     ();
    void        LoadCurrentGT                           (savvy::variant &ThisGenotype);
    void        LoadLooVariant                          (const std::vector<std::int8_t>& gt, const std::vector<float>& lds,int loonumReadRecords, int StartSamId, int EndSamId);
    bool        LoadSampleNames                         (std::string empDoseFN, std::string doseFN);
    bool        doesExistFile                           (string filename);

    void        LoadData                                (int VariantId, savvy::variant &ThisGenotype);
    void        GetData                                 (int VariantId);
    void        ClearBuffer                             ();
};



#endif //METAM_HAPLOTYPESET_H
