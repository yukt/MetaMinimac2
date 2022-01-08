#include "HaplotypeSet.h"
#include "assert.h"

#include <savvy/reader.hpp>

bool HaplotypeSet::LoadSampleNames(std::string empDoseFN, std::string doseFN)
{
    if(!doesExistFile(empDoseFN)) return false;
    if(!doesExistFile(doseFN)) return false;
    EmpDoseFileName = empDoseFN;
    DoseFileName = doseFN;

    if (!GetSampleInformationfromHDS(DoseFileName))
        return false;
    int tempNoSamples=numSamples;
    vector<string> tempindividualName=individualName;
    vector<int> tempSampleNoHaplotypes=SampleNoHaplotypes;
    GetSampleInformation(EmpDoseFileName);
    if(!CheckSampleConsistency(tempNoSamples,tempindividualName,tempSampleNoHaplotypes,DoseFileName,EmpDoseFileName)) return false;

    return true;

}


bool HaplotypeSet::GetSampleInformationfromHDS(string filename)
{
    savvy::reader inFile(filename);
    savvy::variant record;
    individualName.clear();

    if (!inFile)
    {
        cout << "\n Program could NOT open file : " << filename << endl<<endl;
        return false;
    }
    numSamples = inFile.samples().size();
    if(numSamples==0)
    {
        std::cout << "\n Number of Samples read from VCF File    : " << numSamples << endl;
        std::cout << "\n ERROR !!! "<<endl;
        cout << "\n NO samples found in VCF File !! \n Please Check Input File !!!  "<< endl;
        return false;
    }
    individualName.resize(numSamples);
    for (int i = 0; i < numSamples; i++)
    {
        individualName[i]=inFile.samples()[i];
    }

    if (!inFile.read(record))
        return cout << "\nERROR: Program could NOT read record from file : " << filename << endl<<endl, false;

    std::vector<float> hds_vec;
    if (!record.get_format("HDS", hds_vec))
        return cout << "\nERROR: HDS FORMAT field not found in record from file : " << filename << endl<<endl, false;

    std::size_t max_ploidy = hds_vec.size() / numSamples;
    CummulativeSampleNoHaplotypes.resize(numSamples);
    SampleNoHaplotypes.resize(numSamples, max_ploidy);

    if (max_ploidy == 2)
    {
        for (std::size_t i = 0; i < numSamples; ++i)
        {
            if (savvy::typed_value::is_end_of_vector(hds_vec[i * max_ploidy + 1]))
                SampleNoHaplotypes[i] = 1;
        }
    }
    else if (max_ploidy != 1)
    {
        cout << "\nERROR: Only haploid and diploid are supported (detected " << max_ploidy << ") : " << filename << endl<<endl;
        return false;
    }

    numActualHaps= 0;
    for (std::size_t i = 0; i < SampleNoHaplotypes.size(); ++i)
    {
        CummulativeSampleNoHaplotypes[i] = numActualHaps;
        numActualHaps += SampleNoHaplotypes[i];
    }

    return true;

}


bool HaplotypeSet::GetSampleInformation(string filename)
{
    savvy::reader inFile(filename);
    savvy::variant record;
    individualName.clear();
    SampleNoHaplotypes.clear();
    CummulativeSampleNoHaplotypes.clear();

    if (!inFile)
    {
        cout << "\n Program could NOT open file : " << filename << endl<<endl;
        return false;
    }
    numSamples = inFile.samples().size();
    if(numSamples==0)
    {
        std::cout << "\n Number of Samples read from VCF File    : " << numSamples << endl;
        std::cout << "\n ERROR !!! "<<endl;
        cout << "\n NO samples found in VCF File !! \n Please Check Input File !!!  "<< endl;
        return false;
    }
    individualName.resize(numSamples);

    for (int i = 0; i < numSamples; i++)
    {
        individualName[i]=inFile.samples()[i];
    }

    if (!inFile.read(record))
        return cout << "\nERROR: Program could NOT read record from file : " << filename << endl<<endl, false;

    std::vector<std::int8_t> gt_vec;
    if (!record.get_format("GT", gt_vec))
        return cout << "\nERROR: GT FORMAT field not found in record from file : " << filename << endl<<endl, false;

    std::size_t max_ploidy = gt_vec.size() / numSamples;
    CummulativeSampleNoHaplotypes.resize(numSamples);
    SampleNoHaplotypes.resize(numSamples, max_ploidy);

    if (max_ploidy == 2)
    {
        for (std::size_t i = 0; i < numSamples; ++i)
        {
            if (savvy::typed_value::is_end_of_vector(gt_vec[i * max_ploidy + 1]))
                SampleNoHaplotypes[i] = 1;
        }
    }
    else if (max_ploidy != 1)
    {
        cout << "\nERROR: Only haploid and diploid are supported (detected " << max_ploidy << ") : " << filename << endl<<endl;
        return false;
    }

    numActualHaps = 0;
    for (std::size_t i = 0; i < SampleNoHaplotypes.size(); ++i)
    {
        CummulativeSampleNoHaplotypes[i] = numActualHaps;
        numActualHaps += SampleNoHaplotypes[i];
    }

    return true;

}


bool HaplotypeSet::CheckSampleConsistency(int tempNoSamples,
                                          vector<string> &tempindividualName,
                                          vector<int> tempSampleNoHaplotypes,
                                          string File1, string File2)
{
    if(tempNoSamples!=numSamples)
    {
        cout<<"\n ERROR !!! "<<endl;
        cout<< " "<<File1<<" has "<<tempNoSamples<<" samples, but "<< File2<<" has "<<numSamples<<" samples !!!"<<endl;
        return false;
    }
    for(int i=0; i<numSamples; i++)
    {
        if(tempindividualName[i]!=individualName[i])
        {
            cout<<"\n  ERROR !!! "<<endl;
            cout<< " "<<File1<<" and "<< File2<<" have different samples orders !!! "<<endl;
            return false;
        }
        if(tempSampleNoHaplotypes[i]!=SampleNoHaplotypes[i])
        {
            cout<<"\n  ERROR !!! "<<endl;
            cout<< " "<<DoseFileName<<" and "<< EmpDoseFileName<<" have different ploidy for sample ID ["<< individualName[i]<<"]  !!! "<<endl;
            return false;
        }
    }
    return true;
}

void HaplotypeSet::LoadEmpVariantList()
{
    savvy::reader inFile(EmpDoseFileName);
    savvy::variant record;
    TypedVariantList.clear();

    int numReadRecords=0;

    while (inFile.read(record))
    {
        ++numReadRecords;
        if(numReadRecords==1)
            finChromosome = record.chrom();

        variant tempVariant;
        tempVariant.chr=record.chrom();
        tempVariant.bp=record.pos();
        tempVariant.altAlleleString = record.alts().empty() ? "" : record.alts()[0];
        tempVariant.refAlleleString = record.ref();
        tempVariant.name=tempVariant.chr+":"+to_string(tempVariant.bp)+":"+ tempVariant.refAlleleString+":"+tempVariant.altAlleleString;
        TypedVariantList.push_back(tempVariant);
    }

    noTypedMarkers=TypedVariantList.size();
}

void HaplotypeSet::ClearEmpVariantList()
{
    TypedVariantList.clear();
}

void HaplotypeSet::ReadBasedOnSortCommonGenotypeList(vector<string> &SortedCommonGenoList, int StartSamId, int EndSamId)

{
    savvy::reader inFile(EmpDoseFileName);
    savvy::variant record;
    int bp,numReadRecords=0;
    int numHapsInBatch = 2*(EndSamId - StartSamId);
    string cno,name,refAllele,altAllele,prevID="",currID;

    LooDosage.clear();
    TypedGT.clear();
    LooDosage.resize(numHapsInBatch);
    TypedGT.resize(numHapsInBatch);
    for(int i=0; i<numHapsInBatch; i++)
    {
        LooDosage[i].resize(SortedCommonGenoList.size());
        TypedGT[i].resize(SortedCommonGenoList.size());
    }
    int SortIndex = 0;
    int numComRecord = 0;
    std::vector<std::int8_t> gt_vec;
    std::vector<float> lds_vec;
    while (inFile.read(record))
    {
        ++numReadRecords;
        cno=record.chrom();
        bp=record.pos();
        altAllele = record.alts().empty() ? "" : record.alts()[0];
        refAllele = record.ref();
        name = cno+":"+to_string(bp)+":"+refAllele+":"+altAllele;

        if(SortedCommonGenoList[SortIndex]==name)
        {
            if (!record.get_format("GT", gt_vec) || !record.get_format("LDS", lds_vec))
            {
                cout << " ERROR: GT or LDS missing from empirical dosage file (" << EmpDoseFileName << ")" <<endl;
                abort();
            }
            LoadLooVariant(gt_vec, lds_vec, numComRecord, StartSamId, EndSamId);
            numComRecord++;
            SortIndex++;
        }
    }
    if(SortedCommonGenoList.size()!=SortIndex)
    {
        cout<<" ERROR CODE 2819: Please contact author with this code to help with bug fixing ..."<<endl;
        abort();
    }
}

void HaplotypeSet::LoadCurrentGT(savvy::variant & rec)
{
    rec.get_format("GT", CurrentHapDosage);
}

void HaplotypeSet::LoadLooVariant(const std::vector<std::int8_t>& gt, const std::vector<float>& lds, int loonumReadRecords, int StartSamId, int EndSamId)
{
    assert(gt.size() == numActualHaps);
    assert(gt.size() == lds.size());
    int ploidy = numActualHaps / numSamples;
    int start_i = StartSamId * ploidy;
    int end_i = EndSamId * ploidy;
    for (int i = start_i; i < end_i; ++i)
        TypedGT[i - start_i][loonumReadRecords] = gt[i];

    for (int i = StartSamId * ploidy; i < end_i; ++i)
        LooDosage[i - start_i][loonumReadRecords] = lds[i];
}

bool HaplotypeSet::doesExistFile(string filename)
{
    IFILE ifs = ifopen(filename.c_str(), "r");
    if (ifs)
    {
        ifclose(ifs);
        return true;
    }
    else
    {
        ifclose(ifs);
        return false;
    }
}

bool HaplotypeSet::CheckSuffixFile(string prefix, const char* suffix, string &FinalName)
{
    if(doesExistFile(prefix+"."+suffix+".vcf"))
        FinalName=prefix+"."+suffix+".vcf";
    else if(doesExistFile(prefix+"."+suffix+".vcf.gz"))
        FinalName=prefix+"."+suffix+".vcf.gz";
    else
    {
        cout<<"\n No VCF file found ("<<prefix<<"."<<suffix<<".vcf or "<<prefix+"."<<suffix<<".vcf.gz) "<<endl;
        cout<<" Please check input file prefix ["<< prefix <<"] properly ... "<<endl;
        return false;
    }
    return true;
}

void HaplotypeSet::LoadData(int VariantId, savvy::variant &rec)
{
    BufferHapDosage.emplace_back();
    rec.get_format("HDS", BufferHapDosage.back());
    VariantId2Buffer[VariantId] = BufferNoVariants++;
}

void HaplotypeSet::GetData(int VariantId)
{
    CurrentHapDosage.clear();
    int index=VariantId2Buffer[VariantId];
    CurrentHapDosage=BufferHapDosage[index];
}

void HaplotypeSet::ClearBuffer()
{
    noMarkers += BufferNoVariants;
    BufferHapDosage.clear();
    BufferNoVariants = 0;
    VariantId2Buffer.clear();
}