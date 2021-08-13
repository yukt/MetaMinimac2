#include "HaplotypeSet.h"
#include "assert.h"

bool HaplotypeSet::LoadSampleNames(string prefix)
{
    InfilePrefix.Copy(prefix.c_str());
    if(!CheckSuffixFile(prefix,"dose", DoseFileName)) return false;
    if(!CheckSuffixFile(prefix,"empiricalDose", EmpDoseFileName)) return false;

    GetSampleInformationfromHDS(DoseFileName);
    int tempNoSamples=numSamples;
    vector<string> tempindividualName=individualName;
    vector<int> tempSampleNoHaplotypes=SampleNoHaplotypes;
    GetSampleInformation(EmpDoseFileName);
    if(!CheckSampleConsistency(tempNoSamples,tempindividualName,tempSampleNoHaplotypes,DoseFileName,EmpDoseFileName)) return false;

    return true;

}


bool HaplotypeSet::GetSampleInformationfromHDS(string filename)
{
    VcfFileReader inFile;
    VcfHeader header;
    VcfRecord record;
    individualName.clear();

    if (!inFile.open(filename.c_str(), header))
    {
        cout << "\n Program could NOT open file : " << filename << endl<<endl;
        return false;
    }
    numSamples = header.getNumSamples();
    if(numSamples==0)
    {
        std::cout << "\n Number of Samples read from VCF File    : " << numSamples << endl;
        std::cout << "\n ERROR !!! "<<endl;
        cout << "\n NO samples found in VCF File !! \n Please Check Input File !!!  "<< endl;
        return false;
    }
    individualName.resize(numSamples);
    CummulativeSampleNoHaplotypes.resize(numSamples);
    SampleNoHaplotypes.resize(numSamples);

    for (int i = 0; i < numSamples; i++)
    {
        individualName[i]=header.getSampleName(i);
    }

    inFile.setSiteOnly(false);
    inFile.readRecord(record);


    VcfRecordGenotype &ThisGenotype=record.getGenotypeInfo();

    int tempHapCount=0;
    for (int i = 0; i<(numSamples); i++)
    {
        string temp=*ThisGenotype.getString("HDS",i);
        char *end_str;

        char *pch = strtok_r((char *) temp.c_str(), ",", &end_str);
        if(pch==NULL)
        {
            std::cout << "\n ERROR !!! \n Empty Value for Individual : " << individualName[i] << " at First Marker  " << endl;
            std::cout << " Most probably a corrupted VCF file. Please check input VCF file !!! " << endl;
            cout << "\n Program Exiting ... \n\n";
            return false;
        }
        char *pch1 = strtok_r(NULL, "\t", &end_str);

        if(pch1==NULL)
        {
            SampleNoHaplotypes[i]=1;
        }
        else
        {
            SampleNoHaplotypes[i]=2;
        }

        CummulativeSampleNoHaplotypes[i]=tempHapCount;
        tempHapCount+=SampleNoHaplotypes[i];
    }

    inFile.close();

    return true;

}


bool HaplotypeSet::GetSampleInformation(string filename)
{
    VcfFileReader inFile;
    VcfHeader header;
    VcfRecord record;
    individualName.clear();

    if (!inFile.open(filename.c_str(), header))
    {
        cout << "\n Program could NOT open file : " << filename << endl<<endl;
        return false;
    }
    numSamples = header.getNumSamples();
    if(numSamples==0)
    {
        std::cout << "\n Number of Samples read from VCF File    : " << numSamples << endl;
        std::cout << "\n ERROR !!! "<<endl;
        cout << "\n NO samples found in VCF File !! \n Please Check Input File !!!  "<< endl;
        return false;
    }
    individualName.resize(numSamples);
    CummulativeSampleNoHaplotypes.resize(numSamples);
    SampleNoHaplotypes.resize(numSamples);

    for (int i = 0; i < numSamples; i++)
    {
        individualName[i]=header.getSampleName(i);
    }

    inFile.setSiteOnly(false);
    inFile.readRecord(record);
    int tempHapCount=0;
    for (int i = 0; i<(numSamples); i++)
    {
        if(record.getNumGTs(i)==0)
        {
            std::cout << "\n ERROR !!! \n Empty Value for Individual : " << individualName[i] << " at First Marker  " << endl;
            std::cout << " Most probably a corrupted VCF file. Please check input VCF file !!! " << endl;
            cout << "\n Program Exiting ... \n\n";
            return false;
        }
        else
        {
            CummulativeSampleNoHaplotypes[i]=tempHapCount;
            SampleNoHaplotypes[i]=(record.getNumGTs(i));
            tempHapCount+=SampleNoHaplotypes[i];
        }
    }
    inFile.close();
    numActualHaps=tempHapCount;

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
    VcfFileReader inFile;
    VcfHeader header;
    VcfRecord record;
    inFile.open(EmpDoseFileName.c_str(), header);
    inFile.setSiteOnly(true);
    TypedVariantList.clear();

    int numReadRecords=0;

    while (inFile.readRecord(record))
    {
        ++numReadRecords;
        if(numReadRecords==1)
            finChromosome = record.getChromStr();

        variant tempVariant;
        tempVariant.chr=record.getChromStr();
        tempVariant.bp=record.get1BasedPosition();
        tempVariant.altAlleleString = record.getAltStr();
        tempVariant.refAlleleString = record.getRefStr();
        tempVariant.name=tempVariant.chr+":"+to_string(tempVariant.bp)+":"+ tempVariant.refAlleleString+":"+tempVariant.altAlleleString;
        TypedVariantList.push_back(tempVariant);
    }

    noTypedMarkers=TypedVariantList.size();
    inFile.close();
}

void HaplotypeSet::ClearEmpVariantList()
{
    TypedVariantList.clear();
}

void HaplotypeSet::ReadBasedOnSortCommonGenotypeList(vector<string> &SortedCommonGenoList, int StartSamId, int EndSamId)

{
    VcfFileReader inFile;
    VcfHeader header;
    VcfRecord record;
    VcfRecordGenotype *recordGeno;
    inFile.open(EmpDoseFileName.c_str(), header);
    inFile.setSiteOnly(false);
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
    while (inFile.readRecord(record))
    {
        ++numReadRecords;
        cno=record.getChromStr();
        bp=record.get1BasedPosition();
        altAllele = record.getAltStr();
        refAllele = record.getRefStr();
        name = cno+":"+to_string(bp)+":"+refAllele+":"+altAllele;

        if(SortedCommonGenoList[SortIndex]==name)
        {
            recordGeno=&record.getGenotypeInfo();
            LoadLooVariant(*recordGeno, numComRecord, StartSamId, EndSamId);
            numComRecord++;
            SortIndex++;
        }
    }
    if(SortedCommonGenoList.size()!=SortIndex)
    {
        cout<<" ERROR CODE 2819: Please contact author with this code to help with bug fixing ..."<<endl;
        abort();
    }

    inFile.close();
}

void HaplotypeSet::LoadCurrentGT(VcfRecordGenotype & ThisGenotype)
{
    CurrentHapDosage.clear();
    CurrentHapDosage.resize(2*numSamples);
    for (int i=0; i<numSamples; i++)
    {
        string temp = *ThisGenotype.getString("GT", i);
        if(SampleNoHaplotypes[i]==2)
        {
            char *end_str1;
            char *pch1 = strtok_r((char *) temp.c_str(), "|", &end_str1);
            CurrentHapDosage[2*i] = atof(pch1);
            pch1 = strtok_r(NULL, "\t", &end_str1);
            CurrentHapDosage[2*i+1] = atof(pch1);
        }
        else
        {
            CurrentHapDosage[2*i] = atof(temp.c_str());
        }
    }
}

void HaplotypeSet::LoadLooVariant(VcfRecordGenotype &ThisGenotype,int loonumReadRecords, int StartSamId, int EndSamId)
{
    int NoHapsLoad = 0;
    for (int i = StartSamId; i<EndSamId; i++)
    {
        string temp=*ThisGenotype.getString("LDS",i);
        if(SampleNoHaplotypes[i]==2)
        {
            char *end_str;
            char *pch = strtok_r((char *) temp.c_str(), "|", &end_str);
            LooDosage[2*NoHapsLoad][loonumReadRecords] = atof(pch);
            pch = strtok_r(NULL, "\t", &end_str);
            LooDosage[2*NoHapsLoad + 1][loonumReadRecords] = atof(pch);


            temp = *ThisGenotype.getString("GT", i);
            char *end_str1;
            char *pch1 = strtok_r((char *) temp.c_str(), "|", &end_str1);
            TypedGT[2*NoHapsLoad][loonumReadRecords] = atof(pch1);
            pch1 = strtok_r(NULL, "\t", &end_str1);
            TypedGT[2*NoHapsLoad+1][loonumReadRecords] = atof(pch1);
        }
        else
        {
            LooDosage[2*NoHapsLoad][loonumReadRecords] = atof(temp.c_str());
            temp = *ThisGenotype.getString("GT", i);
            TypedGT[2*NoHapsLoad][loonumReadRecords] = atof(temp.c_str());

        }
        NoHapsLoad++;
    }
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

void HaplotypeSet::LoadData(int VariantId, VcfRecordGenotype &ThisGenotype, int StartSamId, int EndSamId)
{
    VariantId2Buffer[VariantId] = BufferNoVariants;

    vector<float> tempHapDosage;
    tempHapDosage.resize(2*(EndSamId-StartSamId), 0.0);

    for (int i = StartSamId; i<EndSamId; i++)
    {
        string temp=*ThisGenotype.getString("HDS",i);
        char *end_str;

        if(SampleNoHaplotypes[i]==2) {
            char *pch = strtok_r((char *) temp.c_str(), ",", &end_str);
            tempHapDosage[2*(i-StartSamId)] = atof(pch);

            pch = strtok_r(NULL, "\t", &end_str);
            tempHapDosage[2*(i-StartSamId)+1] = atof(pch);
        }
        else
        {
            tempHapDosage[2*(i-StartSamId)] = atof(temp.c_str());
        }

    }

    BufferHapDosage.push_back(tempHapDosage);
    BufferNoVariants++;

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