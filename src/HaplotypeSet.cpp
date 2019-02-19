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
    numHaplotypes=tempHapCount;

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
    numHaplotypes=2*numSamples;
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
    LoadVariantList(EmpDoseFileName);
    TypedVariantList=VariantList;
    noTypedMarkers=numMarkers;
}

void HaplotypeSet::LoadVariantList(string FileName)
{
    VcfFileReader inFile;
    VcfHeader header;
    VcfRecord record;
    inFile.open(FileName.c_str(), header);
    inFile.setSiteOnly(true);
    VariantList.clear();

    int bp,numReadRecords=0;
    string cno,id,refAllele,altAllele,prevID="",currID;

    while (inFile.readRecord(record))
    {
        ++numReadRecords;

        if(numReadRecords==1)
            finChromosome=record.getChromStr();
        variant tempVariant;
        tempVariant.chr=record.getChromStr();
        tempVariant.bp=record.get1BasedPosition();
        tempVariant.name=record.getIDStr();
        tempVariant.altAlleleString = record.getAltStr();
        tempVariant.refAlleleString = record.getRefStr();
        VariantList.push_back(tempVariant);
    }

    numMarkers=VariantList.size();
    inFile.close();
}

void HaplotypeSet::SortCommonGenotypeList(std::unordered_set<string> &CommonGenotypeVariantNameList,
                                          vector<string> &SortedCommonGenoList,
                                          vector<variant> &CommonTypedVariantList)

{
    VcfFileReader inFile;
    VcfHeader header;
    VcfRecord record;
    VcfRecordGenotype *recordGeno;
    inFile.open(EmpDoseFileName.c_str(), header);
    inFile.setSiteOnly(false);
    SortedCommonGenoList.clear();
    int bp,numReadRecords=0;
    string cno,id,refAllele,altAllele,prevID="",currID;

    CommonTypedVariantList.clear();
    CommonTypedVariantList.resize(CommonGenotypeVariantNameList.size());
    LooDosage.clear();
    TypedGT.clear();
//    CurrentHapDosage.resize(numHaplotypes);
    LooDosage.resize(numHaplotypes);
    TypedGT.resize(numHaplotypes);
    for(int i=0; i<numHaplotypes; i++)
    {
        LooDosage[i].resize(CommonGenotypeVariantNameList.size());
        TypedGT[i].resize(CommonGenotypeVariantNameList.size());
    }

    int numComRecord = 0;
    while (inFile.readRecord(record))
    {
        ++numReadRecords;

        if(CommonGenotypeVariantNameList.find(record.getIDStr())!=CommonGenotypeVariantNameList.end())
        {
            variant &tempVariant=CommonTypedVariantList[numComRecord];
            tempVariant.chr=record.getChromStr();
            tempVariant.bp=record.get1BasedPosition();
            tempVariant.name=record.getIDStr();
            tempVariant.altAlleleString = record.getAltStr();
            tempVariant.refAlleleString = record.getRefStr();
            VariantList.push_back(tempVariant);

            recordGeno=&record.getGenotypeInfo();
            LoadLooVariant(*recordGeno, numComRecord);

            SortedCommonGenoList.push_back(record.getIDStr());

            numComRecord++;
        }
    }

    if(SortedCommonGenoList.size()!=CommonGenotypeVariantNameList.size())
    {
        cout<<" ERROR CODE 4219: Please contact author with this code to help with bug fixing ..."<<endl;
        abort();
    }
    inFile.close();
}

void HaplotypeSet::ReadBasedOnSortCommonGenotypeList(vector<string> &SortedCommonGenoList)

{
    VcfFileReader inFile;
    VcfHeader header;
    VcfRecord record;
    VcfRecordGenotype *recordGeno;
    inFile.open(EmpDoseFileName.c_str(), header);
    inFile.setSiteOnly(false);
    int bp,numReadRecords=0;
    string cno,id,refAllele,altAllele,prevID="",currID;

    LooDosage.clear();
    TypedGT.clear();
//    CurrentHapDosage.resize(numHaplotypes);
    LooDosage.resize(numHaplotypes);
    TypedGT.resize(numHaplotypes);
    for(int i=0; i<numHaplotypes; i++)
    {
        LooDosage[i].resize(SortedCommonGenoList.size());
        TypedGT[i].resize(SortedCommonGenoList.size());
    }
    int SortIndex = 0;
    int numComRecord = 0;
    while (inFile.readRecord(record))
    {
        ++numReadRecords;

        if(SortedCommonGenoList[SortIndex]==record.getIDStr())
        {
            recordGeno=&record.getGenotypeInfo();
            LoadLooVariant(*recordGeno, numComRecord);
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

void HaplotypeSet::LoadHapDoseVariant(VcfRecordGenotype &ThisGenotype)
{
    for (int i = 0; i<(numSamples); i++)
    {
        string temp=*ThisGenotype.getString("HDS",i);
        char *end_str;

        if(SampleNoHaplotypes[i]==2) {
            char *pch = strtok_r((char *) temp.c_str(), ",", &end_str);
            CurrentHapDosage[2*i] = atof(pch);

            pch = strtok_r(NULL, "\t", &end_str);
            CurrentHapDosage[2*i+1] = atof(pch);
        }
        else
        {
            CurrentHapDosage[2*i] = atof(temp.c_str());
        }

    }
}

void HaplotypeSet::LoadHapDoseVariant(VcfRecordGenotype &ThisGenotype, int StartSamId, int EndSamId)
{
    CurrentHapDosage.clear();
    CurrentHapDosage.resize(2*(EndSamId-StartSamId));

    for (int i = StartSamId; i<EndSamId; i++)
    {
        string temp=*ThisGenotype.getString("HDS",i);
        char *end_str;

        if(SampleNoHaplotypes[i]==2) {
            char *pch = strtok_r((char *) temp.c_str(), ",", &end_str);
            CurrentHapDosage[2*(i-StartSamId)] = atof(pch);

            pch = strtok_r(NULL, "\t", &end_str);
            CurrentHapDosage[2*(i-StartSamId)+1] = atof(pch);
        }
        else
        {
            CurrentHapDosage[2*(i-StartSamId)] = atof(temp.c_str());
        }

    }
}

void HaplotypeSet::LoadLooVariant(VcfRecordGenotype &ThisGenotype,int loonumReadRecords)
{
    for (int i = 0; i<(numSamples); i++)
    {
        string temp=*ThisGenotype.getString("LDS",i);
        if(SampleNoHaplotypes[i]==2)
        {
            char *end_str;
            char *pch = strtok_r((char *) temp.c_str(), "|", &end_str);
            LooDosage[2*i][loonumReadRecords] = atof(pch);
            pch = strtok_r(NULL, "\t", &end_str);
            LooDosage[2*i + 1][loonumReadRecords] = atof(pch);


            temp = *ThisGenotype.getString("GT", i);
            char *end_str1;
            char *pch1 = strtok_r((char *) temp.c_str(), "|", &end_str1);
            TypedGT[2*i][loonumReadRecords] = atof(pch1);
            pch1 = strtok_r(NULL, "\t", &end_str1);
            TypedGT[2*i+1][loonumReadRecords] = atof(pch1);
        }
        else
        {
            LooDosage[2*i][loonumReadRecords] = atof(temp.c_str());
            temp = *ThisGenotype.getString("GT", i);
            TypedGT[2*i][loonumReadRecords] = atof(temp.c_str());

        }
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

