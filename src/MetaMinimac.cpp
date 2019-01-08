#include "MarkovModel.h"
#include <iostream>
#include "simplex.h"

using BT::Simplex;
using namespace std;

const int MAXBP = 999999999;

int MetaMinimac::Analyze()
{
    ParseInputVCFFiles();
    CheckSampleNameCompatibility();

    OpenStreamInputDosageFiles(true);
    LoadVariantInfo();
    CloseStreamInputDosageFiles();

    LoadLooDosage();

    return PerformFinalAnalysis();
}

bool MetaMinimac::ParseInputVCFFiles()
{

    InPrefixList.clear();
    size_t pos = 0;
    std::string delimiter(myUserVariables.FileDelimiter) ;
    std::string token;
    int Count=0;
    string tempName=myUserVariables.inputFiles.c_str();
    while ((pos = tempName.find(delimiter)) != std::string::npos)
    {
        token = tempName.substr(0, pos);
        InPrefixList.push_back(token.c_str());
        tempName.erase(0, pos + delimiter.length());
        Count++;
    }
    InPrefixList.push_back(tempName.c_str());


    NoInPrefix=(int)InPrefixList.size();
    InputData.clear();
    InputData.resize(NoInPrefix);

    cout<<endl<<  " Number of Studies                  : "<<NoInPrefix<<endl;
    for(int i=0;i<NoInPrefix;i++)
    {
        cout<<  " -- Study "<<i+1<<" Prefix                  : "<<InPrefixList[i]<<endl;
    }


    if(NoInPrefix<2)
    {
        cout<<"\n ERROR ! Must have at least 2 studies for meta-imputation to work !!! "<<endl;
        cout<<" Program aborting ... "<<endl<<endl;
        return false;
    }
    if(NoInPrefix>4)
    {
        cout<<"\n ERROR ! Must have less than 5 studies for meta-imputation to work !!! "<<endl;
        cout<<" Program aborting ... "<<endl<<endl;
        return false;
    }

    return true;
}

bool MetaMinimac::CheckSampleNameCompatibility()
{
    cout<<"\n Checking Sample Compatibility across files ... "<<endl;

    for(int i=0;i<NoInPrefix;i++)
    {
        if(!InputData[i].LoadSampleNames(InPrefixList[i].c_str()))
            return false;
        if(i>0)
            if(!InputData[i].CheckSampleConsistency(InputData[i-1].numSamples,
                                                    InputData[i-1].individualName,
                                                    InputData[i-1].SampleNoHaplotypes,
                                                    InputData[i-1].DoseFileName,
                                                    InputData[i].DoseFileName))
                return false;

    }

    cout<<" -- Successful !!! "<<endl;
    return true;
}


void MetaMinimac::OpenStreamInputDosageFiles(bool siteOnly)
{
    InputDosageStream.resize(NoInPrefix);
    CurrentRecordFromStudy.resize(NoInPrefix);
    StudiesHasVariant.resize(NoInPrefix);
    for(int i=0; i<NoInPrefix;i++)
    {
        VcfHeader header;
        InputDosageStream[i] = new VcfFileReader();
        CurrentRecordFromStudy[i]= new VcfRecord();
        InputDosageStream[i]->open( (GetDosageFileFullName(InPrefixList[i])).c_str() , header);
        InputDosageStream[i]->setSiteOnly(siteOnly);
        InputDosageStream[i]->readRecord(*CurrentRecordFromStudy[i]);
    }
}

void MetaMinimac::CloseStreamInputDosageFiles()
{
    for (int i = 0; i < NoInPrefix; i++)
    {
        delete InputDosageStream[i];
        delete CurrentRecordFromStudy[i];
    }
}


void MetaMinimac::OpenStreamInputEmpDosageFiles()
{
    InputEmpDosageStream.resize(NoInPrefix);
    CurrentEmpRecordFromStudy.resize(NoInPrefix);
    for(int i=0; i<NoInPrefix;i++)
    {
        VcfHeader header;
        InputEmpDosageStream[i] = new VcfFileReader();
        CurrentEmpRecordFromStudy[i]= new VcfRecord();
        InputEmpDosageStream[i]->open( (GetEmpDosageFileFullName(InPrefixList[i])).c_str() , header);
        InputEmpDosageStream[i]->readRecord(*CurrentEmpRecordFromStudy[i]);
    }
}

void MetaMinimac::CloseStreamInputEmpDosageFiles()
{
    for (int i = 0; i < NoInPrefix; i++)
    {
        delete InputEmpDosageStream[i];
        delete CurrentEmpRecordFromStudy[i];
    }
}


string MetaMinimac::GetDosageFileFullName(String prefix)
{
    if(doesExistFile(prefix+".dose.vcf"))
        return (string)prefix+".dose.vcf";
    else if(doesExistFile(prefix+".dose.vcf.gz"))
        return (string)prefix+".dose.vcf.gz";
    return "";
}

string MetaMinimac::GetEmpDosageFileFullName(String prefix)
{
    if(doesExistFile(prefix+".empiricalDose.vcf"))
        return (string)prefix+".empiricalDose.vcf";
    else if(doesExistFile(prefix+".empiricalDose.vcf.gz"))
        return (string)prefix+".empiricalDose.vcf.gz";
    return "";
}


bool MetaMinimac::doesExistFile(String filename)
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


bool MetaMinimac::LoadVariantInfo()
{
    do
    {
        FindCurrentMinimumPosition();
        if(CurrentFirstVariantBp==MAXBP)
            break;
        ReadCurrentVariantInfo();
        UpdateCurrentRecords();
    }while(true);

    NoVariants = VariantList.size();
    return true;

}

void MetaMinimac::FindCurrentMinimumPosition() {

    if (NoInPrefix == 2) {
        int a = CurrentRecordFromStudy[0]->get1BasedPosition();
        int b = CurrentRecordFromStudy[1]->get1BasedPosition();
        CurrentFirstVariantBp = a;
        NoStudiesHasVariant = 1;
        StudiesHasVariant[0] = 0;

        if (b == a) {
            if (IsVariantEqual(*CurrentRecordFromStudy[0], *CurrentRecordFromStudy[1]) == 1) {
                NoStudiesHasVariant = 2;
                StudiesHasVariant[1] = 1;
            }
        } else if (b < a) {
            StudiesHasVariant[0] = 1;
            CurrentFirstVariantBp = b;
        }

    }
}

int MetaMinimac::IsVariantEqual(VcfRecord &Rec1, VcfRecord &Rec2)
{
    if(strcmp(Rec1.getRefStr(),Rec2.getRefStr())!=0)
        return 0;
    if(strcmp(Rec1.getAltStr(),Rec2.getAltStr())!=0)
        return 0;
    return 1;
}


void MetaMinimac::ReadCurrentVariantInfo()
{
    variant tempVariant;
    tempVariant.bp = CurrentFirstVariantBp;
    tempVariant.NoStudiesHasVariant = NoStudiesHasVariant;
    tempVariant.StudiesHasVariant.resize(NoStudiesHasVariant);

    for(int i=0; i<NoStudiesHasVariant; i++)
    {
        tempVariant.StudiesHasVariant[i] = StudiesHasVariant[i];
    }

    VcfRecord *record = CurrentRecordFromStudy[StudiesHasVariant[0]];
    tempVariant.chr=record->getChromStr();
    tempVariant.name=record->getIDStr();
    tempVariant.altAlleleString = record->getAltStr();
    tempVariant.refAlleleString = record->getRefStr();
    if(record->getInfo().getString("IMPUTED") == NULL & NoStudiesHasVariant==NoInPrefix)
    {
        tempVariant.typed = true;
        CommonTypedVariantList.push_back(tempVariant);
        CommonGenotypeVariantNameList.push_back(tempVariant.name);
    }
    NoCommonTypedVariants = CommonGenotypeVariantNameList.size();
    VariantList.push_back(tempVariant);

}


void MetaMinimac::UpdateCurrentRecords()
{
    for(int i=0; i<NoStudiesHasVariant;i++)
    {
        int index = StudiesHasVariant[i];
        if(!InputDosageStream[index]->readRecord(*CurrentRecordFromStudy[index]))
            CurrentRecordFromStudy[index]->set1BasedPosition(MAXBP);
    }
}

bool MetaMinimac::LoadLooDosage()
{
    for(int i=0; i<NoInPrefix; i++)
        InputData[i].ReadBasedOnSortCommonGenotypeList(CommonGenotypeVariantNameList);
    return true;
}

int MetaMinimac::PerformFinalAnalysis()
{
    Weights.resize(InputData[0].numHaplotypes);
    for (int i = 0; i < InputData[0].numHaplotypes; i++)
    {
        Weights[i].resize(NoVariants);
    }

    for (int j = 0; j < InputData[0].numSamples; j++)
    {
        if (InputData[0].SampleNoHaplotypes[j] == 2)
        {
            GetMetaEstimate(2 * j);
            GetMetaEstimate(2 * j + 1);
        }
        else
            GetMetaEstimate(2 * j);
    }
    return 0;
}

void MetaMinimac::GetMetaEstimate(int Sample)
{


//    LogOddsModel ThisSampleAnalysis;
//    ThisSampleAnalysis.initialize(Sample, this);
//    vector<double> init(NoInPrefix-1, 0.0);
//    vector<double> MiniMizer = Simplex(ThisSampleAnalysis, init);
//
//    MarkovModel MM;
////    MM.initialize(Sample, this);
//    CurrentInitialProb.resize(NoInPrefix, 0.0);
//    logitTransform(MiniMizer, CurrentInitialProb);

}