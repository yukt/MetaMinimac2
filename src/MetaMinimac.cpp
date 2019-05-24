#include "MetaMinimac.h"
#include "MarkovModel.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include "simplex.h"

using BT::Simplex;

using namespace std;

const int MAXBP = 999999999;

int MetaMinimac::Analyze()
{
    if(!myUserVariables.CheckValidity()) return -1;


    cout<<" ------------------------------------------------------------------------------"<<endl;
    cout<<"                             INPUT VCF DOSAGE FILE                             "<<endl;
    cout<<" ------------------------------------------------------------------------------"<<endl;

    if (!ParseInputVCFFiles())
    {
        cout << "\n Program Exiting ... \n\n";
        return -1;
    }

    if (!CheckSampleNameCompatibility())
    {
        cout << "\n Program Exiting ... \n\n";
        return -1;
    }


    OpenStreamInputDosageFiles(true);
    if (!OpenStreamOutputDosageFiles())
    {
        cout <<" Please check your write permissions in the output directory\n OR maybe the output directory does NOT exist ...\n";
        cout << "\n Program Exiting ... \n\n";
        return -1;
    }

    LoadVariantInfo();
    CloseStreamInputDosageFiles();

//    LoadLooDosage();

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

    NoHaplotypes = InputData[0].numActualHaps;
    NoSamples = InputData[0].numSamples;
    if (myUserVariables.VcfBuffer > NoSamples)
        myUserVariables.VcfBuffer = NoSamples;

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
    finChromosome = CurrentRecordFromStudy[0]->getChromStr();
}

void MetaMinimac::CloseStreamInputDosageFiles()
{
    for (int i = 0; i < NoInPrefix; i++)
    {
        delete InputDosageStream[i];
        delete CurrentRecordFromStudy[i];
    }
}

bool MetaMinimac::OpenStreamOutputDosageFiles()
{
    vcfdosepartial = ifopen(myUserVariables.outfile + ".metaDose.vcf.gz", "wb", InputFile::BGZF);
    VcfPrintStringPointer = (char*)malloc(sizeof(char) * (myUserVariables.PrintBuffer));
    if(vcfdosepartial==NULL)
    {
        cout <<"\n\n ERROR !!! \n Could NOT create the following file : "<< myUserVariables.outfile + ".metaDose.vcf.gz" <<endl;
        return false;
    }
    ifprintf(vcfdosepartial,"##fileformat=VCFv4.1\n");
    time_t t = time(0);
    struct tm * now = localtime( & t );
    ifprintf(vcfdosepartial,"##filedate=%d.%d.%d\n",(now->tm_year + 1900),(now->tm_mon + 1) ,now->tm_mday);
    ifprintf(vcfdosepartial,"##source=MetaMinimac.v%s\n",VERSION);
    ifprintf(vcfdosepartial,"##contig=<ID=%s>\n", finChromosome.c_str());
    ifprintf(vcfdosepartial,"##INFO=<ID=AF,Number=1,Type=Float,Description=\"Estimated Alternate Allele Frequency\">\n");
    ifprintf(vcfdosepartial,"##INFO=<ID=MAF,Number=1,Type=Float,Description=\"Estimated Minor Allele Frequency\">\n");
    ifprintf(vcfdosepartial,"##INFO=<ID=R2,Number=1,Type=Float,Description=\"Estimated Imputation Accuracy (R-square)\">\n");
    ifprintf(vcfdosepartial,"##INFO=<ID=TRAINING,Number=0,Type=Flag,Description=\"Marker was used to train meta-imputation weights\">\n");

    if(myUserVariables.infoDetails)
    {
        ifprintf(vcfdosepartial,"##INFO=<ID=NST,Number=1,Type=Integer,Description=\"Number of studies marker was found during meta-imputation\">\n");
        for(int i=0; i<NoInPrefix; i++)
            ifprintf(vcfdosepartial,"##INFO=<ID=S%d,Number=0,Type=Flag,Description=\"Marker was present in Study %d\">\n",i+1,i+1);
    }
    if(myUserVariables.GT)
        ifprintf(vcfdosepartial,"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
    if(myUserVariables.DS)
        ifprintf(vcfdosepartial,"##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Estimated Alternate Allele Dosage : [P(0/1)+2*P(1/1)]\">\n");
    if(myUserVariables.HDS)
        ifprintf(vcfdosepartial,"##FORMAT=<ID=HDS,Number=2,Type=Float,Description=\"Estimated Haploid Alternate Allele Dosage \">\n");
    if(myUserVariables.GP)
        ifprintf(vcfdosepartial,"##FORMAT=<ID=GP,Number=3,Type=Float,Description=\"Estimated Posterior Probabilities for Genotypes 0/0, 0/1 and 1/1 \">\n");
    if(myUserVariables.SD)
        ifprintf(vcfdosepartial,"##FORMAT=<ID=SD,Number=1,Type=Float,Description=\"Variance of Posterior Genotype Probabilities\">\n");

    ifprintf(vcfdosepartial,"##metaMinimac_Command=%s\n", myUserVariables.CommandLine.c_str());

    ifprintf(vcfdosepartial,"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");

    for(int Id=0;Id<InputData[0].numSamples;Id++)
    {
        ifprintf(vcfdosepartial,"\t%s",InputData[0].individualName[Id].c_str());
    }
    ifprintf(vcfdosepartial,"\n");

    ifclose(vcfdosepartial);

    if(myUserVariables.debug)
    {
        metaWeight = ifopen(myUserVariables.outfile + ".metaWeights.vcf.gz", "wb", InputFile::BGZF);
        WeightPrintStringPointer = (char*)malloc(sizeof(char) * (myUserVariables.PrintBuffer));
        if(metaWeight==NULL)
        {
            cout <<"\n\n ERROR !!! \n Could NOT create the following file : "<< myUserVariables.outfile + ".metaWeights.vcf.gz" <<endl;
            return false;
        }
        ifprintf(metaWeight,"##fileformat=NA\n");
        time_t t = time(0);
        struct tm * now = localtime( & t );
        ifprintf(metaWeight,"##filedate=%d.%d.%d\n",(now->tm_year + 1900),(now->tm_mon + 1) ,now->tm_mday);
        ifprintf(metaWeight,"##source=MetaMinimac.v%s\n",VERSION);
        ifprintf(metaWeight,"##contig=<ID=%s>\n", finChromosome.c_str());
        ifprintf(metaWeight,"##INFO=<ID=AF,Number=1,Type=Float,Description=\"Estimated Alternate Allele Frequency\">\n");
        ifprintf(metaWeight,"##INFO=<ID=MAF,Number=1,Type=Float,Description=\"Estimated Minor Allele Frequency\">\n");
        ifprintf(metaWeight,"##metaMinimac_Command=%s\n",myUserVariables.CommandLine.c_str());

        ifprintf(metaWeight,"#SNP");

        for(int Id=0;Id<InputData[0].numSamples;Id++)
        {
            ifprintf(metaWeight,"\t%s",InputData[0].individualName[Id].c_str());
        }
        ifprintf(metaWeight,"\n");
        ifclose(metaWeight);

    }


    return true;
}


string MetaMinimac::GetDosageFileFullName(String prefix)
{
    if(doesExistFile(prefix+".dose.vcf"))
        return (string)prefix+".dose.vcf";
    else if(doesExistFile(prefix+".dose.vcf.gz"))
        return (string)prefix+".dose.vcf.gz";
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
    cout<<"\n Scanning input VCFs for SNPs ... "<<endl;
    VariantList.clear();
    NoOffsetThisBlock.clear();
    NoVariantsProcessed = 0;
    do
    {
        FindCurrentMinimumPosition();
        if(CurrentFirstVariantBp==MAXBP)
            break;
        ReadCurrentVariantInfo();
        UpdateCurrentRecords();
    }while(true);

    NoVariants = VariantList.size();
    NoCommonTypedVariants = CommonGenotypeVariantNameList.size();
    assert(NoVariantsProcessed==NoVariants);
    cout<<" -- Found " << NoSamples <<" samples on " << NoVariants <<" sites ("<< NoCommonTypedVariants <<" common typed)! "<<endl;
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

    else if (NoInPrefix==3)
    {

        CurrentFirstVariantBp=CurrentRecordFromStudy[0]->get1BasedPosition();

        for(int i=1;i<NoInPrefix;i++)
            if(CurrentRecordFromStudy[i]->get1BasedPosition() < CurrentFirstVariantBp)
                CurrentFirstVariantBp=CurrentRecordFromStudy[i]->get1BasedPosition();

        NoStudiesHasVariant=0;
        VcfRecord *minRecord;
        minRecord = new VcfRecord();
        int Begin=0;
        for(int i=0;i<NoInPrefix;i++)
        {
            if(CurrentRecordFromStudy[i]->get1BasedPosition() == CurrentFirstVariantBp)
            {
                if(Begin==0)
                {
                    Begin=1;
                    minRecord=CurrentRecordFromStudy[i];
                    StudiesHasVariant[NoStudiesHasVariant] = i;
                    NoStudiesHasVariant++;
                }
                else if (IsVariantEqual(*minRecord, *CurrentRecordFromStudy[i])==1)
                {
                    StudiesHasVariant[NoStudiesHasVariant] = i;
                    NoStudiesHasVariant++;
                }
            }
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
    NoVariantsProcessed++;
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
        NoOffsetThisBlock.push_back(NoVariantsProcessed);
    }
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

//void MetaMinimac::LoadLooDosage()
//{
//    for(int i=0; i<NoInPrefix; i++)
//        InputData[i].ReadBasedOnSortCommonGenotypeList(CommonGenotypeVariantNameList);
//}

void MetaMinimac::LoadLooDosage()
{
    for(int i=0; i<NoInPrefix; i++)
        InputData[i].ReadBasedOnSortCommonGenotypeList(CommonGenotypeVariantNameList, StartSamId, EndSamId);
}

int MetaMinimac::PerformFinalAnalysis()
{
    cout << endl;
    cout << " ------------------------------------------------------------------------------" << endl;
    cout << "                           META-IMPUTATION ANALYSIS                            " << endl;
    cout << " ------------------------------------------------------------------------------" << endl;


    int maxVcfSample = myUserVariables.VcfBuffer;

    StartSamId = 0;

    batchNo = 0;

    int start_time, time_tot;

    while(true)
    {
        batchNo++;
        EndSamId = StartSamId + (maxVcfSample) < NoSamples ? StartSamId + (maxVcfSample) : NoSamples;
        cout << "  Meta-Imputing Sample " << StartSamId + 1 << "-" << EndSamId << " [" << setprecision(1) << fixed << 100 * (float) EndSamId / NoSamples << "%] ..." << endl;

        start_time = time(0);

        LoadLooDosage();
        InitiateProbs();
        WalkLeft();

        OpenStreamInputDosageFiles(false);

        if(maxVcfSample<NoSamples)
        {
            HapDosageSum.resize(NoVariants, 0.0);
            HapDosageSumSq.resize(NoVariants, 0.0);
            OpenTempOutputFiles();
            FlushPartialVcf();
        }
        else
        {
            FlushAllVcf();
        }


        CloseStreamInputDosageFiles();

        time_tot = time(0) - start_time;
        cout << "      Successful (" << time_tot % 60 << " seconds) !!! " << endl;

        StartSamId = EndSamId;

        if (StartSamId >= NoSamples)
            break;
    }

    if(batchNo > 1)
    {
        AppendtoMainVcf();

        if(myUserVariables.debug)
        {
            AppendtoMainWeightsFile();
        }
    }

    return 0;
}


//void MetaMinimac::GetMetaEstimate(int Sample, int SampleInBatch)
//{
//    MarkovModel MM;
//    MM.initialize(Sample, this);
//    MM.walkLeft(Sample, this);
//    MM.walkRight(Sample, this, SampleInBatch);
//}
//
//void MetaMinimac::GetMetaEstimate(int SampleInBatch)
//{
//    InitLeftProb(SampleInBatch);
//}

void MetaMinimac::InitiateProbs()
{
    int NoSamplesThisBatch = EndSamId-StartSamId;
    LeftProb.clear();
    LeftProb.resize(NoCommonTypedVariants+1);
    for(int i=0; i<NoCommonTypedVariants+1; i++)
    {
        vector<vector<double>> &ThisLeftProb =  LeftProb[i];
        ThisLeftProb.resize(2*NoSamplesThisBatch);
        for(int j=0; j<2*NoSamplesThisBatch; j++)
            ThisLeftProb[j].resize(NoInPrefix);
    }
    PrevLeftProb.clear();
    PrevLeftProb.resize(2*NoSamplesThisBatch);
    PrevRightProb.clear();
    PrevRightProb.resize(2*NoSamplesThisBatch);


    vector<double> init(NoInPrefix-1, 0.0);

    for (int id=0; id<NoSamplesThisBatch; id++)
    {

        int SampleId = StartSamId + id;
        if (InputData[0].SampleNoHaplotypes[SampleId] == 2)
        {
            InitLeftProb(2*id);
            InitLeftProb(2*id+1);
        }
        else
            InitLeftProb(2*id);
    }
}

void logitTransform(vector<double> &From,vector<double> &To)
{

    double sum=1.0;
    int NoDimensions = (int)To.size();
    for(int i=0; i < (NoDimensions-1); i++) sum+=exp(From[i]);
    for(int i=0; i < (NoDimensions-1); i++)  To[i]=exp(From[i])/sum;
    To[NoDimensions-1]=1.0/sum;


    double checkSum=0.0;
    for(int i=0;i<To.size();i++)
        checkSum+=To[i];
    if(checkSum>1.0001)
        abort();
}


void MetaMinimac::InitLeftProb(int HapInBatch)
{
    vector<double> &InitProb = LeftProb[NoCommonTypedVariants][HapInBatch];
    PrevLeftProb[HapInBatch].resize(NoInPrefix);
    PrevRightProb[HapInBatch].resize(NoInPrefix, 1.0);

    LogOddsModel ThisSampleAnalysis;
    ThisSampleAnalysis.reinitialize(HapInBatch, this);
    vector<double> init(NoInPrefix-1, 0.0);
    vector<double> MiniMizer = Simplex(ThisSampleAnalysis, init);
    logitTransform(MiniMizer, InitProb);

    for(int j=0; j<NoInPrefix; j++)
    {
        InitProb[j]+=backgroundError;
        PrevLeftProb[HapInBatch][j] = InitProb[j];
    }
}

void MetaMinimac::WalkLeft()
{
    NoVariantsProcessed = 0;
    NoCommonVariantsProcessed = 0;


    if (VariantList[NoVariants-1].typed)
    {
        NoVariantsProcessed++;
        NoCommonVariantsProcessed++;
        ProcessTypedLeftProb();
    }

    while(NoVariantsProcessed<NoVariants)
    {
        NoVariantsProcessed++;
        WalkOneStepLeft();
        if (VariantList[NoVariants-NoVariantsProcessed].typed)
        {
            NoCommonVariantsProcessed++;
            ProcessTypedLeftProb();
        }
    }
    WalkOneStepLeft();

    assert(NoCommonVariantsProcessed==NoCommonTypedVariants);
}

void MetaMinimac::WalkOneStepLeft()
{
    for(int id=0; id<EndSamId-StartSamId; id++)
    {
        int SampleId = StartSamId + id;
        if (InputData[0].SampleNoHaplotypes[SampleId] == 2)
        {
            WalkOneStepLeft(2*id);
            WalkOneStepLeft(2*id+1);
        }
        else
            WalkOneStepLeft(2*id);
    }
}

void MetaMinimac::WalkOneStepLeft(int HapInBatch)
{
    vector<double> &PrevProb=PrevLeftProb[HapInBatch];
    double r = Recom*1.0/NoInPrefix, complement = 1-Recom;
    vector<double> ThisProb;
    ThisProb.resize(NoInPrefix, 0.0);

//    double sum = 0.0;
    for(int i=0; i<NoInPrefix; i++)
    {
        for(int j=0; j<NoInPrefix; j++)
            ThisProb[i] += PrevProb[j]*r;
        ThisProb[i] += PrevProb[i]*complement;
//        sum += ThisProb[i];
    }

    for(int i=0; i<NoInPrefix; i++)
        PrevProb[i] = ThisProb[i];
}

void MetaMinimac::ProcessTypedLeftProb()
{
    for(int id=0; id<EndSamId-StartSamId; id++)
    {
        int SampleId = StartSamId + id;
        if (InputData[0].SampleNoHaplotypes[SampleId] == 2)
        {
            ProcessTypedLeftProb(2*id);
            ProcessTypedLeftProb(2*id+1);
        }
        else
            ProcessTypedLeftProb(2*id);
    }
}

void MetaMinimac::ProcessTypedLeftProb(int HapInBatch)
{
    float ThisGT = InputData[0].TypedGT[HapInBatch][NoCommonTypedVariants-NoCommonVariantsProcessed];
    vector<double> &ThisPrevLeftProb = PrevLeftProb[HapInBatch];
    vector<double> &ThisLeftProb = LeftProb[NoCommonTypedVariants-NoCommonVariantsProcessed][HapInBatch];

    double sum=0.0;
    for(int i=0; i<NoInPrefix; i++)
    {
        float ThisLooDosage = InputData[i].LooDosage[HapInBatch][NoCommonTypedVariants-NoCommonVariantsProcessed];
        ThisPrevLeftProb[i] *= (ThisGT==1)?(ThisLooDosage+backgroundError):(1-ThisLooDosage+backgroundError);
        ThisLeftProb[i] = ThisPrevLeftProb[i];
        sum += ThisLeftProb[i];
    }

    while(sum < JumpThreshold)
    {
        sum = 0.0;
        for(int i=0; i<NoInPrefix; i++)
        {
            ThisLeftProb[i] *= JumpFix;
            ThisPrevLeftProb[i] *= JumpFix;
            sum += ThisLeftProb[i];
        }
    }

}

void MetaMinimac::OpenTempOutputFiles()
{
    cout << "      Saving in temporary VCF file ... " << endl;
    VcfPrintStringPointerLength=0;
    stringstream ss;
    ss << (batchNo);
    string PartialVcfFileName(myUserVariables.outfile);
    PartialVcfFileName += ".metaDose.part."+(string)(ss.str())+".vcf.gz";
    vcfdosepartial = ifopen(PartialVcfFileName.c_str(), "wb", InputFile::BGZF);

    if(myUserVariables.debug)
    {
        string PartialWeightFileName(myUserVariables.outfile);
        PartialWeightFileName += ".metaWeights.part."+(string)(ss.str())+".vcf.gz";
        vcfweightpartial = ifopen(PartialWeightFileName.c_str(), "wb", InputFile::BGZF);
        WeightPrintStringPointerLength = 0;
    }

}

void MetaMinimac::FlushPartialVcf()
{
    NoVariantsProcessed = 0;
    NoCommonVariantsProcessed = 0;
    for(int i=0; i<NoVariants; i++)
    {
        CurrentVariant = &VariantList[i];
        WalkOneStepRight();
        ReadCurrentDosageData();

        if(CurrentVariant->typed)
        {
            UpdateLeftProb();
            CreateMetaImputedData();
            CalculateStats();
            PrintMetaImputedData();
            if(myUserVariables.debug)
                PrintMetaWeight();
            UpdateRightProb();
            NoCommonVariantsProcessed++;
        }
        else
        {
            CreateMetaImputedData();
            CalculateStats();
            PrintMetaImputedData();
        }

        UpdateCurrentRecords();
        NoVariantsProcessed++;
    }
    if (VcfPrintStringPointerLength > 0)
    {
        ifprintf(vcfdosepartial,"%s",VcfPrintStringPointer);
        VcfPrintStringPointerLength = 0;
    }
    ifclose(vcfdosepartial);

    if(myUserVariables.debug)
    {
        if(WeightPrintStringPointerLength > 0)
            ifprintf(vcfweightpartial, "%s", WeightPrintStringPointer);
        ifclose(vcfweightpartial);
    }
}


void MetaMinimac::FlushAllVcf()
{
    VcfPrintStringPointerLength=0;
    vcfdosepartial = ifopen(myUserVariables.outfile + ".metaDose.vcf.gz", "a", InputFile::BGZF);
    if(myUserVariables.debug)
    {
        WeightPrintStringPointerLength=0;
        vcfweightpartial = ifopen(myUserVariables.outfile + ".metaWeights.vcf.gz", "a", InputFile::BGZF);
    }

    NoVariantsProcessed = 0;
    NoCommonVariantsProcessed = 0;
    for(int i=0; i<NoVariants; i++)
    {
        CurrentVariant = &VariantList[i];
        WalkOneStepRight();
        ReadCurrentDosageData();

        if(CurrentVariant->typed)
        {
            UpdateLeftProb();
            CreateMetaImputedData();
            PrintVariantInfo();
            PrintMetaImputedData();
            if(myUserVariables.debug)
            {
                PrintWeightVariantInfo();
                PrintMetaWeight();
            }
            UpdateRightProb();
            NoCommonVariantsProcessed++;
        }
        else
        {
            CreateMetaImputedData();
            PrintVariantInfo();
            PrintMetaImputedData();
        }

        UpdateCurrentRecords();
        NoVariantsProcessed++;
    }
    if (VcfPrintStringPointerLength > 0)
    {
        ifprintf(vcfdosepartial,"%s",VcfPrintStringPointer);
        VcfPrintStringPointerLength = 0;
    }
    ifclose(vcfdosepartial);

    if(myUserVariables.debug)
    {
        if(WeightPrintStringPointerLength > 0)
            ifprintf(vcfweightpartial, "%s", WeightPrintStringPointer);
        ifclose(vcfweightpartial);
    }
}


void MetaMinimac::WalkOneStepRight()
{
    for(int id=0; id<EndSamId-StartSamId; id++)
    {
        int SampleId = StartSamId + id;
        if (InputData[0].SampleNoHaplotypes[SampleId] == 2)
        {
            WalkOneStepRight(2*id);
            WalkOneStepRight(2*id+1);
        }
        else
            WalkOneStepRight(2*id);
    }
}

void MetaMinimac::WalkOneStepRight(int HapInBatch)
{
    vector<double> &PrevLeft=PrevLeftProb[HapInBatch];
    vector<double> &PrevRight=PrevRightProb[HapInBatch];
    vector<double> ThisLeft, ThisRight;

    double r = Recom*1.0/NoInPrefix, complement = 1-Recom;
    double r_c = r/complement;
    ThisLeft.resize(NoInPrefix);
    ThisRight.resize(NoInPrefix);

    for(int i=0; i<NoInPrefix; i++)
    {
        ThisLeft[i] = PrevLeft[i]*(1+2*r_c);
        ThisRight[i] = PrevRight[i]*complement;
        for(int j=0; j<NoInPrefix; j++)
        {
            ThisRight[i] += PrevRight[j]*r;
            ThisLeft[i] -= PrevLeft[j]*r_c;
        }
    }



    for(int i=0; i<NoInPrefix; i++)
    {
        PrevRight[i] = ThisRight[i];
        PrevLeft[i] = ThisLeft[i];
    }
}

void MetaMinimac::UpdateLeftProb()
{
    vector<vector<double>> &ThisLeftProb = LeftProb[NoCommonVariantsProcessed];
    for(int id=0; id<EndSamId-StartSamId; id++)
    {
        if(InputData[0].SampleNoHaplotypes[StartSamId+id]==2)
        {
            for(int i=0; i<NoInPrefix; i++)
            {
                PrevLeftProb[2*id][i] = ThisLeftProb[2*id][i];
                PrevLeftProb[2*id+1][i] = ThisLeftProb[2*id+1][i];
            }
        }
        else
        {
            for(int i=0; i<NoInPrefix; i++)
            {
                PrevLeftProb[2*id][i] = ThisLeftProb[2*id][i];
            }
        }
    }
}

void MetaMinimac::UpdateRightProb()
{
    for(int id=0; id<EndSamId-StartSamId; id++)
    {
        if(InputData[0].SampleNoHaplotypes[StartSamId+id]==2)
        {
            UpdateRightProb(2*id);
            UpdateRightProb(2*id+1);
        }
        else
            UpdateRightProb(2*id);
    }
}

void MetaMinimac::UpdateRightProb(int HapInBatch)
{
    float ThisGT = InputData[0].TypedGT[HapInBatch][NoCommonVariantsProcessed];
    vector<double> &ThisRightProb = PrevRightProb[HapInBatch];
    vector<double> &ThisLeftProb  = PrevLeftProb[HapInBatch];

    double sum=0.0;
    for(int i=0; i<NoInPrefix; i++)
    {
        float ThisLooDosage = InputData[i].LooDosage[HapInBatch][NoCommonVariantsProcessed];
        ThisRightProb[i] *= (ThisGT==1)?(ThisLooDosage+backgroundError):(1-ThisLooDosage+backgroundError);
        ThisLeftProb[i] /= (ThisGT==1)?(ThisLooDosage+backgroundError):(1-ThisLooDosage+backgroundError);
        sum += ThisRightProb[i];
    }

    while(sum < JumpThreshold)
    {
        sum = 0.0;
        for(int i=0; i<NoInPrefix; i++)
        {
            ThisRightProb[i] *= JumpFix;
            sum += ThisRightProb[i];
        }
    }
}

void MetaMinimac::ReadCurrentDosageData()
{
    NoStudiesHasVariant = CurrentVariant->NoStudiesHasVariant;
    for(int j=0; j<NoStudiesHasVariant; j++)
    {
        int index = CurrentVariant->StudiesHasVariant[j];
        InputData[index].LoadHapDoseVariant(CurrentRecordFromStudy[index]->getGenotypeInfo(), StartSamId, EndSamId);
    }
}

void MetaMinimac::CreateMetaImputedData()
{
    CurrentHapDosageSum = 0;
    CurrentHapDosageSumSq = 0;
    CurrentMetaImputedDosage.clear();
    CurrentPosterior.clear();

    if(NoStudiesHasVariant==1)
    {
        CurrentMetaImputedDosage = InputData[CurrentVariant->StudiesHasVariant[0]].CurrentHapDosage;
        for(int i=0; i<2*(EndSamId-StartSamId); i++)
        {
            CurrentHapDosageSum += CurrentMetaImputedDosage[i];
            CurrentHapDosageSumSq += CurrentMetaImputedDosage[i]*CurrentMetaImputedDosage[i];
        }
    }
    else
    {
        CurrentMetaImputedDosage.resize(2*(EndSamId-StartSamId));
        CurrentPosterior.resize(2*(EndSamId-StartSamId));
        for(int id=0; id<EndSamId-StartSamId; id++)
        {
            int SampleId = StartSamId+id;
            if(InputData[0].SampleNoHaplotypes[SampleId]==2)
            {
                MetaImpute(2*id);
                MetaImpute(2*id + 1);
            }
            else
            {
                MetaImpute(2*id);
            }
        }
    }
}

void MetaMinimac::CalculateStats()
{
    HapDosageSum[NoVariantsProcessed] += CurrentHapDosageSum;
    HapDosageSumSq[NoVariantsProcessed] += CurrentHapDosageSumSq;
}

void MetaMinimac::MetaImpute(int Sample)
{
    double WeightSum = 0.0;
    double Dosage = 0.0;

    vector<double> &ThisLeft = PrevLeftProb[Sample];
    vector<double> &ThisRight = PrevRightProb[Sample];
    vector<float> &ThisPosterior = CurrentPosterior[Sample];
    ThisPosterior.resize(NoInPrefix);

    for (int j=0; j<NoStudiesHasVariant; j++)
    {
        int index = CurrentVariant->StudiesHasVariant[j];
        double Weight = ThisLeft[index]*ThisRight[index];
        ThisPosterior[index] = Weight;
        WeightSum += Weight;
        Dosage += Weight * InputData[index].CurrentHapDosage[Sample];
    }
    Dosage /= WeightSum;

    CurrentMetaImputedDosage[Sample] = Dosage;
    CurrentHapDosageSum += Dosage;
    CurrentHapDosageSumSq += Dosage * Dosage;
}


void MetaMinimac::PrintMetaImputedData()
{
    for(int id=0; id<EndSamId-StartSamId; id++)
    {
        if( InputData[0].SampleNoHaplotypes[StartSamId+id]==2 )
            PrintDiploidDosage((CurrentMetaImputedDosage[2*id]), (CurrentMetaImputedDosage[2*id+1]));
        else
            PrintHaploidDosage((CurrentMetaImputedDosage[2*id]));
    }

    VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,"\n");
    if(VcfPrintStringPointerLength > 0.9 * (float)(myUserVariables.PrintBuffer))
    {
        ifprintf(vcfdosepartial,"%s",VcfPrintStringPointer);
        VcfPrintStringPointerLength=0;
    }

}

void MetaMinimac::PrintMetaWeight()
{
    for(int id=0; id<EndSamId-StartSamId; id++)
    {
        WeightPrintStringPointerLength += sprintf(WeightPrintStringPointer+WeightPrintStringPointerLength,"\t");
        if(InputData[0].SampleNoHaplotypes[StartSamId+id]==2)
        {
            PrintWeightForHaplotype(2*id);
            WeightPrintStringPointerLength += sprintf(WeightPrintStringPointer+WeightPrintStringPointerLength,"|");
            PrintWeightForHaplotype(2*id+1);
        }
        else
            PrintWeightForHaplotype(2*id);
    }

    WeightPrintStringPointerLength+= sprintf(WeightPrintStringPointer+WeightPrintStringPointerLength,"\n");
    if(WeightPrintStringPointerLength > 0.9 * (float)(myUserVariables.PrintBuffer))
    {
        ifprintf(vcfweightpartial, "%s", WeightPrintStringPointer);
        WeightPrintStringPointerLength = 0;
    }
}


void MetaMinimac::PrintDiploidDosage(float &x, float &y)
{

    bool colonIndex=false;
    VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,"\t");

    if(x<0.0005 && y<0.0005)
    {
        if(myUserVariables.GT)
        {

            VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,"0|0");
            colonIndex=true;
        }
        if(myUserVariables.DS)
        {

            if(colonIndex)
                VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,":");
            VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,"0");
            colonIndex=true;
        }
        if(myUserVariables.HDS)
        {

            if(colonIndex)
                VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,":");
            VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,"0,0");
            colonIndex=true;
        }
        if(myUserVariables.GP)
        {

            if(colonIndex)
                VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,":");
            colonIndex=true;
            VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,"1,0,0");
        }
        if(myUserVariables.SD)
        {
            if(colonIndex)
                VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,":");
            colonIndex=true;
            VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,"0");
        }
        return;
    }


    if(myUserVariables.GT)
    {

        VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,"%d|%d",(x>0.5),(y>0.5));
        colonIndex=true;
    }
    if(myUserVariables.DS)
    {

        if(colonIndex)
            VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,":");
        VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,"%.3f",x+ y);
        colonIndex=true;
    }
    if(myUserVariables.HDS)
    {

        if(colonIndex)
            VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,":");
        VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,"%.3f,%.3f",x , y);
        colonIndex=true;
    }
    if(myUserVariables.GP)
    {

        if(colonIndex)
            VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,":");
        colonIndex=true;
        VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,"%.3f,%.3f,%.3f",(1-x)*(1-y),x*(1-y)+y*(1-x),x*y);
    }
    if(myUserVariables.SD)
    {
        if(colonIndex)
            VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,":");
        colonIndex=true;
        VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,"%.3f", x*(1-x) + y*(1-y));
    }



}

void MetaMinimac::PrintHaploidDosage(float &x)
{
    bool colonIndex=false;
    VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,"\t");

    if(x<0.0005)
    {
        if(myUserVariables.GT)
        {
            VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,"0");
            colonIndex=true;
        }
        if(myUserVariables.DS)
        {

            if(colonIndex)
                VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,":");
            VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,"0");
            colonIndex=true;
        }
        if(myUserVariables.HDS)
        {

            if(colonIndex)
                VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,":");
            VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,"0");
            colonIndex=true;
        }
        if(myUserVariables.GP)
        {

            if(colonIndex)
                VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,":");
            colonIndex=true;
            VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,"1,0");
        }
        if(myUserVariables.SD)
        {
            if(colonIndex)
                VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,":");
            colonIndex=true;
            VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,"0");
        }
        return;

    }



    if(myUserVariables.GT)
    {
        VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,"%d",(x>0.5));
        colonIndex=true;
    }
    if(myUserVariables.DS)
    {

        if(colonIndex)
            VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,":");
        VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,"%.3f",x);
        colonIndex=true;
    }
    if(myUserVariables.HDS)
    {

        if(colonIndex)
            VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,":");
        VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,"%.3f",x );
        colonIndex=true;
    }
    if(myUserVariables.GP)
    {

        if(colonIndex)
            VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,":");
        colonIndex=true;
        VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,"%.3f,%.3f",1-x,x);
    }
    if(myUserVariables.SD)
    {
        if(colonIndex)
            VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,":");
        colonIndex=true;
        VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,"%.3f", x*(1-x));
    }
}


void MetaMinimac::PrintWeightForHaplotype(int haploId)
{
    vector<float> &ThisPosterior = CurrentPosterior[haploId];
    float WeightSum = 0.0;
    for(int i=0; i<NoInPrefix; i++)
        WeightSum += ThisPosterior[i];
    WeightPrintStringPointerLength += sprintf(WeightPrintStringPointer+WeightPrintStringPointerLength,"%0.4f", ThisPosterior[0]/WeightSum);
    for(int i=1;i<NoInPrefix;i++)
        WeightPrintStringPointerLength += sprintf(WeightPrintStringPointer+WeightPrintStringPointerLength,",%0.4f", ThisPosterior[i]/WeightSum);
}



void MetaMinimac::AppendtoMainVcf()
{
    cout << endl;
    cout << "  Appending to final output VCF File : " << myUserVariables.outfile + ".metaDose.vcf.gz" <<endl;

    int start_time = time(0);
    VcfPrintStringPointerLength=0;
    vcfdosepartial = ifopen(myUserVariables.outfile + ".metaDose.vcf.gz", "a", InputFile::BGZF);
    vector<IFILE> vcfdosepartialList(batchNo);

    for(int i=1;i<=batchNo;i++)
    {
        stringstream ss;
        ss << (i);
        string PartialVcfFileName(myUserVariables.outfile);
        PartialVcfFileName += ".metaDose.part."+(string)(ss.str())+".vcf.gz";
        vcfdosepartialList[i-1] = ifopen(PartialVcfFileName.c_str(), "r");
    }

    string line;

    for(int i=0; i<NoVariants; i++)
    {
        CurrentVariant = &VariantList[i];
        VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength, "%s\t%d\t%s\t%s\t%s\t.\tPASS\t%s\t%s",
                                             CurrentVariant->chr.c_str(), CurrentVariant->bp, CurrentVariant->name.c_str(),
                                             CurrentVariant->refAlleleString.c_str(), CurrentVariant->altAlleleString.c_str(),
                                             CreateInfo(i).c_str(), myUserVariables.formatStringForVCF.c_str());
        for(int j=1;j<=batchNo;j++)
        {
            line.clear();
            vcfdosepartialList[j-1]->readLine(line);
            VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer + VcfPrintStringPointerLength,"%s",line.c_str());

        }
        VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer + VcfPrintStringPointerLength,"\n");
        if(VcfPrintStringPointerLength > 0.9 * (float)(myUserVariables.PrintBuffer))
        {
            ifprintf(vcfdosepartial,"%s",VcfPrintStringPointer);
            VcfPrintStringPointerLength=0;
        }
    }
    if(VcfPrintStringPointerLength > 0)
    {
        ifprintf(vcfdosepartial,"%s",VcfPrintStringPointer);
        VcfPrintStringPointerLength=0;
    }

    for(int i=1;i<=batchNo;i++)
    {
        ifclose(vcfdosepartialList[i-1]);
        stringstream ss;
        ss << (i);
        string tempFileIndex(myUserVariables.outfile);
        tempFileIndex += ".metaDose.part."+(string)(ss.str())+".vcf.gz";
        remove(tempFileIndex.c_str());
    }
    ifclose(vcfdosepartial);

    int time_tot = time(0) - start_time;
    cout << "  Appending successful (" << time_tot % 60 << " seconds) !!!" << endl;
}

void MetaMinimac::PrintVariantInfo()
{
    VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength, "%s\t%d\t%s\t%s\t%s\t.\tPASS\t%s\t%s",
                                         CurrentVariant->chr.c_str(), CurrentVariant->bp, CurrentVariant->name.c_str(),
                                         CurrentVariant->refAlleleString.c_str(), CurrentVariant->altAlleleString.c_str(),
                                         CreateInfo().c_str(), myUserVariables.formatStringForVCF.c_str());
}

string MetaMinimac::CreateInfo()
{
    double hapSum = CurrentHapDosageSum, hapSumSq = CurrentHapDosageSumSq;
    double freq = hapSum*1.0/NoHaplotypes;
    double maf = (freq > 0.5) ? (1.0 - freq) : freq;
    double rsq = 0.0, evar = freq*(1-freq), ovar = 0.0;
    if (NoHaplotypes > 2 & (hapSumSq - hapSum * hapSum / NoHaplotypes) >0 )
    {
        ovar = (hapSumSq - hapSum * hapSum / NoHaplotypes)/ NoHaplotypes;
        rsq = ovar / (evar + 1e-30);
    }

    stringstream ss;
    ss<< "AF=" << fixed << setprecision(5) << freq <<";MAF=";
    ss<< fixed << setprecision(5) << maf <<";R2=";
    ss<< fixed << setprecision(5) << rsq ;
    ss<<";";

    if(!myUserVariables.infoDetails)
        return ".";
    ss<<"NST=";
    ss<<CurrentVariant->NoStudiesHasVariant;
    for(int i=0; i<CurrentVariant->NoStudiesHasVariant; i++)
    {
        ss<<";S";
        ss<<CurrentVariant->StudiesHasVariant[i]+1;
    }

    if(CurrentVariant->typed & (CurrentVariant->NoStudiesHasVariant == NoInPrefix))
    {
        ss<<";TRAINING";
    }

    return ss.str();
}


string MetaMinimac::CreateInfo(int id)
{
    double hapSum = HapDosageSum[id], hapSumSq = HapDosageSumSq[id];
    double freq = hapSum*1.0/NoHaplotypes;
    double maf = (freq > 0.5) ? (1.0 - freq) : freq;
    double rsq = 0.0, evar = freq*(1-freq), ovar = 0.0;
    if (NoHaplotypes > 2 & (hapSumSq - hapSum * hapSum / NoHaplotypes) >0 )
    {
        ovar = (hapSumSq - hapSum * hapSum / NoHaplotypes)/ NoHaplotypes;
        rsq = ovar / (evar + 1e-30);
    }

    stringstream ss;
    ss<< "AF=" << fixed << setprecision(5) << freq <<";MAF=";
    ss<< fixed << setprecision(5) << maf <<";R2=";
    ss<< fixed << setprecision(5) << rsq ;
    ss<<";";

    if(!myUserVariables.infoDetails)
        return ".";
    ss<<"NST=";
    ss<<CurrentVariant->NoStudiesHasVariant;
    for(int i=0; i<CurrentVariant->NoStudiesHasVariant; i++)
    {
        ss<<";S";
        ss<<CurrentVariant->StudiesHasVariant[i]+1;
    }

    if(CurrentVariant->typed & (CurrentVariant->NoStudiesHasVariant == NoInPrefix))
    {
        ss<<";TRAINING";
    }

    return ss.str();
}

void MetaMinimac::PrintWeightVariantInfo()
{
    WeightPrintStringPointerLength+=sprintf(WeightPrintStringPointer+WeightPrintStringPointerLength,"%s:%d:%s:%s",
                                            CurrentVariant->chr.c_str(), CurrentVariant->bp,
                                            CurrentVariant->refAlleleString.c_str(), CurrentVariant->altAlleleString.c_str());
}

void MetaMinimac::AppendtoMainWeightsFile()
{
    WeightPrintStringPointerLength=0;
    metaWeight = ifopen(myUserVariables.outfile + ".metaWeights.vcf.gz", "a", InputFile::BGZF);
    vector<IFILE> weightpartialList(batchNo);

    for(int i=1; i<=batchNo; i++)
    {
        stringstream ss;
        ss << (i);
        string PartialWeightFileName(myUserVariables.outfile);
        PartialWeightFileName += ".metaWeights.part."+(string)(ss.str())+".vcf.gz";
        weightpartialList[i-1] = ifopen(PartialWeightFileName.c_str(), "r");
    }

    string line;

    for(int i=0; i<NoCommonTypedVariants; i++)
    {
        CurrentVariant = &CommonTypedVariantList[i];
        PrintWeightVariantInfo();
        for(int j=1;j<=batchNo;j++)
        {
            line.clear();
            weightpartialList[j-1]->readLine(line);
            WeightPrintStringPointerLength+=sprintf(WeightPrintStringPointer + WeightPrintStringPointerLength,"%s",line.c_str());
        }
        WeightPrintStringPointerLength+=sprintf(WeightPrintStringPointer + WeightPrintStringPointerLength,"\n");
        if(WeightPrintStringPointerLength > 0.9 * (float)(myUserVariables.PrintBuffer))
        {
            ifprintf(metaWeight, "%s", WeightPrintStringPointer);
            WeightPrintStringPointerLength = 0;
        }
    }
    if(WeightPrintStringPointerLength > 0)
    {
        ifprintf(metaWeight, "%s", WeightPrintStringPointer);
        WeightPrintStringPointerLength = 0;
    }

    for(int i=1;i<=batchNo;i++)
    {
        ifclose(weightpartialList[i-1]);
        stringstream ss;
        ss << (i);
        string tempFileIndex(myUserVariables.outfile);
        tempFileIndex += ".metaWeights.part."+(string)(ss.str())+".vcf.gz";
        remove(tempFileIndex.c_str());
    }
    ifclose(metaWeight);
}
