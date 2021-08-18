#include "MetaMinimac.h"
#include "MarkovModel.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include "simplex.h"
#define RECOM_MIN 1e-04

using BT::Simplex;

using namespace std;

const int MAXBP = 999999999;

String MetaMinimac::Analyze()
{
    if(!myUserVariables.CheckValidity()) return "Command.Line.Error";


    cout<<" ------------------------------------------------------------------------------"<<endl;
    cout<<"                             INPUT VCF DOSAGE FILE                             "<<endl;
    cout<<" ------------------------------------------------------------------------------"<<endl;

    if (!ParseInputVCFFiles())
    {
        cout << "\n Program Exiting ... \n\n";
        return "Command.Line.Error";
    }

    if (!CheckSampleNameCompatibility())
    {
        cout << "\n Program Exiting ... \n\n";
        return "Input.VCF.Dose.Error";
    }

    if (!LoadEmpVariantInfo())
    {
        cout << "\n Program Exiting ... \n\n";
        return "Input.VCF.Dose.Error";
    }

    if(!CheckPhasingConsistency())
    {
        cout << "\n Program Exiting ... \n\n";
        return "Phasing.Inconsistency.Error";
    }

    if (!OpenStreamOutputDosageFiles())
    {
        cout <<" Please check your write permissions in the output directory\n OR maybe the output directory does NOT exist ...\n";
        cout << "\n Program Exiting ... \n\n";
        return "File.Write.Error";
    }

    if (!PerformWeightEstimation())
    {
        cout <<" Please check your write permissions in the output directory\n OR maybe the output directory does NOT exist ...\n";
        cout << "\n Program Exiting ... \n\n";
        return "Weight.Estimation.Error";
    }

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

    cout<<endl<<  " Number of Studies : "<<NoInPrefix<<endl;
    for(int i=0;i<NoInPrefix;i++)
    {
        cout<<  " -- Study "<<i+1<<" Prefix : "<<InPrefixList[i]<<endl;
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
    cout<<"\n Checking sample compatibility across files ... "<<endl;
    for(int i=0;i<NoInPrefix;i++)
    {
        // First, check if dose file and empiricalDose file exist
        // Second, check if sample names are consistent between dose file and empiricalDose file
        if(!InputData[i].LoadSampleNames(InPrefixList[i].c_str()))
            return false;

        // Check if sample names are consistent with InputData[0]
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
    cout<<" -- Found "<< NoSamples << " samples (" << NoHaplotypes << " haplotypes)." <<endl;

    if (myUserVariables.VcfBuffer > NoSamples)
        myUserVariables.VcfBuffer = NoSamples;
    return true;
}


void MetaMinimac::OpenStreamInputDosageFiles()
{
    InputDosageStream.clear();
    CurrentRecordFromStudy.clear();
    StudiesHasVariant.clear();

    InputDosageStream.resize(NoInPrefix);
    CurrentRecordFromStudy.resize(NoInPrefix);
    StudiesHasVariant.resize(NoInPrefix);
    for(int i=0; i<NoInPrefix;i++)
    {
        VcfHeader header;
        InputDosageStream[i] = new VcfFileReader();
        CurrentRecordFromStudy[i]= new VcfRecord();
        InputDosageStream[i]->open( InputData[i].DoseFileName.c_str(), header);
        InputDosageStream[i]->setSiteOnly(false);
        InputDosageStream[i]->readRecord(*CurrentRecordFromStudy[i]);
        InputData[i].noMarkers = 0;
    }
}

void MetaMinimac::OpenStreamInputEmpiricalDoseFiles()
{
    InputDosageStream.clear();
    CurrentRecordFromStudy.clear();

    InputDosageStream.resize(NoInPrefix);
    CurrentRecordFromStudy.resize(NoInPrefix);

    for(int i=0; i<NoInPrefix;i++)
    {
        VcfHeader header;
        InputDosageStream[i] = new VcfFileReader();
        CurrentRecordFromStudy[i]= new VcfRecord();
        InputDosageStream[i]->open( InputData[i].EmpDoseFileName.c_str(), header);
        InputDosageStream[i]->setSiteOnly(false);
    }
}


void MetaMinimac::OpenStreamInputWeightFiles()
{
    VcfHeader header;
    InputWeightStream = new VcfFileReader();
    CurrentRecordFromWeight = new VcfRecord();
    InputWeightStream->open(myUserVariables.outfile + ".metaWeights"+(myUserVariables.gzip ? ".gz" : ""), header);
    InputWeightStream->setSiteOnly(false);
}


void MetaMinimac::CloseStreamInputDosageFiles()
{
    for (int i = 0; i < NoInPrefix; i++)
    {
        delete InputDosageStream[i];
        delete CurrentRecordFromStudy[i];
    }
}

void MetaMinimac::CloseStreamInputWeightFiles()
{
    delete InputWeightStream;
    delete CurrentRecordFromWeight;
    if(!myUserVariables.debug)
    {
        string InputWeightFileName(myUserVariables.outfile + ".metaWeights"+(myUserVariables.gzip ? ".gz" : ""));
        remove(InputWeightFileName.c_str());
    }
}

bool MetaMinimac::OpenStreamOutputDosageFiles()
{
    vcfdosepartial = ifopen(myUserVariables.outfile + ".metaDose.vcf"+(myUserVariables.gzip ? ".gz" : ""), "wb", myUserVariables.gzip ? InputFile::BGZF : InputFile::UNCOMPRESSED);
    VcfPrintStringPointer = (char*)malloc(sizeof(char) * (myUserVariables.PrintBuffer));
    VcfPrintStringPointerLength = 0;
    if(vcfdosepartial==NULL)
    {
        cout <<"\n\n ERROR !!! \n Could NOT create the following file : "<< myUserVariables.outfile + ".metaDose.vcf"+(myUserVariables.gzip ? ".gz" : "") <<endl;
        return false;
    }
    ifprintf(vcfdosepartial,"##fileformat=VCFv4.1\n");
    time_t t = time(0);
    struct tm * now = localtime( & t );
    ifprintf(vcfdosepartial,"##filedate=%d.%d.%d\n",(now->tm_year + 1900),(now->tm_mon + 1) ,now->tm_mday);
    ifprintf(vcfdosepartial,"##source=MetaMinimac2.v%s\n",VERSION);
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

    return true;
}

bool MetaMinimac::LoadEmpVariantInfo()
{
    cout<<"\n Scanning input empirical VCFs for commonly typed SNPs ... "<<endl;
    int time_start = time(0);
    for(int i=0;i<NoInPrefix;i++)
    {
        InputData[i].LoadEmpVariantList();
        cout<<" -- Study "<<i+1<<" #Genotyped Sites = "<<InputData[i].noTypedMarkers<<endl;
    }

    finChromosome = InputData[0].finChromosome;
    FindCommonGenotypedVariants();

    if (NoCommonTypedVariants < 1)
    {
        cout << " -- No commonly genotyped sites found!!!" << endl;
        return false;
    }

    cout<<" -- Found " << NoCommonTypedVariants <<" commonly genotyped! "<<endl;
    cout<<" -- Successful (" << (time(0)-time_start) << " seconds) !!!" << endl;
    return true;

}

void MetaMinimac::FindCommonGenotypedVariants()
{
    std::map<string, int > HashUnionVariantMap;

    for(int i=1;i<NoInPrefix;i++)
    {
        for(int j=0;j<InputData[i].noTypedMarkers;j++)
        {
            variant *thisVariant=&InputData[i].TypedVariantList[j];
            HashUnionVariantMap[thisVariant->name]++;
        }
    }

    TransitionProb.clear();
    int lastbp = 0;
    for(int j=0; j<InputData[0].noTypedMarkers; j++)
    {
        variant *thisVariant = &InputData[0].TypedVariantList[j];
        if(HashUnionVariantMap[thisVariant->name]==NoInPrefix-1)
        {
            CommonGenotypeVariantNameList.push_back(thisVariant->name);
            variant tempVariant;
            tempVariant.chr=thisVariant->chr;
            tempVariant.bp=thisVariant->bp;
            tempVariant.name=thisVariant->name;
            tempVariant.altAlleleString = thisVariant->altAlleleString;
            tempVariant.refAlleleString = thisVariant->refAlleleString;
            CommonTypedVariantList.push_back(tempVariant);
            int distance = tempVariant.bp - lastbp;
            double prob  = max(RECOM_MIN, 1-exp(-lambda*distance));
            TransitionProb.push_back(prob);
            lastbp = tempVariant.bp;
        }
    }

    NoCommonTypedVariants = CommonGenotypeVariantNameList.size();

    for(int i=0;i<NoInPrefix;i++)
    {
        InputData[i].ClearEmpVariantList();
    }
}

bool MetaMinimac::CheckPhasingConsistency()
{
    if (!myUserVariables.hapcheck)
    {
        cout << "\n Warning !!! Phasing Consistency Check is Skipped." << endl;
        cout << "\n User is responsible for making sure the phasing are consistent across input files.";
        cout << "\n If phasing is inconsistent, meta-imputed results will not make sense." << endl;
        return true;
    }

    int start_time = time(0);

    cout<<"\n Checking phasing consistency across input files ... "<<endl;
    OpenStreamInputEmpiricalDoseFiles();

    string chrom, altAllele, refAllele, name;
    int bp;

    for (int k=0; k<NoCommonTypedVariants; k++)
    {
        string common_variant = CommonGenotypeVariantNameList[k];

        NoStudiesHasVariant = 0;

        for(int i=0; i<NoInPrefix; i++)
        {
            while(InputDosageStream[i]->readRecord(*CurrentRecordFromStudy[i])) {
                chrom = CurrentRecordFromStudy[i]->getChromStr();
                bp = CurrentRecordFromStudy[i]->get1BasedPosition();
                altAllele = CurrentRecordFromStudy[i]->getAltStr();
                refAllele = CurrentRecordFromStudy[i]->getRefStr();
                name = chrom + ":" + to_string(bp) + ":" + refAllele + ":" + altAllele;
                if (name == common_variant)
                {
                    InputData[i].LoadCurrentGT(CurrentRecordFromStudy[i]->getGenotypeInfo());
                    NoStudiesHasVariant++;
                    break;
                }
            }
        }

        if (NoStudiesHasVariant < NoInPrefix)
        {
            cout<<" Did not found SNP " << common_variant << " in all the files" << endl;
            cout<<" ERROR CODE 2819: Please contact author with this code to help with bug fixing ..."<<endl;
            return false;
        }

        bool flag_consistent = true;
        for (int j=0; j<NoSamples; j++)
        {
            if (InputData[0].SampleNoHaplotypes[j] == 2)
            {
                for(int i=1; i<NoInPrefix; i++)
                {
                    if (InputData[i].CurrentHapDosage[2*j] != InputData[0].CurrentHapDosage[2*j])
                        flag_consistent = false;
                    if (InputData[i].CurrentHapDosage[2*j+1] != InputData[0].CurrentHapDosage[2*j+1])
                        flag_consistent = false;
                }
            }
            else
            {
                for(int i=1; i<NoInPrefix; i++)
                {
                    if (InputData[i].CurrentHapDosage[2*j] != InputData[0].CurrentHapDosage[2*j])
                        flag_consistent = false;
                }
            }
        }
        if(!flag_consistent)
        {
            cout << "Phasing is not consistent at SNP " << common_variant << "." << endl;
            cout << "\n Please note that phasing consistency is important, ";
            cout << "\n because meta-imputation is performed on individual haplotype.";
            cout << "\n If the phasing is different between input imputed files,";
            cout << "\n then the resulting meta dosages will not make any sense." << endl;
            return false;
        }

    }

    cout<<" -- Consistent among commonly genotyped sites." <<endl;
    cout<<" -- Completed in " << time(0)-start_time << " seconds." << endl;
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

    else
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



void MetaMinimac::UpdateCurrentRecords()
{
    for(int i=0; i<NoStudiesHasVariant;i++)
    {
        int index = StudiesHasVariant[i];
        if(!InputDosageStream[index]->readRecord(*CurrentRecordFromStudy[index]))
            CurrentRecordFromStudy[index]->set1BasedPosition(MAXBP);
    }
}


void MetaMinimac::LoadLooDosage()
{
    printf(" -- Loading Empirical Dosage Data ...\n");
    for(int i=0; i<NoInPrefix; i++)
        InputData[i].ReadBasedOnSortCommonGenotypeList(CommonGenotypeVariantNameList, StartSamId, EndSamId);
}


bool MetaMinimac::PerformWeightEstimation()
{
    cout << endl;
    cout << " ------------------------------------------------------------------------------" << endl;
    cout << "                               WEIGHT ESTIMATION                               " << endl;
    cout << " ------------------------------------------------------------------------------" << endl;

    if(!OpenStreamOutputWeightFiles())
        return false;

    int start_time = time(0);

    int maxVcfSample = myUserVariables.VcfBuffer;

    StartSamId = 0;

    batchNo = 0;

    int start_time_partial, tot_time_partial;

    while(true)
    {
        batchNo++;
        EndSamId = StartSamId + (maxVcfSample) < NoSamples ? StartSamId + (maxVcfSample) : NoSamples;
        cout << "\n Estimate Meta-Weights for Sample " << StartSamId + 1 << "-" << EndSamId << " [" << setprecision(1) << fixed << 100 * (float) EndSamId / NoSamples << "%] ..." << endl;

        start_time_partial = time(0);
        // Read Data From empiricalDose
        LoadLooDosage();

        // Calculate weights
        CalculateWeights();

        // Output weights
        OutputWeights();

        tot_time_partial = time(0) - start_time_partial;
        cout << " -- Successful (" << tot_time_partial << " seconds) !!! " << endl;

        StartSamId = EndSamId;

        if (StartSamId >= NoSamples)
            break;
    }

    if(batchNo > 1)
    {
        AppendtoMainWeightsFile();
    }
    int time_tot = time(0) - start_time;
    printf("\n Weight Estimation Completed in %d hours, %d mins, %d seconds.\n",
           time_tot / 3600, (time_tot % 3600) / 60, time_tot % 60);

    return true;
}

bool MetaMinimac::OpenStreamOutputWeightFiles()
{
    // output file for weights
    metaWeight = ifopen(myUserVariables.outfile + ".metaWeights"+(myUserVariables.gzip ? ".gz" : ""), "wb", myUserVariables.gzip ? InputFile::BGZF : InputFile::UNCOMPRESSED);
    WeightPrintStringPointer = (char*)malloc(sizeof(char) * (myUserVariables.PrintBuffer));
    WeightPrintStringPointerLength = 0;
    if(metaWeight==NULL)
    {
        cout <<"\n\n ERROR !!! \n Could NOT create the following file : "<< myUserVariables.outfile + ".metaWeights"+(myUserVariables.gzip ? ".gz" : "") <<endl;
        return false;
    }
    ifprintf(metaWeight,"##fileformat=VCFv4.1\n");
    time_t t = time(0);
    struct tm * now = localtime( & t );
    ifprintf(metaWeight,"##filedate=%d.%d.%d\n",(now->tm_year + 1900),(now->tm_mon + 1) ,now->tm_mday);
    ifprintf(metaWeight,"##source=MetaMinimac2.v%s\n",VERSION);
    ifprintf(metaWeight,"##contig=<ID=%s>\n", finChromosome.c_str());
    for (int i=0; i<NoInPrefix; i++)
    {
        ifprintf(metaWeight,"##FORMAT=<ID=WT%d,Number=2,Type=Float,Description=\"Estimated Meta Weights on Study %d: [ weight on haplotype 1 , weight on haplotype 2 ]\">\n", i+1, i+1);
    }

    ifprintf(metaWeight,"##metaMinimac_Command=%s\n",myUserVariables.CommandLine.c_str());

    ifprintf(metaWeight,"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");

    for(int Id=0;Id<InputData[0].numSamples;Id++)
    {
        ifprintf(metaWeight,"\t%s",InputData[0].individualName[Id].c_str());
    }
    ifprintf(metaWeight,"\n");
    ifclose(metaWeight);
    return true;
}


void MetaMinimac::CalculateWeights()
{
    cout << " -- Calculating Weights ... " << endl;
    InitiateWeights();
    CalculateLeftProbs();
    CalculatePosterior();
}

void MetaMinimac::InitiateWeights()
{
    int NoSamplesThisBatch = EndSamId-StartSamId;
    Weights.clear();
    Weights.resize(NoCommonTypedVariants);
    for(int i=0; i<NoCommonTypedVariants; i++)
    {
        vector<vector<double> > &ThisWeights = Weights[i];
        ThisWeights.resize(2*NoSamplesThisBatch);
        for(int j=0; j<2*NoSamplesThisBatch; j++)
            ThisWeights[j].resize(NoInPrefix);
    }
}

void MetaMinimac::CalculateLeftProbs()
{
    int NoSamplesThisBatch = EndSamId-StartSamId;
    for (int id=0; id<NoSamplesThisBatch; id++)
    {
        int SampleId = StartSamId + id;
        if (InputData[0].SampleNoHaplotypes[SampleId] == 2)
        {
            InitiateLeftProb(2*id);
            InitiateLeftProb(2*id+1);
        }
        else
            InitiateLeftProb(2*id);
    }

    NoCommonVariantsProcessed = 1;
    while(NoCommonVariantsProcessed < NoCommonTypedVariants)
    {
        Recom = TransitionProb[NoCommonTypedVariants-NoCommonVariantsProcessed];
        NoCommonVariantsProcessed++;
        for (int id=0; id<NoSamplesThisBatch; id++)
        {
            int SampleId = StartSamId + id;
            if (InputData[0].SampleNoHaplotypes[SampleId] == 2)
            {
                UpdateOneStepLeft(2*id);
                UpdateOneStepLeft(2*id+1);
            }
            else
                UpdateOneStepLeft(2*id);
        }
    }

}

void MetaMinimac::InitiateLeftProb(int HapInBatch)
{
    LogOddsModel ThisSampleAnalysis;
    ThisSampleAnalysis.reinitialize(HapInBatch, this);
    vector<double> init(NoInPrefix-1, 0.0);
    vector<double> MiniMizer = Simplex(ThisSampleAnalysis, init);

    vector<double> InitProb;
    InitProb.resize(NoInPrefix);
    logitTransform(MiniMizer, InitProb);

    float ThisGT = InputData[0].TypedGT[HapInBatch][NoCommonTypedVariants-1];
    vector<double> &ThisWeights = Weights[NoCommonTypedVariants-1][HapInBatch];

    for(int i=0; i<NoInPrefix; i++)
    {
        InitProb[i]+=backgroundError;
        float ThisLooDosage = InputData[i].LooDosage[HapInBatch][NoCommonTypedVariants-1];
        InitProb[i] *= (ThisGT==1)?(ThisLooDosage+backgroundError):(1-ThisLooDosage+backgroundError);
        ThisWeights[i] = InitProb[i];
    }
}

void MetaMinimac::UpdateOneStepLeft(int HapInBatch)
{
    double r = Recom*1.0/NoInPrefix, complement = 1-Recom;
    float ThisGT = InputData[0].TypedGT[HapInBatch][NoCommonTypedVariants-NoCommonVariantsProcessed];
    vector<double> &ThisLeftProb = Weights[NoCommonTypedVariants-NoCommonVariantsProcessed][HapInBatch];
    vector<double> &ThisPrevLeftProb = Weights[NoCommonTypedVariants-NoCommonVariantsProcessed+1][HapInBatch];

    double sum = 0.0;
    for(int i=0; i<NoInPrefix; i++)
    {
        for(int j=0; j<NoInPrefix; j++)
            ThisLeftProb[i] += ThisPrevLeftProb[j]*r;
        ThisLeftProb[i] += ThisPrevLeftProb[i]*complement;

        float ThisLooDosage = InputData[i].LooDosage[HapInBatch][NoCommonTypedVariants-NoCommonVariantsProcessed];
        ThisLeftProb[i] *= (ThisGT==1)?(ThisLooDosage+backgroundError):(1-ThisLooDosage+backgroundError);
        sum += ThisLeftProb[i];
    }


    while(sum < JumpThreshold)
    {
        sum = 0.0;
        for(int i=0; i<NoInPrefix; i++)
        {
            ThisLeftProb[i] *= JumpFix;
            sum += ThisLeftProb[i];
        }
    }

}

void MetaMinimac::CalculatePosterior()
{
    int NoSamplesThisBatch = EndSamId-StartSamId;
    PrevRightProb.clear();
    PrevRightProb.resize(2*NoSamplesThisBatch);
    for (int id=0; id<NoSamplesThisBatch; id++)
    {
        int SampleId = StartSamId + id;
        if (InputData[0].SampleNoHaplotypes[SampleId] == 2)
        {
            InitiateRightProb(2*id);
            InitiateRightProb(2*id+1);
        }
        else
            InitiateRightProb(2*id);
    }

    NoCommonVariantsProcessed = 1;
    while(NoCommonVariantsProcessed < NoCommonTypedVariants)
    {
        Recom = TransitionProb[NoCommonVariantsProcessed];
        for (int id=0; id<NoSamplesThisBatch; id++)
        {
            int SampleId = StartSamId + id;
            if (InputData[0].SampleNoHaplotypes[SampleId] == 2)
            {
                UpdateOneStepRight(2*id);
                UpdateOneStepRight(2*id+1);
            }
            else
                UpdateOneStepRight(2*id);
        }
        NoCommonVariantsProcessed++;
    }

}

void MetaMinimac::InitiateRightProb(int HapInBatch)
{
    PrevRightProb[HapInBatch].resize(NoInPrefix, 1.0);
}


void MetaMinimac::UpdateOneStepRight(int HapInBatch)
{

    double r = Recom*1.0/NoInPrefix, complement = 1-Recom;
    float ThisGT = InputData[0].TypedGT[HapInBatch][NoCommonVariantsProcessed-1];
    vector<double> &ThisWeight = Weights[NoCommonVariantsProcessed][HapInBatch];
    vector<double> &ThisPrevRightProb = PrevRightProb[HapInBatch];
    vector<double> ThisRightProb;
    ThisRightProb.resize(NoInPrefix, 0.0);

    for(int i=0; i<NoInPrefix; i++)
    {
        float ThisLooDosage = InputData[i].LooDosage[HapInBatch][NoCommonVariantsProcessed-1];
        ThisPrevRightProb[i] *= (ThisGT==1)?(ThisLooDosage+backgroundError):(1-ThisLooDosage+backgroundError);
    }

    double sum = 0.0;
    for(int i=0; i<NoInPrefix; i++)
    {
        for(int j=0; j<NoInPrefix; j++)
            ThisRightProb[i] += ThisPrevRightProb[j]*r;
        ThisRightProb[i] += ThisPrevRightProb[i]*complement;
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

    for(int i=0; i<NoInPrefix; i++)
    {
        ThisWeight[i] *= ThisRightProb[i];
        ThisPrevRightProb[i] = ThisRightProb[i];
    }
}

void MetaMinimac::OutputWeights()
{
    if(myUserVariables.VcfBuffer<NoSamples)
    {
        OutputPartialWeights();
    }
    else
    {
        OutputAllWeights();
    }
}

void MetaMinimac::OutputPartialWeights()
{
    stringstream ss;
    ss << (batchNo);
    string PartialWeightFileName(myUserVariables.outfile);
    PartialWeightFileName += ".metaWeights.part."+(string)(ss.str()) + (myUserVariables.gzip ? ".gz" : "");
    vcfweightpartial = ifopen(PartialWeightFileName.c_str(), "wb", myUserVariables.gzip ? InputFile::BGZF : InputFile::UNCOMPRESSED);
    WeightPrintStringPointerLength = 0;

    NoCommonVariantsProcessed = 0;

    while ( NoCommonVariantsProcessed  < NoCommonTypedVariants)
    {
        CurrWeights = &Weights[NoCommonVariantsProcessed];
        PrintMetaWeight();
        NoCommonVariantsProcessed++;
    }

    if(WeightPrintStringPointerLength > 0)
    {
        ifprintf(vcfweightpartial, "%s", WeightPrintStringPointer);
        WeightPrintStringPointerLength = 0;
    }
    ifclose(vcfweightpartial);
}

void MetaMinimac::OutputAllWeights()
{
    WeightPrintStringPointerLength=0;
    vcfweightpartial = ifopen(myUserVariables.outfile + ".metaWeights"+ (myUserVariables.gzip ? ".gz" : ""), "a", myUserVariables.gzip ? InputFile::BGZF : InputFile::UNCOMPRESSED);

    NoCommonVariantsProcessed = 0;

    while ( NoCommonVariantsProcessed  < NoCommonTypedVariants)
    {
        CurrWeights = &Weights[NoCommonVariantsProcessed];
        PrintWeightVariantInfo();
        PrintMetaWeight();
        NoCommonVariantsProcessed++;
    }

    if(WeightPrintStringPointerLength > 0)
    {
        ifprintf(vcfweightpartial, "%s", WeightPrintStringPointer);
        WeightPrintStringPointerLength = 0;
    }
    ifclose(vcfweightpartial);

}

String MetaMinimac::PerformFinalAnalysis()
{
    cout << endl;
    cout << " ------------------------------------------------------------------------------" << endl;
    cout << "                               FINAL ANALYSIS                               " << endl;
    cout << " ------------------------------------------------------------------------------" << endl;

    printf("\n Gathering dosege information and saving results for all samples ...\n");

    int start_time = time(0);

    OpenStreamInputDosageFiles();
    OpenStreamInputWeightFiles();

    vcfdosepartial = ifopen(myUserVariables.outfile + ".metaDose.vcf"+ (myUserVariables.gzip ? ".gz" : ""), "a", myUserVariables.gzip ? InputFile::BGZF : InputFile::UNCOMPRESSED);
    VcfPrintStringPointerLength=0;

    NoVariants = 0;
    NoCommonVariantsProcessed = 0;
    int NoRecordProcessed = 0;

    if(!InitiateWeightsFromRecord())
    {
        cout << " Weight File Mismatch !!! \n";
        cout << "\n Program Exiting ... \n\n";
        return "Weight.Mismatch.Error";
    }

    BufferBp = 0;
    BufferNoVariants = 0;

    do
    {
        FindCurrentMinimumPosition();
        if(CurrentFirstVariantBp>BufferBp)
        {
            MetaImputeCurrentBuffer();
            ClearCurrentBuffer();
            if(CurrentFirstVariantBp == MAXBP) break;
            while(CurrentFirstVariantBp == CurrBp)
            {
                UpdateWeights();
            }
        }
        ReadCurrentDosageData();
        UpdateCurrentRecords();
        NoRecordProcessed++;
    }while(true);
    NoRecords = NoRecordProcessed;

    if (VcfPrintStringPointerLength > 0)
    {
        ifprintf(vcfdosepartial,"%s",VcfPrintStringPointer);
        VcfPrintStringPointerLength = 0;
    }
    ifclose(vcfdosepartial);


    CloseStreamInputDosageFiles();
    CloseStreamInputWeightFiles();

    int time_tot = time(0) - start_time;
    printf("\n Final Analysis Completed in %d hours, %d mins, %d seconds.\n",
           time_tot / 3600, (time_tot % 3600) / 60, time_tot % 60);

    return "Success";
}

bool MetaMinimac::InitiateWeightsFromRecord()
{
    PrevBp = 0, CurrBp = CommonTypedVariantList[0].bp;
    WeightsFromCurrentRecord.clear();
    WeightsFromPreviousRecord.clear();
    WeightsFromCurrentRecord.resize(2*NoSamples);
    WeightsFromPreviousRecord.resize(2*NoSamples);
    for (int i=0; i<2*NoSamples; i++)
    {
        WeightsFromCurrentRecord[i].resize(NoInPrefix);
        WeightsFromPreviousRecord[i].resize(NoInPrefix);
    }

    ReadCurrentWeights();

    int bp = CurrentRecordFromWeight->get1BasedPosition();
    if (bp != CurrBp)
    {
        return false;
    }

    CopyCurrentWeightsToPreviousWeights();

    return true;
}

void MetaMinimac::ReadCurrentWeights()
{
    if(InputWeightStream->readRecord(*CurrentRecordFromWeight))
    {
        VcfRecordGenotype & ThisGenotype = CurrentRecordFromWeight->getGenotypeInfo();

        for (int i=0; i<NoSamples; i++)
        {
            for (int j=0; j<NoInPrefix; j++)
            {
                stringstream ss;
                ss << "WT" << j+1;
                string temp = *ThisGenotype.getString(ss.str(), i);
                if(InputData[0].SampleNoHaplotypes[i]==2)
                {
                    char *end_str;
                    char *pch = strtok_r((char *) temp.c_str(), ",", &end_str);
                    WeightsFromCurrentRecord[2*i][j] = atof(pch);
                    pch = strtok_r(NULL, ",", &end_str);
                    WeightsFromCurrentRecord[2*i+1][j] = atof(pch);
                }
                else
                {
                    WeightsFromCurrentRecord[2*i][j] = atof(temp.c_str());
                }
            }
        }
    }
}

void MetaMinimac::CopyCurrentWeightsToPreviousWeights()
{
    for (int i=0; i<2*NoSamples; i++)
    {
        vector<float> &ThisWeight = WeightsFromCurrentRecord[i];
        vector<float> &PrevWeight = WeightsFromPreviousRecord[i];
        for (int j=0; j<NoInPrefix; j++)
        {
            PrevWeight[j] = ThisWeight[j];
        }

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


void MetaMinimac::ReadCurrentDosageData()
{
    VcfRecord* tempRecord = CurrentRecordFromStudy[StudiesHasVariant[0]];
    variant tempVariant;
    tempVariant.chr  = tempRecord->getChromStr();
    tempVariant.bp   = tempRecord->get1BasedPosition();
    tempVariant.refAlleleString = tempRecord->getRefStr();
    tempVariant.altAlleleString = tempRecord->getAltStr();
    tempVariant.name = tempVariant.chr+":"+to_string(tempVariant.bp)+":"+ tempVariant.refAlleleString+":"+tempVariant.altAlleleString;

    int VariantId = 0;
    for(vector<variant>::iterator thisVariant = BufferVariantList.begin(); thisVariant < BufferVariantList.end(); ++thisVariant)
    {
        if(tempVariant.name == thisVariant->name)
        {
            int count = thisVariant->NoStudiesHasVariant;
            for(int i=0; i<NoStudiesHasVariant; i++)
            {
                thisVariant->StudiesHasVariant[count]=StudiesHasVariant[i];
                count++;
            }
            thisVariant->NoStudiesHasVariant = count;
            break;
        }
        VariantId++;
    }
    if(VariantId==BufferNoVariants)
    {
        tempVariant.NoStudiesHasVariant = NoStudiesHasVariant;
        tempVariant.StudiesHasVariant = StudiesHasVariant;
        BufferVariantList.push_back(tempVariant);
        BufferNoVariants++;
    }

    for(int j=0; j<NoStudiesHasVariant; j++)
    {
        int index = StudiesHasVariant[j];
//        InputData[index].LoadData(VariantId,CurrentRecordFromStudy[index]->getGenotypeInfo(), StartSamId, EndSamId);
        InputData[index].LoadData(VariantId,CurrentRecordFromStudy[index]->getGenotypeInfo(), 0, NoSamples);
    }
}


void MetaMinimac::MetaImputeCurrentBuffer()
{
    for(int VariantId=0; VariantId<BufferNoVariants; VariantId++)
    {
        CurrentVariant = &BufferVariantList[VariantId];
        CreateMetaImputedData(VariantId);
        PrintVariantInfo();
        PrintMetaImputedData();
    }
}

void MetaMinimac::CreateMetaImputedData(int VariantId)
{
    CurrentHapDosageSum = 0;
    CurrentHapDosageSumSq = 0;
    CurrentMetaImputedDosage.clear();
    if(CurrentVariant->NoStudiesHasVariant==1)
    {
        CurrentMetaImputedDosage=InputData[CurrentVariant->StudiesHasVariant[0]].BufferHapDosage[0];
        for(int i=0; i<2*NoSamples; i++)
        {
            CurrentHapDosageSum += CurrentMetaImputedDosage[i];
            CurrentHapDosageSumSq += CurrentMetaImputedDosage[i]*CurrentMetaImputedDosage[i];
        }
    }
    else
    {
        CurrentMetaImputedDosage.resize(2*NoSamples);
        for (int j=0; j<CurrentVariant->NoStudiesHasVariant; j++)
        {
            int index = CurrentVariant->StudiesHasVariant[j];
            InputData[index].GetData(VariantId);
        }
        for(int id=0; id<NoSamples; id++)
        {
            if(InputData[0].SampleNoHaplotypes[id]==2)
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



void MetaMinimac::UpdateWeights()
{
    NoCommonVariantsProcessed++;
    CopyCurrentWeightsToPreviousWeights();
    PrevBp = CurrBp;
    if(NoCommonVariantsProcessed < NoCommonTypedVariants)
    {
        ReadCurrentWeights();
        CurrBp = CommonTypedVariantList[NoCommonVariantsProcessed].bp;
    }
    else
    {
        CurrBp = MAXBP;
    }
}



void MetaMinimac::MetaImpute(int Sample)
{
    vector<float>& ThisPrevWeights = WeightsFromPreviousRecord[Sample];
    vector<float>& ThisCurrWeights = WeightsFromCurrentRecord[Sample];

    float WeightSum = 0.0;
    float Dosage = 0.0;

    for (int j=0; j<CurrentVariant->NoStudiesHasVariant; j++)
    {
        int index = CurrentVariant->StudiesHasVariant[j];
        float Weight = (ThisPrevWeights[index]*(CurrBp-BufferBp)+ThisCurrWeights[index]*(BufferBp-PrevBp))*1.0/(CurrBp-PrevBp);
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
    for(int id=0; id<NoSamples; id++)
    {
        if( InputData[0].SampleNoHaplotypes[id]==2 )
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
            PrintDiploidWeightForSample(id);
//            PrintWeightForHaplotype(2*id);
//            WeightPrintStringPointerLength += sprintf(WeightPrintStringPointer+WeightPrintStringPointerLength,"|");
//            PrintWeightForHaplotype(2*id+1);
        }
        else
        {
            PrintHaploidWeightForSample(id);
//            PrintWeightForHaplotype(2*id);
        }
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


//void MetaMinimac::PrintWeightForHaplotype(int haploId)
//{
//    vector<double>& ThisCurrWeights = (*CurrWeights)[haploId];
//    double WeightSum = 0.0;
//    for(int i=0; i<NoInPrefix; i++)
//        WeightSum += ThisCurrWeights[i];
//    WeightPrintStringPointerLength += sprintf(WeightPrintStringPointer+WeightPrintStringPointerLength,"%0.4f", ThisCurrWeights[0]/WeightSum);
//    for(int i=1;i<NoInPrefix;i++)
//        WeightPrintStringPointerLength += sprintf(WeightPrintStringPointer+WeightPrintStringPointerLength,",%0.4f", ThisCurrWeights[i]/WeightSum);
//}

void MetaMinimac::PrintDiploidWeightForSample(int sampleId)
{
    vector<double>& WeightsHap1 = (*CurrWeights)[2*sampleId];
    vector<double>& WeightsHap2 = (*CurrWeights)[2*sampleId+1];
    double WeightSumHap1 = 0.0, WeightSumHap2 = 0.0;
    for(int i=0; i<NoInPrefix; i++)
    {
        WeightSumHap1 += WeightsHap1[i];
        WeightSumHap2 += WeightsHap2[i];
    }
    WeightPrintStringPointerLength += sprintf(WeightPrintStringPointer+WeightPrintStringPointerLength,"%0.4f,%0.4f", WeightsHap1[0]/WeightSumHap1, WeightsHap2[0]/WeightSumHap2);
    for(int i=1;i<NoInPrefix;i++)
        WeightPrintStringPointerLength += sprintf(WeightPrintStringPointer+WeightPrintStringPointerLength,":%0.4f,%0.4f", WeightsHap1[i]/WeightSumHap1, WeightsHap2[i]/WeightSumHap2);
}

void MetaMinimac::PrintHaploidWeightForSample(int sampleId)
{
    vector<double>& ThisCurrWeights = (*CurrWeights)[2*sampleId];
    double WeightSum = 0.0;
    for(int i=0; i<NoInPrefix; i++)
        WeightSum += ThisCurrWeights[i];
    WeightPrintStringPointerLength += sprintf(WeightPrintStringPointer+WeightPrintStringPointerLength,"%0.4f", ThisCurrWeights[0]/WeightSum);
    for(int i=1;i<NoInPrefix;i++)
        WeightPrintStringPointerLength += sprintf(WeightPrintStringPointer+WeightPrintStringPointerLength,":%0.4f", ThisCurrWeights[i]/WeightSum);

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
    if(!myUserVariables.infoDetails)
        return ".";

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

    ss<<"NST=";
    ss<<CurrentVariant->NoStudiesHasVariant;
    for(int i=0; i<CurrentVariant->NoStudiesHasVariant; i++)
    {
        ss<<";S";
        ss<<CurrentVariant->StudiesHasVariant[i]+1;
    }

    if(CurrentVariant->name == CommonGenotypeVariantNameList[NoCommonVariantsProcessed-1])
    {
        ss<<";TRAINING";
    }

    ss<< ";AF=" << fixed << setprecision(5) << freq <<";MAF=";
    ss<< fixed << setprecision(5) << maf <<";R2=";
    ss<< fixed << setprecision(5) << rsq ;

    return ss.str();
}


void MetaMinimac::PrintWeightVariantInfo()
{
    // "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
    variant& tempVariant = CommonTypedVariantList[NoCommonVariantsProcessed];
    WeightPrintStringPointerLength+=sprintf(WeightPrintStringPointer+WeightPrintStringPointerLength,"%s\t%d\t%s\t%s\t%s\t.\tPASS\t.\tWT1",
                                            tempVariant.chr.c_str(), tempVariant.bp, tempVariant.name.c_str(),
                                            tempVariant.refAlleleString.c_str(), tempVariant.altAlleleString.c_str());
    for(int i=1; i<NoInPrefix; i++)
    {
        WeightPrintStringPointerLength+=sprintf(WeightPrintStringPointer+WeightPrintStringPointerLength,":WT%d", i+1);
    }
    if(WeightPrintStringPointerLength > 0.9 * (float)(myUserVariables.PrintBuffer))
    {
        ifprintf(vcfweightpartial, "%s", WeightPrintStringPointer);
        WeightPrintStringPointerLength = 0;
    }
}

void MetaMinimac::AppendtoMainWeightsFile()
{
    cout << "\n Appending to final output weight file : " << myUserVariables.outfile + ".metaWeights" + (myUserVariables.gzip ? ".gz" : "") <<endl;
    int start_time = time(0);
    WeightPrintStringPointerLength=0;
    metaWeight = ifopen(myUserVariables.outfile + ".metaWeights" + (myUserVariables.gzip ? ".gz" : ""), "a", myUserVariables.gzip ? InputFile::BGZF : InputFile::UNCOMPRESSED);
    vector<IFILE> weightpartialList(batchNo);

    for(int i=1; i<=batchNo; i++)
    {
        stringstream ss;
        ss << (i);
        string PartialWeightFileName(myUserVariables.outfile);
        PartialWeightFileName += ".metaWeights.part."+(string)(ss.str())+(myUserVariables.gzip ? ".gz" : "");
        weightpartialList[i-1] = ifopen(PartialWeightFileName.c_str(), "r");
    }

    string line;

    NoCommonVariantsProcessed = 0;
    for(int i=0; i<NoCommonTypedVariants; i++)
    {
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
        NoCommonVariantsProcessed++;
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
        tempFileIndex += ".metaWeights.part."+(string)(ss.str())+(myUserVariables.gzip ? ".gz" : "");
        remove(tempFileIndex.c_str());
    }
    ifclose(metaWeight);
    int tot_time = time(0) - start_time;
    cout << " -- Successful (" << tot_time << " seconds) !!!" << endl;
}

void MetaMinimac::ClearCurrentBuffer()
{
    NoVariants += BufferNoVariants;
    for(int i=0; i<NoInPrefix; i++)
        InputData[i].ClearBuffer();
    BufferVariantList.clear();
    BufferNoVariants = 0;
    BufferBp = CurrentFirstVariantBp;
}