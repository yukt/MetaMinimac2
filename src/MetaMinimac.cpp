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

std::string MetaMinimac::Analyze()
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

    if (myUserVariables.inputFiles.size())
    {
        size_t pos = 0;
        std::string delimiter(myUserVariables.FileDelimiter) ;
        std::string token;
        int Count=0;
        string tempName=myUserVariables.inputFiles.c_str();
        while ((pos = tempName.find(delimiter)) != std::string::npos)
        {
            token = tempName.substr(0, pos);
            myUserVariables.empInputFiles.emplace_back(token + ".empiricalDose.vcf.gz");
            myUserVariables.doseInputFiles.emplace_back(token + ".dose.vcf.gz");
            tempName.erase(0, pos + delimiter.length());
            Count++;
        }

        myUserVariables.empInputFiles.emplace_back(tempName + ".empiricalDose.vcf.gz");
        myUserVariables.doseInputFiles.emplace_back(tempName + ".dose.vcf.gz");
    }


    NoInPrefix = myUserVariables.empInputFiles.size();
    InputData.clear();
    InputData.resize(NoInPrefix);

    cout<<endl<<  " Number of Studies : "<<NoInPrefix<<endl;
    for(int i=0;i<NoInPrefix;i++)
    {
        cout<<  " -- Study "<<i+1<<" Prefix : "<<myUserVariables.empInputFiles[i]<<endl;
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
        if(!InputData[i].LoadSampleNames(myUserVariables.empInputFiles[i], myUserVariables.doseInputFiles[i]))
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
        InputDosageStream[i] = new savvy::reader(InputData[i].DoseFileName);
        CurrentRecordFromStudy[i]= new savvy::variant("", MAXBP, "", {""});
        InputDosageStream[i]->read(*CurrentRecordFromStudy[i]);
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
        InputDosageStream[i] = new savvy::reader(InputData[i].EmpDoseFileName);
        CurrentRecordFromStudy[i]= new savvy::variant("", MAXBP, "", {""});
    }
}


void MetaMinimac::OpenStreamInputWeightFiles()
{
    InputWeightStream = new savvy::reader(std::string(myUserVariables.outfile.c_str()) + ".metaWeights" + (myUserVariables.gzip ? ".gz" : ""));
    CurrentRecordFromWeight = new savvy::variant();
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
    std::time_t t = std::time(nullptr);
    char datestr[11];
    std::string filedate(datestr, std::strftime(datestr, sizeof(datestr), "%Y%m%d", std::localtime(&t)));
    assert(filedate.size());

    std::vector<std::pair<std::string, std::string>> headers = {
        {"fileformat", "VCFv4.1"},
        {"filedate", filedate},
        {"source", std::string("MetaMinimac2.") + VERSION},
        {"phasing", "full"},
        {"contig", "<ID=" + finChromosome + ">"},
        {"INFO","<ID=AF,Number=1,Type=Float,Description=\"Estimated Alternate Allele Frequency\">"},
        {"INFO","<ID=MAF,Number=1,Type=Float,Description=\"Estimated Minor Allele Frequency\">"},
        {"INFO","<ID=R2,Number=1,Type=Float,Description=\"Estimated Imputation Accuracy (R-square)\">"},
        {"INFO","<ID=TRAINING,Number=0,Type=Flag,Description=\"Marker was used to train meta-imputation weights\">"}};

    headers.reserve(headers.size() + NoInPrefix + 7);

    if(myUserVariables.infoDetails)
    {
        headers.emplace_back("INFO","<ID=NST,Number=1,Type=Integer,Description=\"Number of studies marker was found during meta-imputation\">");
        for(int i=1; i<=NoInPrefix; i++)
            headers.emplace_back("INFO","<ID=S" + std::to_string(i) + ",Number=0,Type=Flag,Description=\"Marker was present in Study " + std::to_string(i) + "\">");
    }
    if(myUserVariables.GT)
        headers.emplace_back("FORMAT","<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    if(myUserVariables.DS)
        headers.emplace_back("FORMAT","<ID=DS,Number=1,Type=Float,Description=\"Estimated Alternate Allele Dosage : [P(0/1)+2*P(1/1)]\">");
    if(myUserVariables.HDS)
        headers.emplace_back("FORMAT","<ID=HDS,Number=2,Type=Float,Description=\"Estimated Haploid Alternate Allele Dosage \">");
    if(myUserVariables.GP)
        headers.emplace_back("FORMAT","<ID=GP,Number=3,Type=Float,Description=\"Estimated Posterior Probabilities for Genotypes 0/0, 0/1 and 1/1 \">");
    if(myUserVariables.SD)
        headers.emplace_back("FORMAT","<ID=SD,Number=1,Type=Float,Description=\"Variance of Posterior Genotype Probabilities\">");

    headers.emplace_back("metaMinimac_Command", myUserVariables.CommandLine.c_str());


    metaDoseOut.reset(new savvy::writer(std::string(myUserVariables.outfile.c_str()) + ".metaDose.vcf"+(myUserVariables.gzip ? ".gz" : ""), savvy::file::format::vcf, headers, InputData[0].individualName, myUserVariables.gzip ? 6 : 0));
    if(!metaDoseOut->good())
    {
        cout <<"\n\n ERROR !!! \n Could NOT create the following file : "<< myUserVariables.outfile + ".metaDose.vcf"+(myUserVariables.gzip ? ".gz" : "") <<endl;
        return false;
    }

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
            while(InputDosageStream[i]->read(*CurrentRecordFromStudy[i])) {
                chrom = CurrentRecordFromStudy[i]->chrom();
                bp = CurrentRecordFromStudy[i]->pos();
                altAllele = CurrentRecordFromStudy[i]->alts().empty() ? "" : CurrentRecordFromStudy[i]->alts()[0];
                refAllele = CurrentRecordFromStudy[i]->ref();
                name = chrom + ":" + to_string(bp) + ":" + refAllele + ":" + altAllele;
                if (name == common_variant)
                {
                    InputData[i].LoadCurrentGT(*CurrentRecordFromStudy[i]);
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
        int a = CurrentRecordFromStudy[0]->pos();
        int b = CurrentRecordFromStudy[1]->pos();
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

        CurrentFirstVariantBp=CurrentRecordFromStudy[0]->pos();

        for(int i=1;i<NoInPrefix;i++)
            if(CurrentRecordFromStudy[i]->pos() < CurrentFirstVariantBp)
                CurrentFirstVariantBp=CurrentRecordFromStudy[i]->pos();

        NoStudiesHasVariant=0;
        savvy::variant *minRecord = nullptr;

        int Begin=0;
        for(int i=0;i<NoInPrefix;i++)
        {
            if(CurrentRecordFromStudy[i]->pos() == CurrentFirstVariantBp)
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

int MetaMinimac::IsVariantEqual(savvy::variant &Rec1, savvy::variant &Rec2)
{
    if(Rec1.ref() != Rec2.ref())
        return 0;
    if(Rec1.alts() != Rec2.alts())
        return 0;
    return 1;
}



void MetaMinimac::UpdateCurrentRecords()
{
    for(int i=0; i<NoStudiesHasVariant;i++)
    {
        int index = StudiesHasVariant[i];
        if(!InputDosageStream[index]->read(*CurrentRecordFromStudy[index]))
            *CurrentRecordFromStudy[index] = savvy::site_info("", MAXBP, "", {""});
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
        if (!AppendtoMainWeightsFile())
          return false;
    }
    int time_tot = time(0) - start_time;
    printf("\n Weight Estimation Completed in %d hours, %d mins, %d seconds.\n",
           time_tot / 3600, (time_tot % 3600) / 60, time_tot % 60);

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
        OutputAllWeights(std::string(myUserVariables.outfile.c_str()) + ".metaWeights.part." + std::to_string(batchNo) + (myUserVariables.gzip ? ".gz" : ""));
    }
    else
    {
        OutputAllWeights(std::string(myUserVariables.outfile.c_str()) + ".metaWeights"+ (myUserVariables.gzip ? ".gz" : ""));
    }
}

void MetaMinimac::OutputAllWeights(const std::string& out_file_path)
{
    std::time_t t = std::time(nullptr);
    char datestr[11];
    std::string filedate(datestr, std::strftime(datestr, sizeof(datestr), "%Y%m%d", std::localtime(&t)));
    assert(filedate.size());

    std::vector<std::pair<std::string, std::string>> headers = {
        {"fileformat", "VCFv4.1"},
        {"filedate", filedate},
        {"source", std::string("MetaMinimac2.") + VERSION},
        {"phasing", "full"},
        {"contig", "<ID=" + finChromosome + ">"}};

    for (int i=1; i<=NoInPrefix; i++)
    {
        headers.emplace_back("FORMAT","<ID=WT" + std::to_string(i) + ",Number=2,Type=Float,Description=\"Estimated Meta Weights on Study " + std::to_string(i) + ": [ weight on haplotype 1 , weight on haplotype 2 ]\">");
    }

    headers.emplace_back("metaMinimac_Command", myUserVariables.CommandLine.c_str());

    savvy::writer out_weight_file(out_file_path, savvy::file::format::vcf, headers, {InputData[0].individualName.begin() + StartSamId,  InputData[0].individualName.begin() + EndSamId}, myUserVariables.gzip ? 6 : 0);
    savvy::variant out_record;
    std::vector<std::vector<float>> wt_vecs(NoInPrefix, std::vector<float>(Weights[0].size(), savvy::typed_value::end_of_vector_value<float>()));
    NoCommonVariantsProcessed = 0;
    while ( NoCommonVariantsProcessed  < NoCommonTypedVariants)
    {
        CurrWeights = &Weights[NoCommonVariantsProcessed];
        SetWeightVariantInfo(out_record);
        SetMetaWeightVecs(wt_vecs);
        for (int i=1; i<=NoInPrefix; ++i)
            out_record.set_format("WT" + std::to_string(i), wt_vecs[i-1]);
        out_weight_file << out_record;
        NoCommonVariantsProcessed++;
    }

}

std::string MetaMinimac::PerformFinalAnalysis()
{
    cout << endl;
    cout << " ------------------------------------------------------------------------------" << endl;
    cout << "                               FINAL ANALYSIS                               " << endl;
    cout << " ------------------------------------------------------------------------------" << endl;

    printf("\n Gathering dosege information and saving results for all samples ...\n");

    int start_time = time(0);

    OpenStreamInputDosageFiles();
    OpenStreamInputWeightFiles();

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

    int bp = CurrentRecordFromWeight->pos();
    if (bp != CurrBp)
    {
        return false;
    }

    CopyCurrentWeightsToPreviousWeights();

    return true;
}

void MetaMinimac::ReadCurrentWeights()
{
    if(InputWeightStream->read(*CurrentRecordFromWeight))
    {
        std::vector<float> wt_vec;
        for (std::size_t i = 0; i < NoInPrefix; ++i)
        {
            CurrentRecordFromWeight->get_format("WT" + std::to_string(i + 1), wt_vec);
            assert(WeightsFromCurrentRecord.size() == wt_vec.size());
            for (std::size_t j = 0; j < wt_vec.size(); ++j)
                WeightsFromCurrentRecord[j][i] = wt_vec[j];
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
    savvy::variant* tempRecord = CurrentRecordFromStudy[StudiesHasVariant[0]];
    variant tempVariant;
    tempVariant.chr  = tempRecord->chrom();
    tempVariant.bp   = tempRecord->pos();
    tempVariant.refAlleleString = tempRecord->ref();
    tempVariant.altAlleleString = tempRecord->alts().empty() ? "" : tempRecord->alts()[0];
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
        InputData[index].LoadData(VariantId,*CurrentRecordFromStudy[index]);
    }
}


void MetaMinimac::MetaImputeCurrentBuffer()
{
    savvy::variant out_record;
    for(int VariantId=0; VariantId<BufferNoVariants; VariantId++)
    {
        CurrentVariant = &BufferVariantList[VariantId];
        CreateMetaImputedData(VariantId);

        out_record = savvy::variant(
            CurrentVariant->chr,
            CurrentVariant->bp,
            CurrentVariant->refAlleleString,
            {CurrentVariant->altAlleleString},
            CurrentVariant->name,
            savvy::typed_value::missing_value<float>(),
            {"PASS"});

        CreateInfo(out_record);
        SetMetaImputedData(out_record);
        if (!metaDoseOut->write(out_record))
        {
            std::cout << "Error: failed to write meta dosage record\n";
            exit(1);
        }
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
    Dosage = float(std::int16_t(Dosage * BinScalar + 0.5f)) / BinScalar; // bin dosages instead of manual VCF printing using %.3f

    CurrentMetaImputedDosage[Sample] = Dosage;
    CurrentHapDosageSum += Dosage;
    CurrentHapDosageSumSq += Dosage * Dosage;
}


void MetaMinimac::SetMetaImputedData(savvy::variant& out_var)
{
    std::size_t stride = CurrentMetaImputedDosage.size() / NoSamples;

    if (myUserVariables.GT || myUserVariables.HDS || myUserVariables.DS)
    {
        sparse_dosages_.assign(CurrentMetaImputedDosage.begin(), CurrentMetaImputedDosage.end());

        if (myUserVariables.GT)
        {
            sparse_gt_.assign(sparse_dosages_.value_data(), sparse_dosages_.value_data() + sparse_dosages_.non_zero_size(), sparse_dosages_.index_data(), sparse_dosages_.size(), [](float v)
            {
                if (savvy::typed_value::is_end_of_vector(v))
                    return savvy::typed_value::end_of_vector_value<std::int8_t>();
                return std::int8_t(v < 0.5f ? 0 : 1);
            });
            out_var.set_format("GT", sparse_gt_);
        }

        if (myUserVariables.HDS)
        {
            out_var.set_format("HDS", sparse_dosages_);
        }

        if (myUserVariables.DS)
        {
            savvy::stride_reduce(sparse_dosages_, sparse_dosages_.size() / NoSamples, savvy::plus_eov<float>());
            out_var.set_format("DS", sparse_dosages_);
        }
    }

    if (myUserVariables.GP)
    {
        if (stride == 1)
        {
            // All samples are haploid
            dense_float_vec_.resize(NoSamples * 2);
            for (std::size_t i = 0; i < NoSamples; ++i)
            {
                std::size_t dest_idx = i * 2;
                dense_float_vec_[dest_idx] = 1.f - CurrentMetaImputedDosage[i];
                dense_float_vec_[dest_idx + 1] = CurrentMetaImputedDosage[i];
            }
        }
        else if (stride == 2)
        {
            dense_float_vec_.resize(NoSamples * 3);
            for (std::size_t i = 0; i < NoSamples; ++i)
            {
                std::size_t src_idx = i * 2;
                std::size_t dest_idx = i * 3;
                float x = CurrentMetaImputedDosage[src_idx];
                float y = CurrentMetaImputedDosage[src_idx + 1];
                if (savvy::typed_value::is_end_of_vector(y))
                {
                    // haploid
                    dense_float_vec_[dest_idx] = 1.f - x;
                    dense_float_vec_[dest_idx + 1] = x;
                    dense_float_vec_[dest_idx + 2] = y;
                }
                else
                {
                    // diploid
                    dense_float_vec_[dest_idx] = (1.f - x) * (1.f - y);
                    dense_float_vec_[dest_idx + 1] = x * (1.f - y) + y * (1.f - x);
                    dense_float_vec_[dest_idx + 2] = x * y;
                }
            }
        }

        out_var.set_format("GP", dense_float_vec_);
    }

    if (myUserVariables.SD)
    {
        dense_float_vec_.resize(NoSamples);
        if (stride == 1)
        {
            // All samples are haploid
            for (std::size_t i = 0; i < NoSamples; ++i)
            {
                dense_float_vec_[i] = CurrentMetaImputedDosage[i] * (1.f - CurrentMetaImputedDosage[i]);
            }

            out_var.set_format("SD", dense_float_vec_);
        }
        else if (stride == 2)
        {
            for (std::size_t i = 0; i < CurrentMetaImputedDosage.size(); i += 2)
            {
                float x = CurrentMetaImputedDosage[i];
                float y = CurrentMetaImputedDosage[i + 1];
                if (savvy::typed_value::is_end_of_vector(y)) // haploid
                    dense_float_vec_[i / 2] = x * (1.f - x);
                else // diploid
                    dense_float_vec_[i / 2] = x * (1.f - x) + y * (1.f - y);
            }

            out_var.set_format("SD", dense_float_vec_);
        }
        else
        {
            // TODO: suppress error excessive error messages
            std::cerr << "Error: only haploid and diploid samples are supported when generating SD\n";
        }
    }
}


void MetaMinimac::SetMetaWeightVecs(std::vector<std::vector<float>>& wt_vecs)
{
    for(int id=0; id<EndSamId-StartSamId; id++)
    {
        if(InputData[0].SampleNoHaplotypes[StartSamId+id]==2)
        {
            vector<double>& WeightsHap1 = (*CurrWeights)[2*id];
            vector<double>& WeightsHap2 = (*CurrWeights)[2*id+1];
            double WeightSumHap1 = 0.0, WeightSumHap2 = 0.0;
            for(int i=0; i<NoInPrefix; i++)
            {
                WeightSumHap1 += WeightsHap1[i];
                WeightSumHap2 += WeightsHap2[i];
            }

            for(int i=0; i<NoInPrefix; i++)
            {
                wt_vecs[i][2*id] = WeightsHap1[i]/WeightSumHap1;
                wt_vecs[i][2*id+1] = WeightsHap2[i]/WeightSumHap2;
            }
        }
        else
        {
            vector<double>& ThisCurrWeights = (*CurrWeights)[2*id];
            double WeightSum = 0.0;
            for(int i=0; i<NoInPrefix; i++)
                WeightSum += ThisCurrWeights[i];

            for(int i=0; i<NoInPrefix; i++)
                wt_vecs[i][2*id] = ThisCurrWeights[i]/WeightSum;
        }
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

void MetaMinimac::CreateInfo(savvy::variant& rec)
{
    if(!myUserVariables.infoDetails)
        return;

    double hapSum = CurrentHapDosageSum, hapSumSq = CurrentHapDosageSumSq;
    double freq = hapSum*1.0/NoHaplotypes;
    double maf = (freq > 0.5) ? (1.0 - freq) : freq;
    double rsq = 0.0, evar = freq*(1-freq), ovar = 0.0;
    if (NoHaplotypes > 2 & (hapSumSq - hapSum * hapSum / NoHaplotypes) >0 )
    {
        ovar = (hapSumSq - hapSum * hapSum / NoHaplotypes)/ NoHaplotypes;
        rsq = ovar / (evar + 1e-30);
    }

    rec.set_info("NST", CurrentVariant->NoStudiesHasVariant);
    for(int i=0; i<CurrentVariant->NoStudiesHasVariant; i++)
    {
        rec.set_info("S" + std::to_string(CurrentVariant->StudiesHasVariant[i]+1), std::vector<std::int8_t>());
    }

    if(CurrentVariant->name == CommonGenotypeVariantNameList[NoCommonVariantsProcessed-1])
    {
        rec.set_info("TRAINING", std::vector<std::int8_t>());
    }

    rec.set_info("AF", float(freq));
    rec.set_info("MAF", float(maf));
    rec.set_info("R2", float(rsq));
}


void MetaMinimac::SetWeightVariantInfo(savvy::variant& dest)
{
    variant& tempVariant = CommonTypedVariantList[NoCommonVariantsProcessed];
    dest = savvy::site_info(
        tempVariant.chr,
        tempVariant.bp,
        tempVariant.refAlleleString,
        {tempVariant.altAlleleString},
        tempVariant.name,
        savvy::typed_value::missing_value<float>(),
        {"PASS"});
}

bool MetaMinimac::AppendtoMainWeightsFile()
{
    if (batchNo < 1)
        return std::cout << "\nError: not enough partial weight files\n", false;

    cout << "\n Appending to final output weight file : " << myUserVariables.outfile + ".metaWeights" + (myUserVariables.gzip ? ".gz" : "") <<endl;
    int start_time = time(0);

    std::vector<std::string> sample_ids;
    std::list<savvy::reader> partial_weight_files;
    for(int i=1; i<=batchNo; i++)
    {
        stringstream ss;
        ss << (i);
        string PartialWeightFileName(myUserVariables.outfile);
        PartialWeightFileName += ".metaWeights.part."+(string)(ss.str())+(myUserVariables.gzip ? ".gz" : "");
        partial_weight_files.emplace_back(PartialWeightFileName);
        sample_ids.insert(sample_ids.end(), partial_weight_files.back().samples().begin(), partial_weight_files.back().samples().end());
    }

    savvy::writer merged_out_file(std::string(myUserVariables.outfile.c_str()) + ".metaWeights" + (myUserVariables.gzip ? ".gz" : ""), savvy::file::format::vcf, partial_weight_files.front().headers(), sample_ids, myUserVariables.gzip ? 6 : 0);

    std::vector<std::vector<float>> wt_vecs(NoInPrefix);
    std::vector<float> tmp_vec;
    savvy::variant record;

    NoCommonVariantsProcessed = 0;
    for(int i=0; i<NoCommonTypedVariants; i++)
    {
        for (auto jt = wt_vecs.begin(); jt != wt_vecs.end(); ++jt)
            jt->clear();

        for(auto jt = partial_weight_files.begin(); jt != partial_weight_files.end(); ++jt)
        {
            if (!jt->read(record))
                return std::cout << "\nError: not enough partial weight records\n", false;

            for (int k = 0; k < wt_vecs.size(); ++k)
            {
                record.get_format("WT" + std::to_string(k + 1), tmp_vec);
                wt_vecs[k].insert(wt_vecs[k].end(), tmp_vec.begin(), tmp_vec.end());
            }
        }

        for (int j = 0; j < wt_vecs.size(); ++j)
            record.set_format("WT" + std::to_string(j + 1), wt_vecs[j]);

        merged_out_file << record;
        NoCommonVariantsProcessed++;
    }

    for(int i=1;i<=batchNo;i++)
    {
        stringstream ss;
        ss << (i);
        string tempFileIndex(myUserVariables.outfile);
        tempFileIndex += ".metaWeights.part."+(string)(ss.str())+(myUserVariables.gzip ? ".gz" : "");
        remove(tempFileIndex.c_str());
    }

    int tot_time = time(0) - start_time;
    cout << " -- Successful (" << tot_time << " seconds) !!!" << endl;
    return true;
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