# MetaMinimac2
MetaMinimac2 is an efficient tool to combine genotype data imputed against multiple reference panels.

## Installation
Prerequisite: [cget](https://cget.readthedocs.io/en/latest/index.html)>=0.1, [cmake](https://cmake.org)>=3.2
```bash
git clone https://github.com/yukt/MetaMinimac2.git
cd MetaMinimac2
bash install.sh
```

## Usages
The first step is to impute _**pre-phased**_ target haplotypes against reference panels separately using [Minimac4](https://github.com/statgen/Minimac4) with `--meta` option on, which will generate a `.empiricalDose.vcf.gz` file (which is required for MetaMinimac2) in addition to the `.dose.vcf.gz` file. Please see http://genome.sph.umich.edu/wiki/Minimac4 for detailed documentation for Minimac4.

```bash
minimac4 --refHaps refPanelA.m3vcf \
         --haps targetStudy.vcf \
         --prefix PanelA.imputed \
         --meta
         
minimac4 --refHaps refPanelB.m3vcf \
         --haps targetStudy.vcf \
         --prefix PanelB.imputed \
         --meta
```

The second step is to integrate the imputed results using MetaMinimac2.

```bash
MetaMinimac2 -i PanelA.imputed:PanelB.imputed -o A_B.meta.testrun
```

## Options
```
-i, --input  <prefix1:prefix2 ...>  (Required) Colon-separated prefixes of input data to meta-impute
-o, --output <prefix>               (Required) Output prefix [MetaMinimac.Output.Prefix]
-f, --format <string>               Comma-separated output FORMAT tags [GT,DS,HDS]
-p, --skipPhasingCheck              OFF by default. If ON, program will skip phasing consistency check before analysis. 
-s, --skipInfo                      OFF by default. If ON, the INFO fields are removed from the output file
-n, --nobgzip                       OFF by default. If ON, output files will NOT be bgzipped
-w, --weight                        OFF by default. If ON, weights will be saved in [MetaMinimac.Output.Prefix].metaWeights(.gz)
-l, --log                           OFF by default. If ON, log will be written to $prefix.logfile
-h, --help                          OFF by default. If ON, detailed documentation on options and usage will be displayed
```

## Important Notes
_**Phasing must be consistent across input files !!!**_ If the phasing is different between input imputed files, the resulting meta dosage (which is supposed to be a weighted average of imputed dosages on the same haplotype) will be messed up. Therefore, we highly recommend that users always keep `--skipPhasingCheck` OFF  to avoid any risk.

The best practice is to do the phasing first and use the same pre-phased vcf file for imputation against different reference panels. 

When the phasing step is performed by imputation server (e.g. [Michigan Imputation Server](imputationserver.sph.umich.edu), [TOPMed Imputation Server](https://imputation.biodatacatalyst.nhlbi.nih.gov)) where the phased vcf file is not provided as output, please consider using the output `.empiricalDose.vcf.gz` file for imputation against other reference panels. Note that `.empiricalDose.vcf.gz` file only contains genotyped sites included in the corresponding reference panel, so the imputation quality using empiricalDose.vcf might be worse than that using a pre-phased vcf including all genotyped sites.

## Output Files
The meta-imputed result will be saved in `[MetaMinimac.Output.Prefix].metaDose.vcf.gz`.

When `--weight` is ON, the weights for meta-imputation will be saved in `[MetaMinimac.Output.Prefix].metaWeights.gz`. The weight file is also in VCF format, which is good for individual filtering by [vcftools](https://vcftools.github.io) or [bcftools](http://samtools.github.io/bcftools/bcftools.html). It contains one special format `WT` only in the genotype fields, and the values are in the format as follows:

`[weight on PanelA hap1], [weight on PanelB hap1] | [weight on PanelA hap2], [weight on PanelB hap2]`
