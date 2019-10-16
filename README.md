# MetaMinimac2
MetaMinimac2 is an efficient tool to combine genotype data imputed against multiple reference panels.

## Installation
```bash
git clone https://github.com/yukt/MetaMinimac2.git
cd MetaMinimac2
bash install.sh
```

## Usages
The first step is to impute against reference panels separately using Minimac4 with `--meta` option on. Please see http://genome.sph.umich.edu/wiki/Minimac4 for detailed documentation for Minimac4.

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

The second step is to combine the imputed results using MetaMinimac2.

```bash
MetaMinimac2 -i PanelA.imputed:PanelB.imputed -o A_B.meta.testrun
```

## Options
```
-i, --input  <prefix1:prefix2 ...>  Colon-separated prefixes of input data to meta-impute
-o, --output <prefix>               Output prefix [MetaMinimac.Output]
-f, --format <string>               Comma-separated FORMAT tags [GT,DS,HDS]
-s, --skipInfo                      If ON, the INFO fields are removed from the output file
-n, --nobgzip                       If ON, output files will NOT be bgzipped
-w, --weight                        If ON, weights will be saved in $prefix.metaWeights(.gz)
-l, --log                           If ON, log will be written to $prefix.logfile
-h, --help                          If ON, detailed help on options and usage
```

