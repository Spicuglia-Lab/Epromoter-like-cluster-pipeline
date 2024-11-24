# Epromoter-like-cluster-pipeline
Pipeline to identify Epromoter like clusters using RNA-seq and ChIP-seq data

### Dependent envioroment
System envioroment: Linux, python>3.7

python libraries:
```
datetime, pandas, plotnine, numpy, os, subprocess, seaborn, itertools, matplotlib,scipy, pybedtools
```
### Required Input files:
All the required input files must be avilable in folder: 'input_data'

1: RefSeq GTF file: 'hg38.refGene.txt.gz'

   Reference gene annotation files were downloaded from the UCSC Genome Browser: 
   
   human (hg19: https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz, hg38: http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz)
   
   mouse (mm9: https://hgdownload.soe.ucsc.edu/goldenPath/mm9/database/refGene.txt.gz, mm10: https://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/refGene.txt.gz)

   Format:
```
[bin,name,chrom,strand,txStart,txEnd,cdsStart,cdsEnd,exonCount,exonStarts,exonEnds,score,name2,cdsStartStat,cdsEndStat,exonFrames]
```

2: Expression Data file required: induced genes expression (gene names and log2 fold change between stimualated and control)
   
   A tab separated file with at least two columns with below header line.
```
Example:
        gene\tlog2FC
        NUDT4\t2.01
        AHSA1\t1.98
        EPB42\t0.31\n
```

3: PeakFile: TF_peak.bed (transcrition factor binding peaks from ChIP-seq data)
   
   Tab separated file with 3 columns: ['chr','start','end']. No header line in peak file.
   ```
   Example:
        chr1\t5639\t5705
        chr2\t3048\t3049
        chr1\t1150\t1509\n
```

4: TAD File: TAD.bed (OPTIONAL)
   
   Tab separated file with 4 columns: ['chr','start','end']. No header line in peak file.
   ```
   Example:
        chr3\t5639\t5705
        chr1\t3048\t3049
        chr2\t1150\t1509\n
```


### Output
Output folder name is: 
```
<TF Name>_<Induced Condition>_Promoter_<bps upstream from TSS>_<bps downstream from TSS>_dist<clustering distance_criteria>
```
Output files: gene_promoter_cluster_tf_binding.txt, Epromoter_cluster.txt

Each output file inlcude 8 columns: 
```
cluster number, number of genes per cluster, number of promoters per cluster, genes in cluster, number of promoter binding TFs in cluster, Epromoter_gene, chr, TSS
```

### reference

