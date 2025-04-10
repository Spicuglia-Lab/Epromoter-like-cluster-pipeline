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

2: Expression Data file required: induced genes list (gene names and log2 fold change between stimualated and control)
   
   A tab separated file with at least two columns with below header line.
```
Example:
        gene	log2FC
        NUDT4	2.01
        AHSA1	1.98
        EPB42	0.31
```

3: PeakFile: TF_peak.bed (transcrition factor binding peaks from ChIP-seq data)
   
   Tab separated file with 3 columns: ['chr','start','end']. No header line in peak file.
   ```
   Example:
        chr1	5639	5705
        chr2	3048	3049
        chr1	1150	1509
```

4: TAD File: TAD.bed (OPTIONAL)
   
   Tab separated file with 4 columns: ['chr','start','end']. No header line in peak file.
   ```
   Example:
        chr3	5639	5705
        chr1	3048	3049
        chr2	1150	1509
```


### Output
Output folder name is: 
```
<TF Name>_<Induced Condition>_Promoter_<bps upstream from TSS>_<bps downstream from TSS>_dist<clustering distance_criteria>
```
Output files: Cluster_summary_stats.txt, gene_promoter_cluster_tf_binding.txt, Epromoter_cluster.txt

The information in Cluster_summary_stats.txt like:
```
Number_of_induced_genes: 675
Number_of_clusters: 75
Number_of_Epromoter_clusters: 7
```
The gene_promoter_cluster_tf_binding.txt is induced gene clusters.

The Epromoter_cluster.txt is Epromoter clusters.

The cluster file inlcude 9 columns, like: 
```
cluster number, number of genes per cluster, number of promoters per cluster, genes in cluster, number of promoter binding TFs in cluster, TF_binding_gene, num_of_TF_binding_genes, chr_promoter, promoter_start, promoter_end
9	2	4	PMVK,CKS1B	1	PMVK	1	chr1	154908195	154910195
14	2	6	STIP1,GPR137	1	STIP1	1	chr11	63952642	63954642
18	2	7	C11orf52,CRYAB	1	CRYAB	1	chr11	111781494	111783494
21	2	5	BAZ2A,PTGES3	1	PTGES3	1	chr12	57081192	57083192
24	2	3	SCARB1,UBC	1	UBC	1	chr12	125398196	125400196
27	2	3	SNAP23,LRRC57	1	SNAP23	1	chr15	42786831	42788831
28	2	3	EPB42,TGM7	1	EPB42	1	chr15	43512323	43514323
```

### reference

