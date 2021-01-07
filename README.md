# Epromoter-like-cluster-pipeline
Pipeline to identify Epromoter like clusters using RNA-seq and ChIP-seq data


Input File:
  input_data/
  
    1. Name of genome annotation file: mm10.refGene.txt.gz
    
        wget -c -O mm9.refGene.txt.gz http://hgdownload.soe.ucsc.edu/goldenPath/mm9/database/refGene.txt.gz
        wget -c -O mm10.refGene.txt.gz http://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/refGene.txt.gz
        wget -c -O hg19.refGene.txt.gz http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz
        wget -c -O hg38.refGene.txt.gz http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz

    2: Expression Data file required: 'induced_genes.expression'
    A tab separated file with at least two columns with below header line.    
    Example:
        gene\tlog2FC
        NUDT4\t2.01
        AHSA1\t1.98
        EPB42\t0.31\n
        
    3: PeakFile: peak.narrowPeak
    Tab separated file with 4 columns: ['chr','start','end','macs score']. No header line in peak file.    
    Example:
        chr1\t5639\t5705\t310
        chr2\t3048\t3049\t572.8
        chr1\t1150\t1509\t620.32\n
        
    4: TAD File: tad.bed (OPTIONAL)    
    Tab separated file with 4 columns: ['chr','start','end']. No header line in peak file.    
    Example:
        chr3\t5639\t5705
        chr1\t3048\t3049
        chr2\t1150\t1509\n

    NOTE:
    Macs Cutoff Filter Score must be defined by user
    No filter score = 0\n

    OUTPUT FILES:
    Output folder name is: <TF Name>_<Induced Condition>_Promoter_<bps upstream from TSS>_<bps downstream from TSS>_dist<clustering distance_criteria>

    exp_bedgraph/
       expression.mm10.bedgraph: Expression Bedgraph file
       expression.mm10.induced_genes: all the induced genes in refSeq
       induced_genes_not_in_refseq_data.exp: genes not found in the refSeq annotation files
       
    tf_analysis/
       bubble_chart_gene_promoter_cluster_tf_binding_frequency.table
       cluster_gene_frequency_at_100kb_data.txt
       cluster_genes_frequency_at_100kb_plot.png
       gene_least_distant_barplot.png
       gene_promoter_cluster_tf_binding.table
       least_distant_induced_genes.distance
       observed_vs_randomly_selected_genes_distribution.png
       randomly_selected_1_times.184_induced_genes.least_distance
       TF_binding_gene_promoter_frequency_bubble_plot.png
