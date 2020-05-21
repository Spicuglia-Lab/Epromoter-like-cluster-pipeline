# Epromoter-like-cluster-pipeline
Pipeline to identify Epromoter like clusters using RNA-seq and ChIP-seq data


Input File:
  input_data/
  
    1. Name of genome annotation file: mm10.refGene.txt.gz
    
        wget -c -O mm9.refGene.txt.gz http://hgdownload.soe.ucsc.edu/goldenPath/mm9/database/refGene.txt.gz
        wget -c -O mm10.refGene.txt.gz http://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/refGene.txt.gz
        wget -c -O hg19.refGene.txt.gz http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz
        wget -c -O hg38.refGene.txt.gz http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz

    2. Name of RNA expression file: induced_genes.expression
        <tab separated file: with header line>
        <columns: Gene	log2FC>
        Gene	log2FC	
        Xkr4		-1.1665	
        Mrpl15	 0.3167	
        Lypla1	 0.8388	

    3. Name of TF: Fos
    4. Name of TF chipseq bed file: peak_file.narrowPeak
        <tab separated file; No header line>
        <columns: chr; start; end; macs_score>
        chr1	1265	1562	168	
        chr1	9235	9467	178
        chrX	6785	6956	98
        chrY	3641	3758	517
    5. MACS Score filterscore: 10 <filter score can be defined by user>

Output Files:

    output/exp_bedgraph/
       Vierbuchen_2017_RefSeq_4hr_vs_ctl.mm10.bedgraph: Expression Bedgraph file
       Vierbuchen_2017_RefSeq_4hr_vs_ctl.mm10.induced_genes
       output/tf_analysis/
       combined_gene_5p_tf_peak_calls_data.bed
       gene_5p_tf_peak_calls_bed_tool_clust.bed
       gene_promoter_cluster.table
       induced_genes_5p_3p_coordinates.bed
       induced_genes_5p_coordinates.bed
       induced_genes_bed_tool_clust.bed
       induced_gene_tss

       induced_genes_final_cluster.bed
       cluster_gene_frequency_at_100kb_data.txt
       cluster_genes_frequency_at_100kb_plot.png
       dist_all_induced_genes.distance
       least_distant_induced_genes.distance
       gene_least_distant_barplot.png
       gene_tss_peak.table
       gene_promoter_cluster_tf_binding.table
       bubble_chart_gene_promoter_cluster_tf_binding_frequency.table
       TF_binding_gene_promoter_frequency_bubble_plot.png
