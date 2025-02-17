# coding: utf-8
import datetime
import pandas as pd
from pandas.api.types import CategoricalDtype
from plotnine import *
import numpy as np
import os, sys
import subprocess
import seaborn as sns
from itertools import combinations
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
from scipy.stats import ks_2samp
import random
import pybedtools
import warnings
warnings.filterwarnings("ignore")

os.chdir(os.getcwd())#moving to current working directory
cwd = os.getcwd()+'/'
start_time = datetime.datetime.now()
print("\nThe job started at: ",datetime.datetime.now())

input_files ="""
=============================================INSTRUCTIONS=============================================\n
Required python libraries:
                datetime, pandas, plotnine, numpy, os,
                subprocess, seaborn, itertools, matplotlib,
                scipy, pybedtools\n
All the required input files must be avilable in folder: 'input_data'\n
Required Input files:
1: RefSeq GTF file: 'hg38.refGene.txt.gz'
   [bin,name,chrom,strand,txStart,txEnd,cdsStart,cdsEnd,exonCount,exonStarts,exonEnds,score,name2,cdsStartStat,cdsEndStat,exonFrames]\n
2: Expression Data file required: 'induced_genes.expression'
   A tab separated file with at least two columns with below header line.
   Example:
        gene\tlog2FC
        NUDT4\t2.01
        AHSA1\t1.98
        EPB42\t0.31\n
3: PeakFile: peak.bed
   Tab separated file with 3 columns: ['chr','start','end']. No header line in peak file.
   Example:
        chr1\t5639\t5705
        chr2\t3048\t3049
        chr1\t1150\t1509\n
4: TAD File: tad.bed (OPTIONAL)
   Tab separated file with 4 columns: ['chr','start','end']. No header line in peak file.
   Example:
        chr3\t5639\t5705
        chr1\t3048\t3049
        chr2\t1150\t1509\n


======================================================================================================\n
"""
print(input_files)

class InstallPythonLibraries:
    def __init__(self):
        '''Check the Required Python Libraries'''
        list = ['datetime','pandas','plotnine','numpy','os','subprocess','seaborn','itertools','matplotlib','scipy','pybedtools']

        print('Status of Python Libraries: \n')
        for lib in list:
            if lib in sys.modules:
                print('\t{} is already installed...'.format(lib))
            else:
                print('\n\t{} is being installed...'.format(lib))
                try:
                    subprocess.call('pip install '+lib, shell=True)
                except:
                    subprocess.call('pip3 install '+lib, shell=True)

class Expression_bedgraph:
    '''Expression_bedgraph'''
    def __init__(self):
        '''Read File'''
    def induced_genes(self,expression_file,gtf):
        os.chdir(cwd)
        print("##### Generation of bedgraph from the induced genes...")
        self.expression_file = expression_file
        self.gtf = gtf

        expression = pd.read_csv('input_data/'+self.expression_file,header=0,sep='\t',usecols=['gene','log2FC'])
        expression.columns = ['name2','log2FoldChange']
        induced_gene_exp_data = expression[(expression.log2FoldChange > induced_gene_log2FC)]

        #writing induced_genes
        induced_gene_exp_data.to_csv(output_folder_name+'/exp_bedgraph/'+self.expression_file.split('.')[0]+'.'+self.gtf.split('.')[0]+'.induced_genes',sep="\t",header=False,index=False)
        # print('  Files generated:\n\t'+output_folder_name+'/exp_bedgraph/'+self.expression_file.split('.')[0]+'.'+self.gtf.split('.')[0]+'.induced_genes')

        refgene_col = ["bin","name","chrom","strand","txStart","txEnd","cdsStart","cdsEnd","exonCount","exonStarts","exonEnds","score","name2","cdsStartStat","cdsEndStat","exonFrames"]
        refgene_data = pd.read_csv('input_data/'+self.gtf,sep="\t",header=None,names=refgene_col,usecols=["name","chrom","strand","txStart","txEnd","name2"])
        refgene_data = refgene_data[~refgene_data.chrom.str.contains('_')]
        refgene_data['txn_length'] = abs(refgene_data['txStart']-refgene_data['txEnd'])
        refgene_data_txn_length = refgene_data.groupby('name2', as_index=False).apply(pd.DataFrame.sort_values, 'txn_length', ascending = False).reset_index(drop=True)

        """ 
        clean_refgene_data = refgene_data_txn_length[['chrom','strand','txStart','txEnd','name2','txn_length']].copy()
        clean_refgene_data['tss'] = clean_refgene_data.groupby(['name2']).cumcount()+1#add tss number
        strand = clean_refgene_data['strand'] == '-'; clean_refgene_data.loc[strand, ['txStart','txEnd']] = clean_refgene_data.loc[strand,['txEnd','txStart']].values
        #WRITE ALL txSTART/TSS IN THE GENOME
        clean_refgene_data.to_csv(output_folder_name+'/tf_analysis/'+os.path.splitext(self.gtf)[0]+'.all.tss.txt',sep="\t",header=True,index=False)
        # print('\t'+output_folder_name+'/tf_analysis/'+os.path.splitext(self.gtf)[0]+'.all.tss.txt')

        #LEAST DISTANT GENES IN THE GENOME
        clean_refgene_data = clean_refgene_data[clean_refgene_data.tss == 1]#extract genes with tss 1
        data = clean_refgene_data[['chrom','txStart','name2','strand']].copy()
        data.rename(columns={'chrom':'chr','name2': 'gene'}, inplace=True)

        all_distance_combination_genes = pd.DataFrame()
        for name, group in data.groupby('chr'):
            gene_pairs,dist_pairs = list(combinations(group.gene,2)),list(combinations(group['txStart'],2))
            gene_combination,dist_combination = pd.DataFrame(gene_pairs,columns=['gene1','gene2']),pd.DataFrame(dist_pairs,columns=['gene1_dist','gene2_dist'])
            df = pd.concat([gene_combination,dist_combination],axis=1)
            df['distance'] = abs(df.gene1_dist-df.gene2_dist)
            all_distance_combination_genes = all_distance_combination_genes.append(df[['gene1','gene2','distance']],ignore_index=True,sort=False)

        all_distance_combination_genes['distance'] = all_distance_combination_genes['distance'].astype('int64')
        all_distance_combination_genes = all_distance_combination_genes[all_distance_combination_genes['gene1'] != all_distance_combination_genes['gene2']]
        gene_least_distant_genes_indices = all_distance_combination_genes.groupby('gene1')['distance'].idxmin()
        gene_least_distant_table = all_distance_combination_genes.loc[gene_least_distant_genes_indices]
        #print(gene_least_distant_table.head(4),"\n",all_distance_combination_genes.head(4));exit(0)
        all_distance_combination_genes.to_csv(output_folder_name+"/tf_analysis/genome.all_distance_combination_genes.distance", sep='\t',index=False)
        gene_least_distant_table.to_csv(output_folder_name+"/tf_analysis/genome.least_distant_genes.distance", sep='\t',index=False)
        # print("\t"+output_folder_name+"/tf_analysis/"+self.gtf[:-11]+"least_distant_genes.distance\n\t"+output_folder_name+"/tf_analysis/"+self.gtf[:-11]+"all_combination_distant_genes.distance")
         """
        
        induced_gene_not_in_refgene_data = induced_gene_exp_data[~induced_gene_exp_data['name2'].isin(refgene_data_txn_length['name2'])].copy()
        induced_gene_not_in_refgene_data.columns = ['gene','log2FoldChange']
        induced_gene_not_in_refgene_data.drop_duplicates(subset='gene', keep="first",inplace=True)
        induced_gene_not_in_refgene_data.to_csv(output_folder_name+'/exp_bedgraph/induced_genes_not_in_refseq_data.exp',sep="\t",header=True,index=False)

        induced_refgene_data = refgene_data_txn_length[refgene_data_txn_length['name2'].isin(induced_gene_exp_data['name2'])].copy()
        induced_refgene_data['tss'] = induced_refgene_data.groupby(['name2']).cumcount()+1
        exp_with_coordinate = induced_gene_exp_data.merge(induced_refgene_data, on='name2')
        exp_with_coordinate_all_tss = exp_with_coordinate[["name2", "log2FoldChange", "chrom","strand", "txStart","txEnd","tss"]].copy()

        exp_bedgraph = pd.DataFrame()
        for item,row in exp_with_coordinate_all_tss.iterrows():
            if(row['strand'] == '+'):
                exp_bedgraph = exp_bedgraph.append(pd.Series([row['chrom'],row['txStart'],row['txStart']+500,round(row['log2FoldChange'],4),row['name2']]),ignore_index=True)
            if(row['strand'] == '-'):
                exp_bedgraph = exp_bedgraph.append(pd.Series([row['chrom'],row['txEnd'],row['txEnd']-500,round(row['log2FoldChange'],4),row['name2']]),ignore_index=True)

        exp_bedgraph.columns = ['chr','sTSS','TSSto500bps','log2FC','geneName']
        exp_bedgraph = exp_bedgraph.astype({"sTSS": int, "TSSto500bps": int})
        exp_bedgraph.to_csv(output_folder_name+'/exp_bedgraph/'+self.expression_file.split('.')[0]+'.'+self.gtf.split('.')[0]+'.bedgraph',sep="\t",header=False,index=False)
        # print('\t'+output_folder_name+'/exp_bedgraph/'+self.expression_file.split('.')[0]+'.'+self.gtf.split('.')[0]+'.bedgraph')
        # print('\t'+output_folder_name+'/exp_bedgraph/induced_genes_not_in_refseq_data.exp\n')

class Extract_induced_genes_coordinates:
    def __init__(self):
        '''Read File'''
    def induced_genes(self,genelist,coordinatefile):
        print("##### Extraction induced genes coordinates from RefSeq Database file: 'RefGene.txt' ...")
        self.coordinatefile = coordinatefile
        self.genelist = genelist

        induced_gene_data = pd.read_csv(output_folder_name+'/exp_bedgraph/'+self.genelist.split('.')[0]+'.'+self.coordinatefile.split('.')[0]+'.induced_genes',sep="\t",header=0)
        induced_gene_data.columns = ['gene','log2FC']

        refgene_col = ["bin","name","chrom","strand","txStart","txEnd","cdsStart","cdsEnd","exonCount","exonStarts","exonEnds","score","name2","cdsStartStat","cdsEndStat","exonFrames"]
        refgene_data = pd.read_csv('input_data/'+self.coordinatefile, sep="\t", header=None,names=refgene_col)
        refgene_data = refgene_data[~refgene_data.chrom.str.contains('_')]
        refgene_data['txn_length'] = abs(refgene_data['txStart']-refgene_data['txEnd'])
        refgene_data_txn_length = refgene_data.groupby('name2', as_index=False).apply(pd.DataFrame.sort_values, 'txn_length', ascending = False).reset_index(drop=True)

        induced_refgene_data = refgene_data_txn_length[refgene_data_txn_length['name2'].isin(induced_gene_data['gene'])].copy()
        induced_refgene_data['tss'] = induced_refgene_data.groupby(['name2']).cumcount()+1
        induced_gene_tss = pd.DataFrame(induced_refgene_data, columns=['name','chrom','strand','txStart','txEnd','name2','tss'])
        induced_gene_tss.to_csv(output_folder_name+'/tf_analysis/induced_gene_tss',sep="\t",index=False)
        # print('  Induced genes tss coordinates data file is generated: \n\t'+output_folder_name+'/tf_analysis/induced_gene_tss\n')
        return induced_gene_tss


class Inducedgenesclustering:
    def __init__(self,induced_gene):
        '''Clustering on the basis of 5' coordinate of induced genes'''
        #input: induced_gene_tss
        self.induced_gene = induced_gene

    def induced_genes_clusteringByDistance(self):
        print("##### Induced genes clustering by distance...")

        self.induced_gene = pd.read_csv(output_folder_name+'/tf_analysis/'+self.induced_gene,sep="\t",header=0,usecols = ['chrom','strand','txStart','txEnd','name2','tss'])
        os.chdir(cwd)
        five_prime_coordinate_gene = pd.DataFrame()

        dtypes = {'chr':'str','strand':'str','5p_gene1':'int','5p_gene2':'int','3p_gene':'int','gene':'str'}
        five_prime_coordinate_gene = pd.DataFrame(columns=['chr','5p_gene1','5p_gene2','3p_gene','gene'])#gene coordinate

        for index, row in self.induced_gene.iterrows():
            if(row['strand'] == '+' and row['tss'] == 1):
                five_prime_coordinate_gene=five_prime_coordinate_gene.append(pd.Series([row['chrom'],row['txStart'],
                row['txStart'],row['txEnd'],row['name2']],index=['chr','5p_gene1','5p_gene2','3p_gene','gene']),ignore_index=True)

            elif(row['strand'] == '-' and row['tss']==1):
                five_prime_coordinate_gene=five_prime_coordinate_gene.append(pd.Series([row['chrom'],row['txEnd'],row['txEnd'],
                row['txStart'],row['name2']],index=['chr','5p_gene1','5p_gene2','3p_gene','gene']),ignore_index=True)

        sort_five_prime_coordinate_gene = five_prime_coordinate_gene.sort_values(['chr','5p_gene1','5p_gene2'])
        five_prime_coordinate = sort_five_prime_coordinate_gene[['chr','5p_gene1','5p_gene2','gene']]

        five_prime_coordinate.to_csv(output_folder_name+'/tf_analysis/induced_genes_5p_coordinates.bed',sep="\t",index=False,header=None)
        sort_five_prime_coordinate_gene.to_csv(output_folder_name+'/tf_analysis/induced_genes_5p_3p_coordinates.bed',sep='\t',index=False,header=None)
        pybedtools.bedtool.BedTool.cluster(pybedtools.bedtool.BedTool(output_folder_name+'/tf_analysis/induced_genes_5p_3p_coordinates.bed'),d=induced_genes_cluster_distance).saveas(output_folder_name+'/tf_analysis/induced_genes_bed_tool_clust.bed')

        final_cluster = pd.read_csv(output_folder_name+'/tf_analysis/induced_genes_bed_tool_clust.bed',sep="\t",header=None,names=['chr','start1','start2','end','gene','clust1'])
        final_cluster = final_cluster[final_cluster.duplicated(subset='clust1',keep=False)]
        order_final_cluster = pd.DataFrame(final_cluster['clust1'].unique(),columns=['clust1'])
        order_final_cluster['clust2'] = order_final_cluster.index+1

        final_cluster = pd.merge(final_cluster,order_final_cluster,on='clust1')
        final_cluster.to_csv(output_folder_name+'/tf_analysis/induced_genes_final_cluster.bed',sep="\t",header=None,index=False)
        # print('  Induced Gene Cluster file is generated:\n\t'+output_folder_name+'/tf_analysis/induced_genes_final_cluster.bed\n')

    def induced_genes_clusteringByTAD(self,tad_file):
        print("##### Induced genes clustering by TAD...")

        self.tad_file = tad_file
        tad = pd.read_table('input_data/'+self.tad_file,sep="\s+",header=None,names=['chr','start','end'])
        induced_gene_data = pd.read_csv(output_folder_name+'/tf_analysis/'+self.induced_gene,sep="\t",header=0,usecols = ['chrom','strand','txStart','txEnd','name2','tss'])
        induced_gene_data['txEnd'],induced_gene_data['txStart'] = np.where(induced_gene_data.strand == '-',[induced_gene_data.txStart,induced_gene_data.txEnd],[induced_gene_data.txEnd,induced_gene_data.txStart])
        induced_gene_data = induced_gene_data[induced_gene_data.tss == 1]
        induced_gene_data = induced_gene_data[['chrom','name2','txStart']]
        induced_gene_data[['chrom','txStart','txStart','name2']].to_csv(output_folder_name+'/tf_analysis/induced_genes_5p_coordinates.bed',sep="\t",header=None,index=False)
        cluster = pd.DataFrame()
        for index, row in tad.iterrows():
            induced_gene_data_temp = induced_gene_data[induced_gene_data.chrom == row.chr]
            induced_gene_data_temp = induced_gene_data_temp.drop_duplicates(subset='txStart',keep='first')
            temp = induced_gene_data_temp[induced_gene_data_temp.txStart.between(row.start,row.end)]
            temp['cluster no'] = index+1
            cluster = cluster.append(temp)

        cluster = cluster[cluster.duplicated(['cluster no'], keep = False)]
        cluster[['chrom','txStart','txStart','name2','cluster no']].to_csv(output_folder_name+'/tf_analysis/induced_genes_bed_tool_clust.bed',sep="\t",header=False,index=False)

        final_cluster = pd.read_csv(output_folder_name+'/tf_analysis/induced_genes_bed_tool_clust.bed',sep="\t",header=None,names=['chr','start','end','gene','clust1'])
        final_cluster = final_cluster[final_cluster.duplicated(subset='clust1',keep=False)]
        order_final_cluster = pd.DataFrame(final_cluster['clust1'].unique(),columns=['clust1'])
        order_final_cluster['clust2'] = order_final_cluster.index+1
        print(final_cluster)
        final_cluster = pd.merge(final_cluster,order_final_cluster,on='clust1')
        print(final_cluster)
        final_cluster.to_csv(output_folder_name+'/tf_analysis/induced_genes_final_cluster.bed',sep="\t",header=None,index=False)
        # print('  Induced Gene Cluster file is generated:\n\t'+output_folder_name+'/tf_analysis/induced_genes_final_cluster.bed\n')

class Genesclusterbarplot:
    def __init__(self):
        '''Save Plot'''
    def genes_clust_genes_frequency(self, filename):
        # print("##### Genes cluster bar plot...")
        self.filename = filename
        cluster_col_names = ['chr','5p_gene1','5p_gene2','gene','clust_no','cluster']

        hist_table_data = pd.read_csv(output_folder_name+'/tf_analysis/'+self.filename, sep="\t", header=None,names=cluster_col_names)
        hist_table = pd.DataFrame(hist_table_data.groupby(['cluster']).size().value_counts().reset_index(name='number of clusters'))

        ax=sns.barplot(x="index",y="number of clusters",data=hist_table,color='blue')
        ax.set_xlabel("genes in cluster",fontsize=14)   #ax.axes.set_title("Title",fontsize=16)
        ax.set_ylabel("frequency of clusters",fontsize=14)
        plt.savefig(output_folder_name+'/tf_analysis/cluster_genes_frequency_at_100kb_plot.png')

        hist_table.to_csv(output_folder_name+'/tf_analysis/cluster_gene_frequency_at_100kb_data.txt',sep="\t",header=['genes in cluster','frequency of clusters'],index=False)
        # print('  Frequency of Induced Genes per Cluster data file is generated:\n\t'+output_folder_name+'/tf_analysis/cluster_gene_frequency_at_100kb_data.txt')
        # print('  Frequency of Induced Genes per Cluster Barplot is generated:\n\t'+output_folder_name+'/tf_analysis/cluster_gene_frequency_at_100kb_barplot.png\n')

class Genesleastdistant:
    def __init__(self):
        '''Save Plot'''

    def genes_least_distnce(self, filename):
        print("##### Induced genes distance from 5p coordinate...")
        self.filename = filename
        data = pd.read_csv(output_folder_name+'/tf_analysis/'+self.filename, sep="\t",header=None,names=['chr','5p_gene','5p_gene1','gene'])

        all_distance_combination_genes = pd.DataFrame()
        for name, group in data.groupby('chr'):
            #print(name)
            gene_pairs = list(combinations(group.gene,2))
            dist_pairs = list(combinations(group['5p_gene'],2))
            gene_combination = pd.DataFrame(gene_pairs,columns=['gene1','gene2'])
            dist_combination = pd.DataFrame(dist_pairs,columns=['gene1_dist','gene2_dist'])

            df = pd.concat([gene_combination,dist_combination],axis=1)
            df['distance'] = abs(df.gene1_dist-df.gene2_dist)
            all_distance_combination_genes = all_distance_combination_genes.append(df[['gene1','gene2','distance']],ignore_index=True,sort=False)

        all_distance_combination_genes['distance'] = all_distance_combination_genes['distance'].astype('int64')
        gene_least_distant_genes_indices = all_distance_combination_genes.groupby('gene1')['distance'].idxmin()
        gene_least_distant_table = all_distance_combination_genes.loc[gene_least_distant_genes_indices]

        title_font = {'fontname':'sans-serif', 'size':'16', 'color':'black', 'weight':'normal','verticalalignment':'bottom'} # Bottom vertical align
        axis_font = {'fontname':'sans-serif', 'size':'14'}
        bar_plot, ax = plt.subplots(nrows=1,ncols=1)
        x_obj = ["≤100 kb","101-200 kb", "201-500 kb", ">500 kb"]
        y_obj = [sum(i < 100001 for i in gene_least_distant_table['distance']),
                sum(100000 < i < 200001 for i in gene_least_distant_table['distance']),
                sum(200000 < i < 500001 for i in gene_least_distant_table['distance']),
                sum(i > 500000 for i in gene_least_distant_table['distance'])]

        plt.bar(x_obj,y_obj, align = 'center')#BAR PLOT
        #plt.plot(x_obj,y_obj)#, align = 'center') #LINE PLOT
        plt.xlabel('Least distance between two genes from 5\' coordinate', **axis_font)
        plt.ylabel('Frequency of genes', **axis_font)
        bar_plot.savefig(output_folder_name+'/tf_analysis/gene_least_distant_barplot.png')
        plt.close(bar_plot)

        all_distance_combination_genes.to_csv(output_folder_name+'/tf_analysis/dist_all_induced_genes.distance', sep='\t',index=False)
        gene_least_distant_table.to_csv(output_folder_name+'/tf_analysis/least_distant_induced_genes.distance', sep='\t',index=False)
        # print('  Files generated: \n\t'+output_folder_name+'/tf_analysis/least_distant_induced_genes.distance\n\t'+output_folder_name+'/tf_analysis/dist_all_induced_genes.distance\n\t'+output_folder_name+'/tf_analysis/gene_least_distant_barplot.png\n')

class Randomization:
    def __init__(self):
        '''Read File'''

    def unique_genes_in_genome(self,coordinatefile):
        print("##### Genes coordinates from the GTF file...")
        self.coordinatefile = coordinatefile

        refgene_col = ["bin","name","chrom","strand","txStart","txEnd","cdsStart","cdsEnd","exonCount","exonStarts","exonEnds","score","name2","cdsStartStat","cdsEndStat","exonFrames"]
        refgene_data = pd.read_csv('input_data/'+self.coordinatefile, sep="\t", header=None,names=refgene_col,usecols=['chrom','strand','txStart','txStart','txEnd','name2'])#read hg38 refseq coordinate file without alt chr lines
        refgene_data = refgene_data[~refgene_data.chrom.str.contains('_') == True]#remove all alternative chromosomes from refGene
        refgene_data['txn_length'] = abs(refgene_data['txStart']-refgene_data['txEnd'])
        refgene_data_txn_length = refgene_data.groupby('name2', as_index=False).apply(pd.DataFrame.sort_values, 'txn_length', ascending = False).reset_index(drop=True)
        refgene_data['tss'] = refgene_data.groupby(['name2']).cumcount()+1
        refgene_data.loc[refgene_data['strand'] == '-',['txStart','txEnd']] = refgene_data.loc[refgene_data['strand'] == '-',['txEnd','txStart']].values

        genes_5p_3p_coord = refgene_data[['chrom','strand','txStart','txStart','txEnd','name2','tss']].copy()
        genes_5p_3p_coord1 = genes_5p_3p_coord.loc[genes_5p_3p_coord['tss'] == 1]#only tss1 is extracted

        genes_5p_coord = genes_5p_3p_coord1[['chrom','txStart','name2']].copy()
        genes_5p_coord.columns = ['chrom','txStart','txStart1','name2']

        sort_genes_5p_coord = genes_5p_coord.sort_values(['chrom','txStart'])#sort_genes_5p_coord.to_csv(cwd+'hg38.genes_5p_coord.bed',sep="\t",index=False,header=None)#print("  hg38 genes coordinate: 'hg38.genes_5p_coord.bed'\n")
        genes_in_genome = list(set(refgene_data.name2))#unique gene list from hg38 geome coordinate file
        #print("\t1) mm10_genes: unique genes list from hg38 RefSeq file\n\t2) sort_genes_5p_coord: tss1 (largest transcript of gene) coordinates of all genes\n")
        return genes_in_genome, sort_genes_5p_coord#two variables will be input for the function: genes_random_selection

    def random_selection_genes(self,rna_expression_file,genes_in_genome,genes_coordinates):
        print('##### Choosing random genes from unique gene list from the genome')
        self.rna_expression_file = rna_expression_file
        self.genes_in_genome = genes_in_genome
        self.genes_coordinates = genes_coordinates
        no_of_induced_genes = len(open('input_data/'+self.rna_expression_file).readlines())-1

        for i in range(1):
            selected_genes = random.sample(self.genes_in_genome,k=no_of_induced_genes)
            selected_genes_coord = self.genes_coordinates[self.genes_coordinates['name2'].isin(selected_genes)]#extraction of genes tss1 coordinates on the basis of selected gene list values
        data = selected_genes_coord.copy(); data.columns = ['chr','5p_gene','5p_gene1','gene']

        all_distance_combination_genes = pd.DataFrame()
        for name, group in data.groupby('chr'):
            gene_pairs = list(combinations(group.gene,2))
            dist_pairs = list(combinations(group['5p_gene'],2))
            gene_combination = pd.DataFrame(gene_pairs,columns=['gene1','gene2'])
            dist_combination = pd.DataFrame(dist_pairs,columns=['gene1_dist','gene2_dist'])
            df = pd.concat([gene_combination,dist_combination],axis=1)
            df['distance'] = abs(df.gene1_dist-df.gene2_dist)
            all_distance_combination_genes = all_distance_combination_genes.append(df[['gene1','gene2','distance']],ignore_index=True,sort=False)
        all_distance_combination_genes['distance'] = all_distance_combination_genes['distance'].astype('int64')
        gene_least_distant_genes_indices = all_distance_combination_genes.groupby('gene1')['distance'].idxmin()
        gene_least_distant_table = all_distance_combination_genes.loc[gene_least_distant_genes_indices]

        gene_least_distant_table.to_csv(output_folder_name+"/tf_analysis/randomly_selected_1_times."+str(no_of_induced_genes)+"_induced_genes.least_distance", sep='\t',index=False)
        # print("  Files generated: \n\t"+output_folder_name+"/tf_analysis/randomly_selected_1_times."+str(no_of_induced_genes)+"_induced_genes.least_distance\n")

        #PLOT OBSERVED VS RANDOMLY SELECTED GENES
        def obs_vs_rand(obs,rand):
            obs = pd.read_csv(obs,sep="\t",header=0,usecols=['distance']);obs['resource'] = 'Observed'
            rand = pd.read_csv(rand,sep="\t",header=0,usecols=['distance']);rand['resource'] = 'Random'
            data = pd.concat([obs[['distance','resource']],rand[['distance','resource']]],ignore_index=True)

            p = (ggplot(data=data)+
                    stat_density(mapping=aes(x='distance', y='..scaled..',colour='resource'),geom='line',position='identity')+
                    xlim(0,10000000)+scale_y_continuous()+
                    labs(y="Density", x="Least distance between TSS of two genes")+
                    ggtitle('pvalue = '+str(ks_2samp(obs['distance'],rand['distance'])[1]))+
                    theme_classic())
            ggsave(plot=p, filename=output_folder_name+'/tf_analysis/observed_vs_randomly_selected_genes_distribution.png', dpi=150)
            # print('\t'+output_folder_name+'/tf_analysis/observed_vs_randomly_selected_genes_distribution.png')

        obs_vs_rand(output_folder_name+'/tf_analysis/least_distant_induced_genes.distance',
                    output_folder_name+'/tf_analysis/randomly_selected_1_times.'+str(no_of_induced_genes)+'_induced_genes.least_distance')

class Tf_peak_tss_binding:
    def __init__(self):
        '''TF peak tss binding'''
    def tf_peak_tss_cluster(self,filename,tf,induced_gene_tss,filterscore):
        os.chdir(cwd)
        print("##### TF Peak Tss Binding...")
        self.filename = filename
        self.tf = tf
        self.induced_gene_tss = induced_gene_tss
        self.filterscore = filterscore

        #======Jing 2024-11-12========
        #solve TF peak input score problem
        #tf_peak_file = pd.read_csv('input_data/'+self.filename,sep="\t",header=0,usecols=[0,1,2,6],names=['chr','start','end','score'])        #macs2 $1,$2,$3,$5
        # tf_peak_file = pd.read_csv('input_data/'+self.filename,sep="\t",header=0,names=['chr','start','end','score'])        #macs2 $1,$2,$3,$5
        # tf_peak_file['score'] = tf_peak_file['score'].apply(pd.to_numeric, downcast='float', errors='coerce')
        # filter = tf_peak_file['score'] > self.filterscore
        # tf_peak_file = pd.read_csv('input_data/'+self.filename,sep="\t",header=0,names=['chr','start','end','score'])
        tf_peak_file = pd.read_csv('input_data/'+self.filename,sep="\t",header=0,names=['chr','start','end'])
        tf_peak_file['score'] = 1
        #convert score into float type
        # tf_peak_file['score'] = pd.to_numeric(tf_peak_file['score'], errors='coerce')
        # filterscore = float(self.filterscore)
        #filter
        filter = tf_peak_file['score'] > filterscore
        

        tf_peak_data = pd.DataFrame(tf_peak_file[filter], columns=['chr','start','end', 'score'])
        tf_peak_data['tf'] = self.tf
        tf_peak_data.loc[tf_peak_data['end'] < tf_peak_data['start'],['start','end']] = tf_peak_data.loc[tf_peak_data['end'] < tf_peak_data['start'],['end','start']].values
        gene_tss=pd.read_csv(output_folder_name+'/tf_analysis/induced_gene_tss',sep="\t",header=0,usecols = ['chrom','strand','txStart','txEnd','name2','tss'])
        #Add promoter upstream and downstream regions defined by users
        gene_tss['promoter_downstream'] = np.where(gene_tss['strand']== '-',gene_tss['txEnd']-promoter_downstream,gene_tss['txStart']+promoter_downstream)
        gene_tss['promoter_upstream'] = np.where(gene_tss['strand']== '-',gene_tss['txEnd']+promoter_upstream,gene_tss['txStart']-promoter_upstream)

        #swap column of negative strand
        strand = gene_tss['strand'] == '-'; gene_tss.loc[strand, ['promoter_downstream','promoter_upstream']] = gene_tss.loc[strand,['promoter_upstream','promoter_downstream']].values
        gene_promoter = gene_tss[['chrom','promoter_upstream','promoter_downstream','name2','tss']].copy()
        gene_promoter.columns = tf_peak_data.columns

        merge_tf_tss = pd.concat([tf_peak_data,gene_promoter])
        merge_tf_gene = merge_tf_tss.sort_values(['chr','start','end'])

        merge_tf_gene.loc[merge_tf_gene['end'] < merge_tf_gene['start'],['start','end']] = merge_tf_gene.loc[merge_tf_gene['end'] < merge_tf_gene['start'],['end','start']].values
        sort_merge_tf_gene = merge_tf_gene.sort_values(['chr','start'])

        sort_merge_tf_gene['start'] = sort_merge_tf_gene['start'].astype(int)
        sort_merge_tf_gene['end'] = sort_merge_tf_gene['end'].astype(int)

        sort_merge_tf_gene.to_csv(output_folder_name+'/tf_analysis/combined_gene_5p_tf_peak_calls_data.bed',sep="\t",index=False,header=None)

        pybedtools.bedtool.BedTool.cluster(pybedtools.bedtool.BedTool(output_folder_name+'/tf_analysis/combined_gene_5p_tf_peak_calls_data.bed'),d=0).saveas(output_folder_name+'/tf_analysis/gene_5p_tf_peak_calls_bed_tool_clust.bed')
        # print('  Gene Tss and TF Peak Cluster file is generated:\n\t'+output_folder_name+'/tf_analysis/gene_5p_tf_peak_calls_bed_tool_clust.bed')

        cluster_file_columns = ['chr', 'sTSS', 'eTSS', 'gene', 'tf', 'clust']
        cluster_file = pd.read_csv(output_folder_name+'/tf_analysis/gene_5p_tf_peak_calls_bed_tool_clust.bed', sep='\t', names=cluster_file_columns, dtype={'clust_no':int, 'tf_gene_tss_distance': int}, low_memory=False)

        tf_df = cluster_file[cluster_file['tf'] == self.tf].reset_index(drop=True)
        non_tf_df = cluster_file[cluster_file['tf'] != self.tf].reset_index(drop=True)

        temp_peak_table = pd.merge(non_tf_df, tf_df, on='clust')
        temp_peak_table.columns=['chr_gene','sTSS_gene','eTSS_gene','gene','tss_no','clust_no','chr_tf','sTSS_tf','eTSS_tf','score','tf']

        peak_table=temp_peak_table[['clust_no','gene','tss_no','score']].copy()
        peak_table['distance'] = np.minimum(abs(temp_peak_table['sTSS_gene']-temp_peak_table['eTSS_tf']),abs(temp_peak_table['sTSS_gene']-temp_peak_table['sTSS_tf']))
        peak_table.to_csv(output_folder_name+'/tf_analysis/gene_tss_peak.table', sep='\t',index=False)
        # print('  The gene tss peak table is generated: \n\t'+output_folder_name+'/tf_analysis/gene_tss_peak.table\n')
        

class Genespromotertfbindingcluster:
    def __init__(self):
        '''Gene Promoter TF Binding'''
    def tf_promoter_binding_cluster(self,induced_gene_cluster,gene_tss_peak_table):
        print("##### Genes promoter cluster...")
        #input: induced_genes_final_cluster.bed
        self.induced_gene_cluster = induced_gene_cluster
        #input: gene_tss_peak.table
        self.gene_tss_peak_table = gene_tss_peak_table

        #CREATE GENE PROMOTER CLUSTERS TABLE
        columns = ['chr', 'gene_5p', 'gene_5p2', 'gene_name', 'clust1', 'clust2']
        gene_cluster_data = pd.read_csv(output_folder_name+'/tf_analysis/'+self.induced_gene_cluster, sep="\t", names=columns)
        cluster_number = pd.Series(gene_cluster_data['clust2'].unique())
        gene_promoter_cluster_table_columns = ['cluster number', 'number of genes', 'number of promoters', 'flag','genes of cluster']
        gene_promoter_cluster_table = pd.DataFrame()

        for i in cluster_number:
            number_of_genes = []
            cluster_genes = []
            for j in range(0, len(gene_cluster_data)-1):
                if(gene_cluster_data.iloc[j, 5] == i):
                    k = j + 1
                    if(gene_cluster_data.iloc[k,5] == i):
                        dist = abs(gene_cluster_data.iloc[j, 1]-gene_cluster_data.iloc[k, 1])
                        number_of_genes.append(dist)
                    cluster_genes.append(gene_cluster_data.iloc[k-1, 3])
                    if(k == len(gene_cluster_data)-1):
                        cluster_genes.append(gene_cluster_data.iloc[k, 3])
            number_of_genes_per_cluster = len(number_of_genes)+1
            
            #count number_of_promoters (gene_5p distance>1000)
            
            number_of_promoters = sum(i > 1000 for i in number_of_genes)+1
            #number_of_promoters > 1, flag = 'U'
            flag = 'U' if number_of_promoters > 1 else 'F'
            gene_promoter_cluster_table=gene_promoter_cluster_table.append(pd.Series([i,number_of_genes_per_cluster, number_of_promoters, flag, ','.join(cluster_genes)], index=gene_promoter_cluster_table_columns), ignore_index=True)

        gene_promoter_cluster_table[['cluster number','number of genes','number of promoters']] = gene_promoter_cluster_table[['cluster number','number of genes','number of promoters']].astype(int)
        gene_promoter_cluster_table[['flag','genes of cluster']] = gene_promoter_cluster_table[['flag','genes of cluster']].astype(str)
        gene_promoter_cluster_table = gene_promoter_cluster_table.reindex(columns=gene_promoter_cluster_table_columns)

        #gene_promoter_cluster_table.to_csv(output_folder_name+'/tf_analysis/gene_promoter_cluster.table', sep='\t',index=False)
        #print('  The Gene Promoter Cluster Table is generated:\n\t'+output_folder_name+'/tf_analysis/gene_promoter_cluster.table\n')

        gene_peak_data = pd.read_csv(output_folder_name+'/tf_analysis/'+self.gene_tss_peak_table,sep="\t",header=0)
        tf_bound_genes = gene_peak_data['gene'].unique()

        gene_promoter_cluster_tf_binding_table_columns = ['cluster number','number of genes per cluster','number of promoters per cluster','flag','genes in cluster','number of promoter binding TFs in cluster']
        gene_promoter_cluster_tf_binding_table = pd.DataFrame(columns=gene_promoter_cluster_tf_binding_table_columns)


        # number of promoter binding TFs in cluster = count
        for index, gene_promoter_cluster_data_row in gene_promoter_cluster_table.iterrows():
            count = 0

            genes_in_cluster = [x.strip() for x in gene_promoter_cluster_data_row[4].split(',')]

            for element in genes_in_cluster:
                if element in tf_bound_genes:
                    count += 1
            
            #modified: consider distance>1000 between two genes (solve divergent issue)
            # for element in genes_in_cluster:
            #   if count == 0 and element in tf_bound_genes:
            #       count += 1
            #   elif count > 0 and element in tf_bound_genes and gene_peak_data[gene_peak_data['gene'] == element]['distance'].min() > 1000:
            #       count += 1
            #add information, count as 6th column: number of promoter binding TFs in cluster
            gene_promoter_cluster_tf_binding_table = gene_promoter_cluster_tf_binding_table.append(pd.Series([gene_promoter_cluster_data_row[0],gene_promoter_cluster_data_row[1],gene_promoter_cluster_data_row[2],gene_promoter_cluster_data_row[3],gene_promoter_cluster_data_row[4],count],index=gene_promoter_cluster_tf_binding_table_columns), ignore_index = True)

        
        gene_promoter_cluster_tf_binding_table['number of promoter binding TFs in cluster'] = np.where(gene_promoter_cluster_tf_binding_table['number of promoters per cluster']>gene_promoter_cluster_tf_binding_table['number of promoter binding TFs in cluster'],gene_promoter_cluster_tf_binding_table['number of promoter binding TFs in cluster'],gene_promoter_cluster_tf_binding_table['number of promoters per cluster'])
        gene_promoter_cluster_tf_binding_table.to_csv(output_folder_name+'/tf_analysis/gene_promoter_cluster_tf_binding.table', sep='\t',index=False)
        gene_promoter_cluster_tf_binding_frequency_table = gene_promoter_cluster_tf_binding_table[gene_promoter_cluster_tf_binding_table.apply(lambda x: x['flag'] == 'U', axis=1)].groupby(['number of promoters per cluster', 'number of promoter binding TFs in cluster']).size().reset_index(name="Frequency")

        bubble_chart_frequency_table = gene_promoter_cluster_tf_binding_frequency_table.groupby(['number of promoters per cluster', 'number of promoter binding TFs in cluster'])['Frequency'].sum().reset_index()
        promoter_per_cluster = bubble_chart_frequency_table.groupby(['number of promoters per cluster'])['Frequency'].sum().reset_index()
        bubble_chart_frequency_table = bubble_chart_frequency_table.merge(promoter_per_cluster, left_on='number of promoters per cluster', right_on='number of promoters per cluster', how='left')
        bubble_chart_frequency_table['Frequency Percentage'] = 100*(bubble_chart_frequency_table['Frequency_x']/bubble_chart_frequency_table['Frequency_y'])
        bubble_chart_frequency_table.columns = ['number of promoters per cluster','number of promoter binding TFs in cluster','frequency of number of promoter binding TFs in cluster','frequency of number of promoters per cluster','Frequency Percentage']
        bubble_chart_frequency_table.to_csv(output_folder_name+'/tf_analysis/bubble_chart_gene_promoter_cluster_tf_binding_frequency.table',sep="\t",index=None)

        # print("  The files generated:\n\t"+output_folder_name+"/tf_analysis/gene_promoter_cluster_tf_binding.table\n\t"+output_folder_name+"/tf_analysis/bubble_chart_gene_promoter_cluster_tf_binding_frequency.table\n")



#=====2025-02-10======

class EnhancedTfBindingAnalysis:
    def __init__(self, refseq_gtf, peak_bed, promoter_upstream, promoter_downstream, output_folder):
        
        self.refseq_gtf = refseq_gtf
        self.peak_bed = peak_bed
        self.promoter_upstream = promoter_upstream
        self.promoter_downstream = promoter_downstream
        self.output_folder = output_folder

    def process_all_tss_promoters(self):
        
        # read RefSeq annotation file
        refgene_col = ["bin", "name", "chrom", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd",
                       "exonCount", "exonStarts", "exonEnds", "score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames"]
        refgene_data = pd.read_csv('input_data/'+self.refseq_gtf, sep="\t", header=None, names=refgene_col)

        # Filter out invalid chromosomes.
        refgene_data = refgene_data[~refgene_data.chrom.str.contains('_')]

        # Generate all TSS promoter regions
        refgene_data['tss'] = refgene_data.groupby('name2').cumcount() + 1
        refgene_data['promoter_start'] = np.where(
            refgene_data['strand'] == '+',
            refgene_data['txStart'] - self.promoter_upstream,
            refgene_data['txEnd'] - self.promoter_downstream
        )
        refgene_data['promoter_end'] = np.where(
            refgene_data['strand'] == '+',
            refgene_data['txStart'] + self.promoter_downstream,
            refgene_data['txEnd'] + self.promoter_upstream
        )

        # filter invalid coordinates
        refgene_data = refgene_data[refgene_data['promoter_start'] >= 0]
        
        #keep inudced_gene all tss
        gene_tss=pd.read_csv(output_folder_name+'/tf_analysis/induced_gene_tss',sep="\t",header=0,usecols = ['chrom','strand','txStart','txEnd','name2','tss'])
        induced_tss_promoter = refgene_data[refgene_data['name2'].isin(gene_tss['name2'])]
        
        # save TSS promoter
        output_path = f"{self.output_folder}/tf_analysis/induced_alt_tss_promoters.bed"
        induced_tss_promoter[['chrom', 'promoter_start', 'promoter_end', 'name2']].to_csv(
            output_path, sep="\t", header=False, index=False
        )
        return output_path

    def enhanced_tf_binding_analysis(self, gene_promoter_cluster_tf_binding):
        
        all_tss_promoter_bed = self.process_all_tss_promoters()

        # BedTool overlap 
        tf_peaks = pybedtools.BedTool(os.path.join('input_data', self.peak_bed))
        all_promoters = pybedtools.BedTool(all_tss_promoter_bed)

        # overlap tf_peaks and promoters
        overlaps = tf_peaks.intersect(all_promoters, wa=True, wb=True)
        
        #ouput TF_binding_promoters
        output_file_path = os.path.join(output_folder_name, 'tf_analysis/TF_binding_promoters.txt')
        with open(output_file_path, 'w') as f:
            f.write(str(overlaps))


#=====Status of Python Libraries======
InstallPythonLibraries()



#======TAKING INPUT FROM USER======
print('\n\nUSER INPUTS:')
refseq_gtf = input('\n   Enter the name of genome annotation RefSeq GTF file: ')#'hg19.refGene.txt'
print('\n   Definition of the promoter region:')
promoter_upstream = int(input('   \tUpstream from TSS (in basepairs): '))#1000
promoter_downstream = int(input('   \tDownstream from TSS (in basepairs): '))#1000
rna_expression_file = input('\n   Enter the name of Expression file: ')#inducedgenes.tsv
# induced_gene_log2FC = int(input('\n   Induced gene log2FC minimum cutoff: '))#2
induced_gene_log2FC = 0
tf_chiqseq_bedfile = input('\n   Enter the name of TF chipseq bed file: ')#peaks.bed
tf_name = input('\n   Enter the name of TF: ')#'TF'
# fc = int(input('\n   MACS score minimum cutoff to filter chipseq peaks: '))#100
fc = 0
clustering_method = int(input('\n   Clustering methods:\n\t "1"   for clustering by Distance\n\t "2"   for clustering by TAD\n\n\t Please enter the option:   '))

if clustering_method == 1:
    induced_genes_cluster_distance = int(input('\n\t Induced gene clustering distance maximum cutoff (in basepairs): '))
    output_folder_name = tf_name+'_'+str(input('\n   Enter the induced condition: '))+'_'+os.path.splitext(rna_expression_file)[0]+'_Promoter_'+str(promoter_upstream)+'_'+str(promoter_downstream)+'_dist'+str(induced_genes_cluster_distance)

else:
    tad_file = input('\n\t Enter TAD file name: ')#'tad.bed'
    output_folder_name = tf_name+'_'+str(input('\n   Enter the induced condition: '))+'_'+os.path.splitext(rna_expression_file)[0]+'_Promoter_'+str(promoter_upstream)+'_'+str(promoter_downstream)+'_TAD'

print('\n')



#print(output_folder_name)
if not os.path.exists(output_folder_name):#'./output/tf_analysis'):#generate output folders: exp_bedgraph and tf_analysis
    for sub_folders in [output_folder_name+'/exp_bedgraph',output_folder_name+'/tf_analysis']:
        os.makedirs(sub_folders)

expression_bedgraph=Expression_bedgraph()
expression_bedgraph.induced_genes(rna_expression_file,refseq_gtf)

extract_induced_genes=Extract_induced_genes_coordinates()
induced_genes_tss = extract_induced_genes.induced_genes(rna_expression_file,refseq_gtf)

induced_genes_clust = Inducedgenesclustering('induced_gene_tss')
if (clustering_method == 1):
    induced_genes_cluster=induced_genes_clust.induced_genes_clusteringByDistance()
else:
    induced_genes_cluster=induced_genes_clust.induced_genes_clusteringByTAD(tad_file)

# induced_genes_cluster_barplot=Genesclusterbarplot()
# induced_genes_cluster_barplot.genes_clust_genes_frequency('induced_genes_final_cluster.bed')

# least_distant_induced_genes= Genesleastdistant()
# least_distant_induced_genes.genes_least_distnce('induced_genes_5p_coordinates.bed')

# random_selection_least_distant_genes = Randomization()
# genes_in_genome, refgene_data = random_selection_least_distant_genes.unique_genes_in_genome(refseq_gtf)
# random_selection_least_distant_genes.random_selection_genes(rna_expression_file,genes_in_genome, refgene_data)

tf_peak_tss_binding = Tf_peak_tss_binding()
tf_peak_tss_binding.tf_peak_tss_cluster(tf_chiqseq_bedfile,tf_name,'induced_genes_tss',filterscore=fc)

tss_tf_clustering = Genespromotertfbindingcluster()
tss_tf_cluster = tss_tf_clustering.tf_promoter_binding_cluster('induced_genes_final_cluster.bed','gene_tss_peak.table')

# freq_bubble_plot = TfBindingGenePromoterFrequencyBubblePlot()
# freq_bubble_plot.bubble_plot('bubble_chart_gene_promoter_cluster_tf_binding_frequency.table')




#============2024-07-31==========
#improve ouptut
#Generate gene_promoter_cluster_tf_binding_Epromoter
# read gene_promoter_cluster_tf_binding gene_tss_peak
gene_promoter_cluster_tf_binding = pd.read_csv(output_folder_name+'/tf_analysis/'+'gene_promoter_cluster_tf_binding.table', sep='\t')
gene_tss_peak = pd.read_csv(output_folder_name+'/tf_analysis/'+'gene_tss_peak.table', sep='\t')

# get gene name from gene_tss_peak
tss_genes = gene_tss_peak.iloc[:, 1].unique()

# add Epromoter_gene
def get_epromoter_genes(row, tss_genes):
    genes_in_cluster = [gene.strip() for gene in row['genes in cluster'].split(',')]
    epromoter_genes = [gene for gene in genes_in_cluster if gene in tss_genes]
    return ','.join(epromoter_genes) if epromoter_genes else 'None'

gene_promoter_cluster_tf_binding['Epromoter_gene'] = gene_promoter_cluster_tf_binding.apply(
    get_epromoter_genes, axis=1, tss_genes=tss_genes
)


#====add gene TSS coordinate======
induced_genes_5p = pd.read_csv(output_folder_name+'/tf_analysis/'+'induced_genes_5p_coordinates.bed', sep='\t', names=['chr', 'start', 'end', 'gene'])

# initiate chr,TSS
gene_promoter_cluster_tf_binding["chr"] = "NA"
gene_promoter_cluster_tf_binding["TSS"] = "NA"

# update gene_promoter_cluster_tf_binding DataFrame
for index, row in induced_genes_5p.iterrows():
    gene = row["gene"]
    chr_val = row["chr"]
    tss_val = row["start"]
    
    for i, epromoter_genes in gene_promoter_cluster_tf_binding["Epromoter_gene"].iteritems():
        if epromoter_genes != "None" and gene in epromoter_genes.split(","):
            if gene_promoter_cluster_tf_binding.at[i, "chr"] == "NA":
                gene_promoter_cluster_tf_binding.at[i, "chr"] = chr_val
                gene_promoter_cluster_tf_binding.at[i, "TSS"] = str(tss_val)
            else:
                gene_promoter_cluster_tf_binding.at[i, "chr"] += f",{chr_val}"
                gene_promoter_cluster_tf_binding.at[i, "TSS"] += f",{tss_val}"




#========consider divergent promoter or close genes(two genes<2kb)=========
# for each row, update number of promoter binding TFs in cluster
for index, row in gene_promoter_cluster_tf_binding.iterrows():
    if row['number of promoter binding TFs in cluster'] == 2 and row['TSS'] != 'NA':
        tss_values = [int(x) for x in row['TSS'].split(',') if x.strip().isdigit()]
        if len(tss_values) >= 2:  # 确保至少有两个TSS值
            dist = abs(tss_values[0] - tss_values[1])
            if dist < 2000:
                gene_promoter_cluster_tf_binding.at[index, 'number of promoter binding TFs in cluster'] = 1



#====2025-02-10======
# initiate
enhanced_analysis = EnhancedTfBindingAnalysis(
    refseq_gtf,  
    tf_chiqseq_bedfile,              
    promoter_upstream,                      
    promoter_downstream,                    
    output_folder_name             
)

# get TF_binding_promoters.txt
enhanced_analysis.enhanced_tf_binding_analysis(gene_promoter_cluster_tf_binding)


#====add gene promoter coordinate======
# TF_binding_promoters.txt
tf_binding_promoters_file = output_folder_name +'/tf_analysis'+ '/TF_binding_promoters.txt'

if os.path.exists(tf_binding_promoters_file) and os.path.getsize(tf_binding_promoters_file) > 0:
    tf_promoters = pd.read_csv(tf_binding_promoters_file, sep='\t', header=None)
    tf_promoters.columns = ['chr_peak', 'peak_start', 'peak_end', 'chr_promoter', 'promoter_start', 'promoter_end', 'gene']

    # create gene to promoter mapping
    gene_to_promoter = dict(zip(tf_promoters['gene'], tf_promoters[['chr_promoter', 'promoter_start', 'promoter_end']].apply(tuple, axis=1)))

    #remove 'chr', 'TSS'
    gene_promoter_cluster_tf_binding = gene_promoter_cluster_tf_binding.drop(columns=['chr', 'TSS'])

    # initiate chr,promoter coordinates
    gene_promoter_cluster_tf_binding["chr_promoter"] = "None"
    gene_promoter_cluster_tf_binding["promoter_start"] = "None"
    gene_promoter_cluster_tf_binding["promoter_end"] = "None"

    # update gene_promoter_cluster_tf_binding DataFrame
    for index, row in gene_promoter_cluster_tf_binding.iterrows():
        genes = row['Epromoter_gene'].split(',')
        chr_promoter_vals = []
        promoter_start_vals = []
        promoter_end_vals = []
        
        for gene in genes:
            gene = gene.strip()
            if gene in gene_to_promoter:
                chr_promoter, start, end = gene_to_promoter[gene]
                chr_promoter_vals.append(chr_promoter)
                promoter_start_vals.append(str(start))
                promoter_end_vals.append(str(end))
        
        # rm dup
        unique_chr = list(set(chr_promoter_vals))
        unique_start = list(set(promoter_start_vals))
        unique_end = list(set(promoter_end_vals))
        
        # update column
        gene_promoter_cluster_tf_binding.at[index, 'chr_promoter'] = ','.join(unique_chr)
        gene_promoter_cluster_tf_binding.at[index, 'promoter_start'] = ','.join(unique_start)
        gene_promoter_cluster_tf_binding.at[index, 'promoter_end'] = ','.join(unique_end)

else:
    gene_promoter_cluster_tf_binding = gene_promoter_cluster_tf_binding.drop(columns=['chr', 'TSS'])




#only keep flag=U for output (clusters)
gene_promoter_cluster_tf_binding = gene_promoter_cluster_tf_binding[(gene_promoter_cluster_tf_binding['flag'] == 'U')]

#set new cluster number
gene_promoter_cluster_tf_binding["cluster number"] = range(1, len(gene_promoter_cluster_tf_binding) + 1)
#colnames Epromoter_gene change to TF_binding_gene
gene_promoter_cluster_tf_binding = gene_promoter_cluster_tf_binding.rename(columns={'Epromoter_gene': 'TF_binding_gene'})

# flag=U, number of promoter binding TFs in cluster == 1
epromoter_cluster_res = gene_promoter_cluster_tf_binding[(gene_promoter_cluster_tf_binding['flag'] == 'U') & 
                                       			(gene_promoter_cluster_tf_binding['number of promoter binding TFs in cluster'] == 1)]

#remove flag column
gene_promoter_cluster_tf_binding = gene_promoter_cluster_tf_binding.drop('flag', axis=1)
epromoter_cluster_res = epromoter_cluster_res.drop('flag', axis=1)


# # save into new file
# gene_promoter_cluster_tf_binding.to_csv(output_folder_name+'/gene_promoter_cluster_tf_binding_Epromoter.table', sep='\t', index=False)
# epromoter_cluster_res.to_csv(output_folder_name+'/Epromoter_cluster', sep='\t', index=False)

# save into new file
gene_promoter_cluster_tf_binding.to_csv(
    os.path.join(output_folder_name, 'gene_promoter_cluster_tf_binding.txt'), 
    sep='\t', 
    index=False
)

epromoter_cluster_res.to_csv(
    os.path.join(output_folder_name, 'Epromoter_cluster.txt'), 
    sep='\t', 
    index=False
)

print("gene_promoter_cluster_tf_binding.txt has been generated.")
print("Epromoter_cluster.txt has been generated.")




#======alternative promoter======
# TF_binding_promoters.txt
tf_binding_promoters_file = output_folder_name +'/tf_analysis'+ '/TF_binding_promoters.txt'

if os.path.exists(tf_binding_promoters_file) and os.path.getsize(tf_binding_promoters_file) > 0:

    tf_promoters = pd.read_csv(tf_binding_promoters_file, sep='\t', header=None)
    tf_promoters.columns = ['chr_peak', 'peak_start', 'peak_end', 'chr_promoter', 'promoter_start', 'promoter_end', 'gene']

    # Ensure that the genes in TF_binding_promoters.txt are included in the TF_binding_gene column of the gene_promoter_cluster_tf_binding dataset.
    final_output = gene_promoter_cluster_tf_binding.copy()
    final_output['TF_binding_gene'] = final_output['genes in cluster'].apply(
        lambda x: ','.join([g for g in x.split(',') if g in set(tf_promoters['gene'])]) or 'None'
    )
    #When the number of promoter-binding TFs in the cluster is 0, recalculate the number of genes in TF_binding_gene.
    final_output.loc[final_output['number of promoter binding TFs in cluster'] == 0, 'number of promoter binding TFs in cluster'] = \
    final_output['TF_binding_gene'].apply(lambda x: len(x.split(',')) if x != 'None' else 0)

    # save res
    final_output.to_csv(
        f"{output_folder_name}/gene_promoter_cluster_tf_binding.txt",
        sep="\t",
        index=False
    )

    # Epromoter
    epromoter_res = final_output[final_output['number of promoter binding TFs in cluster'] == 1]
    epromoter_res.to_csv(
        f"{output_folder_name}/Epromoter_cluster.txt",
        sep="\t",
        index=False
    )



#======2025-02-10========
# Calculate Statistical Information
# Obtain statistical information from the gene_promoter_cluster_tf_binding data
number_of_induced_genes = len(induced_genes_tss)  
number_of_clusters = len(gene_promoter_cluster_tf_binding)  
number_of_epromoter_clusters = len(epromoter_cluster_res)  

# Write the statistical information to the file Cluster_summary_stats.txt
summary_stats_file = os.path.join(output_folder_name, 'Cluster_summary_stats.txt')
with open(summary_stats_file, 'w') as f:
    f.write(f"Number_of_induced_genes: {number_of_induced_genes}\n")
    f.write(f"Number_of_clusters: {number_of_clusters}\n")
    f.write(f"Number_of_Epromoter_clusters: {number_of_epromoter_clusters}\n")

print(f"Cluster_summary_stats.txt has been generated.")



#====2024-11-12======
#tf_peaks=Tfpeakdistribution()
#tf_peaks.tf_peak_distribution(rna_expression_file,tf_chiqseq_bedfile,'mm9.refGene.all.tss.txt')
# list = ["combined_gene_5p_tf_peak_calls_data.bed",'dist_all_induced_genes.distance','gene_5p_tf_peak_calls_bed_tool_clust.bed',
#         'gene_tss_peak.table','genome.all_distance_combination_genes.distance','genome.least_distant_genes.distance',
#         'induced_genes_5p_3p_coordinates.bed','induced_genes_5p_coordinates.bed','induced_genes_bed_tool_clust.bed',
#         'induced_genes_final_cluster.bed','induced_gene_tss','mm9.refGene.all.tss.txt']
#print('##### Temporary files are deleted...')
#[os.remove(os.path.join(output_folder_name+'/tf_analysis/',f)) for f in list]

try:
   # delete tf_analysis
   for f in os.listdir(os.path.join(output_folder_name, 'tf_analysis')):
       os.remove(os.path.join(output_folder_name, 'tf_analysis', f))
   os.rmdir(os.path.join(output_folder_name, 'tf_analysis'))
   
   # delete exp_bedgraph
   for f in os.listdir(os.path.join(output_folder_name, 'exp_bedgraph')):
       os.remove(os.path.join(output_folder_name, 'exp_bedgraph', f))
   os.rmdir(os.path.join(output_folder_name, 'exp_bedgraph'))

except FileNotFoundError:
   print("Folder does not exist")

print("The job is finished by: ",datetime.datetime.now())
print("Total time cost: ",datetime.datetime.now()-start_time,"\n")



