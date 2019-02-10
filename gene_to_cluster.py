# Adds a column to intersect files that has the name of the gene cluster that the gene belongs to

import pandas as pd
from glob import iglob
import os


def make_genedict(gpap, strain_names):
    genedict = {}
    for strain_name in strain_names:
        gene_column = gpap[strain_name]
        for clustername, genes in gene_column.iteritems():
            split_genes = genes.split('\t')
            if len(split_genes) >= 1:
                for gene in split_genes:
                   genedict[gene] = clustername
    return genedict


def columnsplit(table, column, character_to_split):
    isplit = table[column].str.split(character_to_split, expand=True)
    return isplit


def in_data(filepath):
    data = pd.read_csv(filepath, sep='\t', header=None, comment='#')
    colnames = ['seqid1', 'source1', 'type', 'start1', 'end1', 'score1', 'strand1', 'phase1',
                'attributes1', 'seqid2', 'source2', 'feature', 'start2', 'end2', 'score2', 'strand2',
                'phase2', 'attributes2', 'misc']
    data.rename(columns=dict(zip(data.columns[:], colnames)), inplace=True)
    attributes2_split = columnsplit(data, 'attributes2', ';')
    attributes2_split = attributes2_split.drop(attributes2_split.columns[1:], axis=1)
    attributes2_split = attributes2_split.rename(columns={0: 'Gene'})
    gene_split = columnsplit(attributes2_split, 'Gene', '=')
    gene_split = gene_split.drop(gene_split.columns[0], axis=1).rename(columns={1: 'gene_id'})
    # print(gene_split)
    data = data.join(gene_split)
    return data


def intercds_in_data(intercds_filepath):
    intercds_data = pd.read_csv(intercds_filepath, sep='\t', header=None)
    colnames = ['seqid1', 'source1', 'type', 'start1', 'end1', 'score1', 'strand1', 'phase1',
                'attributes1', 'seqid2', 'start2', 'end2', 'flanking_genes', 'score2', 'strand2',
                'intersecting_bases']
    intercds_data.rename(columns=dict(zip(intercds_data.columns[:], colnames)), inplace=True)
    flanking_split = columnsplit(intercds_data, 'flanking_genes', '_____')
    ID1_split = columnsplit(flanking_split, 0, '=')
    ID1_split = ID1_split.drop(ID1_split.columns[0], axis=1).rename(columns={1: 'ID1'})
    ID2_split = columnsplit(flanking_split, 1, '=')
    ID2_split = ID2_split.drop(ID2_split.columns[0], axis=1).rename(columns={1: 'ID2'})
    intercds_data = intercds_data.join(ID1_split).join(ID2_split)
    return intercds_data


def add_column(file, dictionary):
    filename = file.split('/')[-1]
    gff_data = in_data(file)
    gff_data = add_clustername(gff_data, dictionary, 'gene_id', 'roary_genecluster')
    output_file_dir = '/Users/tanayajadhav/drexel_internship/clusterfiles/'
    output_filepath = output_file_dir + filename + '.with_clusters'
    # print(output_filepath)
    gff_data.to_csv(output_filepath, sep=',', header=0)


def add_clustername(df, dictionary, ID_column, columnname):
    df[columnname] = ""
    indexlist = list(df.index)
    for i in indexlist:
        cell_val = df.loc[i, ID_column]
        if cell_val in dictionary:
            df.loc[i, columnname] = dictionary[cell_val]
    return df


def main():
    file_path = '/Users/tanayajadhav/drexel_internship/gene_presence_absence_paralogs_merged.csv'
    gpap = pd.read_csv(file_path, sep=',', header=0, index_col='Gene')
    gpap.index = gpap.index.str.replace('\t', '__')
    gpap = gpap.fillna('.')
    strain_names = list(gpap.columns[10:].values)
    genedict = make_genedict(gpap, strain_names)
    print(genedict)
    return
    intersect_gff_dir = '/Users/tanayajadhav/drexel_internship/intersects'
    for file in iglob(intersect_gff_dir + '/*.intersected_repeats'):
        add_column(file, genedict)

    intercds_dir = '/Users/tanayajadhav/drexel_internship/intercds_intersects'
    for file in iglob(intercds_dir + '/*repeats'):
        if os.stat(file).st_size == 0:
            pass
            # print(file + ' is empty')
        else:
            intercds_data = intercds_in_data(file)
            intercds_data = add_clustername(intercds_data, genedict, 'ID1', 'roary_genecluster1')
            intercds_data = add_clustername(intercds_data, genedict, 'ID2', 'roary_genecluster2')
            output_dir = '/Users/tanayajadhav/drexel_internship/clusterfiles'
            filename = file.split('/')[-1]
            outfile = output_dir + '/' + filename + '.with_clusters'
            intercds_data.to_csv(outfile, sep='\t', header=0)
        break


if __name__ == '__main__':
    main()
