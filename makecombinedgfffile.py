import pandas as pd
from glob import iglob


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


def add_clustername(df, dictionary, id_column, columnname):
    df[columnname] = ""
    indexlist = list(df.index)
    for i in indexlist:
        cell_val = df.loc[i, id_column]
        if cell_val in dictionary:
            df.loc[i, columnname] = dictionary[cell_val]
        else:
            df.loc[i, columnname] = '.'
    return df


def main():
    gpap_path = '/Users/tanayajadhav/drexel_internship/gene_presence_absence_paralogs_merged.csv'
    gpap = pd.read_csv(gpap_path, sep=',', header=0, index_col='Gene')
    # print(gpap.shape)
    # print(gpap.index[27])
    gpap.index = gpap.index.str.replace('\t', '__')
    # print(gpap.index[26])
    # print(gpap.index[27])
    gpap = gpap.fillna('.')
    # print(gpap.shape)
    # return

    strain_names = list(gpap.columns[10:].values)
    genedict = make_genedict(gpap, strain_names)

    #read in gff file with annotations for all strains
    gff_path = '/Users/tanayajadhav/drexel_internship/from_rachel/fixed_gffs/combined.gff'
    gff_data = pd.read_csv(gff_path, sep='\t', header=None, comment='#')

    #rename columns and filter based on type so only CDS intervals are present
    colnames = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase',
                    'attributes']
    gff_data.rename(columns=dict(zip(gff_data.columns[:], colnames)), inplace=True)
    gff_data = gff_data[gff_data['type'] == 'CDS']

    # add column with gene_id by parsing attributes
    attributes_split = columnsplit(gff_data, 'attributes', ';')
    attributes_split = attributes_split.drop(attributes_split.columns[1:], axis=1).rename(columns={0:'Gene'})
    gene_split = columnsplit(attributes_split, 'Gene', '=')
    gene_split = gene_split.drop(gene_split.columns[0], axis=1).rename(columns={1: 'gene_id'})

    gff_data = gff_data.join(gene_split)


    #add a column with clusternames that correspond to the gene id
    gff_data = add_clustername(gff_data, genedict, 'gene_id', 'cluster')
    filtered_gff = gff_data[gff_data.cluster != '.']

    filtered_gff_file = '/Users/tanayajadhav/drexel_internship/combined_gff_withclusters.gff'
    filtered_gff.to_csv(filtered_gff_file, sep='\t', header=0, index=0)




if __name__ == '__main__':
    main()
