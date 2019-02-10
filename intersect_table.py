##using pandas to read the intersected repeats file and analyze the Info column

import pandas as pd
import csv
from glob import iglob, glob


def columnsplit(table, column, character_to_split):
    isplit = table[column].str.split(character_to_split, expand=True)
    return isplit

def input_data(filepath):

    intersect = pd.read_csv(filepath, sep='\t', header=None)
    intersect.rename(columns={8: 'Info', 0: 'Name', 3: 'Repeat start', 4: 'Repeat stop', 11: 'Feature'}, inplace=True)

    ###splitting Info column
    splitcol = ['repeat', 'positions', 'us', 'unit size', 'rn', 'repeat number', 'p', 'perfection', 'u', 'unit']
    infosplit = columnsplit(intersect, 'Info', ' ')
    infosplit.rename(columns=dict(zip(infosplit.columns[:], splitcol)), inplace=True)

    ###splitting positions column to get start and stop positions of genes
    position_split = columnsplit(infosplit, 'positions', '-')
    position_split.rename(columns={0: 'Gene start', 1: 'Gene stop'}, inplace=True)

    infosplit = infosplit.join(position_split)
    infosplit = infosplit.drop('positions', axis=1)
    intersect = intersect.join(infosplit)

    id_split = columnsplit(intersect, 17, ';')
    id_split = columnsplit(id_split, 0, '=')
    id_split = id_split.drop(0, axis=1)
    id_split = id_split.rename(columns={1: 'ID'})
    intersect = intersect.join(id_split)

    col_to_drop = ['Info', 'us', 'rn', 'p', 'u']
    intersect = intersect.drop(col_to_drop, axis=1)


    intersect_cds = intersect[intersect['Feature'] == 'CDS']
    intersect_cds = intersect_cds.sort_values(['Name', 'Gene start'], axis=0)
    return intersect_cds


def main():
    intersect_dir = '/Users/tanayajadhav/drexel_internship/intersects'
    for intersect_path in iglob(intersect_dir + '/*_repeats'):
        strain_name = intersect_path.split('/')[-1].split('.')[0]
        intersect_cds = input_data(intersect_path)

        rows = []
        for i in range(len(intersect_cds)-1):
            chrom_name = intersect_cds.iloc[i]['Name']
            prev_cds_stop = intersect_cds.iloc[i]['Gene stop']
            prev_cds_stop = int(prev_cds_stop)-1
            next_cds_start = intersect_cds.iloc[i+1]['Gene start']
            cds1 = intersect_cds.iloc[i]['ID']
            cds2 = intersect_cds.iloc[i+1]['ID']
            parent_ids = cds1 + '_____' + cds2
            strand_cds1 = intersect_cds.iloc[i][15]
            strand_cds2 = intersect_cds.iloc[i+1][15]
            strands = strand_cds1 + '_' + strand_cds2
            col_values = pd.Series({'Chr':chrom_name, 'Start':prev_cds_stop, 'Stop':next_cds_start,
                                    'Flanking CDS IDs':parent_ids, 'Score':'.', 'Strands':strands})
            rows.append(col_values)

        new_file = pd.concat(rows, axis=1).T
        new_file = new_file.reindex(columns=['Chr', 'Start', 'Stop', 'Flanking CDS IDs', 'Score', 'Strands'])

        inverse_file_path = '/Users/tanayajadhav/drexel_internship/inverse_bed_cds/' + strain_name + '.inverse.bed'
        new_file.to_csv(inverse_file_path, sep='\t', index=False)


if __name__ == '__main__':
    main()