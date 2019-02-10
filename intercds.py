###Create bed files for intercds regions using original annotated gffs from prokka
###Intersect these (using bedtools intersect) with the phobos output files to find repeats that occur in these regions

import pandas as pd
import csv
from glob import iglob, glob
from subprocess import call

def columnsplit(table, column, character_to_split):
    isplit = table[column].str.split(character_to_split, expand=True)
    return isplit

def in_data(filepath):
    data = pd.read_csv(filepath, sep='\t', header=None, comment='#')
    colnames = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase',
                'attributes']
    data.rename(columns=dict(zip(data.columns[:], colnames)), inplace=True)
    attributes_split = columnsplit(data, 'attributes', ';')
    attributes_split = attributes_split.drop(attributes_split.columns[1:], axis=1)
    attributes_split = attributes_split.rename(columns={0:'Gene'})
    data = data.join(attributes_split)
    data_cds = data[data['type'] == 'CDS']
    data_cds = data_cds.sort_values(['seqid', 'start'], axis=0)
    # print(data_cds)
    return data_cds

def intersect_phobos_intercds(phobos_gff_dir, intercds_bed_dir, intersect_dir):
    for phobos_gff in iglob(phobos_gff_dir + '/*.repeats.gff'):
        strain_name = phobos_gff.split('/')[-1].replace('.repeats.gff', '')
        intercds_bed = intercds_bed_dir + '/' + strain_name + '_intercds.bed'
        cmd = ['bedtools', 'intersect',
               '-wao',
               '-a', phobos_gff,
               '-b', intercds_bed]
        intersect_path = intersect_dir + '/' + strain_name + '.intercds_repeats'
        with open(intersect_path, 'w') as f:
            call(cmd, stdout=f)

def main():
    gff_dir = '/Users/tanayajadhav/drexel_internship/from_rachel/fixed_gffs'
    for gff_file in iglob(gff_dir + '/*.gff'):
        strain_name = gff_file.split('/')[-1].split('.')[0]
        # gff_file = '/Users/tanayajadhav/drexel_internship/from_rachel/fixed_gffs/Rd.gff'
        gff_data = in_data(gff_file)
        # print(gff_data)
        rows = []
        for i in range(len(gff_data) - 1):
            chrom_name = gff_data.iloc[i]['seqid']
            prev_cds_stop = gff_data.iloc[i]['end']
            prev_cds_stop = int(prev_cds_stop) - 1
            next_cds_start = gff_data.iloc[i + 1]['start']
            cds1 = gff_data.iloc[i]['Gene']
            cds2 = gff_data.iloc[i + 1]['Gene']
            parent_ids = cds1 + '_____' + cds2
            strand_cds1 = gff_data.iloc[i]['strand']
            strand_cds2 = gff_data.iloc[i + 1]['strand']
            strands = strand_cds1 + '_' + strand_cds2
            col_values = pd.Series({'Chr': chrom_name, 'Start': prev_cds_stop, 'End': next_cds_start,
                                    'Flanking CDS IDs': parent_ids, 'Score': '.', 'Strands': strands})
            rows.append(col_values)

        new_file = pd.concat(rows, axis=1).T
        new_file = new_file.ix[new_file['End'] > new_file['Start']]
        new_file = new_file.reindex(columns=['Chr', 'Start', 'End', 'Flanking CDS IDs', 'Score', 'Strands'])
        intercds_filepath = '/Users/tanayajadhav/drexel_internship/intercds/' + strain_name + '_intercds.bed'
        new_file.to_csv(intercds_filepath, sep='\t', index=False, header=False)



        phobos_gff_dir = '/Users/tanayajadhav/drexel_internship/phobos_repeats'
        intercds_bed_dir = '/Users/tanayajadhav/drexel_internship/intercds/'
        intercds_intersect_dir = '/Users/tanayajadhav/drexel_internship/intercds_intersects'
        intersect_phobos_intercds(phobos_gff_dir, intercds_bed_dir, intercds_intersect_dir)

if __name__ == '__main__':
    main()
