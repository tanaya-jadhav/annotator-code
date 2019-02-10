import pandas as pd


blastresultfile = '/Users/tanayajadhav/drexel_internship/blastresults/znuC.blastresult.txt'
cluster_bedfile = '/Users/tanayajadhav/drexel_internship/cluster_bedfiles/znuC.bed'

brf_data = pd.read_csv(blastresultfile, sep='\t', header=None)
brf_colnames = ['query', 'subject', 'pident', 'length', 'mismatch', 'gapopen',
                'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
brf_data.rename(columns=dict(zip(brf_data.columns[:], brf_colnames)), inplace=True)
# print(brf_data)

brf_data['subject_range'] = list(zip(brf_data['sstart'], brf_data['send']))
# print(brf_data)

cbf_data = pd.read_csv(cluster_bedfile, sep='\t', header=None)
cbf_colnames = ['strain', 'start', 'end']
cbf_data.rename(columns=dict(zip(cbf_data.columns[:3], cbf_colnames)), inplace=True)
# print(cbf_data)

cbf_data['gene_range'] = list(zip(cbf_data['start'], cbf_data['end']))
# print(cbf_data)

for index, row in cbf_data.iterrows():
    strain = row['strain']
    gene_range = range(row['gene_range'])
    print(gene_range)
    # for column, hit in brf_data.iterrows():
    #     if strain == hit['subject']:
    #         print(hit)
    break