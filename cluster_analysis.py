##prints all strains where a gene is missing and reports which gene is missing.

### Open gene absence presence file
##find genes in strains that belong to a gene cluster
##analyze with CDS from intersect files

import pandas as pd



def create_dict(df):
    col_list = list(df.columns[10:].values)
    index_list = list(df.index)
    key_list = []
    for i in index_list:
        for j in col_list:
            tup = (i, j)
            # print(tup)
            key_list.append(tup)
    newdict = {}
    for t in key_list:
        x = df.loc[t[0], t[1]]
        newdict[t] = x
    return newdict


file_path = '/Users/tanayajadhav/drexel_internship/gene_presence_absence_paralogs_merged.csv'
gpap = pd.read_csv(file_path, sep=',', header=0, index_col='Gene')
gpap = gpap.fillna('.')
# print(gpap[:3])


gpap_dict = create_dict(gpap)
# print(gpap_dict)

for i in gpap_dict:
    if gpap_dict[i] == '.':
        print(i)




