import pandas as pd
from glob import iglob
from subprocess import call


file_path = '/Users/tanayajadhav/drexel_internship/blastresults/sdaA_3	group_876	sdaA.blastresult.txt'
blast_data = pd.read_csv(file_path, sep='\t', header=None)
print(blast_data)