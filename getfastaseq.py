import pandas as pd
from glob import iglob
from subprocess import call
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def run_getfasta(bedfile, input_fastafile, outputfile_dir):
    clustername = bedfile.split('/')[-1].split('.')[0]
    outputfile_path = outputfile_dir + clustername + '.sequences.fa'
    cmd = ['bedtools', 'getfasta', '-s',
           '-fi', input_fastafile,
           '-bed', bedfile,
           '-fo', outputfile_path]
    call(cmd)


def findquery_and_blast(file):
    genename = file.split('/')[-1].split('.')[0]
    translatedseq_list = []
    with open(file, 'rU') as f:
        for record in SeqIO.parse(f, 'fasta'):
            translatedseq = record.seq.translate(table='11')
            translatedseq_list.append(translatedseq)
    unique_list = list(set(translatedseq_list))
    unique_list.sort(key=len)
    query = unique_list[-1]
    query = SeqRecord(query, id='longest_sequence')
    query_dir = '/Users/tanayajadhav/drexel_internship/blastquery/'
    query_path = query_dir + genename + '.query.fa'
    with open(query_path, 'w') as fh:
        SeqIO.write(query, fh, 'fasta')
    runblast(genename, query_path)
    # print(genename, 'done')


def runblast(genename, query_path):
    db_path = '/Users/tanayajadhav/drexel_internship/blastdatabase/multidb.fa'
    br_path = '/Users/tanayajadhav/drexel_internship/blastresults/' + genename + '.blastresult.txt'
    cmd = ['tblastn',
           '-query', query_path,
           '-db', db_path,
           '-outfmt', '6',
           '-out', br_path]
    call(cmd)


def main():
    infile = '/Users/tanayajadhav/drexel_internship/combined_gff_withclusters.gff'
    filtered_gff = pd.read_csv(infile, sep='\t', header=None, comment='#')
    colnames = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase',
                'attributes', 'geneid', 'cluster']
    filtered_gff.rename(columns=dict(zip(filtered_gff.columns[:], colnames)), inplace=True)
    bed = filtered_gff[['seqid', 'start', 'end', 'strand', 'cluster']].copy()

    bed['start'] = bed['start'] - 1
    bed['end'] = bed['end']
    bed['name'] = '.'
    bed['score'] = '.'
    bed = bed[['seqid', 'start', 'end', 'name', 'score', 'strand', 'cluster']]
    grouped_bed = bed.groupby('cluster')
    bedfile_dir = '/Users/tanayajadhav/drexel_internship/cluster_bedfiles/'
    for name, group in grouped_bed:
        bedfile_path = bedfile_dir + name + '.bed'
        columns_to_write = ['seqid', 'start', 'end', 'name', 'score', 'strand']
        group.to_csv(bedfile_path, sep='\t', header=0, index=0, columns=columns_to_write)

    cluster_fasta_dir = '/Users/tanayajadhav/drexel_internship/cluster_fastafiles/'
    input_fastafile = '//Users/tanayajadhav/drexel_internship/blastdatabase/multidb.fa'
    for bedfile in iglob(bedfile_dir + '*.bed'):
        run_getfasta(bedfile, input_fastafile, cluster_fasta_dir)

    for file in iglob(cluster_fasta_dir + '*.sequences.fa'):
        findquery_and_blast(file)








if __name__ == '__main__':
    main()
