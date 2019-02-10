#!/usr/bin/env python3
# Authors:  Rachel Ehrlich and Tanaya Jadhav

from glob import iglob, glob
from subprocess import Popen, PIPE, call, check_output
import os
import sys
import re
from collections import defaultdict, Counter, namedtuple
import argparse

import pandas as pd
import numpy as np


def run_phobos(fasta_dir, phobos_gff_dir):
    # print('in run_phobos')
    for fasta_path in iglob(fasta_dir + '/*fasta'):
        strain_name = fasta_path.split('/')[-1].split('.fasta')[0]
        phobos_gff = phobos_gff_dir + '/' + strain_name + '.repeats.gff'
        cmd = ['phobos', '-M', 'exact',
               '-U', '15', '-u', '2',
               '--minLength_b', '5',
               '--outputFormat', '2',
               fasta_path,
               phobos_gff
               ]
        call(cmd)


def intersect_phobos_prokka(phobos_gff_dir, prokka_gff_dir, intersect_dir):
    for phobos_gff in iglob(phobos_gff_dir + '/*.repeats.gff'):
        strain_name = phobos_gff.split('/')[-1].replace('.repeats.gff', '')
        prokka_gff = prokka_gff_dir + '/' + strain_name + '.gff'
        cmd = ['bedtools', 'intersect',
               '-wao',
               '-a', phobos_gff,
               '-b', prokka_gff]
        intersect_path = intersect_dir + '/' + strain_name + '.intersected_repeats'
        with open(intersect_path, 'w') as f:
            call(cmd, stdout=f)


# def parse_intersect(intersect_dir):
#     for intersect_path in iglob(intersect_dir + '/*.intersected_repeats'):
#         strain_name = intersect_path.split('/')[-1].replace('.intersected_repeats', '')
#         intersect = pd.read_csv(intersect_path, sep='\t', header=None)
#         intersect.rename(columns={'old_name': 'new_name'}, inplace=True)
#         print(intersect.head(5))


def main():
    # print('in main')
    fasta_dir = '/Users/tanayajadhav/drexel_internship/fastas'
    phobos_gff_dir = '/Users/tanayajadhav/drexel_internship/phobos_repeats'
    run_phobos(fasta_dir, phobos_gff_dir)
    # print('back in main')
    prokka_gff_dir = '/Users/tanayajadhav/drexel_internship/from_rachel/fixed_gffs'
    intersect_dir = '/Users/tanayajadhav/drexel_internship/intersects'
    intersect_phobos_prokka(phobos_gff_dir, prokka_gff_dir, intersect_dir)
    #parse_intersect(intersect_dir)

if __name__ == '__main__':
    print('if name == main')
    main()