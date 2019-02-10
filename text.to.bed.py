#!/usr/bin/env python3
# Author: Tanaya Jadhav

import pandas as pd


#def main():
repeats_file = '/Users/tanayajadhav/drexel_internship/phobos_repeats/10810.repeats.gff'
with open(repeats_file, 'r+') as f:
    f_lines = f.readlines()
    f.seek(0)
    for line in f_lines:
        if not line.startswith('#'):
            f.write(line)
        f.truncate()



# if __name__ == '__main__':
#     print('if name == main')
#     main()