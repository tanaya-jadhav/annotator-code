#!/bin/bash

repeatfile=/drexel_internship/phobos_pyscript/1080_repeats.gff
annotatedfile=/drexel_internship/FASTA/fixed_gffs/10810.gff
./bin/bedtools2/bin/intersectBed -a $repeatfile \
								 -b $annotatedfile