# Organization

# 2017-12-01

# Already installed bedtools2 and phobos in drexel_intership/apps and linked to bin
# Moved data originally from Rachel to "from_rachel"

cd ~/drexel_internship

mkdir scripts tmp fastas features repeats intersects

mv FASTA from_rachel
cp from_rachel/fixed_gffs/*.gff features/
cp from_rachel/fixed_gffs/*.fasta fasta/

mkdir repeats/original
mv repeats/*_repeats.gff repeats/original


# Dismantled phobos_pyscript by moving *_intersect_repeats.gff to intersects,
# and *.repeats.gff to repeats
# OH NO!  Josh deleted all Tanaya's python stuff! Jerk! 
# Had someone moved *.py to scripts
# Found out about .idea, then messed something up, losing original python script
# Renamed starting intersect table thing script as intersect_table.py

######## PLAN ############
# 1. Make repeats file -- are defaults okay?
# 2. Filter repeats file -- removing what?
# 3. Filter features file -- just CDS?
# 4. Intersect filtered repeats and features -- -wao with -a repeats.gff -b cds.gff
# 5. Tidy it up somehow


######## 1. Done, but lost program.  Redo!

phobos --help > scripts/phobos.help.txt


# Make three versions of the repeats files, just because.
for fasta in $(ls fastas/*.fasta)
do
	prefix=$(basename -s .fasta $fasta)
	phobos -M exact -U 15 -u 2 --minLength_b 5 --outputFormat 2 $fasta repeats/$prefix.repeats.gff
	phobos -M exact -U 15 -u 2 --minLength_b 5 --outputFormat 3 $fasta repeats/$prefix.repeats.txt
	phobos -M exact -U 15 -u 2 --minLength_b 5 --outputFormat 3 $fasta /dev/stdout \
	| sed '/seq-name/,$!d' \
	| grep -v "#" \
	| grep -v seq-name \
	| awk '{OFS="\t"; print $1, $4-1, $5, $0}' \
	>repeats/$prefix.repeats.bed
done

# normal basename: basename file suffix_to_strip, mac osx: basename -s suffix_to_strip file

# Bed version column 5 = repeat unit size, 9 = length, 10 = repeat count

######### 2. Filter these somehow...

# Want to keep high repeat number and not homopolymer
# Figure this stuff out in pandas or R
# Need good filtering syntax

#phobos -M exact -U 15 -u 2 --minLength_b 3 --outputFormat 3 fastas/Rd.fasta /dev/stdout | less -S

# Doing this as part of phobos, because of minLength_b
# Still want to think about this after intersects

########## 3. Filter features file -- just CDS

for gff in $(ls features)
do
	prefix=$(basename -s .gff $gff)
	awk '$3 == "CDS" ' features/$gff > features/$prefix.cds.gff
done


########### 4. Intersect filtered repeats and features

bedtools intersect -wao -a repeats/Rd.repeats.bed -b features/Rd.cds.gff \
| column -t -s $'\t' | less -S
# normal tab character is "\t", mac os x is $'\t' 

bedtools intersect -wao -a repeats/Hi375.repeats.bed -b features/Hi375.cds.gff \
| column -t -s $'\t' | less -S


##### Exercise:  Make an "inverse" bed version of the *.cds.gff

# column 1: chr
# column 2: previous features end coordinate
# column 3: this feature's start coordinate MINUS 1
# column 4: name: BOTH parent IDs (excluding _gene) separated by ___
# column 5: "score", use .
# column 6: strand": BOTH strands, e.g. +_+, -_-, +_-, -_+