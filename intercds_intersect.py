

from subprocess import call

filepath = '/Users/tanayajadhav/drexel_internship/intercds/Rd_new_intercds.bed'
phobos_gff = '/Users/tanayajadhav/drexel_internship/phobos_repeats/Rd.repeats.gff'
intersect_path = '/Users/tanayajadhav/drexel_internship/intercds_intersects/Rd.intercds_repeats'

cmd = ['bedtools', 'intersect',
       '-wao',
       '-a', phobos_gff,
       '-b', filepath]
with open(intersect_path, 'w') as f:
    call(cmd, stdout=f)
    # print(intersect_path)