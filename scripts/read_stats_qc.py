'''
# @ Author: Zhang Tong
# @ Create Time: 2021-01-06 18:16:46
# @ Modified by: Zhang Tong
# @ Modified time: 2021-01-06 18:19:13
# @ Description:
'''
#!/usr/bin/env python 

import sys

## sys.argv[1]: collapsed_read_count_file (2 rows with raw read count and collapsed read count)
## sys.argv[2]: trimmed_read_count_file (1 rows with trimmed read counts)
## sys.argv[3]: mapped_read_count_file (1 rows with raw read count, mapped reads, unique mapping reads)
## sys.argv[4]: sample name 

with(open(sys.argv[1], "r")) as f1:
    lines= f1.readlines()
    line_counter= 0
    for line in lines:
        line_counter += 1
        p=line.strip().split("\t")
        if line_counter ==1:
            raw_reads = float(p[1])
        if line_counter ==2: 
            collapsed_reads = float(p[1])

with(open(sys.argv[2], "r")) as f1:
    lines= f1.readlines()
    line_counter= 0
    for line in lines:
        line_counter += 1
        p = line.strip().split()
        if line_counter ==1:
            trimmed_reads = float(p[1])


with(open(sys.argv[3], "r")) as f1:
    lines= f1.readlines()
    line_counter= 0
    for line in lines:
        line_counter += 1
        p = line.strip().split("\t")
        if line_counter ==1:
            map_raw = float(p[0])
            mapped_reads = float(p[1])
            unique_mapped_reads  = float(p[2])

print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
    sys.argv[4], raw_reads, collapsed_reads, round(collapsed_reads/raw_reads,2),  
    trimmed_reads,  round(trimmed_reads/collapsed_reads, 2),
    mapped_reads, round(mapped_reads/trimmed_reads, 2), 
    unique_mapped_reads, round(unique_mapped_reads/mapped_reads, 2)
))

if trimmed_reads != map_raw:
    print("## WARNING: reads used for mapping are not trimmed reads!!!")

