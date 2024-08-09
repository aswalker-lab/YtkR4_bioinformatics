# -*- coding: utf-8 -*-
"""
Created on Wed Jun 21 13:45:13 2023

@author: Allison Walker
"""
import os

#input directory where genomes are stored
input_dir = "genomes/"
#output file to write circular classifications
outfile = open("genomes_circular_status.txt",'w')
for file in os.listdir(input_dir):
    if ".gbff" not in file:
        continue
    line_count = 0
    cicrular = ""
    complete = False
    length = 0
    for line in open(input_dir + file):
        if line_count == 0:
            split_line = line.split()
            circular = split_line[5]
            if "CP" in line:
                complete = True
                
            length = split_line[2]
        if line_count == 1:
            outfile.write(file + "," + circular + "," + str(complete) + "," + length + "\n")
            break
        line_count += 1