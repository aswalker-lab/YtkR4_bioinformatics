# -*- coding: utf-8 -*-
"""
Created on Thu Mar  2 11:22:01 2023

@author: Allison Walker
"""
import xml.etree.ElementTree as ET
import urllib.request
import subprocess
import json
import os
import time


#make an output folder to save the genomes in
if not os.path.isdir("genomes"):
    os.mkdir("genomes")
    
#text file containinb blastp results    
blast_results = open("Ytk_output/blastp_results.txt")


gene_locs = {}
gene_locs_out = open("gene_locations.txt",'w')
downloaded_genomes = 0
failed_downloads = 0
for line in blast_results:
    if "#" in line:
        continue
    split_line = line.split("\t")
    protein_accession = split_line[1]
    print(protein_accession)
    #download
    download_process = subprocess.run(["datasets","download","gene","accession",protein_accession],stdout=subprocess.PIPE )
    unzip_process = subprocess.run(["unzip","ncbi_dataset.zip"],stdout=subprocess.PIPE)
    if not os.path.isfile("ncbi_dataset/data/annotation_report.jsonl"):
        failed_downloads += 1
        delete_process = subprocess.run(["rm", "ncbi_dataset.zip"])
        delete_process = subprocess.run(["rm", "README.md"])
        delete_process = subprocess.run(["rm","-r", "ncbi_dataset"])
        continue
    data_report = open("ncbi_dataset/data/annotation_report.jsonl")
    for data_text in data_report:
        genome_accession = data_text[data_text.find("assemblyAccession")+len("assemblyAccession")+3:len(data_text)]
        genome_accession = genome_accession[0:genome_accession.find(",")-1]
        gene_begin = data_text[data_text.find("begin")+len("begin")+3:len(data_text)]
        gene_begin = int(gene_begin[0:gene_begin.find(",")-1])
        gene_end = data_text[data_text.find("end")+len("end")+3:len(data_text)]
        gene_end = int(gene_end[0:gene_end.find(",")-1])
        gene_locs[genome_accession] = (gene_begin, gene_end)
    print(protein_accession)
    print(genome_accession)
    print(gene_locs[genome_accession])
    delete_process = subprocess.run(["rm", "ncbi_dataset.zip"])
    delete_process = subprocess.run(["rm", "README.md"])
    delete_process = subprocess.run(["rm","-r", "ncbi_dataset"])
    print(genome_accession)
    print(["datasets","download","genome","accession",genome_accession,"--include","gbff","--dehydrated"])
    download_process2 = subprocess.run(["datasets","download","genome","accession",genome_accession,"--include","gbff"],stdout=subprocess.PIPE )
    unzip_process = subprocess.run(["unzip","ncbi_dataset.zip"],stdout=subprocess.PIPE)
    if os.path.isfile("ncbi_dataset/data/"+genome_accession+"/genomic.gbff"):
        downloaded_genomes += 1
        gene_locs_out.write(genome_accession + ",")
        gene_locs_out.write(str(gene_locs[genome_accession][0]) + ",")
        gene_locs_out.write(str(gene_locs[genome_accession][1]) + "\n")
        #move and rename files
        move_file = subprocess.run(["mv","ncbi_dataset/data/"+genome_accession+"/genomic.gbff", "genomes/"+genome_accession+".gbff"])
    else:
        failed_downloads += 1
    #delete other files
    delete_process = subprocess.run(["rm", "ncbi_dataset.zip"])
    delete_process = subprocess.run(["rm", "README.md"])
    delete_process = subprocess.run(["rm","-r", "ncbi_dataset"])
    

print("successful downloads: " + str(downloaded_genomes))
print("failed downloads: " + str(failed_downloads))
gene_locs_out.close()
    
