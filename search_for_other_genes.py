# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 13:58:19 2023

@author: Allison Walker
"""
import subprocess
from Bio import SeqIO
import os
#directory where the fasta files for genomes are located
genome_dir = "fasta_genomes/"
#output directory for windows of 40,000 bp in each direction for downstream antiSMASH analysis
region_dir = "genomes_gene_region_40000/"
other_genes = ["YtkR5.fasta","YtkR2.fasta"]
outfile = open("blast_other_genes.txt","w")
id_threshold = [25.1, 20]

len_threshold = 100
hits = {}
contigs = {}
gene_locs = {}
circular_genome = {}
complete_genome = {}
genome_size = {}
for line in open("genomes_circular_status.txt"):
    split_line = line.split(",")
    genome = split_line[0].replace(".gbff","")
    is_circular = False
    if split_line[1] == "circular":
        is_circular = True
    else:
        is_circular = False
    circular_genome[genome] = is_circular
    complete_genome[genome] = bool(split_line[2])
    genome_size[genome] = int(split_line[3])

count = 0
for g in other_genes:
    hits[g[0:g.find(".")]] = {}
gene_loc_file = open("gene_locations.txt")
for line in gene_loc_file:
    split_line = line.split(",")
    genome = split_line[0]
    print(genome_dir + genome +".fasta")
    
    #get cotigs for YtkR4
    if not os.path.isfile(region_dir + genome + "_region.gbff"):
        continue
    region_in = SeqIO.parse(open(region_dir + genome + "_region.gbff", 'r'),"genbank")
    for rec in region_in:
        contigs[genome] = rec.id
        print(rec.id)
    
    start  = int(split_line[1])
    end = int(split_line[2])
    gene_locs[genome] = (start, end)
    i = 0
    for other_gene in other_genes:
        print(other_gene)
        print(i)
        print(id_threshold[i])
        hits[other_gene[0:other_gene.find(".")]][genome] = []
        blastp_search = subprocess.run(["tblastn", "-subject", genome_dir + genome +".fasta",  "-query", other_gene, \
                                "-outfmt",'7 qacc sacc sgi evalue pident qstart qend sstart send', "-evalue", "1e-5"],stdout=subprocess.PIPE )
       # blastp_search = subprocess.run(["tblastn", "-subject", genome_dir + genome +".fasta",  "-query", other_gene, "-outfmt", "7"],stdout=subprocess.PIPE ) 
        #print("The exit code was: %d" % blastp_search.returncode)
        #TODO: deal with errors?
        output = blastp_search.stdout.decode('utf-8','ignore')
        blast_out = open("other_gene_blast_results.txt", "w")
        for line in output:
            if "FASTA-Reader" in line:
                continue
            
            blast_out.write(line)
        blast_out.close()
        in_blast = open("other_gene_blast_results.txt")
        for line in in_blast:
            if "#" not in line and len(line.split()) > 0:
                split_line = line.split()
                print(line)
                hit = split_line[1]
                percent_id = float(split_line[4])
                start = int(split_line[7])
                end = int(split_line[8])
                if percent_id > id_threshold[i] and abs(start-end) > len_threshold:
                    hits[other_gene[0:other_gene.find(".")]][genome].append((hit, percent_id, start, end))
                    
        i += 1
    count += 1


print("DONE WITH BLAST!!")
print(hits)
#find distance to hit
for gene in hits:
    outfile = open(gene + "_blastp_other_hits.txt", 'w')
    outfile_loc = open(gene + "_blastp_hit_locs.txt", 'w')
    for genome in hits[gene]:
        gene_loc = gene_locs[genome]
        contig = contigs[genome]
        if len(hits[gene][genome]) == 0:
            outfile.write(genome + ",no_hit\n")
            continue
        distances = []
        min_dist = 0
        count = 0
        hit_index = 0
        for hit in hits[gene][genome]:
            hit_contig = hit[0]
            hit_contig = hit_contig[0:hit_contig.find("|")]
            
            
            if hit_contig == contig:
                dist = min(abs(hit[2]-gene_loc[1]), abs(hit[3]- gene_loc[0]))
                distances.append(dist)
                if count == 0:
                    min_dist = dist
                elif dist < min_dist:
                    min_dist = dist
                    hit_index = count
                if circular_genome[genome] == True and complete_genome[genome] == True:
                    gene1_dist_to_start = gene_loc[0]
                    gene1_dist_to_end = genome_size[genome] - gene_loc[1]
                    gene2_dist_to_start = hit[2]
                    gene2_dist_to_end = genome_size[genome] - hit[3]
                    circ_dist = genome_size[genome]
                    if gene1_dist_to_start< gene2_dist_to_start:
                        circ_dist = gene2_dist_to_end + gene1_dist_to_start
                    else:
                        circ_dist = gene1_dist_to_end + gene2_dist_to_start
                    
                    print("genome size: " + str(genome_size[genome]))
                    print("circular distance: " + str(circ_dist))
                    print("linear dist: " + str(dist))
                    print("gene1 start: " + str(gene_loc[0]))
                    print("gene2 start: " + str(hit[2]))
                    if circ_dist < min_dist:
                        min_dist = dist
                count += 1
        if len(distances) == 0:
            outfile.write(genome + ",no_hit_in_contig\n")
            continue
        print(distances)
        outfile.write(genome + "," + str(min_dist) + "\n")
        hit = hits[gene][genome][hit_index]
        hit_contig = hit[0]
        hit_contig = hit_contig[0:hit_contig.find("|")]
        outfile_loc.write(genome + ","+ hit_contig + "," + str(hit[1]) + "," + str(hit[2]) + "," + str(hit[3]) + "\n")
        
    outfile.close()
    outfile_loc.close()