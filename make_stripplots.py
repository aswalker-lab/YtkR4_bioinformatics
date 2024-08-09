# -*- coding: utf-8 -*-
"""
Created on Fri Aug  9 18:45:06 2024

@author: Allison Walker
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

protein0_name = "YtkR4"
protein1_name = "YtkR2"
protein2_name = "YtkR5"
protein1_results_file = open("YtkR2_blastp_other_hits.txt")
protein2_results_file = open("YtkR5_blastp_other_hits.txt")
SSN_cluster_file = open("blastp_YktR5_SSN_clusters.txt")

def convertListItemsInt(list_to_convert):
    return [int(i) for i in list_to_convert]

protein1_df = pd.read_csv(protein1_results_file, header=None, names=["genome","hit_dist_p1"])
protein2_df = pd.read_csv(protein2_results_file, header=None, names=["genome2","hit_dist_p2"])
combined_df = pd.concat([protein1_df, protein2_df], axis=1)

no_protein1_hits = combined_df.iloc[np.where(combined_df['hit_dist_p1'] == "no_hit")]#.value_counts()["no_hit"]
no_protein2_hits = combined_df.iloc[np.where(combined_df['hit_dist_p2'] == "no_hit")]#.value_counts()["no_hit"]
no_hits = no_protein1_hits.iloc[np.where(no_protein1_hits['hit_dist_p2'] == "no_hit")]

protein1_hits = combined_df.iloc[np.where(combined_df['hit_dist_p1'] != "no_hit")]
protein2_hits = combined_df.iloc[np.where(combined_df['hit_dist_p2'] != "no_hit")]
both_hits=protein1_hits.iloc[np.where(protein1_hits['hit_dist_p1'] != "no_hit")]

nh_protein0_only = len(no_hits)
nh_protein0and1 = len(protein1_hits) - len(both_hits)
nh_protein0and2 = len(protein2_hits) - len(both_hits)
nh_protein0and1and2 = len(both_hits)
print("protein 1 hits")
print(len(protein1_hits))
print("protein 2 hits")
print(len(protein2_hits))
print("no hits")
print(nh_protein0_only)
print(nh_protein0and1)
print(nh_protein0and2)
print(nh_protein0and1and2)

#hit in contig analysis
no_protein1_contig_hits = combined_df.iloc[np.where((combined_df['hit_dist_p1'] == "no_hit") | (combined_df['hit_dist_p1'] == "no_hit_in_contig"))]
no_conting_hits = no_protein1_contig_hits.iloc[np.where((no_protein1_contig_hits['hit_dist_p2'] == "no_hit") | (no_protein1_contig_hits['hit_dist_p2'] == "no_hit_in_contig"))]

protein1_contig_hits = combined_df.iloc[np.where((combined_df['hit_dist_p1'] != "no_hit") & (combined_df['hit_dist_p1'] != "no_hit_in_contig"))]
protein2_contig_hits = combined_df.iloc[np.where((combined_df['hit_dist_p2'] != "no_hit") & (combined_df['hit_dist_p2'] != "no_hit_in_contig"))]
both_contig_hits = protein1_contig_hits.iloc[np.where((protein1_contig_hits['hit_dist_p2'] != "no_hit") & (protein1_contig_hits['hit_dist_p2'] != "no_hit_in_contig"))]
nc_protein0_only = len(no_conting_hits)
nc_protein0and1 = len(protein1_contig_hits) - len(both_contig_hits)
nc_protein0and2 = len(protein2_contig_hits) - len(both_contig_hits)
nc_protein0and1and2 = len(both_contig_hits)
print("same contig results")
print(nc_protein0_only)
print(nc_protein0and1)
print(nc_protein0and2)
print(nc_protein0and1and2)

#box plots
SSN_clusters = {}
for line in SSN_cluster_file:
    SSN_clusters[line.split(",")[0]]= int(line.split(",")[1].replace("\n",""))

print(SSN_clusters['GCA_023887685.1'])    
cluster_label = []
for accession in combined_df["genome"]:
    if accession in SSN_clusters:
        cluster_label.append(SSN_clusters[accession])   
    else:
        cluster_label.append(-2)
 

combined_df['YtkR5_cluster'] = cluster_label
fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(8, 6))
protein1_dists = combined_df.iloc[np.where((combined_df['hit_dist_p1'] != "no_hit") & (combined_df['hit_dist_p1'] != "no_hit_in_contig"))]["hit_dist_p1"]
protein2_dists = combined_df.iloc[np.where((combined_df['hit_dist_p2'] != "no_hit") & (combined_df['hit_dist_p2'] != "no_hit_in_contig"))]["hit_dist_p2"]
protein2_cluster1_dists =combined_df.iloc[np.where((combined_df['hit_dist_p2'] != "no_hit") & (combined_df['hit_dist_p2'] != "no_hit_in_contig") & (combined_df['YtkR5_cluster'] == 1))]["hit_dist_p2"]
protein2_cluster2_dists =combined_df.iloc[np.where((combined_df['hit_dist_p2'] != "no_hit") & (combined_df['hit_dist_p2'] != "no_hit_in_contig") & (combined_df['YtkR5_cluster'] == 2))]["hit_dist_p2"]
protein2_cluster3_dists =combined_df.iloc[np.where((combined_df['hit_dist_p2'] != "no_hit") & (combined_df['hit_dist_p2'] != "no_hit_in_contig") & (combined_df['YtkR5_cluster'] == 3))]["hit_dist_p2"]
protein2_cluster4_dists =combined_df.iloc[np.where((combined_df['hit_dist_p2'] != "no_hit") & (combined_df['hit_dist_p2'] != "no_hit_in_contig") & (combined_df['YtkR5_cluster'] == 4))]["hit_dist_p2"]
#not clusters 1-4 (or -2)
protein2_rare_cluster_dists =combined_df.iloc[np.where((combined_df['hit_dist_p2'] != "no_hit") & (combined_df['hit_dist_p2'] != "no_hit_in_contig") & (combined_df['YtkR5_cluster'] != 1) \
                                                       &  (combined_df['YtkR5_cluster'] != 2) &  (combined_df['YtkR5_cluster'] != 3) &  (combined_df['YtkR5_cluster'] != 4) &  (combined_df['YtkR5_cluster'] != -2))]["hit_dist_p2"]

protein1_dists =convertListItemsInt(protein1_dists)
protein2_dists =convertListItemsInt(protein2_dists)
protein2_cluster1_dists = convertListItemsInt(protein2_cluster1_dists)
protein2_cluster2_dists = convertListItemsInt(protein2_cluster2_dists)
protein2_cluster3_dists = convertListItemsInt(protein2_cluster3_dists)
protein2_cluster4_dists = convertListItemsInt(protein2_cluster4_dists)
protein2_rare_cluster_dists = convertListItemsInt(protein2_rare_cluster_dists)
data1 = [protein1_dists, protein2_dists, protein2_cluster1_dists, protein2_cluster2_dists, protein2_cluster3_dists, protein2_cluster4_dists, protein2_rare_cluster_dists]

all_dist_data = []
for i in protein1_dists:
    all_dist_data.append(['YtkR2', i])
for i in protein2_cluster1_dists:
    all_dist_data.append(['YtkR5_1', i])
for i in protein2_cluster2_dists:
    all_dist_data.append(['YtkR5_2', i])
for i in protein2_cluster3_dists:
    all_dist_data.append(['YtkR5_3', i])
for i in protein2_cluster4_dists:
    all_dist_data.append(['YtkR5_4', i])
for i in protein2_rare_cluster_dists:
    all_dist_data.append(['YtkR5_minor', i])
    
all_dist_df = pd.DataFrame(all_dist_data, columns=['protein_type','distance'])
print(all_dist_df)

pos1 = [1, 2, 3, 4, 5, 6, 7]
data2 = [ protein2_cluster1_dists, protein2_cluster2_dists, protein2_cluster3_dists, protein2_cluster4_dists]
pos2 = [1, 2, 3, 4]
break_loc = 40000
print(all_dist_df['distance'])
filtered_df = all_dist_df[all_dist_df['distance'].astype(float) >=break_loc]
colors = [sns.xkcd_rgb["black"],sns.xkcd_rgb["red"],sns.xkcd_rgb["blue"],sns.xkcd_rgb["orange"],sns.xkcd_rgb["green"],sns.xkcd_rgb["sky blue"]]
sns.stripplot(data=all_dist_df, x='protein_type', y='distance', ax=axs,palette=colors,s=7)
d = .015  # how big to make the diagonal lines in axes coordinates



axs.set_ylabel('distance to YtkR4',fontsize=16)
axs.set_yscale('log')
axs.set_xlabel('protein',fontsize=16)
plt.xticks(fontsize=16,rotation=45,ha = 'right')
axs.tick_params(axis='y',labelsize=16)

fig.tight_layout()
fig.savefig("distance_to_YtkR4_strip_plot.pdf")
