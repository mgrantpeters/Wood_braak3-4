#!/usr/bin/env python

import scanpy as sc
import muon as mu
import pandas as pd
mudata = mu.read('wnn_scaled_clustered.h5mu')
print(mudata)
mudata.obs['sample_id'] = mudata['rna'].obs['sample_id']
mudata.obs['Group'] = mudata['rna'].obs['Group']
mudata.obs['region'] = mudata['rna'].obs['region']
mudata.obs['batch'] = mudata['rna'].obs['batch']
mudata.obs['RIN'] = mudata['rna'].obs['RIN']
mudata.obs['Sex'] = mudata['rna'].obs['Sex']
mudata['rna'].obsm['X_umap'] = mudata.obsm['X_umap']

sc.pl.umap(mudata, color = ['sample_id', 'Group', 'region', 'batch', 'RIN', 'Sex'], save = '_50PCs_30knn_intmultiome_WNN.png')
genes = ["PLP1",          #Oligodendrocyte markers
	    "MBP",
        "PDGFRA",     #OPCs
	    "BCAS1",
        "CLDN5",      #CLDN5
        "PTPRC",      #Leukocytes
        "CD74",       #Phagocytes
	    "TREM2",
	    "APOE",
        "AQP4",       #Astrocytes
	    "GFAP",    
        "GABRB2",     #Neurons
        "GAD2"]       #Inhibitory neurons]
sc.pl.umap(mudata['rna'], color = genes, save = '_50PCs_30knn_intmultiome_WNN_markers.png')
mudata['rna'].obs.to_csv('rna_obs.csv')
print('Done')
