#!/usr/bin/env python

import scanpy as sc
import muon as mu
import pandas as pd
import numpy as np
import seaborn as sns
print('Reading data')
madata = mu.read('wnn_scaled_clustered.h5mu')
ann = pd.read_csv('split_subclustering/all_subclustered_annotations/all_cells_annotations.csv', index_col = 0)
mm_umap = pd.read_csv('multimodal/md0.5_umap.txt.gz', sep = '\t')

sns.set(rc={'figure.figsize':(8,8)}, font_scale=2)
sns.set_style('white')

print('Adding annotation')
merged_ann = madata['rna'].obs.merge(ann, right_index = True, left_index = True, how = 'left')

madata['rna'].obs['annotation_level_0'] = merged_ann['annotation_level_0']
madata['rna'].obs['annotation_level_1'] = merged_ann['annotation_level_1']
madata['rna'].obs['annotation_level_2'] = merged_ann['annotation_level_2']
madata['rna'].obs['annotation_level_3'] = merged_ann['annotation_level_3']

madata['atac'].obs['annotation_level_0'] = merged_ann['annotation_level_0']
madata['atac'].obs['annotation_level_1'] = merged_ann['annotation_level_1']
madata['atac'].obs['annotation_level_2'] = merged_ann['annotation_level_2']
madata['atac'].obs['annotation_level_3'] = merged_ann['annotation_level_3']

madata.obs['annotation_level_0'] = merged_ann['annotation_level_0'].astype(str)
madata.obs['annotation_level_1'] = merged_ann['annotation_level_1'].astype(str)
madata.obs['annotation_level_2'] = merged_ann['annotation_level_2'].astype(str)
madata.obs['annotation_level_3'] = merged_ann['annotation_level_3'].astype(str)

madata['rna'].obsm['X_umap'] = np.array(mm_umap[['0', '1']])

print('Dropping doublets')
madata = madata[madata['rna'].obs['annotation_level_2']!='drop']
madata = madata[madata['rna'].obs['annotation_level_0'].isnull()==False]

#madata.write("split_subclustering/wnn_scaled_subclustered.h5mu")
print(len(madata['rna'].obs['annotation_level_0'].unique()))
sc.pl.umap(madata['rna'], color = ['annotation_level_0'], size =5, palette=sns.color_palette("husl", len(madata['rna'].obs['annotation_level_0'].unique())), save = '_annotation_level_0.png')
sc.pl.umap(madata['rna'], color = ['annotation_level_1'], size=5, palette=sns.color_palette("husl", len(madata['rna'].obs['annotation_level_1'].unique())), save = '_annotation_level_1.png')
sc.pl.umap(madata['rna'], color = ['annotation_level_2'], size=5,palette=sns.color_palette("husl", len(madata['rna'].obs['annotation_level_2'].unique())), save = '_annotation_level_2.png')
sc.pl.umap(madata['rna'], color = ['annotation_level_3'], size=5,palette=sns.color_palette("husl", len(madata['rna'].obs['annotation_level_3'].unique())), save = '_annotation_level_3.png')
print('Done')


