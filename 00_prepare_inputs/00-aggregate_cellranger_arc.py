

import pandas as pd
import numpy as np

df = pd.read_csv('/nemo/lab/gandhis/home/shared/ASAP_singleCell/multiome_cellranger_arc/Samples_NewGenome/Multiome_finished/scripts_cellranger/data/cellranger_arc_aggr.csv', index_col = 0)
paths = list(df['per_barcode_metrics'])

for i in range (0,len(paths)):
    sample=i+1
    if i==0:
        df = pd.read_csv(paths[i])
        print(np.shape(df))
    else:
        df2 = pd.read_csv(paths[i])
        df2['barcode'] = df2['barcode'].str.split('-', expand=True)[0]+'-%i'%sample
        df2['atac_barcode'] = df2['atac_barcode'].str.split('-', expand=True)[0]+'-%i'%sample
        df2['gex_barcode'] = df2['barcode']
        df=pd.concat([df,df2])
        print(np.shape(df))
print('final shape is ', np.shape(df))
        
df.to_csv('data/per_barcode_metrics_aggr.csv', index=False)

print('All done.')