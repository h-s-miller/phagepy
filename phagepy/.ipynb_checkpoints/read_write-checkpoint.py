import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad

def create_anndata(counts_file, metadata_file, transpose=True):
    #read in counts
    counts_df=pd.read_csv(counts_file, index_col=0, header=0) # should have index names and column names

    #read in metadata
    meta_df=pd.read_csv(metadata_file, index_col=0, header=0) #index of metadata should be the same syntax as the counts df

    #if the data is in var x obs format instead of obs x var (this is default for download from PhageDB), need to transpose
    if transpose:
        counts_df=counts_df.T

    # create observation matrix (obs), will later combine with meta and filter. create variables matrix (var) of all peptides
    obs=pd.DataFrame(data=None, index=counts_df.index)
    var=pd.DataFrame(data=None, index=counts_df.columns)
    
    #create anndata
    adata=sc.AnnData(X=counts_df, obs=obs, var=var)

    # filter on observations that are present in the metdata
    adata=adata[adata.obs.index.isin(meta_df.index)]

    # merge metdata with obs
    adata.obs=adata.obs.merge(meta_df, how='left',
                              left_index=True, right_index=True) #merge on keys present in the adata (how='left')

    del meta_df, counts_df

    return adata


    
