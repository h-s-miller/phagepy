import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
import scipy.sparse as sparse

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
    adata=sc.AnnData(X=counts_df, obs=obs, var=var, dtype=np.float32)

    # filter on observations that are present in the metdata
    adata=adata[adata.obs.index.isin(meta_df.index)]

    # merge metdata with obs
    adata.obs=adata.obs.merge(meta_df, how='left',
                              left_index=True, right_index=True) #merge on keys present in the adata (how='left')

    del meta_df, counts_df

    return adata


def load_csv_as_sparse(csv_path, transpose, sparse_fmt):
    """
    Load counts as a sparse matrix 

    Parameters
    ----------
    csv_path: str
        path to count data csv

    Returns
    -------
    var_values: list
        list of peptides from the rows of CSV file 
    obs_values: list
        list of samples from the columns  of CSV file
    X: scipy.sparse.csr 
        sparse count matrix in compressed row format (may change to column later, will see)
    
    """ 
    # get peptides and samples
    f=open(csv_path, 'r')
    if transpose: 
        obs_values=f.readline().strip('\n').split(',')[1:] ##just have to get the first line for this, skip first bc its "peptide"
        var_values=pd.read_csv(csv_path, usecols=['peptide']).tolist() #loading only one column isnt super fast, but limits memory usage
    else: 
        var_values=f.readline().strip('\n').split(',')[1:] ##just have to get the first line for this, skip first bc its "peptide"
        obs_values=pd.read_csv(csv_path, usecols=['peptide']).tolist() #loading only one column isnt super fast, but limits memory usage

    f.close()

    # load count data 
    keeps=list(range(len(obs_values))) # need number of columns to keep to load only counts in numpy
    X=np.loadtxt(csv_path, delimiter=',', skiprows=1, usecols=keeps[1:]) #load matrix into numpy

    #convert to sparse
    if sparse_fmt=='csr':
        X=sparse.csr_matrix(X)
    elif sparse_fmt=='csc':
        X=sparse.csc_matrix(X)
    else: 
        raise ValueError('Input for sparse_fmt is invalid, please use either "csr" for compressed row or "csc" for compressed column')

    return var_values, obs_values, X

def create_sparse_anndata(counts_file, metadata_file, transpose=True, sparse_fmt='csr'):
    """
    Load counts data into AnnData with counts in sparse format for speed and memory efficiency 

    Parameters
    ----------
    counts_file: str
        path to count data csv
    metadata_file: str
        path to metadata file. NOTE: the index of metadata file MUST match the samples in count data
    transpose: bool, default=True
        if the data is in var x obs format instead of obs x var (this is default for download from PhageDB), set Transpose=True
    sparse_fmt: str, default='csr'
        set sparse_fmt to 'csr' for scipy.sparse.csr_matrix format or to 'csc' for scipy.sparse.csc_matrix format 

    Returns
    -------
    adata: sc.AnnData
        anndata objects with sparse counts in X, samples in adata.obs and peptides in adata.var
    
    """ 

    #read in metadata
    meta_df=pd.read_csv(metadata_file, index_col=0, header=0) #index of metadata should be the same syntax as the counts df

    #load counts data and convert to anndata
    var_values, obs_values, X_sparse = load_csv_as_sparse(counts_file, transpose, sparse_fmt)
    adata=sc.AnnData(X=X_sparse, obs=obs_values, var=var_values, dtype=np.float32)

    # filter on observations that are present in the metdata
    adata=adata[adata.obs.index.isin(meta_df.index)]

    # merge metdata with obs
    adata.obs=adata.obs.merge(meta_df, how='left',
                              left_index=True, right_index=True) #merge on keys present in the adata (how='left')
    return adata