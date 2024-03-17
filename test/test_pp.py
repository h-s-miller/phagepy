import pytest
import phagepy as pp
import numpy as np
import scanpy as sc

def data_for_testing():
    counts='./test/test_counts.csv' # a 50 peptide x 5 obs count matrix (from mouseome)
    meta='./test/test_meta.csv' # contains 3/5 of the samples in the count matrix
    adata=pp.create_anndata(counts, meta)

    return adata

def test_read_data():
    counts='./test/test_counts.csv' # a 50 peptide x 5 obs count matrix (from mouseome)
    meta='./test/test_meta.csv' # contains 3/5 of the samples in the count matrix 

    ## read data ##
    adata=pp.create_anndata(counts, meta)

    ## check if shape is correct ##
    assert(adata.shape[0]==8 and adata.shape[1]==6), "Adata not read in correctly, shape is wrong"

    ## check if the metadata read in correctly ##
    assert(len(adata.obs.columns)==3), "Metadata not added to adata object"
    
    ## check if the transpose argument works correctly ##
    counts='./test/test_counts_transpose.csv'
    adata=pp.create_anndata(counts,meta, transpose=False)
    assert(adata.shape[0]==8 and adata.shape[1]==6), "Adata not read in correctly when transpose=False, shape is wrong"

def test_read_sparse(): 
    counts='./test/test_counts_sparse.csv' # sparse matrix
    meta='./test/test_meta.csv' # contains 3/5 of the samples in the count matrix 

    ## read data ## 
    adata=pp.create_sparse_anndata(counts, meta)

    ## check if shape is correct ##
    assert(adata.shape[0]==8 and adata.shape[1]==6), "Adata not read in correctly, shape is wrong"

    ## check if the metadata read in correctly ##
    assert(len(adata.obs.columns)==3), "Metadata not added to adata object"


def test_filter_out_ctrl():
    ad=data_for_testing()
    ad=pp.define_ctrl_set_locs(ad, obs_key='group', obs_value='MockIP')
    ad=pp.define_ctrl_set(ad)
    assert(ad.shape[0]==4 and ad.shape[1]==6, "Not filtering out controls properly, shape is wrong")


def test_fold_change():
    ad=data_for_testing()
    ad=pp.define_ctrl_set_locs(ad, obs_key='group', obs_value='MockIP')
    ad=pp.define_ctrl_set(ad)

    ##mean and median should give the same result, ie for each peptide median=mean in the MockIP
    fc=pp.peptide_fold_change(ad, metric='mean')
    assert(np.array_equal(fc.sum(axis=0),[2., 2., 8., 2., 8., 1.]))

    fc=pp.peptide_fold_change(ad, metric='median')
    assert(np.array_equal(fc.sum(axis=0),[2., 2., 8., 2., 8., 1.]))
    
def test_scaling_factor():
    ad=data_for_testing()
    ad=pp.define_ctrl_set_locs(ad, obs_key='group', obs_value='MockIP')
    ad=pp.define_ctrl_set(ad)
    
    scaling=pp.compute_scaling_factor(ad, n=2)
    assert(np.array_equal(scaling,[3., 3., 3., 3.]))
