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
    # since there are only 3/6 of the samples of the count matrix in the metadata
    # the adata object should only have 3 observations
    assert(adata.shape[0]==3 and adata.shape[1]==50), "Adata not read in correctly, shape is wrong"

    ## check if the metadata read in correctly ##
    assert(len(adata.obs.columns)==11), "Metadata not added to adata object"
    
    ## check if the transpose argument works correctly ##
    counts='./test/test_counts_transpose.csv'
    adata=pp.create_anndata(counts,meta, transpose=False)
    assert(adata.shape[0]==3 and adata.shape[1]==50), "Adata not read in correctly when transpose=False, shape is wrong"

    
def test_compute_expected_rpk():
    ad=data_for_testing()
    ctrl_ids=['AG9_R3_Demux2_S248_R1_001','B6_10C_R3_Demux2_S282_R1_001']

    # pseudocount
    ad.X=ad.X+1
    # need to define control first
    ad.uns['control_ids']=ctrl_ids # normally define_ctrl_set() fxn would do this, but i picked arbitrary obs to be ctrl for test
    ad=pp.filter_out_ctrl_set(ad)

    print(ad)
    exp_rpk=pp.compute_expected_rpk(ad)

    print(exp_rpk.shape)
