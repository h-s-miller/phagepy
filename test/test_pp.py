import pytest
from phagepy import *
#from read_write import create_anndata
import numpy as np
import scanpy as sc

def data_for_testing() -> adata:
    counts='test_counts.csv' # a 50 peptide x 5 obs count matrix (from mouseome)
    meta='test_meta.csv' # contains 3/5 of the samples in the count matrix
    adata=create_anndata(counts, meta)

    return adata

def test_read_data():
    counts='test_counts.csv' # a 50 peptide x 5 obs count matrix (from mouseome)
    meta='test_meta.csv' # contains 3/5 of the samples in the count matrix 

    ## read data ##
    adata=create_anndata(counts, meta)

    ## check if shape is correct ##
    # since there are only 3/6 of the samples of the count matrix in the metadata
    # the adata object should only have 3 observations
    assert(adata.shape[0]==3 and adata.shape[1]==50), "Adata not read in correctly, shape is wrong"

    ## check if the metadata read in correctly ##
    assert(len(adata.obs.columns)==11), "Metadata not added to adata object"
    
    ## check if the transpose argument works correctly ##
    counts='test_counts_transpose.csv'
    adata=create_anndata(counts,meta, transpose)
    assert(adata.shape[0]==3 and adata.shape[1]==50), "Adata not read in correctly when transpose=False, shape is wrong"

    
def test_compute_expected_rpk():
    ad=data_for_testing()

    ctrl_ids=['AG9_R3_Demux2_S248_R1_001','B6_10C_R3_Demux2_S282_R1_001']

    exp_rpk=compute_expected_rpk(ad,ctrl_ids)

    print(exp_rpk.shape)
