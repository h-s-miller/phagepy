# import dependencies
import numpy as np
import pandas as pd

def compute_expected_rpk(ad, ctrl_ids) -> np.array:
    """
    Computes expected rpk value for each peptide wrt the user-defined control group. 
    
    For each peptide, p_i, the expected rpk value is defined as,
    p_i =  ( avg(rpk(J,i)) where J is the control group ids if the rpk is non-zero in at least one of the control samples
           {
           ( median(rpk(J,L)) where J is the control group ids and L is the set of all peptides that belong to the same gene as peptide i 

    """

    X=ad.raw.X[ctrl_ids,:]

    exp_rpk=np.empty((X.shape[1],0),float)
    
    for i in X.shape[1]:
        if np.sum(X[:,i]) > 0: #expected rpk for peptides is the mean expression in ctrl group if expression is nonzero
            np.append(exp_rpk,np.mean(X[:,i]))

        elif np.sum(X[:,i]) == 0: 
            # find the other peptides that belong to same gene
            pep=ad.var_names[i]
            g_locs=ad.var.index.get_loc(g2p[p2g[pep]])

            np.append(exp_rpk,np.median(X[:,locs]))

        else:
            raise ValueError('rpk values should not be negative')

    return exp_rpk
            

def peptide_fold_change(ad, ctrl_ids) -> np.array:
    """
    Returns fold-change enrichment over user-defined control group by calculating expected rpk for each peptide in the control group. 

    eg, calculating the fold change of each peptide over AG bead controls 
    """
    exp_rpk=compute_expected_rpk(ad,ctrl_ids)

    return np.divide(ad.X, exp_rpk)


def compute_scaling_factor(ad, ctrl_ids,n) -> np.array:
    
    
    #compute internal control locs 
    internal_ctrl_locs=compute_internal_ctrl_set(ad, ctrl_ids, n)

    # divide median(internal control set of ctrl)/median(interal control set of cohort)
    scaling_factors=np.divide(np.median(ad.raw.X[ctrl_ids,internal_control_locs],axis=0),
                              np.median(ad.X[:,interal_control_locs],axis=0))

    return scaling_factors
    
def compute_internal_ctrl_set(ad, ctrl_ids, n) -> np.array:
    #compute the mean for each peptide in the ctrl set
    means=np.mean(ad.raw.X[ctrl_ids,:], axis=0)

    #and compute the max n peptides
    locs=np.argpartition(means, -n)[-n:]

    return locs
    
def define_ctrl_set(ad, obs_key, obs_value):
    # find obs names of control observations
    control_locs=ad.obs.index[ad.obs[obs_key]==obs_value]

    # add as unstructured annotation (uns) of adata object
    # accessible as adata.uns['control_keys']=control_locs
    adata.uns['control_keys']=control_locs

def filter_out_ctrl_set(ad, obs_key, obs_value):
    try:
        adata.uns['control_keys']
    except NameError:
        raise ValueError('Need to define control set first')


    adata.raw=adata.copy()

    adata=adata[~adata.obs.isin(adata.uns['control_keys'])]

