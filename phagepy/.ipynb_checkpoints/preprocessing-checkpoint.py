# import dependencies
import numpy as np
import pandas as pd

def compute_expected_rpk(ad, metric, X_ctrl_key='X_control') -> np.array:
    """
    Computes expected rpk value for each peptide wrt the user-defined control group. 
    
    For each peptide, p_i, the expected rpk value is defined as,
    p_i =  ( median(rpk(J,i)) where J is the control group ids if the rpk is non-zero in at least one of the control samples
           {
           ( median(rpk(J,L)) where J is the control group ids and L is the set of all peptides that belong to the same gene as peptide i 

    """
    X=ad.varm[X_ctrl_key].T

    #check that there are no all zero columns or negatives
    if np.any(np.sum(X,axis=0)==0):
        raise ValueError('Zero expression peptides. Pseudocount +1 before running.')
        
    if np.any(X<0):
        raise ValueError('Check count data, there should not be negative rpk values.')
              
    
    #then can compute mean or median (default=median)
    if metric=='mean':
        exp_rpk=np.mean(X, axis=0)
    elif metric=='median':
        exp_rpk=np.median(X, axis=0)
    else: 
        raise ValueError("Metric value not valid. Please use metric='median' or metric='mean'")
    return exp_rpk
            

def peptide_fold_change(ad, X_ctrl_key='X_control', metric='median') -> np.array:
    """
    Returns fold-change enrichment over user-defined control group 
    by calculating expected rpk for each peptide in the control group. 

    eg, calculating the fold change of each peptide over AG bead controls 
    """
    exp_rpk=compute_expected_rpk(ad, metric, X_ctrl_key)

    return ad.X/ exp_rpk[None,:]


def compute_scaling_factor(ad, n, X_ctrl_key='X_control') -> np.array:
    """
    Computes scaling factor which is the fraction of the median of the internal control set in controls over cases.
    So, if rpk values are inflated in the sample they will be scaled down and vice versa.
    
    The internal control set is a set of peptides that has high mean in the control group and nonzero stdev.
    """
    
    #compute internal control locs 
    internal_ctrl_locs=compute_internal_ctrl_set(ad, n, X_ctrl_key)
    
    #take the median of the sum of internal controls in the control samples
    median_internal_ctrl_exprs=np.median(np.sum(ad.varm[X_ctrl_key].T[:,internal_ctrl_locs],axis=1))
    

    # divide median(sum of internal control set of ctrl)/sum(interal control set of individual sample)
    scaling_factors=np.divide(median_internal_ctrl_exprs,
                              np.sum(ad.X[:,internal_ctrl_locs],axis=1))

    return scaling_factors
    
def compute_internal_ctrl_set(ad, n, X_ctrl_key='X_control') -> np.array:
    X=ad.varm[X_ctrl_key].T
    
    #filter out locs with stdev less than mean
    stdev=np.std(X, axis=0)
    means=np.mean(X, axis=0)
    X=X[:,stdev<=means]
    
    #compute the mean for each peptide in the ctrl set
    means=np.mean(X, axis=0)

    #and compute the max n peptides
    locs=np.argpartition(means, -n)[-n:]

    return locs
    
def define_ctrl_set_locs(ad, obs_key, obs_value, key_ids='control_ids'):
    """
    Defines the control set based on obs_key and obs_value and saves key ids to adata in uns. 
    
    obs_key is the row name of the obs object that contains the obs_value which defines control group. 
    
    Eg, if our obs looks like 
    
    seq_id | group  |  round | patient_id
    ------------------------------------
    seq 1    Disease   2       patient 1
    seq 2    Disease   3       patient 1
    seq 3    Disease   2       patient 2
    seq 4    AG        2       AG1
    seq 5    AG        3       AG1
    
    obs_key='group' and obs_value='AG'
    
    can store it in adata.uns however you like using key_ids param,
    but default is adata.uns['control_ids']
    """
    # find obs names of control observations
    control_locs=ad.obs.index[ad.obs[obs_key]==obs_value]

    # add as unstructured annotation (uns) of adata object
    # accessible as adata.uns['control_keys']=control_locs
    ad.uns[key_ids]=control_locs
    return ad

def define_ctrl_set(ad, key_ids='control_ids', key_X='X_control', filter_out=True):
    """
    Removes control set after it is defined in define_ctrl_set()
    
    Control set data is removed from main adata object, but stored in adat.varm
    """
    
    try:
        ad.uns[key_ids]
    except NameError:
        raise ValueError('Need to define control set first')

    ad.varm[key_X]=ad[ad.uns[key_ids],:].X.T
    
    if filter_out:
        ad=ad[~ad.obs.index.isin(ad.uns[key_ids])]
    return ad
