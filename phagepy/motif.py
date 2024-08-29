import skbio.alignment as align
import pandas as pd
import numpy as np
import scipy.stats as ss
import random

dayhoff_matrix={'a':['C'],
                'b':['A','G','P','S','T'], 
                'c':['D','E','N','Q','B','Z'],
                'd':['H','K','R'],
                'e':['I','L','M','V','J'],
                'f':['F','W','Y'],
               'other':['X','*']}

def make_dayhoff_sub_matrix():
    dayhoff_map={}
    for x in dayhoff_matrix.keys():
        for y in dayhoff_matrix[x]:
            dayhoff_map[y]=x
    
    aa=['C','A','G','P','S','T','D','E','N','Q','H','K','R','I','L','M','V','F','W','Y','B','Z','J','X','*']
    sub_matrix={}
    for a in aa: 
        mini_dict={}
        for b in aa:
            if b==a:
                mini_dict[b]=25
            else:
                if dayhoff_map[a]==dayhoff_map[b]:
                    mini_dict[b]=10
                else:
                    mini_dict[b]=-25
        sub_matrix[a]=mini_dict
    return sub_matrix

def waterman(motif_, substitution_):
    """
    
    """
    sw=align.StripedSmithWaterman(motif_, protein=True,
                              substitution_matrix=substitution_, gap_open_penalty=1, gap_extend_penalty=4, )
    
    return sw

def sliding_sw(seq_, motif_):
    scores=[]
    for i in range(len(seq_)-3):
        xx=seq_[i:i+4]
        if ''.join(xx[:3])=='FGD' or ''.join(xx[1:])=='FGD':
            scores.append(500)
        else:
            scores.append(sw(xx)['optimal_alignment_score'])
    
    return(np.sum(scores))

def boot_compare_scores(n_boots_, library_file_, scores_):
    """
    Bootstrap comparsion between a given set of epitope scores and epitope scores from bootstrap replicates of random peptides from a given library

    Parameters
    ----------
    n_boots_: int
        number of bootstrap replicates to run
    library_file_: str
        path to fastq file of peptide sequences to use as background distribution of epitope scoring 
    scores_: list
        list of epitope scores that you will compare background distribution to 

    Returns
    -------
    results: dict
        dictionary of results with 'boot_scores', 't-test' and 'p-val' keys
    
    """ 

    boot_means=[]
    boot_t=[]
    boot_p=[]

    lines=open(library_file_,'r').readlines()
    for i in range(n_boots_):
        #get random sample index
        idx=random.sample(range(482672), 91)
        idx=2*np.array(idx)
        pep_boot=[lines[x] for x in idx]

        alignment_score_boots=[sliding_sw(x) for x in pep_boot]
        boot_means.append(np.mean(alignment_score_boots))
        T_boot=ss.ttest_ind(scores_, alignment_score_boots)
        boot_t.append(T_boot[0])
        boot_p.append(T_boot[1])

    results={'boot_scores':boot_means, 't-test':boot_t, 'p-val':boot_p}
    return results

