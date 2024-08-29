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

def waterman(motif_, substitution_, kwargs={}):
    """
    sets the Smith-Waterman alignment, implementation from skbio, to a certain motif

    Parameters
    ----------
    motif_: str
        motif of interest
    substitution_: 2D dict
        Provides the score for each possible substitution of sequence characters. 
        This may be used for protein or nucleotide sequences. The entire set of possible 
        combinations for the relevant sequence type MUST be enumerated in the dict of dicts. 
        This will override match_score and mismatch_score. Required when protein is True. Default is None. 
    kwargs: dict, default=None
        params to pass through `skbio.alignment.StripedSmithWaterman`

    """
    sw=align.StripedSmithWaterman(motif_, protein=True,
                              substitution_matrix=substitution_, gap_open_penalty=1, gap_extend_penalty=4, **kwargs)
    
    return sw

def sliding_sw(sw_, seq_, motif_, window_width_=None, boost_exact_match_=None):
    """
    scores a peptide for a known motif with sliding windows across the sequence. The score is defined as the sum of either (i) the smith-waterman alignment of the sliding window to the motif  
     or (ii) the `boost_exact_match` parameter, when it is set, if there is an exact match between the sliding window and the epitope. 

    Parameters
    ----------
    sw_: `skbio.alignment.StripedSmithWaterman` object 
        defined smith-waterman alignment function from `waterman()` function output
    seq_: str
        peptide string to score
    motif_: str
        motif of interest
    window_width_: int, default = len(motif_)+1 
        width of sliding window to compare to motif  
    boost_exact_match_: int, default=None  
        if exact match, boost_exact_match is added to the sum instead of the SW alignment score

    Returns
    -------
    results: dict
        dictionary of results with 'boot_scores', 't-test' and 'p-val' keys
    """
    if not window_width_:
        window_width_=len(motif_)+1

    scores=[]
    for i in range(len(seq_)-(window_width_-1)):
        xx=seq_[i:i+window_width_]
        if boost_exact_match_: 
            if ''.join(xx[:len(motif_)])==motif_ or ''.join(xx[1:])==motif_:
                scores.append(boost_exact_match_)
            else:
                scores.append(sw_(xx)['optimal_alignment_score'])
        else: 
            scores.append(sw_(xx)['optimal_alignment_score'])
    
    return(np.sum(scores))

def boot_compare_scores(n_boots_, library_file_, scores_, n_peptides_=None):
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
    n_peptides_: int, default=len(scores_)
        number of peptides in each bootstrap 
    Returns
    -------
    results: dict
        dictionary of results with 'boot_scores', 't-test' and 'p-val' keys
    
    """ 

    boot_means=[]
    boot_t=[]
    boot_p=[]

    if not n_peptides_:
        n_peptides_=len(scores_)

    lines=open(library_file_,'r').readlines()
    for i in range(n_boots_):
        #get random sample index
        idx=random.sample(range(len(lines/2)), n_peptides_)
        idx=2*np.array(idx)
        pep_boot=[lines[x] for x in idx]

        alignment_score_boots=[sliding_sw(x) for x in pep_boot]
        boot_means.append(np.mean(alignment_score_boots))
        T_boot=ss.ttest_ind(scores_, alignment_score_boots)
        boot_t.append(T_boot[0])
        boot_p.append(T_boot[1])

    results={'boot_scores':boot_means, 't-test':boot_t, 'p-val':boot_p}
    return results

