import pandas as pd
import numpy as np
import random

AA2NA = {
    "A": list("GCT,GCC,GCA,GCG".split(",")),
    "R": list("CGT,CGC,CGA,CGG,AGA,AGG".split(",")),
    "N": list("AAT,AAC".split(",")),
    "D": list("GAT,GAC".split(",")),
    "C": list("TGT,TGC".split(",")),
    "Q": list("CAA,CAG".split(",")),
    "E": list("GAA,GAG".split(",")),
    "G": list("GGT,GGC,GGA,GGG".split(",")),
    "H": list("CAT,CAC".split(",")),
    "I": list("ATT,ATC,ATA".split(",")),
    "L": list("TTA,TTG,CTT,CTC,CTA,CTG".split(",")),
    "K": list("AAA,AAG".split(",")),
    "M": list("ATG".split(",")),
    "F": list("TTT,TTC".split(",")),
    "P": list("CCT,CCC,CCA,CCG".split(",")),
    "S": list("TCT,TCC,TCA,TCG,AGT,AGC".split(",")),
    "T": list("ACT,ACC,ACA,ACG".split(",")),
    "W": list("TGG".split(",")),
    "Y": list("TAT,TAC".split(",")),
    "V": list("GTT,GTC,GTA,GTG".split(",")),
    "*": list("TAA,TGA,TAG".split(","))
}

def aa2na(seq):
    na_seq = [random.choice(AA2NA.get(c, ["---"])) for c in seq]
    return "".join(na_seq)

def generate_SS_hits(ad, layer, ctrl_key, ctrl_value, z_cutoff=2):
    # take mean and standard deviation of the control group
    mu=np.mean(ad[ad.obs[ctrl_key]==ctrl_value].layers[layer],axis=0)
    sd=np.std(ad[ad.obs[ctrl_key]==ctrl_value].layers[layer],axis=0)

    #take and save z-score
    ad.layers['Z_{}'.format(ctrl_value)]=ad.layers[layer]-mu/sd

    #take and save SS hits
    ad.obsm['SS_{}_Z-{}'.format(ctrl_value,z_cutoff)]=pd.DataFrame(data=ad.layers['Z_{}'.format(ctrl_value)]>z_cutoff,
                                             index=ad.obs.index,
                                             columns=ad.var.index)
    return None

def generate_peptide_table(adata, map_file):
    #load peptide info 
    pep2gene=pd.read_csv(map_file, index_col=0, header=0)
    p2g=pep2gene.gene.to_dict()
    p2seq=pep2gene.sequence.to_dict()
    del pep2gene
    
    pep_table=pd.DataFrame(index=adata.var_names)
    pep_table['gene']=pep_table.index.map(p2g)
    pep_table['fragment']=[x.split('fragment_')[1] for x in pep_table.index]
    pep_table['seq']=pep_table.index.map(p2seq)
    
    return pep_table


def generate_alanine_lib_fastq(pep_table, out_file):    
    hibit_seqs={} # make a dictionary-- fasta header:fasta seq
    
    for p in pep_table.index:
        #add original seq to the dictionary 
        hibit_seqs['{}_frag{}_ALAscanNULL'.format(pep_table.loc[p,'gene'],pep_table.loc[p,'fragment'])]=a2na(pep_table.loc[p,'seq'])
        
        #do alanine scan of the peptide
        for i in range(49): 
            hibit_seqs['{}_frag{}_ALAscan{}'.format(pep_table.loc[p,'gene'],pep_table.loc[p,'fragment'],i)]=aa2na(pep_table.loc[p,'seq'][: i] + A + pep_table.loc[p,'seq'][: i]) #stops at i-th AA and replaces with alanine
    
    #write dictionary to fasta
    with open(out_file,'w') as f: 
        for header in hibit_seqs.keys():
            f.write('>'+header+'\n')
            f.write(hibit_seqs[header]+'\n')