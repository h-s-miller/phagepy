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
"*": list("TAA,TGA,TAG".split(","))}

NA2AA = {'GCT': 'A',
 'GCC': 'A',
 'GCA': 'A',
 'GCG': 'A',
 'CGT': 'R',
 'CGC': 'R',
 'CGA': 'R',
 'CGG': 'R',
 'AGA': 'R',
 'AGG': 'R',
 'AAT': 'N',
 'AAC': 'N',
 'GAT': 'D',
 'GAC': 'D',
 'TGT': 'C',
 'TGC': 'C',
 'CAA': 'Q',
 'CAG': 'Q',
 'GAA': 'E',
 'GAG': 'E',
 'GGT': 'G',
 'GGC': 'G',
 'GGA': 'G',
 'GGG': 'G',
 'CAT': 'H',
 'CAC': 'H',
 'ATT': 'I',
 'ATC': 'I',
 'ATA': 'I',
 'TTA': 'L',
 'TTG': 'L',
 'CTT': 'L',
 'CTC': 'L',
 'CTA': 'L',
 'CTG': 'L',
 'AAA': 'K',
 'AAG': 'K',
 'ATG': 'M',
 'TTT': 'F',
 'TTC': 'F',
 'CCT': 'P',
 'CCC': 'P',
 'CCA': 'P',
 'CCG': 'P',
 'TCT': 'S',
 'TCC': 'S',
 'TCA': 'S',
 'TCG': 'S',
 'AGT': 'S',
 'AGC': 'S',
 'ACT': 'T',
 'ACC': 'T',
 'ACA': 'T',
 'ACG': 'T',
 'TGG': 'W',
 'TAT': 'Y',
 'TAC': 'Y',
 'GTT': 'V',
 'GTC': 'V',
 'GTA': 'V',
 'GTG': 'V',
 'TAA': '*',
 'TGA': '*',
 'TAG': '*'}

def translate(seq):
    aa=[]
    for i in range(int(len(seq)/3)):
        aa.append(NA2AA[seq[i*3:(i*3)+3]])
    return ''.join(aa)

def aa2na(seq):
    na_seq = [random.choice(AA2NA.get(c, ["---"])) for c in seq]
    return "".join(na_seq)
    hibit_seqs={} # make a dictionary-- fasta header:fasta seq

def generate_SS_hits(ad, layer, ctrl_key, ctrl_value, z_cutoff=3):
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


def generate_alanine_lib_fastq(pep_table, out_file, n_, len_peps=49):
    hibit_seqs={} # make a dictionary-- fasta header:fasta seq'
    linker_5p='AGCCATCCGCAGTTCGAGAAA'
    linker_3p='GACTACAAGGACGACGATGAT'
    
    for p in pep_table.index:
        #add original seq to the dictionary 
        hibit_seqs['{}_frag{}_ALAscanNULL'.format(pep_table.loc[p,'gene'],pep_table.loc[p,'fragment'])]=aa2na(pep_table.loc[p,'seq'])
        
        #do alanine scan of the peptide
        for i in range(len_peps-(n_-1)): #only go up to the last window of length n_
            if i==0: 
                seq=aa2na('A'*(n_) +pep_table.loc[p,'seq'][n_:])
                #check for restriction sites and replace with equivalent codon 
                if replace_restriction_sites(seq):
                    seq=replace_restriction_sites(seq)
                    
                hibit_seqs['{}_frag{}_ALAscan{}'.format(pep_table.loc[p,'gene'],pep_table.loc[p,'fragment'],i)]=seq
                
            else:
                seq=aa2na(pep_table.loc[p,'seq'][: i] + 'A'*(n_) + pep_table.loc[p,'seq'][i+n_:]) #stops at i-th AA and replaces with alanine
                #check for restriction sites and replace with equivalent codon 
                if replace_restriction_sites(seq):
                    seq=replace_restriction_sites(seq)
                
                hibit_seqs['{}_frag{}_ALAscan{}'.format(pep_table.loc[p,'gene'],pep_table.loc[p,'fragment'],i)]=seq
                
    #write dictionary to fasta
    with open(out_file,'w') as f: 
        for header in hibit_seqs.keys():
            f.write('>'+header+'\n')
            f.write(linker_5p+hibit_seqs[header]+linker_3p+'\n') #add linkers to final sequence
            
def replace_restriction_sites(seq) -> int:
    """
    restriction sites:
    ecoRI='GAATTC'
    hindIII='AAGCTT'
    bamHI='GGATCC'
    XhoI='CTCGAG'
    """
    restriction_sites=['GAATTC','AAGCTT','GGATCC','CTCGAG']
    
    new_seq=None
    for r in restriction_sites: 
        if seq.find(r) != -1:
            ##find codon at restriction site
            x=seq.find(r)
            n=x%3 #find start of codon
            codon=seq[x-n:x-n+3]
            
            ##replace with synonomous mutation
            tmp=AA2NA[NA2AA[codon]][:] #get copy of codon list
            tmp.remove(codon) #remove the one causing restriction site
            new_codon=tmp[0] #replace 
            new_seq=seq[:x-n]+new_codon+seq[x-n+3:]
            break
        
        else:
            continue 
    
    if new_seq:
        return new_seq
    else:
        return None
            

