import matplotlib.pyplot as plt
import pandas as pd
from sklearn.decomposition import PCA
import matplotlib as mpl
import seaborn as sns
import scipy.stats as ss
import numpy as np

def PCA_samples_2D(adata, layer_='FC_over_AG', raw=False):
    """
    runs PCA for visualization on samples
    
    Parameters
    ----------
    adata: anndata
    layer_: str, default: 'FC_over_AG'
        the layer (anndata.layers[`layer_`]) of rpk counts to use for dim. reduction.
    raw: bool, default: False
    
    Returns
    -------
    principalDF: pandas.DataFrame
        dataframe w/ shape samplesx2, with coordinates of each sample in principal component space
    """
    
    pca = PCA(n_components=2)
    if raw is True:
        principalComponents = pca.fit_transform(adata.X)
    else:
        principalComponents = pca.fit_transform(adata.layers[layer_])

    principalDf = pd.DataFrame(data = principalComponents
             , columns = ['principal component 1', 'principal component 2'], index = adata.obs.index)
    return principalDf

def plot_PCA_samples_2D(adata, layer_='FC_over_AG', raw=False, obs_key=None, print_title=None,
                       obs_colors=None, save=False, save_title=None, save_dir=None, 
                       zoom=False, zoom_xlim=None, zoom_ylim=None, kwds={}):
    """
    Runs dimensionality reduction with PCA and visualizes samples in 2D PCA space.
    
    Parameters
    ----------
    adata: anndata.AnnData
    layer_: str, default: 'FC_over_AG'
        the layer (anndata.layers[`layer_`]) of rpk counts to use for dim. reduction.
    raw: bool, default: False
    obs_key: str, optional
        if you want to color samples by a category, obs_key is the category variable in the adata object
        (eg, obs_key='disease' w/ adata.obs['disease']=['healthy','disease'])
    obs_colors: list, optional
        if you want to color samples by a category, this provides the colors. needs to be a list or 
        dict of length equal to the num of category values in the obs_key comman (eg, obs_key='disease',
        obs_colors=['blue','red'])
     save: bool, default=False
         directs matplotlib to save plot
     save_title: str, optional
         title for plot
     save_dir: str, optional
         path where plot is saved
    kwds
        Are passed to :func:`matplotlib.pyplot.matshow`.
         
    Returns
    -------
    2D plot of PCA, which is optionally saved
    
    """
    ## run the PCA ##
    principalDf=PCA_samples_2D(adata, layer_, raw)
    
    ## plot the PCA ## 
    fig = plt.figure(figsize = (8,8))
    ax = fig.add_subplot(1,1,1) 
    ax.set_xlabel('Principal Component 1', fontsize = 15)
    ax.set_ylabel('Principal Component 2', fontsize = 15)
    
    # if grouping
    if obs_key:
        targets = adata.obs[obs_key].unique()
        
        #make color dict if not passed
        if not obs_colors: 
            n=len(targets)
            colors=mpl.cm.Set1.colors[:n] # Set1 is default palette
        else:
            colors=obs_colors
            
        for target, color in zip(targets,colors):
            indicesToKeep = adata.obs.index[adata.obs[obs_key] == target]
            ax.scatter(principalDf.loc[indicesToKeep, 'principal component 1'],
                       principalDf.loc[indicesToKeep, 'principal component 2'],
                       c = color,
                       s = 50)
        
        ax.legend(targets)

    #if not grouping
    else:
        ax.scatter(principalDf.loc[:,'principal component 1'],
                   principalDf.loc[:,'principal component 2'], 
                  c = 'r',
                  s= 50, **kwds)
    
    # if zoom
    if zoom:
        ax.set_xlim(zoom_xlim)
        ax.set_ylim(zoom_ylim)
    
    if print_title:
        ax.set_title(print_title,fontsize=25)
    
    if save:
        plt.savefig(save_dir+save_title)
    
    plt.show()
    
def make_corr_matrix(ad, layer_, metric_, save_table): 
    """
    makes a correlation matrix between all samples on raw counts

    Parameters
    ----------
    ad: anndata.AnnData
        input anndata object

    Returns
    -------
    corr_matrix: pd.DataFrame
        obs x obs matrix of pearson correlation values 
    """
    # define layer 
    if layer_:
        XX=ad.layers[layer_]
    else:
        XX=ad.X

    # fill dataframe of correlations
    corr_matrix=pd.DataFrame(data=None, index=ad.obs.index, columns=ad.obs.index)
    for x in range(len(ad.obs.index)):
        for y in range(len(ad.obs.index)):
            x_loc=ad.obs.index[x]
            y_loc=ad.obs.index[y]
            if metric_=='pearson':
                corr_matrix.loc[x_loc,y_loc]=ss.pearsonr(XX[x,:], XX[y,:])[0]
            elif metric_=='spearman':
                corr_matrix.loc[x_loc,y_loc]=ss.spearmanr(XX[x,:], XX[y,:])[0]  
            elif metric_=='cosine':
                corr_matrix.loc[x_loc,y_loc]=np.dot(XX[x,:],XX[y,:])/(np.linalg.norm(XX[x,:])*np.linalg.norm(XX[y,:]))
            else:
                raise ValueError('Invalid "metric_" value. Please enter either "pearson", "spearman" or "cosine" for metric arg.')
    
    #optional save
    if save_table:
        corr_matrix.to_csv(save_table)

    return corr_matrix


def plot_correlation(ad, metric_='pearson', layer_=None, save_title=None, save_dir=None, save_table=None, kwds={}):
    """
    makes a correlation plot between all samples on counts or another layer of your choosing using either pearson(default), spearman or cosine similarity metrics.

    Parameters
    ----------
    ad: anndata.AnnData
        input anndata object
    metric_: str, ['pearson', 'spearman', 'cosine'], default='pearson'
    layer_: str, default: None
        another adata layer to use instead of X
    save_title: str, optional
        title to save figure
    save_dir: str, optional
        directory to save figure, if save_title is used and no save_dir, it will save in CWD
    save_table: str, optional
        file path to save the correlation values 
    kwds
        Are passed to :func:`matplotlib.pyplot.matshow`.

    Returns
    -------
    None 
    """
    corr_matrix=make_corr_matrix(ad, layer_, metric_, save_table)

    fig,ax=plt.subplots(figsize=(10,10))

    #index_values=[x.split('_')[1] for x in corr_matrix.index]
    index_values=corr_matrix.index


    mat=ax.matshow(corr_matrix.astype(float),**kwds)
    ax.set_title('{} correlation of all samples'.format(metric_),
                fontsize=18)

    ax.set_yticks(range(len(index_values)))
    ax.set_yticklabels(index_values)
    ax.set_xticks(range(len(index_values)))
    ax.set_xticklabels(index_values, rotation=90)


    fig.colorbar(mat, ax=ax, shrink=0.82)
    if save_title: 
        if save_dir:
            plt.savefig(save_dir+save_title)
        else:
            plt.savefig(save_title)
    plt.show()
        
def plot_violin(ad, peptides, g_, obs_key, counts='False', layer='FC_over_AG',obs_colors=None, save=False,save_title=None, save_dir=None, kwds={}):
    """
    Runs dimensionality reduction with PCA and visualizes samples in 2D PCA space.
    
    Parameters
    ----------
    ad: anndata.AnnData
    peptides: pd.DataFrame
        petide info table 
    obs_key: str, optional
        if you want to color samples by a category, obs_key is the category variable in the adata object
        (eg, obs_key='disease' w/ adata.obs['disease']=['healthy','disease'])
    counts: bool, default=False
        if True, plots values from adata.X
    layer: str, default='FC_over_AG'
        which adata.layers data you want to use for the plot
    obs_colors: list, optional
        if you want to color samples by a category, this provides the colors. needs to be a list or 
        dict of length equal to the num of category values in the obs_key comman (eg, obs_key='disease',
        obs_colors=['blue','red'])
    save_title: str, optional
        title for plot
    save_dir: str, optional 
        directory to save plot. if save_dir==None and save_title!=None, then plot is saved in CWD.
    kwds
        Are passed to :func:`matplotlib.pyplot.violinplot`.

         
    Returns
    -------
    violin plot, which is optionally saved
    
    """

    title_=peptides.loc[g_,'gene']+'_fragment-'+str(peptides.loc[g_,'fragment'])
    
    if counts:
        df=pd.DataFrame(data=ad[:,ad.var.index==g_].X, 
              index=ad.obs.index, columns=['rpk'])
    else:
        df=pd.DataFrame(data=ad[:,ad.var.index==g_].layers[layer], index=ad.obs.index, columns=['rpk'])
    df[obs_key]=ad.obs[obs_key]

    if not obs_colors: 
        c1='#664DB2'
        c2='#99B24D'
    else:
        c1,c2=obs_colors

    c3='darkgrey'

    print(df.head())
    fig, ax =plt.subplots(figsize=(8,6))
    sns.violinplot(data=df, x=obs_key, y='rpk', hue=obs_key, width=0.8, palette=[c1,c2],ax=ax, density_norm='count', **kwds)
    sns.stripplot(data=df, x=obs_key, y='rpk', hue=obs_key, size=10, ax=ax,palette=[c3,c3], jitter=True)
    
    
    ax.set_xlabel(obs_key, fontsize=15)
    if counts:
        ax.set_ylabel('rpk', fontsize=15)
    else:
        ax.set_ylabel(layer, fontsize=15)
    ax.set_title(title_, fontsize=15)
    

    if save_title: 
        if save_dir:
            plt.savefig(save_dir+save_title)
        else:
            plt.savefig(save_title)

    plt.show()

def plot_jitter(ad, peptides, g_, obs_key, counts=False, layer='FC_over_AG',obs_colors=None, save_title=None, save_dir=None, kwds={}):
    """
    Makes a jitter plot of a peptide of interest which you can split between an observation category
    
    Parameters
    ----------
    ad: anndata.AnnData
    peptides: pd.DataFrame
        petide info table 
    g_: str
        peptide to plot
    obs_key: str, optional
        if you want to color samples by a category, obs_key is the category variable in the adata object
        (eg, obs_key='disease' w/ adata.obs['disease']=['healthy','disease'])
    counts: bool, default=False
        if True, plot ad.X values, not a layer
    layer: str, default='FC_over_AG'
        which adata.layers data you want to use for the plot
    obs_colors: list, optional
        if you want to color samples by a category, this provides the colors. needs to be a list or 
        dict of length equal to the num of category values in the obs_key comman (eg, obs_key='disease',
        obs_colors=['blue','red'])
     save_title: str, optional
         title for plot
    save_dir: str, optional 
        directory to save plot. if save_dir==None and save_title!=None, then plot is saved in CWD.
    kwds
        Are passed to :func:`seaborn.stripplot`.

         
    Returns
    -------
    jitter plot, which is optionally saved
    
    """

    title_=peptides.loc[g_,'gene']+'_fragment-'+str(peptides.loc[g_,'fragment'])
    
    if counts:
        df=pd.DataFrame(data=ad[:,ad.var.index==g_].X, 
              index=ad.obs.index, columns=['rpk'])
    else:
        df=pd.DataFrame(data=ad[:,ad.var.index==g_].layers[layer], index=ad.obs.index, columns=['rpk'])
    df[obs_key]=ad.obs[obs_key]

    if not obs_colors: 
        c1='#33FFB3'
        c2='#FF4D59'
    else:
        c1,c2=obs_colors

    fig, ax =plt.subplots(figsize=(8,6))
    sns.stripplot(data=df, x=obs_key, hue=obs_key, y='rpk', size=10, ax=ax, palette=[c1,c2], jitter=True, **kwds)
    
    
    ax.set_xlabel(obs_key, fontsize=15)
    if counts:
        ax.set_ylabel('rpk', fontsize=15)
    else:
        ax.set_ylabel(layer, fontsize=15)
    ax.set_title(title_, fontsize=15)
    

    if save_title: 
        if save_dir:
            plt.savefig(save_dir+save_title)
        else:
            plt.savefig(save_title)

    plt.show()
