import matplotlib.pyplot as plt
import pandas as pd
from sklearn.decomposition import PCA
import matplotlib as mpl
import seaborn as sns

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
                       zoom=False, zoom_xlim=None, zoom_ylim=None):
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
                  s= 50)
    
    # if zoom
    if zoom:
        ax.set_xlim(zoom_xlim)
        ax.set_ylim(zoom_ylim)
    
    if print_title:
        ax.set_title(print_title,fontsize=25)
    
    if save:
        plt.savefig(save_dir+save_title)
    
    plt.show()
    
    
def plot_violin(ad, peptides, g_, obs_key, layer='FC_over_AG',obs_colors=None, save=False,save_title=None, save_dir=None):
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
    layer: str, default='FC_over_AG'
        which adata.layers data you want to use for the plot
    obs_colors: list, optional
        if you want to color samples by a category, this provides the colors. needs to be a list or 
        dict of length equal to the num of category values in the obs_key comman (eg, obs_key='disease',
        obs_colors=['blue','red'])
     save: bool, default=False
         directs matplotlib to save plot
     save_title: str, optional
         title for plot
    

         
    Returns
    -------
    violin plot, which is optionally saved
    
    """

    title_=peptides.loc[g_,'gene']+'_fragment-'+str(peptides.loc[g_,'fragment'])
    
    df=pd.DataFrame(data=ad[:,ad.var.index==g_].layers[layer], 
              index=ad.obs.index, columns=['rpk'])
    df[obs_key]=ad.obs[obs_key]

    if not obs_colors: 
        c1='#33FFB3'
        c2='#FF4D59'

    fig, ax =plt.subplots(figsize=(8,6))
    sns.violinplot(data=df, x='cohort', y='rpk',palette=['#33FFB3','#FF4D59'], width=0.8, ax=ax, scale='count')
    sns.stripplot(data=df, x='cohort', y='rpk', size=10, ax=ax, palette=['#0057D9'], jitter=True)
    
    
    ax.set_xlabel(obs_key, fontsize=15)
    ax.set_ylabel(layer, fontsize=15)
    ax.set_title(title_, fontsize=15)
    

    if save:
        plt.savefig(save_dir+save_title)

