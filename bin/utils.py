import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import scipy

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler as SS
from statsmodels.stats.multitest import fdrcorrection
from matplotlib import rc

rc('text', usetex=False)
rc('text.latex', preamble=r'\usepackage{cmbright}')
rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
rc = {'lines.linewidth': 2,
      'axes.labelsize': 18, 
      'axes.titlesize': 18, 
      'axes.facecolor': 'DFDFE5'}
sns.set_context('notebook', rc=rc)
sns.set_style("dark")

mpl.rcParams['xtick.labelsize'] = 16 
mpl.rcParams['ytick.labelsize'] = 16 
mpl.rcParams['legend.fontsize'] = 14

def tidy_up(df):
    """
    Casts data into final format, and also computes the regularized log1TPM and skewness
    
    Parameters:
    -----------
    df : DataFrame
        The data to tidy up
        
    Returns:
    --------
    tidy : DataFrame
        The tidy DataFrame
    """
    # cast the dataframe from matrix shape into tidy format
    tidy = df.reset_index()\
            .melt(id_vars='gene', var_name='Sample', value_name='log1TPM') 
    # annotate metadata:
    tidy['treatment'] = tidy.Sample.str[:3]
    tidy['Cell_ID'] = tidy.Sample.str[:-1]
    tidy['Polarity'] = tidy.Sample.str[-1]
    # compute mean log1TPM per sample
    tidy['sample_mean_log1TPM'] = tidy.groupby('Sample').log1TPM.transform(np.mean).values

    # Correct for bias in log1TPM, which is equal to log1TPM + mean log1TPM for that sample
    # calling the variable RPM 
    
    tidy['log1RPM'] = tidy.log1TPM + tidy.sample_mean_log1TPM
    # split tidy dataframe into two, A and P dataframes
    atidy = tidy[tidy.Polarity == 'A'].set_index(['gene', 'Cell_ID'])
    ptidy = tidy[tidy.Polarity == 'P'].set_index(['gene', 'Cell_ID'])

    # stack the two tidy dataframes so that A and P from the same cell are side by side
    stack = atidy.join(ptidy[['log1RPM', 'sample_mean_log1TPM']], rsuffix=('_posterior'))
    
    # compute skew as regularized Delta / regularized Total counts
    stack['Skew'] = ((stack.log1RPM - stack.log1RPM_posterior) /
                    (stack.log1RPM + stack.log1RPM_posterior))
    
    return stack


def get_mean_skew(df):
    """
    Groups data by treatment and gene, and computes the mean of log1RPM, log1RPM_posterior, and Skew
    
    Parameters:
    -----------
    df : DataFrame
        The stacky DataFrame
        
    """
    mean_skews = df\
                .groupby(['treatment', 'gene'])[['log1RPM', 'log1RPM_posterior', 'Skew']]\
                .mean()\
                .reset_index()
    mean_skews['mean_log1RPM'] = (mean_skews.log1RPM + mean_skews.log1RPM_posterior) / 2
    return mean_skews


def get_pc(df, n_pcs=2, which='pca'):
    """
    Performs PCA on the stacky data, returns the PCA object or the principal components as needed
    
    Parameters:
    -----------
    df : DataFrame
        The stacky DataFrame
    component : int
        The component to return
    n_pcs : int
        The number of principal components to return
    which : str
        Whether to return the PCA object + standard_scaler object or a dataframe with the PC coordinates
    Returns:
    --------
    pca : PCA object
        The PCA object
    pc : DataFrame
        The principal components
    """
    pca = PCA(n_pcs)
    scaler = SS()
    coords = pca.fit_transform(scaler.fit_transform(df.T))
    pc = pd.DataFrame(coords, index=df.columns, columns=[f'PC{i+1}' for i in range(n_pcs)])
    pc['treatment'] = pd.Series(pc.index).astype(str).str[:3].values
    pc['polarity'] = pd.Series(pc.index).astype(str).str[-1].values
    if which == 'pca':
        return scaler, pca
    elif which == 'pc':
        return pc

 
def filter_genes(df, condition='G04'):
    # Filter to keep only genes that are present in both 'con' and condition
    filtered_genes = df.groupby('gene').filter(lambda x: set(['con', condition]).issubset(x['treatment']))

    # Final filtering of stacky to include only the desired genes
    df_G04 = df[df['gene'].isin(filtered_genes['gene'].unique())]

    # Further filter to keep only rows where treatment is 'con' or 'G04'
    df_G04 = df_G04[df_G04['treatment'].isin(['con', condition])]
    return df_G04


def sig_test(df, condition='G04', qval_threshold=0.1, pc_index=1):
    """
    Tests condition 'condition' vs 'control' for skewness in the stacky DataFrame

    Parameters:
    -----------
    condition : str
        The condition to test against control
    df : DataFrame
        The stacky DataFrame
    qval_threshold : float
        The threshold for determining significance based on q-values
    pc_index : int
        The index of the principal component to use for selecting genes (1-based)
        
    Returns:
    --------
    df : DataFrame
        The stacky DataFrame with p-values and q-values added
    """
    df_G04 = filter_genes(df, condition)
            
    # Display the number of unique genes
    print('Number of genes: ', len(df_G04['gene'].unique()))
                        
    # Use PCA to select genes that are likely to have a significant skew value
    tmp = df_G04.pivot(index='gene', columns='Sample', values='Skew')
    scaler, pca = get_pc(tmp)
    
    # Convert pc_index from 1-based to 0-based for indexing pca.components_
    if pc_index > pca.components_.shape[0]:
        raise ValueError(f"pc_index {pc_index} is greater than the number of PCs computed ({pca.components_.shape[0]})")
    
    wanted = tmp.index[np.argsort(np.abs(pca.components_[pc_index-1]))[-500:]]

    df_G04 = df_G04[df_G04.gene.isin(wanted)]

    PV = []
    G = []
    for gene, group in df_G04.reset_index()[df_G04.reset_index().gene.isin(wanted)].groupby('gene'):
        npolar = group[group.treatment == condition].Skew
        polar = group[group.treatment == 'con'].Skew
        t, p = scipy.stats.ttest_ind(npolar, polar, equal_var=False)
        PV += [p]
        G += [gene]
        
    # Store results in a DataFrame
    pvals = pd.DataFrame([PV, G], index=['p_value', 'gene']).T.set_index('gene')
    pvals['q_value'] = fdrcorrection(pvals.p_value.values)[1]
    pvals['sig'] = (pvals.q_value < qval_threshold).values
    
    df_G04 = df_G04.join(pvals, on='gene')
    return df_G04



def ecdf(x):
    """ Plot the fraction of values that are >= the x-axis value"""
    xs = np.sort(x)
    y = [np.sum(xi > x).sum() / len(x) for xi in xs]
    return xs, y
