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

def plot_pc(data, explained_variance=None, x=1, y=2, colors=None):
    """
    Plots the principal components of the data.
    
    Parameters:
    -----------
    data : DataFrame
        The data to plot. Must include columns for principal components and treatment.
    explained_variance : array, optional
        The explained variance of the PCA. If provided, axis labels will include these values.
    x : int
        The x-axis component (PC number).
    y : int
        The y-axis component (PC number).
    colors : list, optional
        List of colors to use for different hues. Must match the number of unique values in the 'treatment' column.
        
    Returns:
    --------
    A seaborn.scatterplot object
    """

    g = sns.scatterplot(x=f'PC{x}', y=f'PC{y}', hue='treatment', palette=colors,
                        data=data, style='polarity', s=150)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)
    if explained_variance is not None:
        plt.xlabel(f'PC{x}, {explained_variance[x - 1]*100:.2f}%')
        plt.ylabel(f'PC{y}, {explained_variance[y - 1]*100:.2f}%')
    return g

