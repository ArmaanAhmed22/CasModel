from dis import dis
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from math import e

df = pd.read_csv(snakemake.input[0])

COMBINATORIAL_FUNCTION = lambda distance: 1-0.5*e**(-distance/6)

def main():
    pre_df = np.ndarray((df.shape[0],df.shape[0]))
    mask = np.ndarray((df.shape[0],df.shape[0]))
    for i in range(df.shape[0]):
        for j in range(df.shape[0]):
            cur_tolerance = df["Normalized Signal"].iloc[i] * df["Normalized Signal"].iloc[j]
            print(f"i:{i} | j:{j}")
            additional_penalty = 1
            
            nuc_distance = abs(i - j)
            
            
            if nuc_distance == 1:
                additional_penalty*=1-(1-COMBINATORIAL_FUNCTION(nuc_distance)) * 0.8
            elif nuc_distance >= 2 and nuc_distance <= 6:
                additional_penalty*=COMBINATORIAL_FUNCTION(nuc_distance)
            
            cur_tolerance*=additional_penalty
            
            mask[i, j] = i <= j
            
            pre_df[i, j] = cur_tolerance
    
    
    sns.heatmap(pre_df, mask=mask)
    plt.savefig(snakemake.output[0], dpi=600, bbox_inches='tight')        
            

main()