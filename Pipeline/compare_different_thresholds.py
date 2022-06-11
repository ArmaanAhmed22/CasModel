import pandas as pd
from math import e
import numpy as np

import seaborn as sns
import matplotlib.pyplot as plt

from numpy import trapz

df = pd.read_csv(snakemake.input[0])

def get_combinatorial(loc1, loc2):
    distance = abs(loc1-loc2)
    if distance > 6:
        return 1
    penalty = 1-0.5*e**(-distance/6)
    if distance == 0:
        return 1-(1-penalty) * 0.6
    else:
        return penalty
    
    


def get_all_combinations(number_mismatches, length):
    locations = list(range(number_mismatches))
    
    while True:
        yield locations
        
        cur_i = number_mismatches - 1
        while cur_i >= 0:
            if locations[cur_i] == length - 1:
                cur_i-=1
            else:
                locations[cur_i] += 1
                break
        else:
            break


def return_percent_targeted(data, number_mismatches, threshold):
    if number_mismatches == 0:
        return 1
    num_combinations = 0
    pre_array = []
    for cur_combination in get_all_combinations(number_mismatches, data.shape[0]):
        num_combinations += 1
        cur_tolerance = 1
        for current_mismatch_location in cur_combination:
            cur_tolerance*=data["Normalized Signal"].iloc[current_mismatch_location]
        for cur_mismatch_index in range(len(cur_combination)):
            for second_mismatch_index in range(cur_mismatch_index, len(cur_combination)):
                cur_mismatch_location = cur_combination[cur_mismatch_index]
                second_mismatch_location = cur_combination[second_mismatch_index]
                cur_tolerance *=get_combinatorial(cur_mismatch_location, second_mismatch_location)
        
        pre_array.append(cur_tolerance)
    cur_array = np.array(pre_array)
    return (cur_array > threshold).sum() / num_combinations
                

def main():
    dx = 0.01
    NUMBER_MISMATCHES = range(6)
    thresholds = np.arange(0, 1+dx, dx)
    
    _, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 8))
    
    auc_pre_df = {"AUC": [], "Number Mismatches": []}
    
    tolerance_pre_df = {"Percent Tolerance": [], "Threshold": [], "Number Mismatches": []}
    
    for i, cur_mismatch_count in enumerate(NUMBER_MISMATCHES):
        for cur_threshold in thresholds:
            tolerance_pre_df["Percent Tolerance"].append(return_percent_targeted(df, cur_mismatch_count, cur_threshold))
            tolerance_pre_df["Threshold"].append(cur_threshold)
            tolerance_pre_df["Number Mismatches"].append(cur_mismatch_count)
    
    tolerance_df = pd.DataFrame(tolerance_pre_df)
    
    for number_mismatches in tolerance_df["Number Mismatches"].unique():
        percent_tolerances = tolerance_df[tolerance_df["Number Mismatches"] == number_mismatches]["Percent Tolerance"]
        auc_pre_df["AUC"].append(trapz(percent_tolerances, dx=dx))
        auc_pre_df["Number Mismatches"].append(number_mismatches)
    
    
    sns.lineplot(x="Threshold", y="Percent Tolerance", data=tolerance_df, ax=ax, hue="Number Mismatches")
    colors = []
    for mismatch_index in range(len(NUMBER_MISMATCHES)):
        cur_line = ax.get_lines()[mismatch_index]
        x_data = cur_line.get_xydata()[:, 0]
        y_data = cur_line.get_xydata()[:, 1]
        
        print(x_data)
        ax.fill_between(x_data, y_data, color=cur_line._color, alpha=1)
        colors.append(cur_line._color)
        ax.margins(x=0, y=0)
        sns.despine(ax=ax)
    
    plt.savefig(snakemake.output[0], dpi=600, bbox_inches="tight")
    plt.close()
    
    
    auc_df = pd.DataFrame(auc_pre_df)
    sns.barplot(x="Number Mismatches", y="AUC", data=auc_df, palette=colors)
    plt.savefig(snakemake.output[1], dpi=600, bbox_inches="tight")
            


main()