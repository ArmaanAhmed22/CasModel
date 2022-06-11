import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


df = pd.read_csv(snakemake.input[0])


def main():
    df["Normalized Signal"] = df["Signal"] / df["Signal"].max()
    df.drop(["Signal"], axis=1, inplace=True)
    print(df)
    df.drop(df[df["Location"] == 0].index, axis=0, inplace=True)
    
    df.to_csv(snakemake.output[0], index=False)
    

    sns.heatmap(data=df["Normalized Signal"].to_numpy().reshape((1,df.shape[0])))
    plt.savefig(snakemake.output[1], dpi=600, bbox_inches="tight")
    

main()