import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# ---- Load abundance table ----
df = pd.read_csv("3_strains_abundances.csv")  # first col = date, rest = genomes
df["date"] = pd.to_datetime(df["date"])

# ---- Normalize so each row sums to 1 (100%) ----
genomes = df.columns[1:]  # all genome columns
df[genomes] = df[genomes].div(df[genomes].sum(axis=1), axis=0)

# ---- Convert dates to days since first sample ----
df["days"] = (df["date"] - df["date"].min()).dt.days

colors = sns.color_palette("Set2", n_colors=len(genomes))
colors[0] = "#003399"
colors[1] = "#FF9933"
colors[2] = "#669933"

# ---- Make stacked area plot ----
plt.figure(figsize=(11,8))

plt.stackplot(df["days"], df[genomes].T, labels=genomes, colors=colors)

plt.xlim(left=0, right=df["days"].max())
plt.ylim(bottom=0, top=1)

plt.legend(loc="upper left", bbox_to_anchor=(1,1), fontsize=16, title="Strains", title_fontsize=16)
plt.xlabel("Time (Day)", fontsize=18)
plt.ylabel("Relative Abundance", fontsize=18)
plt.title("$\it{Bacteroides\ ovatus}$ strain abundances over time", fontsize=18)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.tight_layout()
plt.savefig("3_strain_abundances.png",dpi=600, bbox_inches="tight")