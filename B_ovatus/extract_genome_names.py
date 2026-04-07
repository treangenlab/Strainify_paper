import re
import pandas as pd

# -------- Step 1: Extract genome names from clade tree files --------
def extract_genomes_from_tree(treefile):
    with open(treefile) as f:
        text = f.read()
    # Match GCA_######### with optional version (.#)
    genomes = re.findall(r"(GCA_\d{9}\.\d+)", text)
    return genomes

clade1_genomes = extract_genomes_from_tree("strain_1_assembly_names.tree")
clade2_genomes = extract_genomes_from_tree("strain_2_assembly_names.tree")
clade3_genomes = extract_genomes_from_tree("strain_3_assembly_names.tree")

# -------- Step 2: Load abundance CSV --------
df = pd.read_csv(
    "/home/Users/rl152/Strainify/longitudinal/Bacteroides_ovatus/results_all/abundance_reordered.csv"
)  # first col = date, others = genome IDs

# -------- Step 3: Map tree accessions to actual CSV columns --------
def map_to_csv_columns(accessions, df_columns):
    mapped = []
    for acc in accessions:
        matches = [col for col in df_columns if col.startswith(acc)]
        mapped.extend(matches)
    return mapped

clade1_cols = map_to_csv_columns(clade1_genomes, df.columns)
clade2_cols = map_to_csv_columns(clade2_genomes, df.columns)
clade3_cols = map_to_csv_columns(clade3_genomes, df.columns)

# -------- Step 4: Sum abundances per clade --------
df_clades = pd.DataFrame()
df_clades["date"] = df.iloc[:, 0]

df_clades["strain_1"] = df[clade1_cols].sum(axis=1)
df_clades["strain_2"] = df[clade2_cols].sum(axis=1)
df_clades["strain_3"] = df[clade3_cols].sum(axis=1)






# -------- Step 5: Save result --------
df_clades.to_csv("/home/Users/rl152/Strainify/longitudinal/Bacteroides_ovatus/strain_grouping/3_strains_abundances.csv", index=False)

print(df_clades.head())


