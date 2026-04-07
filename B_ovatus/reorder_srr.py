import pandas as pd

# Load SRR list with desired column order
with open('/home/Users/rl152/PycharmProjects/longitudinal/srr_list.txt') as f:
    srr_order = [line.strip() for line in f]

# Load the abundance CSV
df = pd.read_csv('/home/Users/rl152/PycharmProjects/longitudinal/Bacteroides_ovatus/results_all/abundance_estimates_combined.csv')

# Ensure the first column (strain names) is kept
first_col = df.columns[0]

# Rearranged column list: keep the first column, then reorder the SRR columns
new_columns = [first_col] + srr_order

# Reorder the dataframe
df_reordered = df[new_columns]

# Save to a new CSV
df_reordered.to_csv('abundance_reordered.csv', index=False)
