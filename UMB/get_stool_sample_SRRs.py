import pandas as pd

# Please change the paths to your local setup
df_1 = pd.read_csv("/home/Users/rl152/Strainify/UMB/umb_samples.csv")
df_1 = df_1[df_1['ID']== 'UMB18']

df_2 = pd.read_csv("/dodo/rl152/UMB/PRJNA400628_UMB18/UMB18_mapping.tsv", sep ='\t')
df_3 = pd.concat([df_2[df_2['Run(SRR)'].isin(df_1['Run'].values)], df_2.iloc[-4:,:]], ignore_index = True)


df_3.to_csv("/home/Users/rl152/Strainify/UMB/UMB18_stool_SRRs.csv", index = False)

