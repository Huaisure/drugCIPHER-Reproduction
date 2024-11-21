from drugCIPHER_TS import drug_target_score_TS
import networkx as nx
import pandas as pd
from collections import Counter
import math
from sklearn.model_selection import train_test_split

# df_atc = df[['dg_id', 'dg_name','dg_atc_codes','dg_atc_levels']]
df_atc = pd.read_csv("data_atc.csv")
# for rows with same dg_id, keep only one row
df_atc = df_atc.drop_duplicates(subset="dg_id")
df_atc = df_atc[df_atc["dg_atc_codes"].apply(lambda x: x != "[]")].reset_index(
    drop=True
)
df_atc["dg_atc_codes"] = df_atc["dg_atc_codes"].apply(lambda x: eval(x))
df_atc["dg_atc_levels"] = df_atc["dg_atc_levels"].apply(lambda x: eval(x))

# train test split for df_atc by 8:2
train, test = train_test_split(df_atc, test_size=0.2, random_state=42)

flatten_train = [j for i in train["dg_atc_codes"] for j in i] + [
    k for i in train["dg_atc_levels"] for j in i for k in j
]
counts_dict = dict(Counter(flatten_train))



drug_ts_test = cal_test_ts(train, test)

print(drug_ts_test)

