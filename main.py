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

def cal_code_ts(code1, code2, level1, level2):
    # find the longest common element of the levels
    common_level = set(level1).intersection(set(level2))
    if common_level:
        p_mica = counts_dict[list(common_level)[0]]
    else:
        p_mica = 0
    p_c1 = 0 if code1 not in counts_dict else counts_dict[code1]
    p_c2 = 0 if code2 not in counts_dict else counts_dict[code2]
    if p_c1 == 0 or p_c2 == 0 or p_mica == 0:
        return 0  # If any probability is 0, return similarity as 0 (no shared information)

    # Compute similarity using the provided formula
    similarity = (2 * math.log(p_mica + 1e-6)) / (
        math.log(p_c1 + 1e-6) + math.log(p_c2 + 1e-6)
    )
    return similarity

def cal_ts(code_l1, code_l2, level_l1, level_l2):
    ts = 0
    for idx1, code1 in enumerate(code_l1):
        for idx2, code2 in enumerate(code_l2):
            ts = max(cal_code_ts(code1, code2, level_l1[idx1], level_l2[idx2]), ts)
    return ts

def cal_test_ts(train, test):
    return [
        [cal_ts(train["dg_atc_codes"].iloc[i], test["dg_atc_codes"].iloc[j], train["dg_atc_levels"].iloc[i], test["dg_atc_levels"].iloc[j]) for i in range(len(train))]
    for j in range(len(test))]

drug_ts_test = cal_test_ts(train, test)

print(drug_ts_test)

