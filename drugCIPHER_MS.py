import os
import pandas as pd
from sklearn.model_selection import train_test_split
import json
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.DataStructs import TanimotoSimilarity
from utils import load_fingerprint, cal_test_ts, cal_test_closeness, shortest_path_length_dict, get_relevance_vector_for_drug
from collections import Counter
import numpy as np
import networkx as nx

df = pd.read_csv("kegg/data_with_atc_kegg.csv")

train, test = train_test_split(df, test_size=0.2, random_state=42)

flatten_train = [j for i in train["dg_atc_codes"] for j in i] + [
    k for i in train["dg_atc_levels"] for j in i for k in j
]
counts_dict = dict(Counter(flatten_train))
ppi_net = nx.read_graphml("ppi_network.graphml")
targets = list(ppi_net.nodes())

drug_ts_test = cal_test_ts(train, test, counts_dict)


concordance_score_ts_list = []
concordance_score_cs_list = []

concordance_score_ms_list = []

for j in range(len(test)):
    drugTS_list = drug_ts_test[j]
    test_dg_id = test["dg_id"].iloc[j]
    targets_closeness = get_relevance_vector_for_drug(test_dg_id)
    ts_score_list = [abs(np.corrcoef(drugTS_list, target_closeness)[0, 1]) for target_closeness in targets_closeness]
    concordance_score_ts_list.append(list(ts_score_list))
    
# Compute fingerprints for train and test drugs
train_fps = {dg_id: load_fingerprint(dg_id) for dg_id in train['dg_id']}
test_fps = {dg_id: load_fingerprint(dg_id) for dg_id in test['dg_id']}

# Calculate Tanimoto similarity for each test drug with all train drugs
drug_ts_test = []

for test_id, test_fp in test_fps.items():
    similarities = [
        TanimotoSimilarity(test_fp, train_fp)
        for train_id, train_fp in train_fps.items()
    ]
    drug_ts_test.append(similarities)


for i in range(len(test)):
    test_ts = drug_ts_test[i]
    test_cs = drug_ts_test[i]
    test_dg_id = test["dg_id"].iloc[i]
    targets_closeness = get_relevance_vector_for_drug(test_dg_id)
    cs_score_list = []
    ms_score_list = []
    for j, target_closeness in enumerate(targets_closeness):
        rho_T_pd = concordance_score_ts_list[i][j]
        concordance_score_cs = abs(np.corrcoef(test_cs, target_closeness)[0, 1])
        cs_score_list.append(concordance_score_cs)
        
        X = np.column_stack((np.array(test_ts), np.array(test_cs), np.ones_like(test_ts)))
        beta = np.linalg.pinv(X.T @ X) @ X.T @ np.array(target_closeness)
        
        a_pd, b_pd, c = beta
        sigma_TS_d =  np.std(test_ts) # Standard deviation of TS_d
        sigma_CS_d = np.std(test_cs)  # Standard deviation of CS_d
        rho_C_pd = concordance_score_cs   # Correlation coefficient for C
        numerator = (sigma_TS_d * rho_C_pd / abs(b_pd)) + (sigma_CS_d * rho_T_pd / abs(a_pd))
        denominator = np.sqrt((sigma_TS_d ** 2 / b_pd ** 2) + (sigma_CS_d ** 2 / a_pd ** 2))

        concordance_score_ms = numerator / denominator
        ms_score_list.append(concordance_score_ms)
    concordance_score_cs_list.append(list(cs_score_list))
    concordance_score_ms_list.append(list(ms_score_list))
    





        