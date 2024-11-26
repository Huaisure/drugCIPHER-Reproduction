import os
import pandas as pd
from sklearn.model_selection import train_test_split
import json
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.DataStructs import TanimotoSimilarity
from utils import load_fingerprint, get_relevance_vector_for_drug

# ------------------------------------
# CS

df = pd.read_csv("kegg/data_with_atc_kegg.csv")

train, test = train_test_split(df, test_size=0.2, random_state=42)

# Compute fingerprints for train and test drugs
train_fps = {dg_id: load_fingerprint(dg_id) for dg_id in train["dg_id"]}
test_fps = {dg_id: load_fingerprint(dg_id) for dg_id in test["dg_id"]}

# Calculate Tanimoto similarity for each test drug with all train drugs
similarity_results = []

for test_id, test_fp in test_fps.items():
    similarities = [
        TanimotoSimilarity(test_fp, train_fp)
        for train_id, train_fp in train_fps.items()
    ]
    similarity_results.append(similarities)

# CS_matrix = np.array(similarity_results)
# np.save("CS_matrix.npy", CS_matrix)

# Save the results to a file
output_file = "tanimoto_similarity_results.json"

with open(output_file, "w") as f:
    json.dump(similarity_results, f, indent=4)

print(f"Similarity results saved to {output_file}")

concordance_score_cs_list = []

for i in range(len(test)):
    drugCS_list = similarity_results[i]
    targets_closeness = get_relevance_vector_for_drug(test["dg_id"].iloc[i])
    concordance_score_cs_list.append(
        [
            abs(np.corrcoef(drugCS_list, target_closeness)[0, 1])
            for target_closeness in targets_closeness
        ]
    )

# ------------------------------------
