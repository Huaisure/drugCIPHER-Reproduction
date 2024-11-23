import networkx as nx
import numpy as np
import pandas as pd
import json
import os
import tqdm
from collections import Counter
import math

if os.path.exists("shortest_path_length_dict.json"):
    shortest_path_length_dict = json.load(open("shortest_path_length_dict.json"))
else:
    shortest_path_length_dict = {}

if os.path.exists("data_target_ENSP.json"):
    drug_targets = json.load(open("data_target_ENSP.json"))
else:
    raise FileNotFoundError("data_target_ENSP.json not found")


def shortest_path_length(ppi_network, ensp1, ensp2):
    """
    Given ppi_network, ensp1, ensp2, return the shortest path length between ensp1 and ensp2 in the network.
    """
    try:
        check1 = "(" + ensp1 + "," + ensp2 + ")"
        check2 = "(" + ensp2 + "," + ensp1 + ")"
        if check1 in shortest_path_length_dict:
            length = shortest_path_length_dict[check1]
        elif check2 in shortest_path_length_dict:
            length = shortest_path_length_dict[check2]
        else:
            length = nx.shortest_path_length(
                ppi_network, source=ensp1, target=ensp2, weight="weight"
            )
            shortest_path_length_dict[check1] = length
    except nx.NetworkXNoPath:
        length = -1
    except nx.NodeNotFound:
        # If the node is not in the network, return -1
        length = -1
    return length


def drug_target_interaction(ppi, drug_targets, target):
    """
    drug_targets: list of proteins that are targets of a drug
    target: protein
    return: a score that represents the interaction between the drug targets and the target protein
    """

    res = 0
    for p in drug_targets:
        length = shortest_path_length(ppi, p, target)
        if length != -1:
            length = length / 1000
            res += np.exp(-(length**2))
    return res


def cal_code_ts(counts_dict, code1, code2, level1, level2):
    # find the longest common element of the levels
    common_level = set(level1).intersection(set(level2))
    if common_level:
        p_mica = counts_dict[list(common_level)[0]]
    else:
        p_mica = 0
    p_c1 = 0 if code1 not in counts_dict else counts_dict[code1]
    p_c2 = 0 if code2 not in counts_dict else counts_dict[code2]
    if p_c1 == 0 or p_c2 == 0 or p_mica == 0:
        return (
            0  # If any probability is 0, return similarity as 0 (no shared information)
        )

    # Compute similarity using the provided formula
    similarity = (2 * math.log(p_mica + 1e-6)) / (
        math.log(p_c1 + 1e-6) + math.log(p_c2 + 1e-6)
    )
    return similarity


def cal_ts(counts_dict, code_l1, code_l2, level_l1, level_l2):
    ts = 0
    for idx1, code1 in enumerate(code_l1):
        for idx2, code2 in enumerate(code_l2):
            ts = max(
                cal_code_ts(counts_dict, code1, code2, level_l1[idx1], level_l2[idx2]),
                ts,
            )
    return ts


def cal_test_ts(train, test, counts_dict):
    return [
        [
            cal_ts(
                counts_dict,
                train["dg_atc_codes"].iloc[i],
                test["dg_atc_codes"].iloc[j],
                train["dg_atc_levels"].iloc[i],
                test["dg_atc_levels"].iloc[j],
            )
            for i in range(len(train))
        ]
        for j in range(len(test))
    ]


def target_to_drugs_interaction(ppi, drug, target):
    targets = None
    for d in drug_targets:
        if d["dg_id"] == drug:
            targets = d["target_ENSP"]
            break
    if targets is None:
        raise ValueError(f"Drug {drug} not found in drug_targets")
    return drug_target_interaction(ppi, targets, target)


def cal_test_closeness(train, target, ppi):
    print(f"Calculating closeness for target {target}")
    res = [
        target_to_drugs_interaction(ppi, train["dg_id"].iloc[i], target)
        for i in tqdm.tqdm(range(len(train)))
    ]
    with open("shortest_path_length_dict.json", "w") as f:
        json.dump(shortest_path_length_dict, f, indent=2)
    return res

# Function to load a molecule from its .mol file and compute fingerprints
def load_fingerprint(dg_id, mol_dir="./kegg/mol"):
    mol_path = os.path.join(mol_dir, f"{dg_id}.mol")
    mol = Chem.MolFromMolFile(mol_path)
    if mol is None:
        raise ValueError(f"Error loading {mol_path}")
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
