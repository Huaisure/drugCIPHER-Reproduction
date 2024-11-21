import networkx as nx
import numpy as np
import pandas as pd
import json
from collections import Counter
import math

def shortest_path_length(ppi_network, ensp1, ensp2):
    """
    Given ppi_network, ensp1, ensp2, return the shortest path length between ensp1 and ensp2 in the network.
    """
    try:
        length = nx.shortest_path_length(
            ppi_network, source=ensp1, target=ensp2, weight="weight"
        )
    except nx.NetworkXNoPath:
        length = -1
    return length


def drug_target_interaction(drug_targets, target):
    """
    drug_targets: list of proteins that are targets of a drug
    target: protein
    return: a score that represents the interaction between the drug targets and the target protein
    """
    ppi = nx.read_graphml("ppi_network.graphml")

    res = 0
    for p in drug_targets:
        length = shortest_path_length(ppi, p, target)
        if length != -1:
            res += np.exp(-(length**2))
    return res


def drugTS(drug_atc_list, drug_level_list):
    """
    return: a score that represents the interaction between the two drugs
    """


    def cal_code_ts(code1, code2, level1, level2):
        # find the longest common element of the levels
        common_level = set(level1).intersection(set(level2))
        if common_level:
            p_mica = counts_dict[list(common_level)[0]]
        else:
            p_mica = 0
        p_c1 = counts_dict[code1]
        p_c2 = counts_dict[code2]
        if p_c1 == 0 or p_c2 == 0 or p_mica == 0:
            return 0  # If any probability is 0, return similarity as 0 (no shared information)

        # Compute similarity using the provided formula
        similarity = (2 * math.log(p_mica + 1e-6)) / (
            math.log(p_c1 + 1e-6) + math.log(p_c2 + 1e-6)
        )
        return similarity

    def cal_ts(code_l1, code_l2, level_l1, level_l2):
        ts = 0
        print(code_l1, code_l2, level_l1, level_l2)
        for idx1, code1 in enumerate(code_l1):
            for idx2, code2 in enumerate(code_l2):
                print(code1, code2, level_l1, level_l2)
                ts = max(cal_code_ts(code1, code2, level_l1[idx1], level_l2[idx2]), ts)

        return ts

    return [
        cal_ts(
            drug_atc_list,
            df_atc.iloc[idx].dg_atc_codes,
            drug_level_list,
            df_atc.iloc[idx].dg_atc_levels,
        )
        for idx in range(len(df_atc))
    ]


def target_to_drugs_interaction(target):
    drug_targets = json.load(open("drug_target_ENSP.json"))
    return [
        drug_target_interaction(drug_targets[drug], target) for drug in drug_targets
    ]
