import pandas as pd
import numpy as np
from utils import drugTS, target_to_drugs_interaction


class drug_target_score_TS:
    """
    output a score that represents the interaction between the drug targets and the target protein
    """

    data_atc = pd.read_csv("data_atc.csv")

    def __init__(self, drug, target):
        self.drug = drug
        self.target = target
        self._convert_drugid()

    def _convert_drugid(self):
        if self.drug in self.data_atc.dg_id.values:
            self.drug_atc = self.data_atc[
                self.data_atc.dg_id == self.drug
            ].dg_atc_codes.values[0]
            self.drug_level = self.data_atc[
                self.data_atc.dg_id == self.drug
            ].dg_atc_levels.values[0]

    def score(self):
        drugTS_list = drugTS(self.drug_atc, self.drug_level)
        target_drug_list = target_to_drugs_interaction(self.target)
        cov = np.cov(drugTS_list, target_drug_list)
        std_drug = np.std(drugTS_list)
        std_target = np.std(target_drug_list)
        score = cov / (std_drug * std_target)
        return score
