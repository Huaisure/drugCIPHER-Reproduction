import pandas as pd
from collections import Counter
import json
import numpy as np
import networkx as nx
from sklearn.model_selection import train_test_split
from utils import cal_test_ts, cal_test_closeness, shortest_path_length_dict

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

# print the length of train and test
print(f"Train length: {len(train)}")
print(f"Test length: {len(test)}")

print("**************************")
print("|     Drug CIPHER-TS     |")
print("**************************")

flatten_train = [j for i in train["dg_atc_codes"] for j in i] + [
    k for i in train["dg_atc_levels"] for j in i for k in j
]
counts_dict = dict(Counter(flatten_train))
ppi_net = nx.read_graphml("ppi_network.graphml")
targets = list(ppi_net.nodes())

drug_ts_test = cal_test_ts(train, test, counts_dict)

ensp_count = json.load(open("target_ENSP_count.json"))
ensp = sorted(ensp_count, key=lambda x: ensp_count[x], reverse=True)
# target = "ENSP00000261707"

# choose the top 10 targets
# 这里选择的target为在数据集中出现次数最多的前10个target，为了方便测试
for i in range(10):
    target = ensp[i]
    target_closeness = cal_test_closeness(train, target, ppi_net)
    score_list = []
    for j in range(len(test)):
        drugTS_list = drug_ts_test[i]
        if np.std(drugTS_list) == 0 or np.std(target_closeness) == 0:
            score = 0
        else:
            score = np.corrcoef(drugTS_list, target_closeness)[0, 1]
            score = abs(score)

        score_list.append((score, j))
    score_list_sorted = sorted(score_list, key=lambda x: x[0], reverse=True)
    # print the top 10 drugs
    print("Target: ", target)
    for i in range(10):
        print("=====================================")
        print(f"Drug: {test['dg_name'].iloc[score_list_sorted[i][1]]}")
        print(f"Score: {score}")

# 如果想要测试其他的target，可以将上面的for循环注释掉，然后取消下面的注释，将target改为想要测试的target
# target = "ENSP00000261707"
# target_closeness = cal_test_closeness(train, target, ppi_net)
# score_list = []
# for j in range(len(test)):
#     drugTS_list = drug_ts_test[j]
#     if np.std(drugTS_list) == 0 or np.std(target_closeness) == 0:
#         score = 0
#     else:
#         score = np.corrcoef(drugTS_list, target_closeness)[0, 1]
#         score = abs(score)
#   score_list.append((score, j))
# score_list_sorted = sorted(score_list, key=lambda x: x[0], reverse=True)
# print("Target: ", target)
# for i in range(10):
#     print("=====================================")
#     print(f"Drug: {test['dg_name'].iloc[score_list_sorted[i][1]]}")
#     print(f"Score: {score}")
# 这里我们纸打印了前10个药物，如果要打印想要的药物，只需对应score_list的索引即可
