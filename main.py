from drugCIPHER_TS import drug_target_score_TS
import networkx as nx

drug = "DB00014"
target = "ENSP00000323308"
dts = drug_target_score_TS(drug, target)
print("score:", dts.score())
