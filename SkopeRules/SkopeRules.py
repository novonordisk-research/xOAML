# Code to run SkopeRules
# Nielsen et al, Nat Comms 2024

# The code below was run in a jupyter notebook, and the final output saved as txt format

################# Load libraries ###############

import pandas as pd
import numpy as np
import six
import sys
sys.modules['sklearn.externals.six'] = six
from skrules import SkopeRules

################# Load data ###############

X = pd.read_table("skoperule_input/skrules_raw_values.tsv")
print(X.head(5))

s_labels = pd.read_table("skoperule_input/Clin_model_cluster_ID.tsv")
print(s_labels.head(5))

s_labels = s_labels["cluster_ID"].values

y = pd.read_table("skoperule_input/skrules_obs.tsv")
print(y.head(5))

y = y["obs"].values

################# Run SkopeRules ###############

# from https://github.com/AidanCooper/shap-clustering/blob/main/analysis.ipynb

def skoperules_models(
    X: pd.DataFrame,
    y: np.array,
    labels: np.array,
    n_estimators: int = 10,
    recall_min: float = 0.2,
    max_depth: int = 5,
    max_depth_duplication: int = 7,
    max_samples: float = 1.0,
    random_state: int = 1234,
):
    X["cluster"] = labels
    sr_models = {}
    for cluster in X["cluster"].value_counts().index.sort_values():
        print(cluster)
        Xc = X.drop("cluster", axis=1)
        yc = (X["cluster"] == cluster) * 1
        sr = SkopeRules(
            feature_names=Xc.columns,
            random_state=random_state,
            n_estimators=n_estimators,
            recall_min=recall_min,
            max_depth=max_depth,
            max_depth_duplication=max_depth_duplication,
            max_samples=max_samples,
            max_features=None
        )
        sr.fit(Xc, yc)
        sr_models[cluster] = sr
    return sr_models


sv_models = skoperules_models(X, y, s_labels)

for cluster in sv_models.keys():
    print("Cluster {cluster}")
    if len(sv_models[cluster].rules_) == 0:
        print("N/A")
    else:
        for i, rule in enumerate(sv_models[cluster].rules_):
            if i == 0:
                for line in rule[0].split(" and "):
                    print(f"\t{line}")
                print(f"Precision: {rule[1][0]:.4f}")
                print(f"Recall   : {rule[1][1]:.4f}")
    print()
    
    
################# Rule output ###############

for cluster in sv_models.keys():
    print(cluster)
    print(sv_models[cluster].rules_)
    
    

