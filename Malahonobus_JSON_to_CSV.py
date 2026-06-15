# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 20:33:12 2026

@author: smhil
"""

import json
import pandas as pd

with open("mahalanobis_clusters.json", "r") as f:
    data = json.load(f)

rows = []

for obs_id, obs_data in data.items():

    for analysis, analysis_data in obs_data.items():

        for n_clusters, cluster_data in analysis_data.items():

            for variable, variable_data in cluster_data.items():

                for cluster_id, stats in variable_data.items():

                    rows.append({
                        "observation": obs_id,
                        "analysis": analysis,
                        "n_clusters": int(n_clusters),
                        "variable": variable,
                        "cluster": int(cluster_id),
                        "mean": float(stats["mean"]),
                        "stddev": float(stats["standard deviation"])
                    })

df = pd.DataFrame(rows)

df.to_csv("mahalanobis_clusters.csv", index=False)