#!/usr/bin/env python
# coding: utf-8
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib

rates = pd.read_csv(r"C:\Users\Maziya\Documents\49vs49_GXcommspecific.csv", index_col=0)
rates = rates.replace(0, np.nan)

# differences between monoculture and cocultures

rates_diff_test = []

for i in range(len(rates)):
    diff = [(x - rates.iloc[i, i]) for x in rates.iloc[i , :] if x != np.nan]
    rates_diff_test.append(diff)
    
rates_diff_test = pd.DataFrame(rates_diff_test, index=rates.columns, columns=rates.columns)

# differences between monoculture and cocultures with cutoff

rates_diff = []

for i in range(len(rates)):
    cutoff_pos = rates.iloc[i, i] + (rates.iloc[i, i] * 0.1)
    cutoff_neg = rates.iloc[i, i] - (rates.iloc[i, i] * 0.1)
    diff = [1 if x > cutoff_pos # +1 if difference greater than +10%
            else -1 if x < cutoff_neg # -1 if difference less than -10%
            else 0 if x > cutoff_neg and x < cutoff_pos # else 0
            else None for x in rates.iloc[i , :]] # null values
    rates_diff.append(diff)
    
rates_diff = pd.DataFrame(rates_diff, index=rates.columns, columns=rates.columns)

for i in range(len(rates_diff)):
    rates_diff.iloc[i, i] = None


rates_diff.head()


"""
change first value of each list to change color of "Decrease"
change second value of each list to change color of "No Significant change"
change third value of each list to change color of "Increase"
"""
fig, ax = plt.subplots(figsize=(10,8))

reds = [65, 200, 225]
greens = [105, 200, 99]
blues = [255, 200, 71]

myColors = []
for i, j, k, in zip(reds, greens, blues):
    colors = [k/255, j/255, i/255, 1]
    myColors.append(colors)
cmap = matplotlib.colors.LinearSegmentedColormap.from_list('Custom', myColors, len(myColors))

ax = sns.heatmap(data=rates_diff, cmap=cmap, linecolor='black', linewidth=0.01, xticklabels=True, yticklabels=True)
plt.yticks(fontstyle='italic', fontsize=8)
plt.xticks(fontstyle='italic', fontsize=8)
plt.xlabel("(in co-culture with)",fontsize=16)
plt.ylabel("Organism",fontsize=16)

colorbar = ax.collections[0].colorbar
colorbar.set_ticks([-0.667, 0.000, 0.667])
colorbar.set_ticklabels(['Decrease','No Significant change','Increase'])
ax.tick_params(labelsize=8)
plt.savefig("heatmap_monovsco.png", dpi=900, bbox_inches='tight')