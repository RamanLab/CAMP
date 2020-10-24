#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib


# In[2]:


rates = pd.read_csv(r"C:\Users\Maziya\Documents\49vs49_GXcommspecific.csv", index_col=0)
rates = rates.replace(0, np.nan)


# In[3]:



fig, ax = plt.subplots(figsize=(12,10)) 
#sns.heatmap(rates,fmt=".2f",annot=True, cmap="YlGnBu",ax=ax,linewidth=0.01, xticklabels=True,yticklabels=True);
sns.heatmap(rates,cmap="YlGnBu",ax=ax,linewidth=0.01, xticklabels=True,yticklabels=True);
plt.yticks(fontstyle='italic', fontsize=8)
plt.xticks(fontstyle='italic', fontsize=8)
plt.xlabel("(in co-culture with)",fontsize=16)
plt.ylabel("Organism",fontsize=16)
colorbar = ax.collections[0].colorbar
params = {'mathtext.default': 'regular' }          
plt.rcParams.update(params)
colorbar.set_label('Growth Rates$(h^{-1})$')                                                                                      
plt.savefig("heatmap_GXcommspecific.png", dpi=900, bbox_inches='tight')


# In[ ]:




