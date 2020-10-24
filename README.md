# **CAMP** 
# **C**o-culture/**C**ommunity **A**nalyses for **M**etabolite **P**roduction
#### CAMP is a computational analyses framework which allows generation and evaluation of numerous microbial communities for optimal production of a metabolite of interest. 

#### Overview
1) Microbial GSMMs are retrieved from databases such as [VMH](https://www.vmh.life/) (Virtual Metabolic Human). Each of these GSMMs is simulated in three different nutrient conditions (Excess, Minimal and Species-specific). Predicted growth rates and product flux are obtained using flux balance analysis (FBA) and flux variability analysis (FVA). The product yield is computed as the maximum product flux obtained per unit flux of substrate uptake. 
2) Two-species communities are created using [SteadyCom](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005539). Community models are again simulated in three nutrient conditions (Excess, Minimal, Community Specific). FBA and FVA are used to predict community growth rates and product yield in the community. Monoculture and co-culture growth rates are compared to identify an increase or decrease in growth when an organism is simulated in the presence of another. 
3) Expected product yield in a community is compared to the observed product yield. Communities which have a 10-fold increase in product yield are regarded as candidate communities for optimal production of the target metabolite. Communities are assessed for their relative abundances, type of interaction behaviour observed and the cross-fed metabolites.
4) In silico community optimisation is performed using [FSEOF](https://aem.asm.org/content/76/10/3097.long), which enables to shortlist potential reaction knock-outs that will increase product flux in the community. Reaction knock-outs can be from either species in the community. 

![CAMP](CAMP.png)

#### Prerequisites
All simulations were performed in MATLAB R2018a (MathWorks Inc., USA) using: 
1. [COBRA Toolbox v3.0](https://opencobra.github.io/cobratoolbox/stable/)
2. [IBM ILOG CPLEX](https://www.ibm.com/in-en/products/ilog-cplex-optimization-studio) 


#### Acknowledgements
[Initiative for Biological Systems Engineering](https://ibse.iitm.ac.in/)


