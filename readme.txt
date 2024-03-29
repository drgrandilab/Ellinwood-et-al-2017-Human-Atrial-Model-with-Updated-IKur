Matlab code for Ellinwood et al. human atrial model with updated IKur formulation.

Ellinwood et al. integrated a Markov model of the IKur current derived from Zhou et al. (PLoS ONE
2012; e42295), modified to fit voltage-clamp data from human atrial myocytes, into the Morotti et al.
human atrial model (J Mol Cell Cardiol. 2016 Jul; 96:63-71, available for download on this webpage).
The IKur current model was extended to reproduce dose-dependent drug binding to various conformational
states of the channel with variable affinity and binding kinetics.


_____________________________________________________________________________________________________
Contents:

readme.txt						this file

Ellinwood_et_al_ham_ikur_drug_model_main.m		loads initial conditions and runs the simulation
Ellinwood_et_al_ham_ikur_drug_model.m			excitation-contraction coupling model
Ellinwood_et_al_ham_ikur_calcCurrents.m			supporting function for simulation output analysis

.mat files			initial conditions (steady-state, drug-free)
- yf_ham_ikur_df_1Hz		nSR condition, 1-Hz pacing
- yf_ham_ikur_df_1Hz		nSR condition, 3-Hz pacing
- yf_ham_ikur_df_1Hz_cAF	cAF condition, 1-Hz pacing
- yf_ham_ikur_df_3Hz_cAF	cAF condition, 3-Hz pacing

_____________________________________________________________________________________________________

References:

Ellinwood N, Dobrev D, Morotti S, Grandi E. (2017).
Revealing kinetics and state-dependent binding properties of IKur-targeting drugs that maximize atrial
fibrillation selectivity.
Chaos 27, 093918 (2017). doi: https://doi.org/10.1063/1.5000226
[See Erratum, Chaos 27, 109902 (2017). doi: https://doi.org/10.1063/1.5007051]

Ellinwood N, Dobrev D, Morotti S, Grandi E. (2017).
In silico assessment of efficacy and safety of IKur inhibitors in chronic atrial fibrillation: role of 
kinetics and state-dependence of drug binding
Front. Pharmacol. 2017; 8: 799. doi: https://doi.org/10.3389/fphar.2017.00799

Morotti S, McCulloch AD, Bers DM, Edwards AG, Grandi E. (2015).
Atrial-selective targeting of arrhythmogenic phase-3 earlyafterdepolarizations in human myocytes.
J Mol Cell Cardiol. 2016 Jul;96:63-71. doi: https://doi.org/10.1016/j.yjmcc.2015.07.030

Grandi E, Pandit SV, Voigt N, Workman AJ, Dobrev D, Jalife J, Bers DM. (2011).
Human Atrial Action Potential and Ca2+ Model: Sinus Rhythm and Chronic Atrial Fibrillation.
Circ Res. 2011 Oct 14;109(9):1055-66. doi: https://doi.org/10.1161/CIRCRESAHA.111.253955

Please cite the above papers when using this model.
