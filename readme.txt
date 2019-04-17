Matlab code for Ellinwood et al. human atrial model.

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
Chaos 27, 093918. doi: http://dx.doi.org/10.1063/1.5000226

Ellinwood N, Dobrev D, Morotti S, Grandi E. (2017).
In silico assessment of efficacy and safety of IKur inhibitors in chronic atrial fibrillation: role of 
kinetics and state-dependence of drug binding
Front. Pharmacol.

Morotti S, McCulloch AD, Bers DM, Edwards AG, Grandi E. (2015).
Atrial-selective targeting of arrhythmogenic phase-3 earlyafterdepolarizations in human myocytes.
J Mol Cell Cardiol. 2016 Jul;96:63-71. doi: 10.1016/j.yjmcc.2015.07.030

Grandi E, Pandit SV, Voigt N, Workman AJ, Dobrev D, Jalife J, Bers DM. (2011).
Human Atrial Action Potential and Ca2+ Model: Sinus Rhythm and Chronic Atrial Fibrillation.
Circ Res. 2011 Oct 14;109(9):1055-66. Epub 2011 Sep 15.

Please cite the above papers when using this model.