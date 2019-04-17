%% Ellinwood et al - Human Atrial Myocyte model

% This file loads initial conditions, defines simulation parameters, calls
% the ode solver, and plots simulation results
% 
% Please cite the following papers when using this model:
% 
% Ellinwood N, Dobrev D, Morotti S, Grandi E. (2017).
% Revealing kinetics and state-dependent binding properties of IKur-targeting
% drugs that maximize atrial fibrillation selectivity.
% Chaos 27, 093918. doi: http://dx.doi.org/10.1063/1.5000226
% 
% Ellinwood N, Dobrev D, Morotti S, Grandi E. (2017).
% In silico assessment of efficacy and safety of IKur inhibitors in chronic
% atrial fibrillation: role of kinetics and state-dependence of drug binding.
% Front. Pharmacol.

clear all
close all
clc

%% Load initial conditions
load yf_ham_ikur_df_1Hz % drug-free, steady-state 1-Hz, nSR condition
% load yf_ham_ikur_df_3Hz % drug-free, steady-state 3-Hz, nSR condition
% load yf_ham_ikur_df_1Hz_cAF % drug-free, steady-state 1-Hz, cAF condition
% load yf_ham_ikur_df_3Hz_cAF % drug-free, steady-state 3-Hz, cAF condition

y0 = yfinal;

%% Choice of nSR or AF condition
AF = 0; % (boolean - 0 for nSR conditions, 1 for cAF conditions)

%% Experimental protocol parameters
% Pacing frequency
prot_rate = 1; % PACING RATE (Hz)

% Isoproterenol administration
exp_ISO = 0; % (boolean - 0 for no ISO, 1 for ISO)

% Acetylcholine concentration
exp_Ach = 0; % (uM) DEFINE [ACh] HERE (0-0.1-1 uM)

% Other experimental conditions
exp_Temp = 310; % temperature (300 or 310 K)
exp_Nao = 140; % extracellular [Na] (mM)

%% I_Kur drug-binding parameters
% Set ikur_drug flag to 1 to simulate administration of IKur blockers
ikur_drug_flag = 0; % (boolean - 0 for no drug, 1 for drug)
% This parameter can be used to scale the conductance of IKur (pore block)
% in drug-free or drug-binding conditions. Should have a value between 0 and 1
fraction_pore_block = 0; % 0 = no pore block, 0.5 = 50% pore block, 1 = complete block of IKur
if ikur_drug_flag == 1,
    % Define drug-binding model
    drug_index_ikur = 1; % (1-6)
    % INDEX:	State-specific drug binding description:
    % (1)       To the open state only
    % (2)       To the inactivated state only
    % (3)       To the closed states only
    % (4)       To the open and inactivated states with foot-in-the-door mechanism
    % (5)       To the open and inactivated states with transition between drug-bound layers
    % (6)       To the open and inactivated states with variable affinity
    % Use this parameter to set the concentration of drug available for drug binding
    drug_conc_kur = 1; % uM
    % On rate of drug binding
    kon = 0.1; % 1/(uM*ms) - given the diffusion rate is in 1/(uM*ms)
    % Offrate of drug-binding
    koff = 0.1;  % 1/ms
    if drug_index_ikur == 6, % Ko = kon_O/koff_O, and Ki = kon_I/koff_I
        drug_conc_kur = 1; % uM
        % On rate of drug binding to the open state of Kv1.5
        kon_O = 0.1; % 1/(uM*ms)
        % Off rate of drug binding to the open state of Kv1.5
        koff_O = 0.1; % 1/ms
        % On rate of drug binding to the open state of Kv1.5
        kon_I = 0.1; % 1/(uM*ms)
        % Off rate of drug binding to the open state of Kv1.5
        koff_I = 0.1; % 1/ms
    else 
        kon_O = 0; koff_O = 0; kon_I = 0; koff_I = 0;
    end
else 
    % set drug-binding information to zero for drug-free conditions
    drug_index_ikur = 0; drug_conc_kur = 0; kon_diffusion_rate = 0;
    kon = 0; kon_O = 0; kon_I = 0;
    koff = 0; koff_O = 0; koff_I = 0;
end

%% Ranolazine parameters
% Set ran_flag to 1 to simulate administration of ranolazine
ran_flag = 0; % (boolean - 0 for no drug, 1 for drug)
if ran_flag == 1
    drug_index = 1; drug_conc = 10 * (1E-6); % (M) DEFINE [RAN] HERE
else
    drug_index = 0; drug_conc = 0; % Drug Free
end

%% Parameter array for passing nondefault conditions
prot_par = [AF prot_rate];                                  % 1 2
drug_par = [drug_index drug_conc];                          % 3 4
exp_par = [exp_Temp exp_Nao exp_ISO exp_Ach];               % 5 6 7 8
ikur_drug_par = [drug_index_ikur drug_conc_kur kon koff kon_O koff_O kon_I koff_I fraction_pore_block]; % 9 10 11 12 13 14 15 16 17
p = [prot_par drug_par exp_par ikur_drug_par];

%% Single Run Simulation
% Simulation duration
duration = 3e3; % (ms) DEFINE SIMULATION DURATION HERE
tspan = [0; duration];

tic
options = odeset('RelTol',1e-6,'MaxStep',1,'Stats','on'); 
[t,y] = ode15s(@Ellinwood_et_al_ham_ikur_drug_model,tspan,y0,options,p);
toc

%% Final conditions
% yfinal contains the final values of all state variables
yfinal = y(end,:);

% save yf_ham_ikur_df_1Hz yfinal % drug-free, steady-state 1-Hz, nSR condition
% save yf_ham_ikur_df_3Hz yfinal % drug-free, steady-state 3-Hz, nSR condition
% save yf_ham_ikur_df_1Hz_cAF yfinal % drug-free, steady-state 1-Hz, cAF condition
% save yf_ham_ikur_df_3Hz_cAF yfinal % drug-free, steady-state 3-Hz, cAF condition

%% Extract non-state variables outputs
currents = Ellinwood_et_al_ham_ikur_calcCurrents(t,y,p);

%% Plot simulation outputs
figure, set(gcf,'color','w')
subplot(3,1,1); hold on, plot(t*1e-3,y(:,39),'-')
set(gca,'box','off','tickdir','out','fontsize',16)
ylabel('E_m (mV)')
subplot(3,1,2); hold on, plot(t*1e-3,currents(:,7),'-')
set(gca,'box','off','tickdir','out','fontsize',16)
ylabel('I_K_u_r (A/F)')
subplot(3,1,3); hold on, plot(t*1e-3,y(:,38)*1e3,'-')
set(gca,'box','off','tickdir','out','fontsize',16)
ylabel('[Ca^2^+]_i (\muM)')
xlabel('Time (s)');