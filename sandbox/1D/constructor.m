% Noah Germolus 26 Jan 2021
% This will define test vars for the ecosystem model.

po = struct; % optimizable params 
% (* marks those that were optimized for EQP.)
% (** marks those optimized for HOT)
ps = struct; % fixed params
y = struct;  % state vars

base = ones(100,1); %Something to make into 100-element vectors

% Tfun
po.AE = 2000;

% PHY
po.mu_PHY = 3.3763; %*
po.alpha_PHY_chl = 1.9044;%*
po.beta_PHY = 0.001; %*
po.vref_PHY_n = 0.3393; %*
po.k_PHY_nh4 = 0.19000; %*
po.k_PHY_no3 = 4.0;
po.vref_PHY_p = 0.0140; %*
po.k_PHY_po4 = 0.1;
po.zeta_no3 = 2.0; % Universal
po.theta =  3.4; % Universal  %*
po.ex_PHY_psv = 0.05;
po.ex_PHY_cho = 0.05;
po.pom_PHY = 0.13;  %*

% TR

po.mu_TR = 3.0;
po.alpha_TR_chl = 0.25;
po.beta_TR = 0.002;
po.vref_TR_n = 0.4;
po.k_TR_nh4 = 0.05;
po.k_TR_no3 = 0.5;
po.vref_TR_p = 0.02;
po.k_TR_po4 = 0.047;
po.zeta_nfix = 2.5; % Universal
po.picktrpo4 = 1.0; % Not the uncapitalized parameter. There's a var in the function with "TR"
po.ex_TR_psv = 0.02;
po.ex_TR_nfix = 0.36;
po.ex_TR_cho = 0.05;
po.pom_TR = 0.1;


% % UN
% 
% po.mu_UN = 1.6; %For UN, only this and alpha are unique
% po.alpha_UN_chl = 0.5;
% po.beta_UN = po.beta_PHY;
% po.vref_UN_n = po.vref_PHY_n;
% po.k_UN_nh4 = po.k_PHY_nh4;
% po.k_UN_no3 = po.k_PHY_no3;
% po.vref_UN_p = po.vref_PHY_p;
% po.k_UN_po4 = po.k_PHY_po4;
% po.ex_UN_psv = po.ex_PHY_psv;
% po.ex_UN_nfix = po.ex_TR_nfix;
% po.ex_UN_cho = po.ex_PHY_cho;
% po.pom_UN = po.pom_PHY;

% BA

po.k_DOM = 0.5500; %*
po.r_SDOM = 0.087; %* 0.0087 Adjusted NPG 15 Mar 2021 to get less SDOMc.
po.mu_BA = 2.3; %*
po.BAresp_B = 0.15;
po.r_BAresp_min = 0.35;
po.r_BAresp_max = 0.8665; %*
po.b_BAresp = 0.01;
po.r_BArefr = 0.02;
po.r_BAremi = 6.0;
po.r_BAadj = 2.0;
po.r_BAmort = 0.02;
po.f_BAslct = 0.25;

% PRT

po.mu_PRT = 1.1477; %*
po.g_PHY = 0.9267; %*
po.g_BA = 1.4; %*
po.ex_PRT = 0.20;
po.f_ex_PRT = 0.75;
po.resp_B_PRT = 0.01;
po.resp_A_PRT = 0.6052; %*
po.ex_adj_PRT = 2.0;
po.remi_PRT = 4.0;
po.pom_PRT = 0.0260; %*

% MZ & HZ

po.mu_MZ = 1.4; %*
po.g_PRT = 2.0831; %*
% po.g_UN = po.g_PHY;
po.g_TR = 0.2007; %**
po.ex_MZ = 0.2;
po.f_ex_MZ = 0.75;
po.resp_B_MZ = 0.03;
po.resp_A_MZ = 0.3683; %*
po.ex_adj_MZ = 2.0;
po.remi_MZ = 4.0;
po.pom_MZ = 0.1;
po.ex_refr_MZ = 0.02;
po.remv_MZ = 0.2; %*
po.f_POM_HZ = 0.35;
po.f_SDOM_HZ = 0.2;

% REFR and POM and inorg nuts

po.ex_REFR_SDOM = 0.0009;
po.q_REFR_n = 0.05;
po.q_REFR_p = 0.0007;
po.q_POM_n = 0.12;
po.q_POM_p = 0.01;
po.diss = 0.1593; %*
po.r_ntrf = 0.1700; %*
po.prf_N = 1.05; %*
po.prf_P = 3.8; %*

% Setup for ps

ps.q_PHY_n = 0.034; 
ps.q_PHY_n_rdf = 0.15;
ps.q_PHY_n_max = 0.17;

ps.q_PHY_p = 0.00188; 
ps.q_PHY_p_rdf = 0.0094;
ps.q_PHY_p_max = 0.0169;

ps.q_TR_n = 0.12; 
ps.q_TR_n_rdf = 0.16;
ps.q_TR_n_max = 0.2;

ps.q_TR_p = 0.001; 
ps.q_TR_p_rdf = 0.0035;
ps.q_TR_p_max = 0.006;

% ps.q_UN_n = 0.12; 
% ps.q_UN_n_rdf = 0.16;
% ps.q_UN_n_max = 0.2;
% 
% ps.q_UN_p = 0.001; 
% ps.q_UN_p_rdf = 0.0035;
% ps.q_UN_p_max = 0.006;

ps.q_BA_n = 0.18;
ps.q_BA_p = 0.02;

ps.q_PRT_n = 0.2;
ps.q_PRT_p = 0.022;

ps.q_MZ_n = 0.2;
ps.q_MZ_p = 0.008;


%y.T = 20; % in deg C
%y.PAR = 0; % in W m^-2

y.iPHYc = 0.876 * 1.056* base; %*
y.iPHYn = 0.477* base; %*
y.iPHYp = 0.00863 * 1.1838* base; %*
y.iPHYchl = 0.2* base; % Guess in mg/m3 based on even split

y.iTRc = 0.312 * 0.269* base; %**
y.iTRn = 0.0494 * 0.309* base;%**
y.iTRp = 0.00236 * 0.1725* base;%**
y.iTRchl = 0.051* base;

% y.iUNc = 0.312 * 0.626;%**
% y.iUNn = 0.0494 * 0.716;%**
% y.iUNp = 0.00236 * 0.4219;%**
% y.iUNchl = 0.1;

y.iBAc = 0.741* base; %*
y.iBAn = 0.122 * 0.376* base; %*
y.iBAp = 0.00863 * 3.3373* base; %*

y.iPRTc = 0.876 * 1.067* base; %*
y.iPRTn = 0.122 * 1.643* base; %*
y.iPRTp = 0.00863 * 2.5191* base; %*

y.iMZc = 0.793* base; %*
y.iMZn = 0.122 * 0.342* base; %*
y.iMZp = 0.00863 * 0.2163* base; %*

y.iDETc = 1* base; %* (Altered NPG)
y.iDETn = 0.486* base; %*
y.iDETp = 0* base; %*

y.iNO3 = 3.41* base; %*
y.iNH4 = y.iNO3(end)/10.* base; % all mmol m^-3. Guessses from EQP site avg
y.iPO4 = 0.559* base; %*

y.iLDOMc = 0.005* base; %* (Altered NPG)
y.iLDOMn = 0.122 * 0.220* base; %*
y.iLDOMp = 0.00863 * 0.3451* base; %*

y.iSDOMc = 9* base; %*
y.iSDOMn = 2.90* base; %*
y.iSDOMp = 0* base; %*

% y.fluxBApo4 = 0; %
% y.fluxBAno3 = 0; %
% y.fluxBAnh4 = 0; %
y.pp = 0* base; %
y.bp = 0* base; %
z = [0.5:1:99.5]';