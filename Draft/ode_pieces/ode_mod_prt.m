%% Noah Germolus, 13 Jan 2021
% This is the system of differential equations representing protozoan
% grazers in Luo et al. (2010)


function [dydtt] = ode_mod_prt(y)
%-----------------------------------------------------------------------
%      Protozoan Grazers
%-----------------------------------------------------------------------

% Gross growth rate

growPRTc = graz_PHY_c + graz_UN_c + graz_BA_c;
growPRTn = graz_PHY_n + graz_UN_n + graz_BA_n;
growPRTp = graz_PHY_p + graz_UN_p + graz_BA_p;

% DOM excretion

excr_PRT_LDOC = ex_PRT * f_ex_PRT * growPRTc;
excr_PRT_LDON = ex_PRT * f_ex_PRT * growPRTn;
excr_PRT_LDOP = ex_PRT * f_ex_PRT * growPRTp;

excr_PRT_SDOC = ex_PRT * (1 - f_ex_PRT) * growPRTc;
excr_PRT_SDON = ex_PRT * (1 - f_ex_PRT) * growPRTn * (y.iPRTn/(y.iPRTc*q_PRT_n));
excr_PRT_SDOP = ex_PRT * (1 - f_ex_PRT) * growPRTp * (y.iPRTp/(y.iPRTc*q_PRT_p));

% Respiration

resp_PRT = (resp_B_PRT * tfun(y) * y.iPRTc) + (resp_A_PRT * growPRTc)

% SDOM release to adjust stoichiometry

excr2_PRT_SDOC = ex_adj_PRT * y.iPRTc * max(0,...
    1 - (y.iPRTn/(y.iPRTc*q_PRT_n)),...
    1 - (y.iPRTp/(y.iPRTc*q_PRT_p)));
excr2_PRT_SDON = 0.5 * excr2_PRT_SDOC * (y.iPRTn/y.iPRTc);
excr2_PRT_SDOP = 0.5 * excr2_PRT_SDOC * (y.iPRTp/y.iPRTc);

% Inorganic nutrient remin to adjust stoichiometry.

remi_PRT_n = remi_prt * max(0, y.iPRTn - (q_PRT_n * y.iPRTc),...
    y.iPRTn - (q_PRT_n * y.iPRTp / q_PRT_p));
remi_PRT_p = remi_prt * max(0, y.iPRTp - (q_PRT_p * y.iPRTc),...
    y.iPRTp - (q_PRT_p * y.iPRTn / q_PRT_n));

% POM production

POM_PRT_c = pom_PRT * growPRTc;
POM_PRT_n = q_POM_n * POM_PRT_c;
POM_PRT_p = q_POM_p * POM_PRT_c;

% Grazing on protozoa by metazoa.

graz_PRT_c = tfun(y) * mu_MZ * y.iMZc * y.iPRTc^2 / ...
    (y.iPRTc^2 + g_prt^2 + (y.iTRc * g_prt / g_tr)^2);
graz_PRT_n = graz_PRT_c * y.iPRTn / y.iPRTc;
graz_PRT_p = graz_PRT_c * y.iPRTp / y.iPRTc;

% Net Growth

dydtt.iPRTc = growPRTc - resp_PRT - excr_PRT_LDOC - excr_PRT_SDOC - excr2_PRT_SDOC - POM_PRT_c - graz_PRT_c;
dydtt.iPRTn = growPRTn - remi_PRT_n - excr_PRT_LDON - excr_PRT_SDON - excr2_PRT_SDON - POM_PRT_n - graz_PRT_n;
dydtt.iPRTp = growPRTp - remi_PRT_p - excr_PRT_LDOP - excr_PRT_SDOP - excr2_PRT_SDOP - POM_PRT_p - graz_PRT_p;

end