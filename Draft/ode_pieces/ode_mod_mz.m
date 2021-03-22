%% Noah Germolus, 13 Jan 2021
% This is the system of differential equations representing metazoan
% grazers in Luo et al. (2010)


function [dydtt] = ode_mod_mz(y)
%-----------------------------------------------------------------------
%      Metazoan Grazers
%-----------------------------------------------------------------------

% Gross growth rate

growMZc = graz_PRT_c + graz_TR_c;
growMZn = graz_PRT_n + graz_TR_n;
growMZp = graz_PRT_p + graz_TR_p;

% DOM Excretion

excr_MZ_LDOC = ex_MZ * f_ex_MZ * growMZc;
excr_MZ_LDON = ex_MZ * f_ex_MZ * growMZn;
excr_MZ_LDOP = ex_MZ * f_ex_MZ * growMZp;

excr_MZ_SDOC = ex_MZ * (1 - f_ex_MZ) * growMZc;
excr_MZ_SDON = ex_MZ * (1 - f_ex_MZ) * growMZn * (y.iMZn/(y.iMZc*q_MZ_n));
excr_MZ_SDOP = ex_MZ * (1 - f_ex_MZ) * growMZp * (y.iMZp/(y.iMZc*q_MZ_p));

% Respiration

resp_MZ = (resp_B_MZ * tfun(y) * y.iMZc) + (resp_A_MZ * growMZc)

% Semilabile DOM stoich adjustment

excr2_MZ_SDOC = ex_adj_MZ * y.iMZc * max(0,...
    1 - (y.iMZn/(y.iMZc*q_MZ_n)),...
    1 - (y.iMZp/(y.iMZc*q_MZ_p)));
excr2_MZ_SDON = 0.5 * excr2_MZ_SDOC * (y.iMZn/y.iMZc);
excr2_MZ_SDOP = 0.5 * excr2_MZ_SDOC * (y.iMZp/y.iMZc);

% Inorganic nutrient remin to adjust stoichiometry.

remi_MZ_n = remi_mz * max(0, y.iMZn - (q_MZ_n * y.iMZc),...
    y.iMZn - (q_MZ_n * y.iMZp / q_MZ_p));
remi_MZ_p = remi_mz * max(0, y.iMZp - (q_MZ_p * y.iMZc),...
    y.iMZp - (q_MZ_p * y.iMZn / q_MZ_n));

% POM production

POM_MZ_c = pom_MZ * growMZc;
POM_MZ_n = q_POM_n * POM_MZ_c;
POM_MZ_p = q_POM_p * POM_MZ_c;

% Refractory DOM release.

REFR_MZ_c = ex_refr_MZ * y.iMZc;
REFR_MZ_n = REFR_MZ_c * q_REFR_n;
REFR_MZ_p = REFR_MZ_c * q_REFR_p;

% Removal by higher level zoops

REMV_MZ_c = remv_MZ * y.iMZc^2;
REMV_MZ_n = (y.iMZn / y.iMZc) * REMV_MZ_c;
REMV_MZ_p = (y.iMZp / y.iMZc) * REMV_MZ_c;

% POM production by higher zoops

POM_HZ_c = f_POM_HZ * REMV_MZ_c;
POM_HZ_n = f_POM_HZ * REMV_MZ_n;
POM_HZ_p = f_POM_HZ * REMV_MZ_p;

% SDOM production by higher zoops

excr_HZ_SDOC = f_SDOM_HZ * REMV_MZ_c;
excr_HZ_SDON = f_SDOM_HZ * REMV_MZ_n;
excr_HZ_SDOP = f_SDOM_HZ * REMV_MZ_p;

% Inorganic nuts remineralized by HZ

remi_HZ_n = REMV_MZ_n - POM_HZ_n - excr_HZ_SDON;
remi_HZ_p = REMV_MZ_p - POM_HZ_p - excr_HZ_SDOP;

% Net growth of metazoa.

dydtt.iMZc = growMZc - resp_MZ - excr_MZ_LDOC - excr_MZ_SDOC - excr2_MZ_SDOC - POM_MZ_c - REFR_MZ_c - REMV_MZ_c;
dydtt.iMZn = growMZn - remi_MZ_n - excr_MZ_LDON - excr_MZ_SDON - excr2_MZ_SDON - POM_MZ_n - REFR_MZ_n - REMV_MZ_n;
dydtt.iMZp = growMZp - remi_MZ_p - excr_MZ_LDOP - excr_MZ_SDOP - excr2_MZ_SDOP - POM_MZ_p - REFR_MZ_p - REMV_MZ_p;


end