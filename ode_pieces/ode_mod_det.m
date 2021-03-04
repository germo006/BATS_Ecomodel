%% Noah Germolus, 13 Jan 2021
% This is the system of differential equations representing detritus
% in Luo et al. (2010)


function [dydtt] = ode_mod_det(y)
%-----------------------------------------------------------------------
%      Detrital Processes
%-----------------------------------------------------------------------

% Dissolution

DISS_DET_c = diss * y.iDETc;
DISS_DET_n = prf_n * diss * y.iDETn;
DISS_DET_p = prf_p * diss * y.iDETp;

% Net Change

dydtt.iDETc = POM_PHY_c + POM_TR_c + POM_UN_c + POM_PRT_c + POM_MZ_c + POM_HZ_c - DISS_DET_c; 
dydtt.iDETn = POM_PHY_n + POM_TR_n + POM_UN_n + POM_PRT_n + POM_MZ_n + POM_HZ_n - DISS_DET_n;
dydtt.iDETp = POM_PHY_p + POM_TR_p + POM_UN_p + POM_PRT_p + POM_MZ_p + POM_HZ_p - DISS_DET_p;

end