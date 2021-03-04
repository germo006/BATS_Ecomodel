%% Noah Germolus, 13 Jan 2021
% This is the system of differential equations representing nutrients
% in Luo et al. (2010)


function [dydtt] = ode_mod_nut(y)
%-----------------------------------------------------------------------
%      Inorganic Nutrients
%-----------------------------------------------------------------------

% Nitrification

NTRF = r_ntrf * y.iNH4;

% Dissolved Nutrient Rates of Change

dydtt.iNH4 = dydtt.fluxBAnh4 + remi_PRT_n + remi_MZ_n + remi_HZ_n + ...
    excr_TR_n_nh4 + excr_UN_n_nh4 -...
    growPHYnh4 - growTRnh4 - growUNnh4 - NTRF;
dydtt.iNO3 = dydtt.fluxBAno3 - growPHYno3 - growTRno3 - growUNno3 + NTRF;
dydtt.iPO4 = dydtt.fluxBApo4 + remi_PRT_p + remi_MZ_p + remi_HZ_p -...
    growPHYpo4 - growTRpo4 - growUNpo4;

end