%% Noah Germolus, 11 Jan 2021
% This is the system of differential equations representing phytoplankton
% processes in Luo et al. (2010)

function [dydtt, fluxPHYno3, fluxPHYnh4, fluxPHYpo4] = ode_mod_phy(y)
%-----------------------------------------------------------------------
%      Phytoplankton Processes
%-----------------------------------------------------------------------
% 1. Gross Grow
% Determining growth rate limits
Nfunc_phy_n = ((y.iPHYn/y.iPHYc)-q_PHY_n)/(q_PHY_n_rdf - q_PHY_n);
Nfunc_phy_p = ((y.iPHYp/y.iPHYc)-q_PHY_p)/(q_PHY_p_rdf - q_PHY_p);

temp = max(min(Nfunc_phy_n, Nfunc_phy_p, 1),0); % Force elemental limits
% on rate to between 0 and 100%

Pmax_phy = mu_PHY * tfun(y) * temp; % Max specific growth rate.

% Light limitation and specific gross PP (as carbon)
growPHYc = y.iPHYc * Pmax_phy *...
    (1 - exp(-alpha_PHY_chl * (y.iPHYchl/y.iPHYc) * y.PAR / Pmax_phy)) *...
    exp(-beta_PHY * y.PAR);

% Nitrogen assimilation

Vmax_phy_n = (q_PHY_n_max - (y.iPHYn/y.iPHYc))/(q_PHY_n_max - q_PHY_n_rdf);
temp = max(min(Vmax_phy_n,1),0); % Yet again forcing this.

growPHYnh4 = y.iPHYc * vref_PHY_n * tfun(y) * temp * ...
    y.iNH4 / (y.iNH4 + k_PHY_nh4 + (y.iNO3*k_PHY_nh4/k_PHY_no3));

growPHYno3 = y.iPHYc * vref_PHY_n * tfun(y) * temp * ...
    y.iNO3 / (y.iNO3 + k_PHY_no3 + (y.iNH4*k_PHY_no3/k_PHY_nh4));

growPHYn = growPHYnh4 + growPHYno3; % Total growth on nitrogen limit.

% Phosphorus limitation.

Vmax_phy_p = (q_PHY_p_max - (y.iPHYp/y.iPHYc))/(q_PHY_p_max - q_PHY_p_rdf);
temp = max(min(Vmax_phy_p,1),0); % Yet again forcing this.

% Assuming inorganic P is the limiter and supply, here's the max growth.

growPHYpo4 = y.iPHYc * vref_PHY_p * tfun(y) * temp * ...
    y.iPO4 / (y.iPO4 + k_PHY_PO4);

% The use of nitrate requires reduction, which requires energy. Here we
% evaluate the respiration necessary for this.

respPHY = zeta_no3 * growPHYno3;

% Chlorophyll Production

growPHYchl = theta * growPHYn * growPHYc / ...
    (alpha_PHY_chl * y.iPHYchl * y.PAR * exp(-beta_PHY * i.PAR));

% DOM excretion!

% Passive (fixed proportion)

excr_PHY_c_psv = ex_PHY_psv * y.iPHYc;
excr_PHY_n_psv = ex_PHY_psv * y.iPHYn;
excr_PHY_p_psv = ex_PHY_psv * y.iPHYp;

% Active carbohydrate excretion

excr_PHY_c_cho = ex_PHY_cho * growPHYc;

% Active excretion to adjust stoichiometry.

excr_PHY_c_act = 0.5 * y.iPHYc * max(0,...
    1-(y.iPHYn/(y.iPHYc * q_PHY_n_rdf)),...
    1-(y.iPHYp/(y.iPHYc * q_PHY_p_rdf)));

if excr_PHY_c_act > 0
    excr_PHY_n_act = 0.5 * 0.25 * y.iPHYn * max(0,...
        1 - (y.iPHYp / y.iPHYn)/(q_PHY_p_rdf / q_PHY_n_rdf));
    
    excr_PHY_p_act = 0.5 * 0.25 * y.iPHYp * max(0,...
        1 - (y.iPHYn / y.iPHYp)/(q_PHY_n_rdf / q_PHY_p_rdf));
    
else
    excr_PHY_n_act = 0;
    excr_PHY_p_act = 0;
end

% Partitioning of labile and semilabile DOM

excr_PHY_LDOC = excr_PHY_c_psv + (0.75 * excr_PHY_c_cho);
excr_PHY_LDON = excr_PHY_n_psv;
excr_PHY_LDOP = excr_PHY_p_psv;

excr_PHY_SDOC = excr_PHY_c_act + (0.25 * excr_PHY_c_cho);
excr_PHY_SDON = excr_PHY_n_act;
excr_PHY_SDOP = excr_PHY_p_act;

% POM production by aggregation

POM_PHY_c = pom_PHY * y.iPHYc^2;
POM_PHY_n = (y.iPHYn/y.iPHYc) * POM_PHY_c;
POM_PHY_p = (y.iPHYp/y.iPHYc) * POM_PHY_c;
POM_PHY_chl = (y.iPHYchl/y.iPHYc) * POM_PHY_c;

% Grazing

graz_PHY_c = tfun(y) * mu_PRT * y.iPRTc * y.iPHYc^2 / ...
    (y.iPHYc^2 + g_PHY^2 +...
    (y.iUNc*g_phy/g_un)^2 +...
    (y.iBAc*g_phy/g_ba)^2);
graz_PHY_n = (y.iPHYn/y.iPHYc) * graz_PHY_c;
graz_PHY_p = (y.iPHYp/y.iPHYc) * graz_PHY_c;
graz_PHY_chl = (y.iPHYchl/y.iPHYc) * graz_PHY_c;

% Net growth rate

dydtt.iPHYc = growPHYc - respPHY - excr_PHY_LDOC - excr_PHY_SDOC - POM_PHY_c - graz_PHY_c;
dydtt.iPHYn = growPHYn           - excr_PHY_LDON - excr_PHY_SDON - POM_PHY_n - graz_PHY_n;
dydtt.iPHYp = growPHYp           - excr_PHY_LDOP - excr_PHY_SDOP - POM_PHY_p - graz_PHY_p;
dydtt.iPHYchl = growPHYchl                                           - POM_PHY_chl - graz_PHY_chl;


end