%% Noah Germolus, 12 Jan 2021
% This is the system of differential equations representing unicellular
% N-fixers  in Luo et al. (2010)
% This process was disabled for the Arabian Sea and the Equatorial Pacific
% in the original paper.

function [dydtt] = ode_mod_un(y)
%-----------------------------------------------------------------------
%      Unicellular N2-fixer Processes
%-----------------------------------------------------------------------

% The ol' nutrient limitation thing.
Nfunc_UN_n = ((y.iUNn/y.iUNc)-q_UN_n)/(q_UN_n_rdf - q_UN_n);
Nfunc_UN_p = ((y.iUNp/y.iUNc)-q_UN_p)/(q_UN_p_rdf - q_UN_p);

temp = max(min(Nfunc_UN_n, Nfunc_UN_p, 1),0);

Pmax_UN = mu_UN * tfun(y) * temp; % Max specific growth rate.

% Light limitation and specific gross PP (as carbon)

if Pmax_UN == 0
    growUNc = 0;
else
    growUNc = y.iUNc * Pmax_UN *...
        (1 - exp(-alpha_UN_chl * (y.iUNchl/y.iUNc) * y.PAR / Pmax_UN)) *...
        exp(-beta_UN * y.PAR);
end

% N assimilation, including fixation.

Vmax_UN_n = (q_UN_n_max - (y.iUNn/y.iUNc))/(q_UN_n_max - q_UN_n_rdf);
temp = max(min(Vmax_UN_n,1),0); % Yet again forcing this.

growUNnh4 = y.iUNc * vref_UN_n * tfun(y) * temp * ...
    y.iNH4 / (y.iNH4 + k_UN_nh4 + (y.iNO3*k_UN_nh4/k_UN_no3));

growUNno3 = y.iUNc * vref_UN_n * tfun(y) * temp * ...
    y.iNO3 / (y.iNO3 + k_UN_no3 + (y.iNH4*k_UN_no3/k_UN_nh4));

% Fixation

if y.T <= 15
    nfixUNmax = 0;
else
    nfixUNmax = tfun(y) * Vmax_UN_n * ...
        (((growUNc - (zeta_no3*growUNno3)) * q_UN_n_max) - ...
        growUNno3 - growUNnh4)/...
        (1 + (zeta_nfix * q_UN_n_max));
    nfixUNmax = max(0, nfixUNmax);
end

% total assimilation of N

growUNn = min(y.iUNc * vref_UN_n * tfun(y) * Vmax_UN_n,...
    growUNnh4 + growUNno3 + nfixUNmax);

% Actual (not max) fixation rate

growUNnfix = growUNn - growUNnh4 - growUNno3;

% Phosphorus assimilation

Vmax_UN_p = (q_UN_p_max - (y.iUNp/y.iUNc))/(q_UN_p_max - q_UN_p_rdf);
temp = max(min(Vmax_UN_p,1),0); % Yet again forcing this.

% Assuming inorganic P is the limiter and supply, here's the max growth.

growUNp = y.iUNc * vref_UN_p * tfun(y) * temp * ...
    y.iPO4 / (y.iPO4 + k_UN_PO4);

% Respiration required for nitrogen metabolism, only this time it also
% includes n fixation.

respUN = (zeta_no3 * growUNno3) + (zeta_nfix * growUNnfix);

% Chlorophyll Production

if y.PAR == 0
    growUNchl=0;
else
    growUNchl = theta * growUNn * growUNc / ...
        (alpha_UN_chl * y.iUNchl * y.PAR * exp(-beta_UN * i.PAR));
end

% DOM excretion!

% Passive (fixed proportion)

excr_UN_c_psv = ex_UN_psv * y.iUNc;
excr_UN_n_psv = ex_UN_psv * y.iUNn;
excr_UN_p_psv = ex_UN_psv * y.iUNp;

% Active carbohydrate excretion

excr_UN_c_cho = ex_UN_cho * growUNc;

% Active excretion to adjust stoichiometry.

excr_UN_c_act = 0.5 * y.iUNc * max(0,...
    1-(y.iUNn/(y.iUNc * q_UN_n_rdf)),...
    1-(y.iUNp/(y.iUNc * q_UN_p_rdf)));

if excr_UN_c_act > 0
    excr_UN_n_act = 0.5 * 0.25 * y.iUNn * max(0,...
        1 - (y.iUNp / y.iUNn)/(q_UN_p_rdf / q_UN_n_rdf));
    
    excr_UN_p_act = 0.5 * 0.25 * y.iUNp * max(0,...
        1 - (y.iUNn / y.iUNp)/(q_UN_n_rdf / q_UN_p_rdf));
    
else
    excr_UN_n_act = 0;
    excr_UN_p_act = 0;
end

% Active release of newly fixed N

excr_UN_nfix_DON = 0.5 * ex_UN_nfix * growUNnfix * Nfunc_UN_n;

% Partitioning of labile and semilabile DOM

excr_UN_LDOC = excr_UN_c_psv + (0.75 * excr_UN_c_cho);
excr_UN_LDON = excr_UN_n_psv + excr_UN_nfix_DON;
excr_UN_LDOP = excr_UN_p_psv;

excr_UN_SDOC = excr_UN_c_act + (0.25 * excr_UN_c_cho);
excr_UN_SDON = excr_UN_n_act;
excr_UN_SDOP = excr_UN_p_act;

% Ammonium excretion

excr_UN_n_nh4 = 0.5 * ex_UN_nfix * growUNnfix * Nfunc_UN_n;

% POM production by aggregation

POM_UN_c = pom_UN * y.iUNc^2;
POM_UN_n = (y.iUNn/y.iUNc) * POM_UN_c;
POM_UN_p = (y.iUNp/y.iUNc) * POM_UN_c;
POM_UN_chl = (y.iUNchl/y.iUNc) * POM_UN_c;

% Grazing

graz_UN_c = tfun(y) * mu_PRT * y.iPRTc * y.iUNc^2 / ...
    (y.iUNc^2 + g_UN^2 +...
    (y.iPHYc*g_un/g_phy)^2 + ...
    (y.iBAc*g_un/g_ba)^2);
graz_UN_n = (y.iUNn/y.iUNc) * graz_UN_c;
graz_UN_p = (y.iUNp/y.iUNc) * graz_UN_c;
graz_UN_chl = (y.iUNchl/y.iUNc) * graz_UN_c;

% Net Growth

dydtt.iUNc = growUNc - respUN - excr_UN_LDOC - excr_UN_SDOC - POM_UN_c - graz_UN_c;
dydtt.iUNn = growUNn - excr_UN_n_nh4 - excr_UN_LDON - excr_UN_SDON - POM_UN_n - graz_UN_n;
dydtt.iUNp = growUNp           - excr_UN_LDOP - excr_UN_SDOP - POM_UN_p - graz_UN_p;
dydtt.iUNchl = growUNchl                                           - POM_UN_chl - graz_UN_chl;

end