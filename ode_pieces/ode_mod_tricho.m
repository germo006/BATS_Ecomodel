%% Noah Germolus, 12 Jan 2021
% This is the system of differential equations representing Trichdesmium
% processes in Luo et al. (2010)
% This process was disabled for the Arabian Sea and the Equatorial Pacific
% in the original paper.

function [dydtt] = ode_mod_tricho(y)
%-----------------------------------------------------------------------
%      Trichdesmium Processes
%-----------------------------------------------------------------------

% The ol' nutrient limitation thing.
Nfunc_TR_n = ((y.iTRn/y.iTRc)-q_TR_n)/(q_TR_n_rdf - q_TR_n);
Nfunc_TR_p = ((y.iTRp/y.iTRc)-q_TR_p)/(q_TR_p_rdf - q_TR_p);

temp = max(min(Nfunc_TR_n, Nfunc_TR_p, 1),0);

Pmax_TR = mu_TR * tfun(y) * temp; % Max specific growth rate.

% Light limitation and specific gross PP (as carbon)

if Pmax_TR == 0
    growTRc = 0;
else
    growTRc = y.iTRc * Pmax_TR *...
        (1 - exp(-alpha_TR_chl * (y.iTRchl/y.iTRc) * y.PAR / Pmax_TR)) *...
        exp(-beta_TR * y.PAR);
end

% N assimilation, including fixation.

Vmax_TR_n = (q_TR_n_max - (y.iTRn/y.iTRc))/(q_TR_n_max - q_TR_n_rdf);
temp = max(min(Vmax_TR_n,1),0); % Yet again forcing this.

growTRnh4 = y.iTRc * vref_TR_n * tfun(y) * temp * ...
    y.iNH4 / (y.iNH4 + k_TR_nh4 + (y.iNO3*k_TR_nh4/k_TR_no3));

growTRno3 = y.iTRc * vref_TR_n * tfun(y) * temp * ...
    y.iNO3 / (y.iNO3 + k_TR_no3 + (y.iNH4*k_TR_no3/k_TR_nh4));

% Fixation

if y.T <= 20
    nfixTRmax = 0;
else
    nfixTRmax = tfun(y) * Vmax_TR_n * ...
        (((growTRc - (zeta_no3*growTRno3)) * q_TR_n_max) - ...
        growTRno3 - growTRnh4)/...
        (1 + (zeta_nfix * q_TR_n_max));
    nfixTRmax = max(0, nfixTRmax);
end

% total assimilation of N

growTRn = min(y.iTRc * vref_TR_n * tfun(y) * Vmax_TR_n,...
    growTRnh4 + growTRno3 + nfixTRmax);

% Actual (not max) fixation rate

growTRnfix = growTRn - growTRnh4 - growTRno3;

% Phosphorus assimilation

Vmax_TR_p = (q_TR_p_max - (y.iTRp/y.iTRc))/(q_TR_p_max - q_TR_p_rdf);
temp = max(min(Vmax_TR_p,1),0); % Yet again forcing this.

% Assuming inorganic P is the limiter and supply, here's the max growth.

growTRpo4 = y.iTRc * vref_TR_p * tfun(y) * temp * ...
    y.iPO4 / (y.iPO4 + k_TR_PO4);

% Picking up PO4 from deep water (kind of unclear how this term makes
% sense)

pickTRpo4 = max(0, picktrpo4 * (y.iTRc * ((q_TR_p_max + q_TR_p_rdf)/2) - ...
    y.iTRp));

% Growth on phosphorus

growTRp = growTRpo4 + pickTRpo4;

% Respiration required for nitrogen metabolism, only this time it also
% includes n fixation.

respTR = (zeta_no3 * growTRno3) + (zeta_nfix * growTRnfix);

% Chlorophyll Production

if y.PAR==0
    growTRchl = 0;
else
    growTRchl = theta * growTRn * growTRc / ...
        (alpha_TR_chl * y.iTRchl * y.PAR * exp(-beta_TR * i.PAR));
end

% DOM excretion!

% Passive (fixed proportion)

excr_TR_c_psv = ex_TR_psv * y.iTRc;
excr_TR_n_psv = ex_TR_psv * y.iTRn;
excr_TR_p_psv = ex_TR_psv * y.iTRp;

% Active carbohydrate excretion

excr_TR_c_cho = ex_TR_cho * growTRc;

% Active excretion to adjust stoichiometry.

excr_TR_c_act = 0.5 * y.iTRc * max(0,...
    1-(y.iTRn/(y.iTRc * q_TR_n_rdf)),...
    1-(y.iTRp/(y.iTRc * q_TR_p_rdf)));

if excr_TR_c_act > 0
    excr_TR_n_act = 0.5 * 0.25 * y.iTRn * max(0,...
        1 - (y.iTRp / y.iTRn)/(q_TR_p_rdf / q_TR_n_rdf));
    
    excr_TR_p_act = 0.5 * 0.25 * y.iTRp * max(0,...
        1 - (y.iTRn / y.iTRp)/(q_TR_n_rdf / q_TR_p_rdf));
    
else
    excr_TR_n_act = 0;
    excr_TR_p_act = 0;
end

% Active release of newly fixed N

excr_TR_nfix_DON = 0.5 * ex_TR_nfix * growTRnfix * Nfunc_TR_n;

% Partitioning of labile and semilabile DOM

excr_TR_LDOC = excr_TR_c_psv + (0.75 * excr_TR_c_cho);
excr_TR_LDON = excr_TR_n_psv + excr_TR_nfix_DON;
excr_TR_LDOP = excr_TR_p_psv;

excr_TR_SDOC = excr_TR_c_act + (0.25 * excr_TR_c_cho);
excr_TR_SDON = excr_TR_n_act;
excr_TR_SDOP = excr_TR_p_act;

% Ammonium excretion

excr_TR_n_nh4 = 0.5 * ex_TR_nfix * growTRnfix * Nfunc_TR_n;

% POM production by aggregation

POM_TR_c = pom_TR * y.iTRc^2;
POM_TR_n = (y.iTRn/y.iTRc) * POM_TR_c;
POM_TR_p = (y.iTRp/y.iTRc) * POM_TR_c;
POM_TR_chl = (y.iTRchl/y.iTRc) * POM_TR_c;

% Grazing

graz_TR_c = tfun(y) * mu_MZ * y.iMZc * y.iTRc^2 / ...
    (y.iTRc^2 + g_TR^2 +...
    (y.iPRTc*g_tr/g_prt)^2);
graz_TR_n = (y.iTRn/y.iTRc) * graz_TR_c;
graz_TR_p = (y.iTRp/y.iTRc) * graz_TR_c;
graz_TR_chl = (y.iTRchl/y.iTRc) * graz_TR_c;

% Net Growth

dydtt.iTRc = growTRc - respTR - excr_TR_LDOC - excr_TR_SDOC - POM_TR_c - graz_TR_c;
dydtt.iTRn = growTRn - excr_TR_n_nh4 - excr_TR_LDON - excr_TR_SDON - POM_TR_n - graz_TR_n;
dydtt.iTRp = growTRp           - excr_TR_LDOP - excr_TR_SDOP - POM_TR_p - graz_TR_p;
dydtt.iTRchl = growTRchl                                           - POM_TR_chl - graz_TR_chl;



end