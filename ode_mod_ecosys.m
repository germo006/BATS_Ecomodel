%% Noah Germolus, 13 Jan 2021
% This is the system of differential equations in Luo et al. (2010). It
% might actually be the play to just combine all 11 sections of the
% equations in one place, actually.


function [dydtt] = ode_mod_ecosys(y)

%%-----------------------------------------------------------------------
%      Temperature Effects
%-----------------------------------------------------------------------

tfun = @(T) exp(-AE*((1/(T+273.15))-(1/(25+273.15))));

%%-----------------------------------------------------------------------
%      Phytoplankton Processes
%-----------------------------------------------------------------------
% 1. Gross Grow
% Determining growth rate limits
Nfunc_phy_n = ((y.iPHYn/y.iPHYc)-q_PHY_n)/(q_PHY_n_rdf - q_PHY_n);
Nfunc_phy_p = ((y.iPHYp/y.iPHYc)-q_PHY_p)/(q_PHY_p_rdf - q_PHY_p);

temp = max(min(Nfunc_phy_n, Nfunc_phy_p, 1),0); % Force elemental limits
% on rate to between 0 and 100%

Pmax_phy = mu_PHY * tfun(y.T) * temp; % Max specific growth rate.

% Light limitation and specific gross PP (as carbon)
growPHYc = y.iPHYc * Pmax_phy *...
    (1 - exp(-alpha_PHY_chl * (y.iPHYchl/y.iPHYc) * y.PAR / Pmax_phy)) *...
    exp(-beta_PHY * y.PAR);

% Nitrogen assimilation

Vmax_phy_n = (q_PHY_n_max - (y.iPHYn/y.iPHYc))/(q_PHY_n_max - q_PHY_n_rdf);
temp = max(min(Vmax_phy_n,1),0); % Yet again forcing this.

growPHYnh4 = y.iPHYc * vref_PHY_n * tfun(y.T) * temp * ...
    y.iNH4 / (y.iNH4 + k_PHY_nh4 + (y.iNO3*k_PHY_nh4/k_PHY_no3));

growPHYno3 = y.iPHYc * vref_PHY_n * tfun(y.T) * temp * ...
    y.iNO3 / (y.iNO3 + k_PHY_no3 + (y.iNH4*k_PHY_no3/k_PHY_nh4));

growPHYn = growPHYnh4 + growPHYno3; % Total growth on nitrogen limit.

% Phosphorus limitation.

Vmax_phy_p = (q_PHY_p_max - (y.iPHYp/y.iPHYc))/(q_PHY_p_max - q_PHY_p_rdf);
temp = max(min(Vmax_phy_p,1),0); % Yet again forcing this.

% Assuming inorganic P is the limiter and supply, here's the max growth.

growPHYpo4 = y.iPHYc * vref_PHY_p * tfun(y.T) * temp * ...
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

graz_PHY_c = tfun(y.T) * mu_PRT * y.iPRTc * y.iPHYc^2 / ...
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

%%-----------------------------------------------------------------------
%      Trichdesmium Processes
%-----------------------------------------------------------------------

% The ol' nutrient limitation thing.
Nfunc_TR_n = ((y.iTRn/y.iTRc)-q_TR_n)/(q_TR_n_rdf - q_TR_n);
Nfunc_TR_p = ((y.iTRp/y.iTRc)-q_TR_p)/(q_TR_p_rdf - q_TR_p);

temp = max(min(Nfunc_TR_n, Nfunc_TR_p, 1),0);

Pmax_TR = mu_TR * tfun(y.T) * temp; % Max specific growth rate.

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

growTRnh4 = y.iTRc * vref_TR_n * tfun(y.T) * temp * ...
    y.iNH4 / (y.iNH4 + k_TR_nh4 + (y.iNO3*k_TR_nh4/k_TR_no3));

growTRno3 = y.iTRc * vref_TR_n * tfun(y.T) * temp * ...
    y.iNO3 / (y.iNO3 + k_TR_no3 + (y.iNH4*k_TR_no3/k_TR_nh4));

% Fixation

if y.T <= 20
    nfixTRmax = 0;
else
    nfixTRmax = tfun(y.T) * Vmax_TR_n * ...
        (((growTRc - (zeta_no3*growTRno3)) * q_TR_n_max) - ...
        growTRno3 - growTRnh4)/...
        (1 + (zeta_nfix * q_TR_n_max));
    nfixTRmax = max(0, nfixTRmax);
end

% total assimilation of N

growTRn = min(y.iTRc * vref_TR_n * tfun(y.T) * Vmax_TR_n,...
    growTRnh4 + growTRno3 + nfixTRmax);

% Actual (not max) fixation rate

growTRnfix = growTRn - growTRnh4 - growTRno3;

% Phosphorus assimilation

Vmax_TR_p = (q_TR_p_max - (y.iTRp/y.iTRc))/(q_TR_p_max - q_TR_p_rdf);
temp = max(min(Vmax_TR_p,1),0); % Yet again forcing this.

% Assuming inorganic P is the limiter and supply, here's the max growth.

growTRpo4 = y.iTRc * vref_TR_p * tfun(y.T) * temp * ...
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

graz_TR_c = tfun(y.T) * mu_MZ * y.iMZc * y.iTRc^2 / ...
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

%%-----------------------------------------------------------------------
%      Unicellular N2-fixer Processes
%-----------------------------------------------------------------------

% The ol' nutrient limitation thing.
Nfunc_UN_n = ((y.iUNn/y.iUNc)-q_UN_n)/(q_UN_n_rdf - q_UN_n);
Nfunc_UN_p = ((y.iUNp/y.iUNc)-q_UN_p)/(q_UN_p_rdf - q_UN_p);

temp = max(min(Nfunc_UN_n, Nfunc_UN_p, 1),0);

Pmax_UN = mu_UN * tfun(y.T) * temp; % Max specific growth rate.

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

growUNnh4 = y.iUNc * vref_UN_n * tfun(y.T) * temp * ...
    y.iNH4 / (y.iNH4 + k_UN_nh4 + (y.iNO3*k_UN_nh4/k_UN_no3));

growUNno3 = y.iUNc * vref_UN_n * tfun(y.T) * temp * ...
    y.iNO3 / (y.iNO3 + k_UN_no3 + (y.iNH4*k_UN_no3/k_UN_nh4));

% Fixation

if y.T <= 15
    nfixUNmax = 0;
else
    nfixUNmax = tfun(y.T) * Vmax_UN_n * ...
        (((growUNc - (zeta_no3*growUNno3)) * q_UN_n_max) - ...
        growUNno3 - growUNnh4)/...
        (1 + (zeta_nfix * q_UN_n_max));
    nfixUNmax = max(0, nfixUNmax);
end

% total assimilation of N

growUNn = min(y.iUNc * vref_UN_n * tfun(y.T) * Vmax_UN_n,...
    growUNnh4 + growUNno3 + nfixUNmax);

% Actual (not max) fixation rate

growUNnfix = growUNn - growUNnh4 - growUNno3;

% Phosphorus assimilation

Vmax_UN_p = (q_UN_p_max - (y.iUNp/y.iUNc))/(q_UN_p_max - q_UN_p_rdf);
temp = max(min(Vmax_UN_p,1),0); % Yet again forcing this.

% Assuming inorganic P is the limiter and supply, here's the max growth.

growUNp = y.iUNc * vref_UN_p * tfun(y.T) * temp * ...
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

graz_UN_c = tfun(y.T) * mu_PRT * y.iPRTc * y.iUNc^2 / ...
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

%%-----------------------------------------------------------------------
%      Bacterial Processes
%-----------------------------------------------------------------------
% 1. Gross Grow
% Maximum possible C amount for bacterial use
ALC = y.iLDOMc;
ASC = y.iSDOMc * r_SDOM;
% Carbon Usage
Nfunc_ba_n = y.iBAn/(y.iBAc*q_BA_n); %bacterial N/C quota /ref quota
Nfunc_ba_p = y.iBAp/(y.iBAc*q_BA_p);
temp = min(Nfunc_ba_n, Nfunc_ba_p);
temp = min(temp, 1); % Make sure the limit is unity (not the case for mtabs?)
growBAldoc = mu_BA * y.iBAc * tfun(y.T) * temp * ALC / (ALC+ k_DOM + ASC);
growBAsdoc = mu_BA * y.iBAc * tfun(y.T) * temp * ASC / (ASC+ k_DOM + ALC);
%^^^I added the tfun back into this version (it wasn't in the one H. Kim
%sent.)
% DON and DOP usage
growBAldon = growBAldoc*y.iLDOMn/y.iLDOMc; % available labile N
growBAldop = growBAldoc*y.iLDOMp/y.iLDOMc; % available labile P
growBAsdon = growBAsdoc * min(q_BA_n,...
    y.iSDOMn/y.iSDOMc + f_BAslct/Nfunc_ba_n*(q_BA_n-(y.iSDOMn/y.iSDOMc)));
growBAsdop = growBAsdoc * min(q_BA_p,...
    y.iSDOMp/y.iSDOMc + f_BAslct/Nfunc_ba_p*(q_BA_p-(y.iSDOMp/y.iSDOMc)));

% inorganic nutrients uptake
growBAnh4 = (growBAldon * y.iNH4 * min(1, 1/Nfunc_ba_n))/ y.iLDOMn;
if Nfunc_ba_n < 1
    growBAno3 = min(0.1 * ((growBAldon + growBAsdon) / (y.iLDOMn + y.iSDOMn)) * y.iNO3 * min(1, 1/Nfunc_ba_n),...
        ((growBAldon + growBAsdon) / (y.iLDOMn + y.iSDOMn)) * (y.iNO3  + y.iNH4) - growBAnh4);
    growBAno3 = max(0, growBAno3);
else
    growBAno3 = 0;
end
growBApo4 = (growBAldop / y.iLDOMp) * y(iPO4) * min(1, 1/Nfunc_ba_p);
% Bacteria gross growth
growBAc = growBAldoc + growBAsdoc;
growBAn = growBAldon + growBAsdon + growBAnh4 + growBAno3;
growBAp = growBAldop + growBAsdop + growBApo4;

% 2. respiration
% zeta is an energy required for NO3 use
respBA = r_BAresp_1 * y.iBAc + zeta * growBAno3 + ...
    (r_BAresp_min + (r_BAresp_max - r_BAresp_min)*exp(-b_BAresp*growBAc) ) *...
    growBAc;

% 3. excreting refractory DOM
refrBAc = r_BArefr * y.iBAc;
refrBAn = q_refrDOM_n * refrBAc;
refrBAp = q_refrDOM_p * refrBAc;

% 4. excreting semi-labile DOM and regenerating DIN
if (y.iBAc < y.iBAn/q_BA_n) && (y.iBAc < y.iBAp/q_BA_p) %Carbon in short
    excrBAc = 0;
    excrBAn = 0;
    excrBAp = 0;
    remiBAn = r_BAremi * (y.iBAn - y.iBAc * q_BA_n);
    remiBAp = r_BAremi * (y.iBAp - y.iBAc * q_BA_p);
elseif (y.iBAc > y.iBAn/q_BA_n) && (y.iBAp/q_BA_p > y.iBAn/q_BA_n) %Nitrogen in short
    excrBAc = r_BAadju * (y.iBAc - (y.iBAn/q_BA_n));
    excrBAn = 0;
    excrBAp = r_BAadju * (y.iBAp - ((y.iBAn/q_BA_n) * q_BA_p));
    remiBAn = 0;
    remiBAp = 0;
else %Phosphorus in short
    excrBAc = r_BAadju * (y.iBAc - (y.iBAp/q_BA_p));
    excrBAn = r_BAadju * (y.iBAn - ((y.iBAp/q_BA_p) * q_BA_n));
    excrBAp = 0;
    remiBAn = 0;
    remiBAp = 0;
end

%6. removal by grazing
grazBAc = mu_PRT * y.iPRTc * y.iBAc^2 /...
    (y.iBAc^2 + g_ba^2 + ...
    (y.iSPc * g_ba/g_sp)^2 + ...
    (y.iUNc*g_ba/g_un)^2);
grazBAn = grazBAc * y.iBAn / y.iBAc;
grazBAp = grazBAc * y.iBAp / y.iBAc;

%6b. Mortality due to viruses
mortBAc = r_BAmort * y.iBAc;
mortBAn = r_BAmort * y.iBAn;
mortBAp = r_BAmort * y.iBAp;

%7. BA Derivs
dydtt.iBAc = (growBAc - refrBAc - excrBAc - grazBAc - respBA - mortBAc)/SecPerDay;
dydtt.iBAn = (growBAn - refrBAn - excrBAn - remiBAn - grazBAn - mortBAn)/SecPerDay;
dydtt.iBAp = (growBAp - refrBAp - excrBAp - remiBAp - grazBAp - mortBAp)/SecPerDay;

%8. Flux of inorganic nutrients through bacteria
dydtt.fluxBAnh4 = growBAnh4 - remiBAn;
dydtt.fluxBAno3 = growBAno3;
dydtt.fluxBApo4 = growBApo4 - remiBAp;

%%-----------------------------------------------------------------------
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

resp_PRT = (resp_B_PRT * tfun(y.T) * y.iPRTc) + (resp_A_PRT * growPRTc)

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

graz_PRT_c = tfun(y.T) * mu_MZ * y.iMZc * y.iPRTc^2 / ...
    (y.iPRTc^2 + g_prt^2 + (y.iTRc * g_prt / g_tr)^2);
graz_PRT_n = graz_PRT_c * y.iPRTn / y.iPRTc;
graz_PRT_p = graz_PRT_c * y.iPRTp / y.iPRTc;

% Net Growth

dydtt.iPRTc = growPRTc - resp_PRT - excr_PRT_LDOC - excr_PRT_SDOC - excr2_PRT_SDOC - POM_PRT_c - graz_PRT_c;
dydtt.iPRTn = growPRTn - remi_PRT_n - excr_PRT_LDON - excr_PRT_SDON - excr2_PRT_SDON - POM_PRT_n - graz_PRT_n;
dydtt.iPRTp = growPRTp - remi_PRT_p - excr_PRT_LDOP - excr_PRT_SDOP - excr2_PRT_SDOP - POM_PRT_p - graz_PRT_p;

%%-----------------------------------------------------------------------
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

resp_MZ = (resp_B_MZ * tfun(y.T) * y.iMZc) + (resp_A_MZ * growMZc)

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

%%-----------------------------------------------------------------------
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

%%-----------------------------------------------------------------------
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

%%-----------------------------------------------------------------------
%      Change in Dissolved Organic Matter
%-----------------------------------------------------------------------

% Conversion of SDOM to RDOM

REFR_SDOC = ex_REFR_SDOM * y.iSDOMc * ...
    exp( 1 - min((y.iSDOMn / (y.iSDOMc*q_REFR_n)),...
    (y.iSDOMp / (y.iSDOMc*q_REFR_p))) ) + ...
    max( 0 , y.iSDOMc - (y.iSDOMn/q_REFR_n), y.iSDOMc - (y.iSDOMp/q_REFR_p));
REFR_SDON = q_REFR_n * REFR_SDOC;
REFR_SDOP = q_REFR_p * REFR_SDOC;

% Change in DOM species

dydtt.iLDOMc = excr_PHY_LDOC + excr_TR_LDOC + excr_UN_LDOC + excr_PRT_LDOC + ...
    excr_MZ_LDOC + mortBAc - growBAldoc;
dydtt.iLDOMn = excr_PHY_LDON + excr_TR_LDON + excr_UN_LDON + excr_PRT_LDON + ...
    excr_MZ_LDON + mortBAn - growBAldon;
dydtt.iLDOMp = excr_PHY_LDOP + excr_TR_LDOP + excr_UN_LDOP + excr_PRT_LDOP + ...
    excr_MZ_LDOP + mortBAp - growBAldop;

dydtt.iSDOMc = excr_PHY_SDOC + excr_TR_SDOC + excr_UN_SDOC + excrBAc +...
    excr_PRT_SDOC + excr2_PRT_SDOC + excr_MZ_SDOC + excr2_MZ_SDOC +...
    excr_HZ_SDOC + DISS_DET_c - REFR_SDOC - growBAsdoc;

dydtt.iSDOMn = excr_PHY_SDON + excr_TR_SDON + excr_UN_SDON + excrBAn +...
    excr_PRT_SDON + excr2_PRT_SDON + excr_MZ_SDON + excr2_MZ_SDON +...
    excr_HZ_SDON + DISS_DET_n - REFR_SDON - growBAsdon;

dydtt.iSDOMp = excr_PHY_SDOP + excr_TR_SDOP + excr_UN_SDOP + excrBAp +...
    excr_PRT_SDOP + excr2_PRT_SDOP + excr_MZ_SDOP + excr2_MZ_SDOP +...
    excr_HZ_SDOP + DISS_DET_p - REFR_SDOP - growBAsdop;

%%-----------------------------------------------------------------------
%      Diagnostic Variables
%-----------------------------------------------------------------------

% Primary Production

dydtt.PrPr = growPHYc + growTRc + growUNc - respPHY - respTR - respUN - ...
    excr_PHY_LDOC - excr_TR_LDOC - excr_UN_LDOC -...
    excr_PHY_SDOC - excr_TR_SDOC - excr_UN_SDOC;

% Bacterial Production

dydtt.BPr = growBAc - respBA;

end