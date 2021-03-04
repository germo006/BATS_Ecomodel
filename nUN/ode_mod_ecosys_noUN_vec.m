%% Noah Germolus, 13 Jan 2021
% This is the system of differential equations in Luo et al. (2010). It
% might actually be the play to just combine all 11 sections of the
% equations in one place, actually.
% replacement for when vectorized
% ode_mod_ecosys_noUN_vec(y_vals, ps_vals, po_vals, ts)
function [dydtt] = ode_mod_ecosys_noUN_vec(y_vals, ps_vals, po_vals, ts, z)
if mod(ts, 10) == 0
    msg = 'Evaluating at t = ' + string(ts) + ' days';
    disp(msg)
end
%%-----------------------------------------------------------------------
%      Temperature Effects and PAR
%-----------------------------------------------------------------------
forceNH4 = 0.01;
forceNO3 = 2;
forcePO4 = 0.3;
T = 21 + randn(1);
MLD = 50;
tfun = @(T) exp(-po_vals(73)*((1/(T+273.15))-(1/(25+273.15))));
pfun = @(t) max(0, -10*cos(2*pi*t));
PAR0 = pfun(ts);
PAR = PAR0*exp(z*(-0.038 + 0.05*(y_vals(15) + y_vals(4))));
Kzfun = @(z, MLD) 1.1e-4*24*3600*exp(-0.01*(z - MLD));
Kz = Kzfun(z, MLD);
%%-----------------------------------------------------------------------
%      Phytoplankton Processes
%-----------------------------------------------------------------------
% 1. Gross Grow
% Determining growth rate limits
Nfunc_phy_n = ((y_vals(14)/y_vals(16))-ps_vals(14))/(ps_vals(12) - ps_vals(14));
Nfunc_phy_p = ((y_vals(13)/y_vals(16))-ps_vals(11))/(ps_vals(9) - ps_vals(11));
temp = max([min([Nfunc_phy_n, Nfunc_phy_p, 1]),0]); % Force elemental limits
% on rate to between 0 and 100%
Pmax_phy = po_vals(36) * tfun(T) * temp; % Max specific growth rate.
% Light limitation and specific gross PP (as carbon)
growPHYc = y_vals(16) * Pmax_phy *...
    (1 - exp(-po_vals(71) * (y_vals(15)/y_vals(16)) * PAR / Pmax_phy)) *...
    exp(-po_vals(68) * PAR);
% Nitrogen assimilation
Vmax_phy_n = (ps_vals(13) - (y_vals(14)/y_vals(16)))/(ps_vals(13) - ps_vals(12));
temp = max([min([Vmax_phy_n,1]),0]); % Yet again forcing this.
growPHYnh4 = y_vals(16) * po_vals(6) * tfun(T) * temp * ...
    y_vals(18) / (y_vals(18) + po_vals(44) + (y_vals(17)*po_vals(44)/po_vals(43)));
growPHYno3 = y_vals(16) * po_vals(6) * tfun(T) * temp * ...
    y_vals(17) / (y_vals(17) + po_vals(43) + (y_vals(18)*po_vals(43)/po_vals(44)));
growPHYn = growPHYnh4 + growPHYno3; % Total growth on nitrogen limit.
% Phosphorus limitation.
Vmax_phy_p = (ps_vals(10) - (y_vals(13)/y_vals(16)))/(ps_vals(10) - ps_vals(9));
temp = max([min([Vmax_phy_p,1]),0]); % Yet again forcing this.
% Assuming inorganic P is the limiter and supply, here's the max growth.
growPHYpo4 = y_vals(16) * po_vals(5) * tfun(T) * temp * ...
    y_vals(12) / (y_vals(12) + po_vals(42));
% The use of nitrate requires reduction, which requires energy. Here we
% evaluate the respiration necessary for this.
respPHY = po_vals(1) * growPHYno3;
% Chlorophyll Production
growPHYchl = po_vals(7) * growPHYn * growPHYc / ...
    (po_vals(71) * y_vals(15) * PAR * exp(-po_vals(68) * PAR));
if isnan(growPHYchl)
    growPHYchl = 0;
end
% DOM excretion!
% Passive (fixed proportion)
excr_PHY_c_psv = po_vals(63) * y_vals(16);
excr_PHY_n_psv = po_vals(63) * y_vals(14);
excr_PHY_p_psv = po_vals(63) * y_vals(13);
% Active carbohydrate excretion
excr_PHY_c_cho = po_vals(64) * growPHYc;
% Active excretion to adjust stoichiometry.
excr_PHY_c_act = 0.5 * y_vals(16) * max([0,...
    1-(y_vals(14)/(y_vals(16) * ps_vals(12))),...
    1-(y_vals(13)/(y_vals(16) * ps_vals(9)))]);
if excr_PHY_c_act > 0
    excr_PHY_n_act = 0.5 * 0.25 * y_vals(14) * max([0,...
        1 - (y_vals(13) / y_vals(14))/(ps_vals(9) / ps_vals(12))]);
    
    excr_PHY_p_act = 0.5 * 0.25 * y_vals(13) * max([0,...
        1 - (y_vals(14) / y_vals(13))/(ps_vals(12) / ps_vals(9))]);
    
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
POM_PHY_c = po_vals(31) * y_vals(16)^2;
POM_PHY_n = (y_vals(14)/y_vals(16)) * POM_PHY_c;
POM_PHY_p = (y_vals(13)/y_vals(16)) * POM_PHY_c;
POM_PHY_chl = (y_vals(15)/y_vals(16)) * POM_PHY_c;
% Grazing
graz_PHY_c = tfun(T) * po_vals(35) * y_vals(11) * y_vals(16)^2 / ...
    (y_vals(16)^2 + po_vals(48)^2 +...%    (y.iUNc*po_vals(48)/po.g_UN)^2 +
    (y_vals(30)*po_vals(48)/po_vals(49))^2);
graz_PHY_n = (y_vals(14)/y_vals(16)) * graz_PHY_c;
graz_PHY_p = (y_vals(13)/y_vals(16)) * graz_PHY_c;
graz_PHY_chl = (y_vals(15)/y_vals(16)) * graz_PHY_c;
% Net growth rate
dydtt(16) = growPHYc - respPHY - excr_PHY_LDOC - excr_PHY_SDOC - POM_PHY_c - graz_PHY_c;
dydtt(14) = growPHYn           - excr_PHY_LDON - excr_PHY_SDON - POM_PHY_n - graz_PHY_n;
dydtt(13) = growPHYpo4           - excr_PHY_LDOP - excr_PHY_SDOP - POM_PHY_p - graz_PHY_p;
dydtt(15) = growPHYchl                                           - POM_PHY_chl - graz_PHY_chl;
%%-----------------------------------------------------------------------
%      Trichdesmium Processes
%-----------------------------------------------------------------------
% The ol' nutrient limitation thing.
Nfunc_TR_n = ((y_vals(3)/y_vals(5))-ps_vals(6))/(ps_vals(4) - ps_vals(6));
Nfunc_TR_p = ((y_vals(2)/y_vals(5))-ps_vals(3))/(ps_vals(1) - ps_vals(3));
temp = max([min([Nfunc_TR_n, Nfunc_TR_p, 1]),0]);
Pmax_TR = po_vals(34) * tfun(T) * temp; % Max specific growth rate.
% Light limitation and specific gross PP (as carbon)
if Pmax_TR == 0
    growTRc = 0;
else
    growTRc = y_vals(5) * Pmax_TR *...
        (1 - exp(-po_vals(70) * (y_vals(4)/y_vals(5)) * PAR / Pmax_TR)) *...
        exp(-po_vals(67) * PAR);
end
% N assimilation, including fixation.
Vmax_TR_n = (ps_vals(5) - (y_vals(3)/y_vals(5)))/(ps_vals(5) - ps_vals(4));
temp = max(min(Vmax_TR_n,1),0); % Yet again forcing this.
growTRnh4 = y_vals(5) * po_vals(4) * tfun(T) * temp * ...
    y_vals(18) / (y_vals(18) + po_vals(41) + (y_vals(17)*po_vals(41)/po_vals(40)));
growTRno3 = y_vals(5) * po_vals(4) * tfun(T) * temp * ...
    y_vals(17) / (y_vals(17) + po_vals(40) + (y_vals(18)*po_vals(40)/po_vals(41)));
% Fixation
if T <= 20
    nfixTRmax = 0;
else
    nfixTRmax = tfun(T) * Vmax_TR_n * ...
        (((growTRc - (po_vals(1)*growTRno3)) * ps_vals(5)) - ...
        growTRno3 - growTRnh4)/...
        (1 + (po_vals(2) * ps_vals(5)));
    nfixTRmax = max([0, nfixTRmax]);
end
% total assimilation of N
growTRn = min(y_vals(5) * po_vals(4) * tfun(T) * Vmax_TR_n,...
    growTRnh4 + growTRno3 + nfixTRmax);
% Actual (not max) fixation rate
growTRnfix = growTRn - growTRnh4 - growTRno3;
% Phosphorus assimilation
Vmax_TR_p = (ps_vals(2) - (y_vals(2)/y_vals(5)))/(ps_vals(2) - ps_vals(1));
temp = max([min([Vmax_TR_p,1]),0]); % Yet again forcing this.
% Assuming inorganic P is the limiter and supply, here's the max growth.
growTRpo4 = y_vals(5) * po_vals(3) * tfun(T) * temp * ...
    y_vals(12) / (y_vals(12) + po_vals(39));
% Picking up PO4 from deep water (kind of unclear how this term makes
% sense)
pickTRpo4 = max([0, po_vals(33) * (y_vals(5) * ((ps_vals(2) + ps_vals(1))/2) - ...
    y_vals(2))]);
% Growth on phosphorus
growTRp = growTRpo4 + pickTRpo4;
% Respiration required for nitrogen metabolism, only this time it also
% includes n fixation.
respTR = (po_vals(1) * growTRno3) + (po_vals(2) * growTRnfix);
% Chlorophyll Production
if PAR==0
    growTRchl = 0;
else
    growTRchl = po_vals(7) * growTRn * growTRc / ...
        (po_vals(70) * y_vals(4) * PAR * exp(-po_vals(67) * PAR));
    if isnan(growTRchl)
        growTRchl = 0;
    end
end
% DOM excretion!
% Passive (fixed proportion)
excr_TR_c_psv = po_vals(58) * y_vals(5);
excr_TR_n_psv = po_vals(58) * y_vals(3);
excr_TR_p_psv = po_vals(58) * y_vals(2);
% Active carbohydrate excretion
excr_TR_c_cho = po_vals(60) * growTRc;
% Active excretion to adjust stoichiometry.
excr_TR_c_act = 0.5 * y_vals(5) * max([0,...
    1-(y_vals(3)/(y_vals(5) * ps_vals(4))),...
    1-(y_vals(2)/(y_vals(5) * ps_vals(1)))]);
if excr_TR_c_act > 0
    excr_TR_n_act = 0.5 * 0.25 * y_vals(3) * max([0,...
        1 - (y_vals(2) / y_vals(3))/(ps_vals(1) / ps_vals(4))]);
    
    excr_TR_p_act = 0.5 * 0.25 * y_vals(2) * max([0,...
        1 - (y_vals(3) / y_vals(2))/(ps_vals(4) / ps_vals(1))]);
    
else
    excr_TR_n_act = 0;
    excr_TR_p_act = 0;
end
% Active release of newly fixed N
excr_TR_nfix_DON = 0.5 * po_vals(59) * growTRnfix * Nfunc_TR_n;
% Partitioning of labile and semilabile DOM
excr_TR_LDOC = excr_TR_c_psv + (0.75 * excr_TR_c_cho);
excr_TR_LDON = excr_TR_n_psv + excr_TR_nfix_DON;
excr_TR_LDOP = excr_TR_p_psv;
excr_TR_SDOC = excr_TR_c_act + (0.25 * excr_TR_c_cho);
excr_TR_SDON = excr_TR_n_act;
excr_TR_SDOP = excr_TR_p_act;
% Ammonium excretion
excr_TR_n_nh4 = 0.5 * po_vals(59) * growTRnfix * Nfunc_TR_n;
% POM production by aggregation
POM_TR_c = po_vals(29) * y_vals(5)^2;
POM_TR_n = (y_vals(3)/y_vals(5)) * POM_TR_c;
POM_TR_p = (y_vals(2)/y_vals(5)) * POM_TR_c;
POM_TR_chl = (y_vals(4)/y_vals(5)) * POM_TR_c;
% Grazing
graz_TR_c = tfun(T) * po_vals(37) * y_vals(21) * y_vals(5)^2 / ...
    (y_vals(5)^2 + po_vals(46)^2 +...
    (y_vals(11)*po_vals(46)/po_vals(47))^2);
graz_TR_n = (y_vals(3)/y_vals(5)) * graz_TR_c;
graz_TR_p = (y_vals(2)/y_vals(5)) * graz_TR_c;
graz_TR_chl = (y_vals(4)/y_vals(5)) * graz_TR_c;
% Net Growth
dydtt(5) = growTRc - respTR - excr_TR_LDOC - excr_TR_SDOC - POM_TR_c - graz_TR_c;
dydtt(3) = growTRn - excr_TR_n_nh4 - excr_TR_LDON - excr_TR_SDON - POM_TR_n - graz_TR_n;
dydtt(2) = growTRp           - excr_TR_LDOP - excr_TR_SDOP - POM_TR_p - graz_TR_p;
dydtt(4) = growTRchl                                           - POM_TR_chl - graz_TR_chl;
% %%-----------------------------------------------------------------------
% %      Unicellular N2-fixer Processes
% %-----------------------------------------------------------------------
% 
% % The ol' nutrient limitation thing.
% Nfunc_UN_n = ((y.iUNn/y.iUNc)-ps.q_UN_n)/(ps.q_UN_n_rdf - ps.q_UN_n);
% Nfunc_UN_p = ((y.iUNp/y.iUNc)-ps.q_UN_p)/(ps.q_UN_p_rdf - ps.q_UN_p);
% 
% temp = max([min([Nfunc_UN_n, Nfunc_UN_p, 1]),0]);
% 
% Pmax_UN = po.mu_UN * tfun(T) * temp; % Max specific growth rate.
% 
% % Light limitation and specific gross PP (as carbon)
% 
% if Pmax_UN == 0
%     growUNc = 0;
% else
%     growUNc = nanmax(y.iUNc * Pmax_UN *...
%         (1 - exp(-po.alpha_UN_chl * (y.iUNchl/y.iUNc) * PAR / Pmax_UN)) *...
%         exp(-po.beta_UN * PAR), 0);
% end
% 
% % N assimilation, including fixation.
% 
% Vmax_UN_n = (ps.q_UN_n_max - (y.iUNn/y.iUNc))/(ps.q_UN_n_max - ps.q_UN_n_rdf);
% temp = max([min([Vmax_UN_n,1]),0]); % Yet again forcing this.
% 
% growUNnh4 = y.iUNc * po.vref_UN_n * tfun(T) * temp * ...
%     y_vals(18) / (y_vals(18) + po.k_UN_nh4 + (y_vals(17)*po.k_UN_nh4/po.k_UN_no3));
% 
% growUNno3 = y.iUNc * po.vref_UN_n * tfun(T) * temp * ...
%     y_vals(17) / (y_vals(17) + po.k_UN_no3 + (y_vals(18)*po.k_UN_no3/po.k_UN_nh4));
% 
% % Fixation
% 
% if T <= 15
%     nfixUNmax = 0;
% else
%     nfixUNmax = tfun(T) * Vmax_UN_n * ...
%         (((growUNc - (po_vals(1)*growUNno3)) * ps.q_UN_n_max) - ...
%         growUNno3 - growUNnh4)/...
%         (1 + (po_vals(2) * ps.q_UN_n_max));
%     nfixUNmax = max([0, nfixUNmax]);
% end
% 
% % total assimilation of N
% 
% growUNn = min([y.iUNc * po.vref_UN_n * tfun(T) * Vmax_UN_n,...
%     growUNnh4 + growUNno3 + nfixUNmax]);
% 
% % Actual (not max) fixation rate
% 
% growUNnfix = growUNn - growUNnh4 - growUNno3;
% 
% % Phosphorus assimilation
% 
% Vmax_UN_p = (ps.q_UN_p_max - (y.iUNp/y.iUNc))/(ps.q_UN_p_max - ps.q_UN_p_rdf);
% temp = max([min([Vmax_UN_p,1]),0]); % Yet again forcing this.
% 
% % Assuming inorganic P is the limiter and supply, here's the max growth.
% 
% growUNp = y.iUNc * po.vref_UN_p * tfun(T) * temp * ...
%     y_vals(12) / (y_vals(12) + po.k_UN_po4);
% 
% % Respiration required for nitrogen metabolism, only this time it also
% % includes n fixation.
% 
% respUN = (po_vals(1) * growUNno3) + (po_vals(2) * growUNnfix);
% 
% % Chlorophyll Production
% 
% if PAR == 0
%     growUNchl=0;
% else
%     growUNchl = po_vals(7) * growUNn * growUNc / ...
%         (po.alpha_UN_chl * y.iUNchl * PAR * exp(-po.beta_UN * PAR));
%     if isnan(growUNchl)
%         growUNchl = 0;
%     end
% end
% 
% % DOM excretion!
% 
% % Passive (fixed proportion)
% 
% excr_UN_c_psv = po.ex_UN_psv * y.iUNc;
% excr_UN_n_psv = po.ex_UN_psv * y.iUNn;
% excr_UN_p_psv = po.ex_UN_psv * y.iUNp;
% 
% % Active carbohydrate excretion
% 
% excr_UN_c_cho = po.ex_UN_cho * growUNc;
% 
% % Active excretion to adjust stoichiometry.
% 
% excr_UN_c_act = 0.5 * y.iUNc * max([0,...
%     1-(y.iUNn/(y.iUNc * ps.q_UN_n_rdf)),...
%     1-(y.iUNp/(y.iUNc * ps.q_UN_p_rdf))]);
% 
% if excr_UN_c_act > 0
%     excr_UN_n_act = 0.5 * 0.25 * y.iUNn * max([0,...
%         1 - (y.iUNp / y.iUNn)/(ps.q_UN_p_rdf / ps.q_UN_n_rdf)]);
%     
%     excr_UN_p_act = 0.5 * 0.25 * y.iUNp * max([0,...
%         1 - (y.iUNn / y.iUNp)/(ps.q_UN_n_rdf / ps.q_UN_p_rdf)]);
%     
% else
%     excr_UN_n_act = 0;
%     excr_UN_p_act = 0;
% end
% 
% % Active release of newly fixed N
% 
% excr_UN_nfix_DON = 0.5 * po.ex_UN_nfix * growUNnfix * Nfunc_UN_n;
% 
% % Partitioning of labile and semilabile DOM
% 
% excr_UN_LDOC = excr_UN_c_psv + (0.75 * excr_UN_c_cho);
% excr_UN_LDON = excr_UN_n_psv + excr_UN_nfix_DON;
% excr_UN_LDOP = excr_UN_p_psv;
% 
% excr_UN_SDOC = excr_UN_c_act + (0.25 * excr_UN_c_cho);
% excr_UN_SDON = excr_UN_n_act;
% excr_UN_SDOP = excr_UN_p_act;
% 
% % Ammonium excretion
% 
% excr_UN_n_nh4 = 0.5 * po.ex_UN_nfix * growUNnfix * Nfunc_UN_n;
% 
% % POM production by aggregation
% 
% POM_UN_c = po.pom_UN * y.iUNc^2;
% POM_UN_n = (y.iUNn/y.iUNc) * POM_UN_c;
% POM_UN_p = (y.iUNp/y.iUNc) * POM_UN_c;
% POM_UN_chl = (y.iUNchl/y.iUNc) * POM_UN_c;
% 
% % Grazing
% 
% graz_UN_c = tfun(T) * po_vals(35) * y_vals(11) * y.iUNc^2 / ...
%     (y.iUNc^2 + po.g_UN^2 +...
%     (y_vals(16)*po.g_UN/po_vals(48))^2 + ...
%     (y_vals(30)*po.g_UN/po_vals(49))^2);
% graz_UN_n = (y.iUNn/y.iUNc) * graz_UN_c;
% graz_UN_p = (y.iUNp/y.iUNc) * graz_UN_c;
% graz_UN_chl = (y.iUNchl/y.iUNc) * graz_UN_c;
% 
% % Net Growth
% 
% dydtt.iUNc = growUNc - respUN - excr_UN_LDOC - excr_UN_SDOC - POM_UN_c - graz_UN_c;
% dydtt.iUNn = growUNn - excr_UN_n_nh4 - excr_UN_LDON - excr_UN_SDON - POM_UN_n - graz_UN_n;
% dydtt.iUNp = growUNp           - excr_UN_LDOP - excr_UN_SDOP - POM_UN_p - graz_UN_p;
% dydtt.iUNchl = growUNchl                                           - POM_UN_chl - graz_UN_chl;
%%-----------------------------------------------------------------------
%      Bacterial Processes
%-----------------------------------------------------------------------
% 1. Gross Grow
% Maximum possible C amount for bacterial use
ALC = y_vals(24);
ASC = y_vals(8) * po_vals(16);
% Carbon Usage
Nfunc_ba_n = y_vals(29)/(y_vals(30)*ps_vals(18)); %bacterial N/C quota /ref quota
Nfunc_ba_p = y_vals(28)/(y_vals(30)*ps_vals(17));
temp = min([Nfunc_ba_n, Nfunc_ba_p]);
temp = min(temp, 1); % Make sure the limit is unity (not the case for mtabs?)
growBAldoc = po_vals(38) * y_vals(30) * tfun(T) * temp * ALC / (ALC+ po_vals(45) + ASC);
growBAsdoc = po_vals(38) * y_vals(30) * tfun(T) * temp * ASC / (ASC+ po_vals(45) + ALC);
%^^^I added the tfun back into this version (it wasn't in the one H. Kim
%sent.)
% DON and DOP usage
growBAldon = growBAldoc*y_vals(23)/y_vals(24); % available labile N
growBAldop = growBAldoc*y_vals(22)/y_vals(24); % available labile P
growBAsdon = growBAsdoc * min([ps_vals(18),...
    y_vals(7)/y_vals(8) + po_vals(54)/Nfunc_ba_n*(ps_vals(18)-(y_vals(7)/y_vals(8)))]);
growBAsdop = growBAsdoc * min([ps_vals(17),...
    y_vals(6)/y_vals(8) + po_vals(54)/Nfunc_ba_p*(ps_vals(17)-(y_vals(6)/y_vals(8)))]);
% inorganic nutrients uptake
growBAnh4 = (growBAldon * y_vals(18) * min(1, 1/Nfunc_ba_n))/ y_vals(23);
if Nfunc_ba_n < 1
    growBAno3 = min([0.1 * ((growBAldon + growBAsdon) / (y_vals(23) + y_vals(7))) * y_vals(17) * min(1, 1/Nfunc_ba_n),...
        ((growBAldon + growBAsdon) / (y_vals(23) + y_vals(7))) * (y_vals(17)  + y_vals(18)) - growBAnh4]);
    growBAno3 = max([0, growBAno3]);
else
    growBAno3 = 0;
end
growBApo4 = (growBAldop / y_vals(22)) * y_vals(12) * min([1, 1/Nfunc_ba_p]);
% Bacteria gross growth
growBAc = growBAldoc + growBAsdoc;
growBAn = growBAldon + growBAsdon + growBAnh4 + growBAno3;
growBAp = growBAldop + growBAsdop + growBApo4;
% 2. respiration
% zeta is an energy required for NO3 use
respBA = po_vals(72) * y_vals(30) + po_vals(1) * growBAno3 + ...
    ((po_vals(17) + (po_vals(18) - po_vals(17))*exp(-po_vals(69)*growBAc) ) *...
    growBAc);
% 3. excreting refractory DOM
refrBAc = po_vals(20) * y_vals(30);
refrBAn = po_vals(24) * refrBAc;
refrBAp = po_vals(23) * refrBAc;
% 4. excreting semi-labile DOM and regenerating DIN
if (y_vals(30) < y_vals(29)/ps_vals(18)) && (y_vals(30) < y_vals(28)/ps_vals(17)) %Carbon in short
    excrBAc = 0;
    excrBAn = 0;
    excrBAp = 0;
    remiBAn = po_vals(19) * (y_vals(29) - y_vals(30) * ps_vals(18));
    remiBAp = po_vals(19) * (y_vals(28) - y_vals(30) * ps_vals(17));
elseif (y_vals(30) > y_vals(29)/ps_vals(18)) && (y_vals(28)/ps_vals(17) > y_vals(29)/ps_vals(18)) %Nitrogen in short
    excrBAc = po_vals(22) * (y_vals(30) - (y_vals(29)/ps_vals(18)));
    excrBAn = 0;
    excrBAp = po_vals(22) * (y_vals(28) - ((y_vals(29)/ps_vals(18)) * ps_vals(17)));
    remiBAn = 0;
    remiBAp = 0;
else %Phosphorus in short
    excrBAc = po_vals(22) * (y_vals(30) - (y_vals(28)/ps_vals(17)));
    excrBAn = po_vals(22) * (y_vals(29) - ((y_vals(28)/ps_vals(17)) * ps_vals(18)));
    excrBAp = 0;
    remiBAn = 0;
    remiBAp = 0;
end
%6. removal by grazing
grazBAc = po_vals(35) * y_vals(11) * y_vals(30)^2 /...
    (y_vals(30)^2 + po_vals(49)^2 + ...%     (y.iUNc*po_vals(49)/po.g_UN)^2) +
    (y_vals(16) * po_vals(49)/po_vals(48))^2);
grazBAn = grazBAc * y_vals(29) / y_vals(30);
grazBAp = grazBAc * y_vals(28) / y_vals(30);
%6b. Mortality due to viruses
mortBAc = po_vals(21) * y_vals(30);
mortBAn = po_vals(21) * y_vals(29);
mortBAp = po_vals(21) * y_vals(28);
%7. BA Derivs
dydtt(30) = (growBAc - refrBAc - excrBAc - grazBAc - respBA - mortBAc); % Originally in seconds and converted by SecPerDay
dydtt(29) = (growBAn - refrBAn - excrBAn - remiBAn - grazBAn - mortBAn);
dydtt(28) = (growBAp - refrBAp - excrBAp - remiBAp - grazBAp - mortBAp);
%8. Flux of inorganic nutrients through bacteria (Formerly part of dydtt
%diagnostics)
fluxBAnh4 = growBAnh4 - remiBAn;
fluxBAno3 = growBAno3;
fluxBApo4 = growBApo4 - remiBAp;
%%-----------------------------------------------------------------------
%      Protozoan Grazers
%-----------------------------------------------------------------------
% Gross growth rate
growPRTc = graz_PHY_c  + grazBAc; %+ graz_UN_c
growPRTn = graz_PHY_n + grazBAn; % + graz_UN_n
growPRTp = graz_PHY_p + grazBAp; % + graz_UN_p
% DOM excretion
excr_PRT_LDOC = po_vals(62) * po_vals(50) * growPRTc;
excr_PRT_LDON = po_vals(62) * po_vals(50) * growPRTn;
excr_PRT_LDOP = po_vals(62) * po_vals(50) * growPRTp;
excr_PRT_SDOC = po_vals(62) * (1 - po_vals(50)) * growPRTc;
excr_PRT_SDON = po_vals(62) * (1 - po_vals(50)) * growPRTn * (y_vals(10)/(y_vals(11)*ps_vals(8)));
excr_PRT_SDOP = po_vals(62) * (1 - po_vals(50)) * growPRTp * (y_vals(9)/(y_vals(11)*ps_vals(7)));
% Respiration
resp_PRT = (po_vals(8) * tfun(T) * y_vals(11)) + (po_vals(10) * growPRTc);
% SDOM release to adjust stoichiometry
excr2_PRT_SDOC = po_vals(56) * y_vals(11) * max([0,...
    1 - (y_vals(10)/(y_vals(11)*ps_vals(8))),...
    1 - (y_vals(9)/(y_vals(11)*ps_vals(7)))]);
excr2_PRT_SDON = 0.5 * excr2_PRT_SDOC * (y_vals(10)/y_vals(11));
excr2_PRT_SDOP = 0.5 * excr2_PRT_SDOC * (y_vals(9)/y_vals(11));
% Inorganic nutrient remin to adjust stoichiometry.
remi_PRT_n = po_vals(13) * max([0, y_vals(10) - (ps_vals(8) * y_vals(11)),...
    y_vals(10) - (ps_vals(8) * y_vals(9) / ps_vals(7))]);
remi_PRT_p = po_vals(13) * max([0, y_vals(9) - (ps_vals(7) * y_vals(11)),...
    y_vals(9) - (ps_vals(7) * y_vals(10) / ps_vals(8))]);
% POM production
POM_PRT_c = po_vals(30) * growPRTc;
POM_PRT_n = po_vals(26) * POM_PRT_c;
POM_PRT_p = po_vals(25) * POM_PRT_c;
% Grazing on protozoa by metazoa.
graz_PRT_c = tfun(T) * po_vals(37) * y_vals(21) * y_vals(11)^2 / ...
    (y_vals(11)^2 + po_vals(47)^2 + (y_vals(5) * po_vals(47) / po_vals(46))^2);
graz_PRT_n = graz_PRT_c * y_vals(10) / y_vals(11);
graz_PRT_p = graz_PRT_c * y_vals(9) / y_vals(11);
% Net Growth
dydtt(11) = growPRTc - resp_PRT - excr_PRT_LDOC - excr_PRT_SDOC - excr2_PRT_SDOC - POM_PRT_c - graz_PRT_c;
dydtt(10) = growPRTn - remi_PRT_n - excr_PRT_LDON - excr_PRT_SDON - excr2_PRT_SDON - POM_PRT_n - graz_PRT_n;
dydtt(9) = growPRTp - remi_PRT_p - excr_PRT_LDOP - excr_PRT_SDOP - excr2_PRT_SDOP - POM_PRT_p - graz_PRT_p;
%%-----------------------------------------------------------------------
%      Metazoan Grazers
%-----------------------------------------------------------------------
% Gross growth rate
growMZc = graz_PRT_c + graz_TR_c;
growMZn = graz_PRT_n + graz_TR_n;
growMZp = graz_PRT_p + graz_TR_p;
% DOM Excretion
excr_MZ_LDOC = po_vals(65) * po_vals(51) * growMZc;
excr_MZ_LDON = po_vals(65) * po_vals(51) * growMZn;
excr_MZ_LDOP = po_vals(65) * po_vals(51) * growMZp;
excr_MZ_SDOC = po_vals(65) * (1 - po_vals(51)) * growMZc;
excr_MZ_SDON = po_vals(65) * (1 - po_vals(51)) * growMZn * (y_vals(20)/(y_vals(21)*ps_vals(16)));
excr_MZ_SDOP = po_vals(65) * (1 - po_vals(51)) * growMZp * (y_vals(19)/(y_vals(21)*ps_vals(15)));
% Respiration
resp_MZ = (po_vals(9) * tfun(T) * y_vals(21)) + (po_vals(11) * growMZc);
% Semilabile DOM stoich adjustment
excr2_MZ_SDOC = po_vals(57) * y_vals(21) * max([0,...
    1 - (y_vals(20)/(y_vals(21)*ps_vals(16))),...
    1 - (y_vals(19)/(y_vals(21)*ps_vals(15)))]);
excr2_MZ_SDON = 0.5 * excr2_MZ_SDOC * (y_vals(20)/y_vals(21));
excr2_MZ_SDOP = 0.5 * excr2_MZ_SDOC * (y_vals(19)/y_vals(21));
% Inorganic nutrient remin to adjust stoichiometry.
remi_MZ_n = po_vals(14) * max([0, y_vals(20) - (ps_vals(16) * y_vals(21)),...
    y_vals(20) - (ps_vals(16) * y_vals(19) / ps_vals(15))]);
remi_MZ_p = po_vals(14) * max([0, y_vals(19) - (ps_vals(15) * y_vals(21)),...
    y_vals(19) - (ps_vals(15) * y_vals(20) / ps_vals(16))]);
% POM production
POM_MZ_c = po_vals(32) * growMZc;
POM_MZ_n = po_vals(26) * POM_MZ_c;
POM_MZ_p = po_vals(25) * POM_MZ_c;
% Refractory DOM release.
REFR_MZ_c = po_vals(55) * y_vals(21);
REFR_MZ_n = REFR_MZ_c * po_vals(24);
REFR_MZ_p = REFR_MZ_c * po_vals(23);
% Removal by higher level zoops
REMV_MZ_c = po_vals(12) * y_vals(21)^2;
REMV_MZ_n = (y_vals(20) / y_vals(21)) * REMV_MZ_c;
REMV_MZ_p = (y_vals(19) / y_vals(21)) * REMV_MZ_c;
% POM production by higher zoops
POM_HZ_c = po_vals(53) * REMV_MZ_c;
POM_HZ_n = po_vals(53) * REMV_MZ_n;
POM_HZ_p = po_vals(53) * REMV_MZ_p;
% SDOM production by higher zoops
excr_HZ_SDOC = po_vals(52) * REMV_MZ_c;
excr_HZ_SDON = po_vals(52) * REMV_MZ_n;
excr_HZ_SDOP = po_vals(52) * REMV_MZ_p;
% Inorganic nuts remineralized by HZ
remi_HZ_n = REMV_MZ_n - POM_HZ_n - excr_HZ_SDON;
remi_HZ_p = REMV_MZ_p - POM_HZ_p - excr_HZ_SDOP;
% Net growth of metazoa.
dydtt(21) = growMZc - resp_MZ - excr_MZ_LDOC - excr_MZ_SDOC - excr2_MZ_SDOC - POM_MZ_c - REFR_MZ_c - REMV_MZ_c;
dydtt(20) = growMZn - remi_MZ_n - excr_MZ_LDON - excr_MZ_SDON - excr2_MZ_SDON - POM_MZ_n - REFR_MZ_n - REMV_MZ_n;
dydtt(19) = growMZp - remi_MZ_p - excr_MZ_LDOP - excr_MZ_SDOP - excr2_MZ_SDOP - POM_MZ_p - REFR_MZ_p - REMV_MZ_p;
%%-----------------------------------------------------------------------
%      Detrital Processes
%-----------------------------------------------------------------------
% Dissolution
DISS_DET_c = po_vals(66) * y_vals(27);
DISS_DET_n = po_vals(28) * po_vals(66) * y_vals(26);
DISS_DET_p = po_vals(27) * po_vals(66) * y_vals(25);
% Net Change
dydtt(27) = POM_PHY_c + POM_TR_c + POM_PRT_c + POM_MZ_c + POM_HZ_c - DISS_DET_c;  %+ POM_UN_c
dydtt(26) = POM_PHY_n + POM_TR_n + POM_PRT_n + POM_MZ_n + POM_HZ_n - DISS_DET_n;  %+ POM_UN_n
dydtt(25) = POM_PHY_p + POM_TR_p + POM_PRT_p + POM_MZ_p + POM_HZ_p - DISS_DET_p; % + POM_UN_p
%%-----------------------------------------------------------------------
%      Inorganic Nutrients
%-----------------------------------------------------------------------
% Nitrification
NTRF = po_vals(15) * y_vals(18);
% Dissolved Nutrient Rates of Change
dydtt(18) = fluxBAnh4 + remi_PRT_n + remi_MZ_n + remi_HZ_n + ...
    excr_TR_n_nh4  -... + excr_UN_n_nh4
    growPHYnh4 - growTRnh4  - NTRF +... %- growUNnh4
    Kz*(forceNH4 - y_vals(18));
dydtt(17) = fluxBAno3 - growPHYno3 - growTRno3 + NTRF +...; % - growUNno3
    Kz*(forceNO3 - y_vals(17));
dydtt(12) = fluxBApo4 + remi_PRT_p + remi_MZ_p + remi_HZ_p -...
    growPHYpo4 - growTRpo4 +...; % - growUNp
    Kz*(forcePO4 - y_vals(12));
%%-----------------------------------------------------------------------
%      Change in Dissolved Organic Matter
%-----------------------------------------------------------------------
% Conversion of SDOM to RDOM
REFR_SDOC = po_vals(61) * y_vals(8) * ...
    exp( 1 - min([(y_vals(7) / (y_vals(8)*po_vals(24))),...
    (y_vals(6) / (y_vals(8)*po_vals(23)))]) ) + ...
    max([ 0 , y_vals(8) - (y_vals(7)/po_vals(24)), y_vals(8) - (y_vals(6)/po_vals(23))]);
REFR_SDON = po_vals(24) * REFR_SDOC;
REFR_SDOP = po_vals(23) * REFR_SDOC;
% Change in DOM species
dydtt(24) = excr_PHY_LDOC + excr_TR_LDOC + excr_PRT_LDOC + ... + excr_UN_LDOC
    excr_MZ_LDOC + mortBAc - growBAldoc;
dydtt(23) = excr_PHY_LDON + excr_TR_LDON + excr_PRT_LDON + ... + excr_UN_LDON
    excr_MZ_LDON + mortBAn - growBAldon;
dydtt(22) = excr_PHY_LDOP + excr_TR_LDOP + excr_PRT_LDOP + ... + excr_UN_LDOP
    excr_MZ_LDOP + mortBAp - growBAldop;
dydtt(8) = excr_PHY_SDOC + excr_TR_SDOC + excrBAc +... + excr_UN_SDOC
    excr_PRT_SDOC + excr2_PRT_SDOC + excr_MZ_SDOC + excr2_MZ_SDOC +...
    excr_HZ_SDOC + DISS_DET_c - REFR_SDOC - growBAsdoc;
dydtt(7) = excr_PHY_SDON + excr_TR_SDON + excrBAn +... + excr_UN_SDON
    excr_PRT_SDON + excr2_PRT_SDON + excr_MZ_SDON + excr2_MZ_SDON +...
    excr_HZ_SDON + DISS_DET_n - REFR_SDON - growBAsdon;
dydtt(6) = excr_PHY_SDOP + excr_TR_SDOP + excrBAp +... + excr_UN_SDOP
    excr_PRT_SDOP + excr2_PRT_SDOP + excr_MZ_SDOP + excr2_MZ_SDOP +...
    excr_HZ_SDOP + DISS_DET_p - REFR_SDOP - growBAsdop;
%%-----------------------------------------------------------------------
%      Diagnostic Variables
%-----------------------------------------------------------------------
diaginput = [growPHYc, growTRc, -respPHY, -respTR, -excr_PHY_LDOC,...
    -excr_TR_LDOC, -excr_PHY_SDOC, -excr_TR_SDOC, growBAc, respBA];
[dydtt(1),dydtt(31)] = ppbp(diaginput);
dydtt = dydtt';
% if sum(isnan(dydtt))~=0
%     disp('NaN at ' + string(find(isnan(dydtt))) +...
%         ', time ' +string(ts)+' days.')
%     dydtt(isnan(dydtt))=0;
%     return
% end
end