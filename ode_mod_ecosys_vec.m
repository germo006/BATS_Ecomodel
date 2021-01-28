%% Noah Germolus, 13 Jan 2021
% This is the system of differential equations in Luo et al. (2010). It
% might actually be the play to just combine all 11 sections of the
% equations in one place, actually.
function [dydtt] = ode_mod_ecosys_vec(y_vals, ps_vals, po_vals, ts)
%global y
%%-----------------------------------------------------------------------
%      Temperature Effects and PAR
%-----------------------------------------------------------------------
tfun = @(T) exp(-po_vals(86)*((1/(T+273.15))-(1/(25+273.15))));
pfun = @(t) max(0, -10*cos(2*pi*t));
y_vals(35) = pfun(ts);
%%-----------------------------------------------------------------------
%      Phytoplankton Processes
%-----------------------------------------------------------------------
% 1. Gross Grow
% Determining growth rate limits
Nfunc_phy_n = ((y_vals(17)/y_vals(19))-ps_vals(20))/(ps_vals(18) - ps_vals(20));
Nfunc_phy_p = ((y_vals(16)/y_vals(19))-ps_vals(17))/(ps_vals(15) - ps_vals(17));
temp = max([min([Nfunc_phy_n, Nfunc_phy_p, 1]),0]); % Force elemental limits
% on rate to between 0 and 100%
Pmax_phy = po_vals(40) * tfun(y_vals(34)) * temp; % Max specific growth rate.
% Light limitation and specific gross PP (as carbon)
growPHYc = y_vals(19) * Pmax_phy *...
    (1 - exp(-po_vals(84) * (y_vals(18)/y_vals(19)) * y_vals(35) / Pmax_phy)) *...
    exp(-po_vals(80) * y_vals(35));
% Nitrogen assimilation
Vmax_phy_n = (ps_vals(19) - (y_vals(17)/y_vals(19)))/(ps_vals(19) - ps_vals(18));
temp = max([min([Vmax_phy_n,1]),0]); % Yet again forcing this.
growPHYnh4 = y_vals(19) * po_vals(8) * tfun(y_vals(34)) * temp * ...
    y_vals(21) / (y_vals(21) + po_vals(51) + (y_vals(20)*po_vals(51)/po_vals(50)));
growPHYno3 = y_vals(19) * po_vals(8) * tfun(y_vals(34)) * temp * ...
    y_vals(20) / (y_vals(20) + po_vals(50) + (y_vals(21)*po_vals(50)/po_vals(51)));
growPHYn = growPHYnh4 + growPHYno3; % Total growth on nitrogen limit.
% Phosphorus limitation.
Vmax_phy_p = (ps_vals(16) - (y_vals(16)/y_vals(19)))/(ps_vals(16) - ps_vals(15));
temp = max([min([Vmax_phy_p,1]),0]); % Yet again forcing this.
% Assuming inorganic P is the limiter and supply, here's the max growth.
growPHYpo4 = y_vals(19) * po_vals(7) * tfun(y_vals(34)) * temp * ...
    y_vals(15) / (y_vals(15) + po_vals(49));
% The use of nitrate requires reduction, which requires energy. Here we
% evaluate the respiration necessary for this.
respPHY = po_vals(1) * growPHYno3;
% Chlorophyll Production
growPHYchl = po_vals(9) * growPHYn * growPHYc / ...
    (po_vals(84) * y_vals(18) * y_vals(35) * exp(-po_vals(80) * y_vals(35)));
% DOM excretion!
% Passive (fixed proportion)
excr_PHY_c_psv = po_vals(74) * y_vals(19);
excr_PHY_n_psv = po_vals(74) * y_vals(17);
excr_PHY_p_psv = po_vals(74) * y_vals(16);
% Active carbohydrate excretion
excr_PHY_c_cho = po_vals(75) * growPHYc;
% Active excretion to adjust stoichiometry.
excr_PHY_c_act = 0.5 * y_vals(19) * max([0,...
    1-(y_vals(17)/(y_vals(19) * ps_vals(18))),...
    1-(y_vals(16)/(y_vals(19) * ps_vals(15)))]);
if excr_PHY_c_act > 0
    excr_PHY_n_act = 0.5 * 0.25 * y_vals(17) * max([0,...
        1 - (y_vals(16) / y_vals(17))/(ps_vals(15) / ps_vals(18))]);
    
    excr_PHY_p_act = 0.5 * 0.25 * y_vals(16) * max([0,...
        1 - (y_vals(17) / y_vals(16))/(ps_vals(18) / ps_vals(15))]);
    
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
POM_PHY_c = po_vals(34) * y_vals(19)^2;
POM_PHY_n = (y_vals(17)/y_vals(19)) * POM_PHY_c;
POM_PHY_p = (y_vals(16)/y_vals(19)) * POM_PHY_c;
POM_PHY_chl = (y_vals(18)/y_vals(19)) * POM_PHY_c;
% Grazing
graz_PHY_c = tfun(y_vals(34)) * po_vals(39) * y_vals(14) * y_vals(19)^2 / ...
    (y_vals(19)^2 + po_vals(56)^2 +...
    (y_vals(4)*po_vals(56)/po_vals(53))^2 +...
    (y_vals(33)*po_vals(56)/po_vals(57))^2);
graz_PHY_n = (y_vals(17)/y_vals(19)) * graz_PHY_c;
graz_PHY_p = (y_vals(16)/y_vals(19)) * graz_PHY_c;
graz_PHY_chl = (y_vals(18)/y_vals(19)) * graz_PHY_c;
% Net growth rate
dydtt(19) = growPHYc - respPHY - excr_PHY_LDOC - excr_PHY_SDOC - POM_PHY_c - graz_PHY_c;
dydtt(17) = growPHYn           - excr_PHY_LDON - excr_PHY_SDON - POM_PHY_n - graz_PHY_n;
dydtt(16) = growPHYpo4           - excr_PHY_LDOP - excr_PHY_SDOP - POM_PHY_p - graz_PHY_p;
dydtt(18) = growPHYchl                                           - POM_PHY_chl - graz_PHY_chl;
%%-----------------------------------------------------------------------
%      Trichdesmium Processes
%-----------------------------------------------------------------------
% The ol' nutrient limitation thing.
Nfunc_TR_n = ((y_vals(6)/y_vals(8))-ps_vals(12))/(ps_vals(10) - ps_vals(12));
Nfunc_TR_p = ((y_vals(5)/y_vals(8))-ps_vals(9))/(ps_vals(7) - ps_vals(9));
temp = max([min([Nfunc_TR_n, Nfunc_TR_p, 1]),0]);
Pmax_TR = po_vals(38) * tfun(y_vals(34)) * temp; % Max specific growth rate.
% Light limitation and specific gross PP (as carbon)
if Pmax_TR == 0
    growTRc = 0;
else
    growTRc = y_vals(8) * Pmax_TR *...
        (1 - exp(-po_vals(83) * (y_vals(7)/y_vals(8)) * y_vals(35) / Pmax_TR)) *...
        exp(-po_vals(79) * y_vals(35));
end
% N assimilation, including fixation.
Vmax_TR_n = (ps_vals(11) - (y_vals(6)/y_vals(8)))/(ps_vals(11) - ps_vals(10));
temp = max(min(Vmax_TR_n,1),0); % Yet again forcing this.
growTRnh4 = y_vals(8) * po_vals(6) * tfun(y_vals(34)) * temp * ...
    y_vals(21) / (y_vals(21) + po_vals(48) + (y_vals(20)*po_vals(48)/po_vals(47)));
growTRno3 = y_vals(8) * po_vals(6) * tfun(y_vals(34)) * temp * ...
    y_vals(20) / (y_vals(20) + po_vals(47) + (y_vals(21)*po_vals(47)/po_vals(48)));
% Fixation
if y_vals(34) <= 20
    nfixTRmax = 0;
else
    nfixTRmax = tfun(y_vals(34)) * Vmax_TR_n * ...
        (((growTRc - (po_vals(1)*growTRno3)) * ps_vals(11)) - ...
        growTRno3 - growTRnh4)/...
        (1 + (po_vals(2) * ps_vals(11)));
    nfixTRmax = max([0, nfixTRmax]);
end
% total assimilation of N
growTRn = min(y_vals(8) * po_vals(6) * tfun(y_vals(34)) * Vmax_TR_n,...
    growTRnh4 + growTRno3 + nfixTRmax);
% Actual (not max) fixation rate
growTRnfix = growTRn - growTRnh4 - growTRno3;
% Phosphorus assimilation
Vmax_TR_p = (ps_vals(8) - (y_vals(5)/y_vals(8)))/(ps_vals(8) - ps_vals(7));
temp = max([min([Vmax_TR_p,1]),0]); % Yet again forcing this.
% Assuming inorganic P is the limiter and supply, here's the max growth.
growTRpo4 = y_vals(8) * po_vals(5) * tfun(y_vals(34)) * temp * ...
    y_vals(15) / (y_vals(15) + po_vals(46));
% Picking up PO4 from deep water (kind of unclear how this term makes
% sense)
pickTRpo4 = max([0, po_vals(36) * (y_vals(8) * ((ps_vals(8) + ps_vals(7))/2) - ...
    y_vals(5))]);
% Growth on phosphorus
growTRp = growTRpo4 + pickTRpo4;
% Respiration required for nitrogen metabolism, only this time it also
% includes n fixation.
respTR = (po_vals(1) * growTRno3) + (po_vals(2) * growTRnfix);
% Chlorophyll Production
if y_vals(35)==0
    growTRchl = 0;
else
    growTRchl = po_vals(9) * growTRn * growTRc / ...
        (po_vals(83) * y_vals(7) * y_vals(35) * exp(-po_vals(79) * y_vals(35)));
end
% DOM excretion!
% Passive (fixed proportion)
excr_TR_c_psv = po_vals(69) * y_vals(8);
excr_TR_n_psv = po_vals(69) * y_vals(6);
excr_TR_p_psv = po_vals(69) * y_vals(5);
% Active carbohydrate excretion
excr_TR_c_cho = po_vals(71) * growTRc;
% Active excretion to adjust stoichiometry.
excr_TR_c_act = 0.5 * y_vals(8) * max([0,...
    1-(y_vals(6)/(y_vals(8) * ps_vals(10))),...
    1-(y_vals(5)/(y_vals(8) * ps_vals(7)))]);
if excr_TR_c_act > 0
    excr_TR_n_act = 0.5 * 0.25 * y_vals(6) * max([0,...
        1 - (y_vals(5) / y_vals(6))/(ps_vals(7) / ps_vals(10))]);
    
    excr_TR_p_act = 0.5 * 0.25 * y_vals(5) * max([0,...
        1 - (y_vals(6) / y_vals(5))/(ps_vals(10) / ps_vals(7))]);
    
else
    excr_TR_n_act = 0;
    excr_TR_p_act = 0;
end
% Active release of newly fixed N
excr_TR_nfix_DON = 0.5 * po_vals(70) * growTRnfix * Nfunc_TR_n;
% Partitioning of labile and semilabile DOM
excr_TR_LDOC = excr_TR_c_psv + (0.75 * excr_TR_c_cho);
excr_TR_LDON = excr_TR_n_psv + excr_TR_nfix_DON;
excr_TR_LDOP = excr_TR_p_psv;
excr_TR_SDOC = excr_TR_c_act + (0.25 * excr_TR_c_cho);
excr_TR_SDON = excr_TR_n_act;
excr_TR_SDOP = excr_TR_p_act;
% Ammonium excretion
excr_TR_n_nh4 = 0.5 * po_vals(70) * growTRnfix * Nfunc_TR_n;
% POM production by aggregation
POM_TR_c = po_vals(32) * y_vals(8)^2;
POM_TR_n = (y_vals(6)/y_vals(8)) * POM_TR_c;
POM_TR_p = (y_vals(5)/y_vals(8)) * POM_TR_c;
POM_TR_chl = (y_vals(7)/y_vals(8)) * POM_TR_c;
% Grazing
graz_TR_c = tfun(y_vals(34)) * po_vals(41) * y_vals(24) * y_vals(8)^2 / ...
    (y_vals(8)^2 + po_vals(54)^2 +...
    (y_vals(14)*po_vals(54)/po_vals(55))^2);
graz_TR_n = (y_vals(6)/y_vals(8)) * graz_TR_c;
graz_TR_p = (y_vals(5)/y_vals(8)) * graz_TR_c;
graz_TR_chl = (y_vals(7)/y_vals(8)) * graz_TR_c;
% Net Growth
dydtt(8) = growTRc - respTR - excr_TR_LDOC - excr_TR_SDOC - POM_TR_c - graz_TR_c;
dydtt(6) = growTRn - excr_TR_n_nh4 - excr_TR_LDON - excr_TR_SDON - POM_TR_n - graz_TR_n;
dydtt(5) = growTRp           - excr_TR_LDOP - excr_TR_SDOP - POM_TR_p - graz_TR_p;
dydtt(7) = growTRchl                                           - POM_TR_chl - graz_TR_chl;
%%-----------------------------------------------------------------------
%      Unicellular N2-fixer Processes
%-----------------------------------------------------------------------
% The ol' nutrient limitation thing.
Nfunc_UN_n = ((y_vals(2)/y_vals(4))-ps_vals(6))/(ps_vals(4) - ps_vals(6));
Nfunc_UN_p = ((y_vals(1)/y_vals(4))-ps_vals(3))/(ps_vals(1) - ps_vals(3));
temp = max([min([Nfunc_UN_n, Nfunc_UN_p, 1]),0]);
Pmax_UN = po_vals(37) * tfun(y_vals(34)) * temp; % Max specific growth rate.
% Light limitation and specific gross PP (as carbon)
if Pmax_UN == 0
    growUNc = 0;
else
    growUNc = y_vals(4) * Pmax_UN *...
        (1 - exp(-po_vals(82) * (y_vals(3)/y_vals(4)) * y_vals(35) / Pmax_UN)) *...
        exp(-po_vals(78) * y_vals(35));
end
% N assimilation, including fixation.
Vmax_UN_n = (ps_vals(5) - (y_vals(2)/y_vals(4)))/(ps_vals(5) - ps_vals(4));
temp = max([min([Vmax_UN_n,1]),0]); % Yet again forcing this.
growUNnh4 = y_vals(4) * po_vals(4) * tfun(y_vals(34)) * temp * ...
    y_vals(21) / (y_vals(21) + po_vals(45) + (y_vals(20)*po_vals(45)/po_vals(44)));
growUNno3 = y_vals(4) * po_vals(4) * tfun(y_vals(34)) * temp * ...
    y_vals(20) / (y_vals(20) + po_vals(44) + (y_vals(21)*po_vals(44)/po_vals(45)));
% Fixation
if y_vals(34) <= 15
    nfixUNmax = 0;
else
    nfixUNmax = tfun(y_vals(34)) * Vmax_UN_n * ...
        (((growUNc - (po_vals(1)*growUNno3)) * ps_vals(5)) - ...
        growUNno3 - growUNnh4)/...
        (1 + (po_vals(2) * ps_vals(5)));
    nfixUNmax = max([0, nfixUNmax]);
end
% total assimilation of N
growUNn = min([y_vals(4) * po_vals(4) * tfun(y_vals(34)) * Vmax_UN_n,...
    growUNnh4 + growUNno3 + nfixUNmax]);
% Actual (not max) fixation rate
growUNnfix = growUNn - growUNnh4 - growUNno3;
% Phosphorus assimilation
Vmax_UN_p = (ps_vals(2) - (y_vals(1)/y_vals(4)))/(ps_vals(2) - ps_vals(1));
temp = max([min([Vmax_UN_p,1]),0]); % Yet again forcing this.
% Assuming inorganic P is the limiter and supply, here's the max growth.
growUNp = y_vals(4) * po_vals(3) * tfun(y_vals(34)) * temp * ...
    y_vals(15) / (y_vals(15) + po_vals(43));
% Respiration required for nitrogen metabolism, only this time it also
% includes n fixation.
respUN = (po_vals(1) * growUNno3) + (po_vals(2) * growUNnfix);
% Chlorophyll Production
if y_vals(35) == 0
    growUNchl=0;
else
    growUNchl = po_vals(9) * growUNn * growUNc / ...
        (po_vals(82) * y_vals(3) * y_vals(35) * exp(-po_vals(78) * y_vals(35)));
end
% DOM excretion!
% Passive (fixed proportion)
excr_UN_c_psv = po_vals(66) * y_vals(4);
excr_UN_n_psv = po_vals(66) * y_vals(2);
excr_UN_p_psv = po_vals(66) * y_vals(1);
% Active carbohydrate excretion
excr_UN_c_cho = po_vals(68) * growUNc;
% Active excretion to adjust stoichiometry.
excr_UN_c_act = 0.5 * y_vals(4) * max([0,...
    1-(y_vals(2)/(y_vals(4) * ps_vals(4))),...
    1-(y_vals(1)/(y_vals(4) * ps_vals(1)))]);
if excr_UN_c_act > 0
    excr_UN_n_act = 0.5 * 0.25 * y_vals(2) * max([0,...
        1 - (y_vals(1) / y_vals(2))/(ps_vals(1) / ps_vals(4))]);
    
    excr_UN_p_act = 0.5 * 0.25 * y_vals(1) * max([0,...
        1 - (y_vals(2) / y_vals(1))/(ps_vals(4) / ps_vals(1))]);
    
else
    excr_UN_n_act = 0;
    excr_UN_p_act = 0;
end
% Active release of newly fixed N
excr_UN_nfix_DON = 0.5 * po_vals(67) * growUNnfix * Nfunc_UN_n;
% Partitioning of labile and semilabile DOM
excr_UN_LDOC = excr_UN_c_psv + (0.75 * excr_UN_c_cho);
excr_UN_LDON = excr_UN_n_psv + excr_UN_nfix_DON;
excr_UN_LDOP = excr_UN_p_psv;
excr_UN_SDOC = excr_UN_c_act + (0.25 * excr_UN_c_cho);
excr_UN_SDON = excr_UN_n_act;
excr_UN_SDOP = excr_UN_p_act;
% Ammonium excretion
excr_UN_n_nh4 = 0.5 * po_vals(67) * growUNnfix * Nfunc_UN_n;
% POM production by aggregation
POM_UN_c = po_vals(31) * y_vals(4)^2;
POM_UN_n = (y_vals(2)/y_vals(4)) * POM_UN_c;
POM_UN_p = (y_vals(1)/y_vals(4)) * POM_UN_c;
POM_UN_chl = (y_vals(3)/y_vals(4)) * POM_UN_c;
% Grazing
graz_UN_c = tfun(y_vals(34)) * po_vals(39) * y_vals(14) * y_vals(4)^2 / ...
    (y_vals(4)^2 + po_vals(53)^2 +...
    (y_vals(19)*po_vals(53)/po_vals(56))^2 + ...
    (y_vals(33)*po_vals(53)/po_vals(57))^2);
graz_UN_n = (y_vals(2)/y_vals(4)) * graz_UN_c;
graz_UN_p = (y_vals(1)/y_vals(4)) * graz_UN_c;
graz_UN_chl = (y_vals(3)/y_vals(4)) * graz_UN_c;
% Net Growth
dydtt(4) = growUNc - respUN - excr_UN_LDOC - excr_UN_SDOC - POM_UN_c - graz_UN_c;
dydtt(2) = growUNn - excr_UN_n_nh4 - excr_UN_LDON - excr_UN_SDON - POM_UN_n - graz_UN_n;
dydtt(1) = growUNp           - excr_UN_LDOP - excr_UN_SDOP - POM_UN_p - graz_UN_p;
dydtt(3) = growUNchl                                           - POM_UN_chl - graz_UN_chl;
%%-----------------------------------------------------------------------
%      Bacterial Processes
%-----------------------------------------------------------------------
% 1. Gross Grow
% Maximum possible C amount for bacterial use
ALC = y_vals(27);
ASC = y_vals(11) * po_vals(18);
% Carbon Usage
Nfunc_ba_n = y_vals(32)/(y_vals(33)*ps_vals(24)); %bacterial N/C quota /ref quota
Nfunc_ba_p = y_vals(31)/(y_vals(33)*ps_vals(23));
temp = min([Nfunc_ba_n, Nfunc_ba_p]);
temp = min(temp, 1); % Make sure the limit is unity (not the case for mtabs?)
growBAldoc = po_vals(42) * y_vals(33) * tfun(y_vals(34)) * temp * ALC / (ALC+ po_vals(52) + ASC);
growBAsdoc = po_vals(42) * y_vals(33) * tfun(y_vals(34)) * temp * ASC / (ASC+ po_vals(52) + ALC);
%^^^I added the tfun back into this version (it wasn't in the one H. Kim
%sent.)
% DON and DOP usage
growBAldon = growBAldoc*y_vals(26)/y_vals(27); % available labile N
growBAldop = growBAldoc*y_vals(25)/y_vals(27); % available labile P
growBAsdon = growBAsdoc * min([ps_vals(24),...
    y_vals(10)/y_vals(11) + po_vals(62)/Nfunc_ba_n*(ps_vals(24)-(y_vals(10)/y_vals(11)))]);
growBAsdop = growBAsdoc * min([ps_vals(23),...
    y_vals(9)/y_vals(11) + po_vals(62)/Nfunc_ba_p*(ps_vals(23)-(y_vals(9)/y_vals(11)))]);
% inorganic nutrients uptake
growBAnh4 = (growBAldon * y_vals(21) * min(1, 1/Nfunc_ba_n))/ y_vals(26);
if Nfunc_ba_n < 1
    growBAno3 = min([0.1 * ((growBAldon + growBAsdon) / (y_vals(26) + y_vals(10))) * y_vals(20) * min(1, 1/Nfunc_ba_n),...
        ((growBAldon + growBAsdon) / (y_vals(26) + y_vals(10))) * (y_vals(20)  + y_vals(21)) - growBAnh4]);
    growBAno3 = max([0, growBAno3]);
else
    growBAno3 = 0;
end
growBApo4 = (growBAldop / y_vals(25)) * y_vals(15) * min([1, 1/Nfunc_ba_p]);
% Bacteria gross growth
growBAc = growBAldoc + growBAsdoc;
growBAn = growBAldon + growBAsdon + growBAnh4 + growBAno3;
growBAp = growBAldop + growBAsdop + growBApo4;
% 2. respiration
% zeta is an energy required for NO3 use
respBA = po_vals(85) * y_vals(33) + po_vals(1) * growBAno3 + ...
    ((po_vals(19) + (po_vals(20) - po_vals(19))*exp(-po_vals(81)*growBAc) ) *...
    growBAc);
% 3. excreting refractory DOM
refrBAc = po_vals(22) * y_vals(33);
refrBAn = po_vals(26) * refrBAc;
refrBAp = po_vals(25) * refrBAc;
% 4. excreting semi-labile DOM and regenerating DIN
if (y_vals(33) < y_vals(32)/ps_vals(24)) && (y_vals(33) < y_vals(31)/ps_vals(23)) %Carbon in short
    excrBAc = 0;
    excrBAn = 0;
    excrBAp = 0;
    remiBAn = po_vals(21) * (y_vals(32) - y_vals(33) * ps_vals(24));
    remiBAp = po_vals(21) * (y_vals(31) - y_vals(33) * ps_vals(23));
elseif (y_vals(33) > y_vals(32)/ps_vals(24)) && (y_vals(31)/ps_vals(23) > y_vals(32)/ps_vals(24)) %Nitrogen in short
    excrBAc = po_vals(24) * (y_vals(33) - (y_vals(32)/ps_vals(24)));
    excrBAn = 0;
    excrBAp = po_vals(24) * (y_vals(31) - ((y_vals(32)/ps_vals(24)) * ps_vals(23)));
    remiBAn = 0;
    remiBAp = 0;
else %Phosphorus in short
    excrBAc = po_vals(24) * (y_vals(33) - (y_vals(31)/ps_vals(23)));
    excrBAn = po_vals(24) * (y_vals(32) - ((y_vals(31)/ps_vals(23)) * ps_vals(24)));
    excrBAp = 0;
    remiBAn = 0;
    remiBAp = 0;
end
%6. removal by grazing
grazBAc = po_vals(39) * y_vals(14) * y_vals(33)^2 /...
    (y_vals(33)^2 + po_vals(57)^2 + ...
    (y_vals(19) * po_vals(57)/po_vals(56))^2 + ...
    (y_vals(4)*po_vals(57)/po_vals(53))^2);
grazBAn = grazBAc * y_vals(32) / y_vals(33);
grazBAp = grazBAc * y_vals(31) / y_vals(33);
%6b. Mortality due to viruses
mortBAc = po_vals(23) * y_vals(33);
mortBAn = po_vals(23) * y_vals(32);
mortBAp = po_vals(23) * y_vals(31);
%7. BA Derivs
dydtt(33) = (growBAc - refrBAc - excrBAc - grazBAc - respBA - mortBAc); % Originally in seconds and converted by SecPerDay
dydtt(32) = (growBAn - refrBAn - excrBAn - remiBAn - grazBAn - mortBAn);
dydtt(31) = (growBAp - refrBAp - excrBAp - remiBAp - grazBAp - mortBAp);
%8. Flux of inorganic nutrients through bacteria
dydtt(36) = growBAnh4 - remiBAn;
dydtt(35) = growBAno3;
dydtt(34) = growBApo4 - remiBAp;
%%-----------------------------------------------------------------------
%      Protozoan Grazers
%-----------------------------------------------------------------------
% Gross growth rate
growPRTc = graz_PHY_c + graz_UN_c + grazBAc;
growPRTn = graz_PHY_n + graz_UN_n + grazBAn;
growPRTp = graz_PHY_p + graz_UN_p + grazBAp;
% DOM excretion
excr_PRT_LDOC = po_vals(73) * po_vals(58) * growPRTc;
excr_PRT_LDON = po_vals(73) * po_vals(58) * growPRTn;
excr_PRT_LDOP = po_vals(73) * po_vals(58) * growPRTp;
excr_PRT_SDOC = po_vals(73) * (1 - po_vals(58)) * growPRTc;
excr_PRT_SDON = po_vals(73) * (1 - po_vals(58)) * growPRTn * (y_vals(13)/(y_vals(14)*ps_vals(14)));
excr_PRT_SDOP = po_vals(73) * (1 - po_vals(58)) * growPRTp * (y_vals(12)/(y_vals(14)*ps_vals(13)));
% Respiration
resp_PRT = (po_vals(10) * tfun(y_vals(34)) * y_vals(14)) + (po_vals(12) * growPRTc);
% SDOM release to adjust stoichiometry
excr2_PRT_SDOC = po_vals(64) * y_vals(14) * max([0,...
    1 - (y_vals(13)/(y_vals(14)*ps_vals(14))),...
    1 - (y_vals(12)/(y_vals(14)*ps_vals(13)))]);
excr2_PRT_SDON = 0.5 * excr2_PRT_SDOC * (y_vals(13)/y_vals(14));
excr2_PRT_SDOP = 0.5 * excr2_PRT_SDOC * (y_vals(12)/y_vals(14));
% Inorganic nutrient remin to adjust stoichiometry.
remi_PRT_n = po_vals(15) * max([0, y_vals(13) - (ps_vals(14) * y_vals(14)),...
    y_vals(13) - (ps_vals(14) * y_vals(12) / ps_vals(13))]);
remi_PRT_p = po_vals(15) * max([0, y_vals(12) - (ps_vals(13) * y_vals(14)),...
    y_vals(12) - (ps_vals(13) * y_vals(13) / ps_vals(14))]);
% POM production
POM_PRT_c = po_vals(33) * growPRTc;
POM_PRT_n = po_vals(28) * POM_PRT_c;
POM_PRT_p = po_vals(27) * POM_PRT_c;
% Grazing on protozoa by metazoa.
graz_PRT_c = tfun(y_vals(34)) * po_vals(41) * y_vals(24) * y_vals(14)^2 / ...
    (y_vals(14)^2 + po_vals(55)^2 + (y_vals(8) * po_vals(55) / po_vals(54))^2);
graz_PRT_n = graz_PRT_c * y_vals(13) / y_vals(14);
graz_PRT_p = graz_PRT_c * y_vals(12) / y_vals(14);
% Net Growth
dydtt(14) = growPRTc - resp_PRT - excr_PRT_LDOC - excr_PRT_SDOC - excr2_PRT_SDOC - POM_PRT_c - graz_PRT_c;
dydtt(13) = growPRTn - remi_PRT_n - excr_PRT_LDON - excr_PRT_SDON - excr2_PRT_SDON - POM_PRT_n - graz_PRT_n;
dydtt(12) = growPRTp - remi_PRT_p - excr_PRT_LDOP - excr_PRT_SDOP - excr2_PRT_SDOP - POM_PRT_p - graz_PRT_p;
%%-----------------------------------------------------------------------
%      Metazoan Grazers
%-----------------------------------------------------------------------
% Gross growth rate
growMZc = graz_PRT_c + graz_TR_c;
growMZn = graz_PRT_n + graz_TR_n;
growMZp = graz_PRT_p + graz_TR_p;
% DOM Excretion
excr_MZ_LDOC = po_vals(76) * po_vals(59) * growMZc;
excr_MZ_LDON = po_vals(76) * po_vals(59) * growMZn;
excr_MZ_LDOP = po_vals(76) * po_vals(59) * growMZp;
excr_MZ_SDOC = po_vals(76) * (1 - po_vals(59)) * growMZc;
excr_MZ_SDON = po_vals(76) * (1 - po_vals(59)) * growMZn * (y_vals(23)/(y_vals(24)*ps_vals(22)));
excr_MZ_SDOP = po_vals(76) * (1 - po_vals(59)) * growMZp * (y_vals(22)/(y_vals(24)*ps_vals(21)));
% Respiration
resp_MZ = (po_vals(11) * tfun(y_vals(34)) * y_vals(24)) + (po_vals(13) * growMZc);
% Semilabile DOM stoich adjustment
excr2_MZ_SDOC = po_vals(65) * y_vals(24) * max([0,...
    1 - (y_vals(23)/(y_vals(24)*ps_vals(22))),...
    1 - (y_vals(22)/(y_vals(24)*ps_vals(21)))]);
excr2_MZ_SDON = 0.5 * excr2_MZ_SDOC * (y_vals(23)/y_vals(24));
excr2_MZ_SDOP = 0.5 * excr2_MZ_SDOC * (y_vals(22)/y_vals(24));
% Inorganic nutrient remin to adjust stoichiometry.
remi_MZ_n = po_vals(16) * max([0, y_vals(23) - (ps_vals(22) * y_vals(24)),...
    y_vals(23) - (ps_vals(22) * y_vals(22) / ps_vals(21))]);
remi_MZ_p = po_vals(16) * max([0, y_vals(22) - (ps_vals(21) * y_vals(24)),...
    y_vals(22) - (ps_vals(21) * y_vals(23) / ps_vals(22))]);
% POM production
POM_MZ_c = po_vals(35) * growMZc;
POM_MZ_n = po_vals(28) * POM_MZ_c;
POM_MZ_p = po_vals(27) * POM_MZ_c;
% Refractory DOM release.
REFR_MZ_c = po_vals(63) * y_vals(24);
REFR_MZ_n = REFR_MZ_c * po_vals(26);
REFR_MZ_p = REFR_MZ_c * po_vals(25);
% Removal by higher level zoops
REMV_MZ_c = po_vals(14) * y_vals(24)^2;
REMV_MZ_n = (y_vals(23) / y_vals(24)) * REMV_MZ_c;
REMV_MZ_p = (y_vals(22) / y_vals(24)) * REMV_MZ_c;
% POM production by higher zoops
POM_HZ_c = po_vals(61) * REMV_MZ_c;
POM_HZ_n = po_vals(61) * REMV_MZ_n;
POM_HZ_p = po_vals(61) * REMV_MZ_p;
% SDOM production by higher zoops
excr_HZ_SDOC = po_vals(60) * REMV_MZ_c;
excr_HZ_SDON = po_vals(60) * REMV_MZ_n;
excr_HZ_SDOP = po_vals(60) * REMV_MZ_p;
% Inorganic nuts remineralized by HZ
remi_HZ_n = REMV_MZ_n - POM_HZ_n - excr_HZ_SDON;
remi_HZ_p = REMV_MZ_p - POM_HZ_p - excr_HZ_SDOP;
% Net growth of metazoa.
dydtt(24) = growMZc - resp_MZ - excr_MZ_LDOC - excr_MZ_SDOC - excr2_MZ_SDOC - POM_MZ_c - REFR_MZ_c - REMV_MZ_c;
dydtt(23) = growMZn - remi_MZ_n - excr_MZ_LDON - excr_MZ_SDON - excr2_MZ_SDON - POM_MZ_n - REFR_MZ_n - REMV_MZ_n;
dydtt(22) = growMZp - remi_MZ_p - excr_MZ_LDOP - excr_MZ_SDOP - excr2_MZ_SDOP - POM_MZ_p - REFR_MZ_p - REMV_MZ_p;
%%-----------------------------------------------------------------------
%      Detrital Processes
%-----------------------------------------------------------------------
% Dissolution
DISS_DET_c = po_vals(77) * y_vals(30);
DISS_DET_n = po_vals(30) * po_vals(77) * y_vals(29);
DISS_DET_p = po_vals(29) * po_vals(77) * y_vals(28);
% Net Change
dydtt(30) = POM_PHY_c + POM_TR_c + POM_UN_c + POM_PRT_c + POM_MZ_c + POM_HZ_c - DISS_DET_c; 
dydtt(29) = POM_PHY_n + POM_TR_n + POM_UN_n + POM_PRT_n + POM_MZ_n + POM_HZ_n - DISS_DET_n;
dydtt(28) = POM_PHY_p + POM_TR_p + POM_UN_p + POM_PRT_p + POM_MZ_p + POM_HZ_p - DISS_DET_p;
%%-----------------------------------------------------------------------
%      Inorganic Nutrients
%-----------------------------------------------------------------------
% Nitrification
NTRF = po_vals(17) * y_vals(21);
% Dissolved Nutrient Rates of Change
dydtt(21) = dydtt(36) + remi_PRT_n + remi_MZ_n + remi_HZ_n + ...
    excr_TR_n_nh4 + excr_UN_n_nh4 -...
    growPHYnh4 - growTRnh4 - growUNnh4 - NTRF;
dydtt(20) = dydtt(35) - growPHYno3 - growTRno3 - growUNno3 + NTRF;
dydtt(15) = dydtt(34) + remi_PRT_p + remi_MZ_p + remi_HZ_p -...
    growPHYpo4 - growTRpo4 - growUNp;
%%-----------------------------------------------------------------------
%      Change in Dissolved Organic Matter
%-----------------------------------------------------------------------
% Conversion of SDOM to RDOM
REFR_SDOC = po_vals(72) * y_vals(11) * ...
    exp( 1 - min([(y_vals(10) / (y_vals(11)*po_vals(26))),...
    (y_vals(9) / (y_vals(11)*po_vals(25)))]) ) + ...
    max([ 0 , y_vals(11) - (y_vals(10)/po_vals(26)), y_vals(11) - (y_vals(9)/po_vals(25))]);
REFR_SDON = po_vals(26) * REFR_SDOC;
REFR_SDOP = po_vals(25) * REFR_SDOC;
% Change in DOM species
dydtt(27) = excr_PHY_LDOC + excr_TR_LDOC + excr_UN_LDOC + excr_PRT_LDOC + ...
    excr_MZ_LDOC + mortBAc - growBAldoc;
dydtt(26) = excr_PHY_LDON + excr_TR_LDON + excr_UN_LDON + excr_PRT_LDON + ...
    excr_MZ_LDON + mortBAn - growBAldon;
dydtt(25) = excr_PHY_LDOP + excr_TR_LDOP + excr_UN_LDOP + excr_PRT_LDOP + ...
    excr_MZ_LDOP + mortBAp - growBAldop;
dydtt(11) = excr_PHY_SDOC + excr_TR_SDOC + excr_UN_SDOC + excrBAc +...
    excr_PRT_SDOC + excr2_PRT_SDOC + excr_MZ_SDOC + excr2_MZ_SDOC +...
    excr_HZ_SDOC + DISS_DET_c - REFR_SDOC - growBAsdoc;
dydtt(10) = excr_PHY_SDON + excr_TR_SDON + excr_UN_SDON + excrBAn +...
    excr_PRT_SDON + excr2_PRT_SDON + excr_MZ_SDON + excr2_MZ_SDON +...
    excr_HZ_SDON + DISS_DET_n - REFR_SDON - growBAsdon;
dydtt(9) = excr_PHY_SDOP + excr_TR_SDOP + excr_UN_SDOP + excrBAp +...
    excr_PRT_SDOP + excr2_PRT_SDOP + excr_MZ_SDOP + excr2_MZ_SDOP +...
    excr_HZ_SDOP + DISS_DET_p - REFR_SDOP - growBAsdop;
%%-----------------------------------------------------------------------
%      Diagnostic Variables
%-----------------------------------------------------------------------
% Primary Production
dydtt(37) = growPHYc + growTRc + growUNc - respPHY - respTR - respUN - ...
    excr_PHY_LDOC - excr_TR_LDOC - excr_UN_LDOC -...
    excr_PHY_SDOC - excr_TR_SDOC - excr_UN_SDOC;
% Bacterial Production
dydtt(38) = growBAc - respBA;
dydtt=dydtt';
end