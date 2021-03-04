%% Noah Germolus, 11 Jan 2021
% This is the system of differential equations representing bacterial
% processes in Luo et al. (2010)

function [dydtt] = ode_mod_bact(y)

%-----------------------------------------------------------------------
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
growBAldoc = mu_BA * y.iBAc * tfun(y) * temp * ALC / (ALC+ k_DOM + ASC);
growBAsdoc = mu_BA * y.iBAc * tfun(y) * temp * ASC / (ASC+ k_DOM + ALC);
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
end