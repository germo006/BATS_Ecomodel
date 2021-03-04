%% Noah Germolus, 13 Jan 2021
% This is the system of differential equations representing DOM
% in Luo et al. (2010)


function [dydtt] = ode_mod_dom(y)
%-----------------------------------------------------------------------
%      Change in Dissolved Organic Matter
%-----------------------------------------------------------------------

% Conversion of SDOM to RDOM

REFR_SDOC = ex_REFR_SDOM * y.iSDOC * ...
    exp( 1 - min((y.iSDON / (y.iSDOC*q_REFR_n)),...
    (y.iSDOP / (y.iSDOC*q_REFR_p))) ) + ...
    max( 0 , y.iSDOC - (y.iSDON/q_REFR_n), y.iSDOC - (y.iSDOP/q_REFR_p));
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

end