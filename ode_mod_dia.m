%% Noah Germolus, 13 Jan 2021
% This is the system of differential equations representing diagnostics
% in Luo et al. (2010)
% It was especially obvious that this should all be in one function file
% when the input to this wasn't used at all.

function [dydtt] = ode_mod_dia(y)

%-----------------------------------------------------------------------
%      Diagnostic Variables
%-----------------------------------------------------------------------

% Primary Production

dydtt.PrPr = growPHYc + growTRc + growUNc - respPHY - respTR - respUN - ...
    excr_PHY_LDOC - excr_TR_LDOC - excr_UN_LDOC -...
    excr_PHY_SDOC - excr_TR_SDOC - excr_UN_SDOC;

% Bacterial Production

dydtt.BPr = growBAc - respBA;

end