%% Noah Germolus, 11 Jan 2021
% This is the simple Tfun for Luo et al. (2010)
% I assume for now that the Arrhenius parameter is one that needs to be
% optimized.

function Tfun = tfun(y)
% Here I also assume that the structure y contains the temperature values.
% AE is to-be-optimized.

Tfun = exp(-AE*((1/(y.T+273.15))-(1/(25+273.15))));

end