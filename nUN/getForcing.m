%% Noah Germolus 12 Mar 2021
% This is a script that pulls forcing values for the ode_mod functions.

function [forceNH4, forceNO3, forcePO4, T, MLD, PAR0] = getForcing(t, concat3PAR, concat3PARsig)

% NUTRIENTS
% These are average nutrient concentrations at the base of the mixed layer
% for HOTS in Luo et al. Check for Bermuda.
forceNH4 = 0.03;
forceNO3 = 4.89;
forcePO4 = 0.671;

% TEMPERATURE
% Only letting this vary by a single degree-ish
T = 21 + randn(1);

% MIXED LAYER DEPTH
% I have simplified this almost as much as I did the mixing equations that
% rely upon it. At Bermuda, this is very shallow in the summer, and Ruth
% Curry has demonstrated several distinct vertical zones with her glider
% data, which I might hope to use.
MLD = 50;

% SUNLIGHT
% To emulate variation in sunlight, I took three days of shipboard PAR data
% from a July in Bermuda, and calculated an hourly average and standard
% deviation  (n = 1800 measurements/hr * 3 days = 5400). 

% I then take whatever time is being used in the differential equation and 
% round it to the nearest hour.
tr = round(24*t);

% Because the PAR input vectors are 24 elements and I need to scan through 
% them again every 24 h, I create an index that cycles from 1-24.
tf = mod(tr, 24) +1; 

% I use that hour code as a seed for the rng so that if ode45 
% evaluates multiple times in that vicinity, the PAR doesn't change.
rng(tr)

% Finally, the value of PAR at the surface is determined. It must be >0.
PAR0 = max(0, normrnd(concat3PAR(tf),concat3PARsig(tf)));

end