function dPdt=phos2(t,P)

% PHOS2 is a function file representing a 2 box global ocean phosphate 
%    model to be used by ODE23 or ODE45 to integrate the ODEs
%    NB: Output must be a column vector for the integrator to work
%
%    Value initialization
%
% Started 2004:08:30 D. Glover, copied from Jenkins' web page (2002)
% Modif'd 2006:10:18 DMG to fix the V(1) typo in second equation
% Modif'd 2006:10:20 DMG to update some of the documentation

V=[3 100] * 1e16;   % volume of reservoirs in m3
dPdt=zeros(2,1);    % Initialize output as a column vector
Fx = 6e14;          % overturning water flux in m3 per year
Fr = 3e13;          % river water flux in m3 per year
Pr = 1.5;           % river water P concentration in mMol per m3
Sed = 0.01*Fx*P(1); % sedimentary loss of P in deep box
%
% difference equations
dPdt(1) = (Fr*Pr - Fx*P(1) + Fx*P(2))/V(1);
dPdt(2) = (Fx*P(1) - Fx*P(2) - Sed)/V(2);
