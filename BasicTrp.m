function dTdt=BasicTrp(t,T,P, z)

% PHOS2 is a function file representing a 2 box global ocean phosphate 
%    model to be used by ODE23 or ODE45 to integrate the ODEs
%    NB: Output must be a column vector for the integrator to work
%
%    Value initialization

Q_trp = 0.9 - 0.5*cos(2*pi*(t));
[~, I] = Sunlight(t, z);

kph = 1e-17;             % Photochemical rate constant, M molph^-1
Qm_trp = 0.7;               % Minimum trp cell quota ng/ug C
quo = (Q_trp./Qm_trp) - 1;       % How satisfied is this requirement?
k_trp = 1e-12;             % How fast is transfer? /ug C /s
k_ex = quo*1e-20;


if T<=1e-12 || quo<0
    dTdt = sum(k_ex.*P)- sum(kph.*I);
else
    dTdt = -sum(k_trp.*P.*T) - sum(kph.*I);

end