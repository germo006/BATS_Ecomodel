%% Noah Germolus 15 Oct 2020
% This sets parameters and runs a solution of BasicTrp.m and requires
% Sunlight.m for light attenuation.

tstep = 1/(24*60); % One minute
tspan = 100;         % Days



t = 0:tstep:tspan;
z = 1;
P = 0.5;

dT = @(t, T) BasicTrp(t, T, P, z);
[ts, T] = ode45(dT, t, 2e-10);
plot(t,T)
% hold on
% yyaxis right
% [~, I] = Sunlight(t, z);
% Q_trp = 0.9 - 0.5*cos(2*pi*(t - 200/(24*60)));
% plot(t, Q_trp/max(Q_trp))
% plot(t, I/max(I))