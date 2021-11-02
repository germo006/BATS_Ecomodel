%% Initializing Vals
% 'constructor' is a script that sets up the set and optimizable parameters
% and the initial state variables. It initializes these as three
% structures, as that was how I initially programmed the model and it gives
% me at least one script where the variable names are clear. 

constructor

% Here I also load two variables containing information about PAR. These
% will be used to generate a semirandom number within realistic estimates
% of surface PAR flux at BATS.  These two variables are required inputs for
% the model function as of 15 Mar 2021.

load('PARdata.mat', 'concat3PAR', 'concat3PARsig')

%% Creating Input Vectors
% This section runs the original code (that can handle a structure as 
% input) and then strips the parameter structures down to three column
% vectors. Each column vector is sorted by the name of the original
% variable, so to keep track of them the `_names` variables are created
% here, where they are in the same order as the rearranged values and can
% be used to identify variables of interest. The reason I have to evaluate
% dydtt here is because it was also written as a structure, whereas ODE
% solvers require it be a column vector of the same size as y (both yo and
% the output).

test = ode_mod_ecosys_noUN(y, ps, po, 0.5, z, concat3PAR, concat3PARsig);

dydtt_names = 'dydtt.' + string(fieldnames(test)); dydtt_vals = struct2array(test)';

po_names = 'po.' + string(fieldnames(po)); po_vals = struct2array(po)';
ps_names = 'ps.' + string(fieldnames(ps)); ps_vals = struct2array(ps)';
y_names =  'y.'+ string(fieldnames(y)); y_vals = struct2array(y)';

[dydtt_names, I] = sort(dydtt_names, 'descend'); dydtt_vals = dydtt_vals(I);
[po_names, I] = sort(po_names, 'descend'); po_vals = po_vals(I);
[ps_names, I] = sort(ps_names, 'descend'); ps_vals = ps_vals(I);
[y_names, I] = sort(y_names, 'descend'); y_vals = y_vals(I);


%% Running the forawrd model with preset params. 
% This runs the forward model, then using the integrated solution pulls out
% the primary productivity and bacterial productivity variables that were
% used by Luo et al. I have separated them by section because the getppbp
% function does take a couple minutes, and you don't have to run it if you
% don't want to. That being said, the pp and bp values produced by the
% solver are *not the real values*, since it is integrating the
% differential equations and pp, bp are the instantaneous rates rather than
% the integrals.

odefun = @(ts, y_vals) ode_mod_ecosys_noUN_vec(y_vals, ps_vals, po_vals, ts, z, concat3PAR, concat3PARsig);
opts = odeset('RelTol' , 1e-6, 'MaxStep', 0.1);
[to, yo] = ode45(odefun, [0:1/24:10000]', y_vals, opts);

%%
tic
[ppo, bpo, ~] = getppbp(odefun, to, yo);
toc
%% Plots.
% Depending on where you need to look and how much you need to zoom, choose
% x1 and x2 accordingly. 
load('AlbumMaps.mat')

x1 = 5000;
x2 = 10000;

subplot(3,1,1)
plot(to, yo(:,24), 'LineWidth', 1.5, 'Color', drift{4})
ylabel('LDOC, mmol C m^{-3}')
xlim([x1, x2])
subplot(3,1,2)
plot(to, yo(:,8), 'LineWidth', 1.5, 'Color', drift{2})
ylabel('SDOC, mmol C m^{-3}')
xlim([x1, x2])
subplot(3,1,3)
plot(to, yo(:,27), 'LineWidth', 1.5, 'Color', drift{1})
ylabel('DetC, mmol C m^{-3}')
xlabel('days')
xlim([x1, x2])


figure
subplot(3,1,1)
plot(to, yo(:, 30), 'LineWidth', 1.5, 'Color', drift{2})
ylabel('BAc, mmol C m^{-3}')
xlim([x1, x2])
subplot(3,1,2)
plot(to, yo(:, 16), 'LineWidth', 1.5, 'Color', drift{4})
ylabel('PHYc, mmol C m^{-3}')
xlim([x1, x2])
subplot(3,1,3)
plot(to, yo(:, 11), 'LineWidth', 1.5, 'Color', drift{1})
ylabel('PRTc, mmol C m^{-3}')
xlabel('days')
xlim([x1, x2])

figure
subplot(5,1,1)
plot(to, yo(:, 17), 'LineWidth', 1.5, 'Color', drift{2})
ylabel('NO_3^-, mmol m^{-3}')
xlim([x1, x2])
subplot(5,1,2)
plot(to, yo(:, 12), 'LineWidth', 1.5, 'Color', drift{4})
ylabel('PO_4^{-3}, mmol m^{-3}')
xlim([x1, x2])
subplot(5,1,3)
plot(to, yo(:, 4)+yo(:, 15), 'LineWidth', 1.5, 'Color', drift{1})
ylabel('Chl_{a}, mg m^{-3}')
xlabel('days')
xlim([x1, x2])
subplot(5,1,4)
plot(to, ppo, 'LineWidth', 1.5, 'Color', drift{3})
ylabel('Primary Prod., mmolC m^{-3} d^{-1}')
xlim([x1, x2])
subplot(5,1,5)
plot(to, bpo, 'LineWidth', 1.5, 'Color', drift{5})
ylabel('Bact Prod., mmolC m^{-3} d^{-1}')
xlabel('days')
xlim([x1, x2])