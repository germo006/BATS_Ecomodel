constructor
%%
test = ode_mod_ecosys_noUN(y, ps, po, 0.5);

dydtt_names = 'dydtt.' + string(fieldnames(test)); dydtt_vals = struct2array(test)';

po_names = 'po.' + string(fieldnames(po)); po_vals = struct2array(po)';
ps_names = 'ps.' + string(fieldnames(ps)); ps_vals = struct2array(ps)';
y_names =  'y.'+ string(fieldnames(y)); y_vals = struct2array(y)';

[dydtt_names, I] = sort(dydtt_names, 'descend'); dydtt_vals = dydtt_vals(I);
[po_names, I] = sort(po_names, 'descend'); po_vals = po_vals(I);
[ps_names, I] = sort(ps_names, 'descend'); ps_vals = ps_vals(I);
[y_names, I] = sort(y_names, 'descend'); y_vals = y_vals(I);


%% Running the forawrd model with preset params. 

odefun = @(ts, y_vals) ode_mod_ecosys_noUN_vec(y_vals, ps_vals, po_vals, ts);
opts = odeset('RelTol' , 1e-6, 'MaxStep', 0.1);
[to, yo] = ode45(odefun, [0:0.1:100000]', y_vals, opts);

%%
load('AlbumMaps.mat')

x1 = 9980;
x2 = 10000;

subplot(3,1,1)
plot(to, yo(:,21), 'LineWidth', 1.5, 'Color', CP1{4})
ylabel('LDOC, mmol C m^{-3}')
xlim([x1, x2])
%set(gca, 'Color', [0.5 0.5 0.5], 'LineWidth', 1)
subplot(3,1,2)
plot(to, yo(:,7), 'LineWidth', 1.5, 'Color', CP1{2})
ylabel('SDOC, mmol C m^{-3}')
xlim([x1, x2])
%set(gca, 'Color', [0.6 0.6 0.6], 'LineWidth', 1)
subplot(3,1,3)
plot(to, yo(:,26), 'LineWidth', 1.5, 'Color', CP1{1})
ylabel('DetC, mmol C m^{-3}')
xlabel('days')
xlim([x1, x2])
%set(gca, 'Color', [0.9 0.9 0.9], 'LineWidth', 1)

figure
subplot(3,1,1)
plot(to, yo(:, 29), 'LineWidth', 1.5, 'Color', CP1{2})
ylabel('BAc, mmol C m^{-3}')
xlim([x1, x2])
%set(gca, 'Color', [0.6 0.6 0.6], 'LineWidth', 1)
subplot(3,1,2)
plot(to, yo(:, 15), 'LineWidth', 1.5, 'Color', CP1{4})
ylabel('PHYc, mmol C m^{-3}')
xlim([x1, x2])
%set(gca, 'Color', [0.5 0.5 0.5], 'LineWidth', 1)
subplot(3,1,3)
plot(to, yo(:, 10), 'LineWidth', 1.5, 'Color', CP1{1})
ylabel('PRTc, mmol C m^{-3}')
xlabel('days')
%set(gca, 'Color', [0.9 0.9 0.9], 'LineWidth', 1)
xlim([x1, x2])