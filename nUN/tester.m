constructor
%%
test = ode_mod_ecosys_noUN(y, ps, po, 0.5, z);

dydtt_names = 'dydtt.' + string(fieldnames(test)); dydtt_vals = struct2array(test)';

po_names = 'po.' + string(fieldnames(po)); po_vals = struct2array(po)';
ps_names = 'ps.' + string(fieldnames(ps)); ps_vals = struct2array(ps)';
y_names =  'y.'+ string(fieldnames(y)); y_vals = struct2array(y)';

[dydtt_names, I] = sort(dydtt_names, 'descend'); dydtt_vals = dydtt_vals(I);
[po_names, I] = sort(po_names, 'descend'); po_vals = po_vals(I);
[ps_names, I] = sort(ps_names, 'descend'); ps_vals = ps_vals(I);
[y_names, I] = sort(y_names, 'descend'); y_vals = y_vals(I);


%% Running the forawrd model with preset params. 

odefun = @(ts, y_vals) ode_mod_ecosys_noUN_vec(y_vals, ps_vals, po_vals, ts, z);
opts = odeset('RelTol' , 1e-6, 'MaxStep', 0.1);
[to, yo] = ode45(odefun, [0:1/24:10000]', y_vals, opts);
[ppo, bpo] = getppbp(odefun, to, yo);

%%
load('AlbumMaps.mat')

x1 = 9900;
x2 = 10000;

subplot(3,1,1)
plot(to, yo(:,24), 'LineWidth', 1.5, 'Color', CP1{4})
ylabel('LDOC, mmol C m^{-3}')
xlim([x1, x2])
subplot(3,1,2)
plot(to, yo(:,8), 'LineWidth', 1.5, 'Color', CP1{2})
ylabel('SDOC, mmol C m^{-3}')
xlim([x1, x2])
subplot(3,1,3)
plot(to, yo(:,27), 'LineWidth', 1.5, 'Color', CP1{1})
ylabel('DetC, mmol C m^{-3}')
xlabel('days')
xlim([x1, x2])


figure
subplot(3,1,1)
plot(to, yo(:, 30), 'LineWidth', 1.5, 'Color', CP1{2})
ylabel('BAc, mmol C m^{-3}')
xlim([x1, x2])
subplot(3,1,2)
plot(to, yo(:, 16), 'LineWidth', 1.5, 'Color', CP1{4})
ylabel('PHYc, mmol C m^{-3}')
xlim([x1, x2])
subplot(3,1,3)
plot(to, yo(:, 11), 'LineWidth', 1.5, 'Color', CP1{1})
ylabel('PRTc, mmol C m^{-3}')
xlabel('days')
xlim([x1, x2])

figure
subplot(5,1,1)
plot(to, yo(:, 17), 'LineWidth', 1.5, 'Color', CP1{2})
ylabel('NO_3^-, mmol m^{-3}')
xlim([x1, x2])
subplot(5,1,2)
plot(to, yo(:, 12), 'LineWidth', 1.5, 'Color', CP1{4})
ylabel('PO_4^{-3}, mmol m^{-3}')
xlim([x1, x2])
subplot(5,1,3)
plot(to, yo(:, 4)+yo(:, 15), 'LineWidth', 1.5, 'Color', CP1{1})
ylabel('Chl_{a}, mg m^{-3}')
xlabel('days')
xlim([x1, x2])
subplot(5,1,4)
plot(to, ppo, 'LineWidth', 1.5, 'Color', CP1{3})
ylabel('Primary Prod., mmolC m^{-3} d^{-1}')
xlim([x1, x2])
subplot(5,1,5)
plot(to, bpo, 'LineWidth', 1.5, 'Color', CP1{1})
ylabel('Bact Prod., mmolC m^{-3} d^{-1}')
xlabel('days')
xlim([x1, x2])