constructor
%%
test = ode_mod_ecosys(y, ps, po, 0.5);

dydtt_names = 'dydtt.' + string(fieldnames(test)); dydtt_vals = struct2array(test)';

po_names = 'po.' + string(fieldnames(po)); po_vals = struct2array(po)';
ps_names = 'ps.' + string(fieldnames(ps)); ps_vals = struct2array(ps)';
y_names =  'y.'+ string(fieldnames(y)); y_vals = struct2array(y)';

[dydtt_names, I] = sort(dydtt_names, 'descend'); dydtt_vals = dydtt_vals(I);
[po_names, I] = sort(po_names, 'descend'); po_vals = po_vals(I);
[ps_names, I] = sort(ps_names, 'descend'); ps_vals = ps_vals(I);
[y_names, I] = sort(y_names, 'descend'); y_vals = y_vals(I);


%% Running the forawrd model with preset params. 

odefun = @(ts, y_vals) ode_mod_ecosys_vec(y_vals, ps_vals, po_vals, ts);
[to, yo] = ode45(odefun, [0:0.001:1000]', y_vals);

plot(to, yo(:,[11,27,30]), 'LineWidth', 1.5)
legend('SDOC', 'LDOC', 'RDOC')

figure
plot(to, yo(:, [33,19,14]), 'LineWidth', 1.5)
legend('BAc', 'PHYc', 'PRTc')

