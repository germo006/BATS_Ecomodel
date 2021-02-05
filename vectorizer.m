%% Noah Germolus 02 Feb 2021
% The purpose of this script is to alter the text of a copy of
% ode_mod_ecosys.m so that the inputs are *not* structures, but vectors.
% Vectorizer is a bit of a misnomer given that the model equations are in
% no way vectorized.

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

%run once
filetxt = string(fileread('Copy_of_ode_mod_ecosys.m'));
po_replace = 'po_vals(' + string(1:length(po_vals))' + ')';
ps_replace = 'ps_vals(' + string(1:length(ps_vals))' + ')';
y_replace = 'y_vals(' + string(1:length(y_vals))' + ')';
dydtt_replace = 'dydtt(' + string(1:length(dydtt_vals))' + ')';
allnames = [po_names;ps_names;y_names;dydtt_names];
allreplace = [po_replace;ps_replace;y_replace;dydtt_replace];

filetxt = string(fileread('Copy_of_ode_mod_ecosys.m'));
for ii = 1:length(allreplace)
    filetxt = strrep(filetxt, allnames{ii}, allreplace{ii});
end
filetxt = regexprep(filetxt, '[\n\r]+', '\n');

fid = fopen('ode_mod_ecosys_vec.m', 'wt');
fwrite(fid, filetxt);
fclose(fid);
%After the fclose, the function intro  must be modded.