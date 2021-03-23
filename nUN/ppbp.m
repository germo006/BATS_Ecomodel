function [pp,bp] = ppbp(diaginput)
%PPBP is designed to pull out the "diagnostic variables" primary production
%and bacterial production from the results of an ODE solve of the model.
%   DEPRECATED 23 MAR 2021

pp = sum(diaginput(1:8)); bp = sum(diaginput(9:10));

end

