function k = pantothenic_khydro(pH)
%PANTOTHENIC_HYDRO Summary of this function goes here
%   Detailed explanation goes here

k = ((7.73e5).*10.^(-pH)) + 9.67e5.*10.^(pH-14);

end

