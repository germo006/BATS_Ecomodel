function [ppo, bpo, dydtto] = getppbp(odefcn, to, yo)

bpo = zeros(length(to),1); ppo = bpo;
dydtto = zeros(size(to,1), size(yo,2));

for ii = 1:length(to)
    dydtt = odefcn(to(ii), yo(ii,:)');
    ppo(ii) = dydtt(1);
    bpo(ii) = dydtt(31);
    dydtto(ii, :) = dydtt';
end

end