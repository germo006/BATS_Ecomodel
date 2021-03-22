function [ppo, bpo] = getppbp(odefcn, to, yo)

bpo = zeros(length(to),1); ppo = bpo;

for ii = 1:length(to)
    dydtt = odefcn(to(ii), yo(ii,:)');
    ppo(ii) = dydtt(1);
    bpo(ii) = dydtt(31);
end

end