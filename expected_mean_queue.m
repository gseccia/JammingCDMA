function [ nQ, p0, pQ, rho] = expected_mean_queue(lambda, mu, k)
%NQ Summary of this function goes here
%   Detailed explanation goes here
rho = lambda/(k*mu);
p_constant = ((k*rho)^k)/(factorial(k)*(1-rho));

p0 = p_constant;
for n = 0:k-1
    p0 = p0 + ((k*rho)^n)/factorial(n);
end
p0 = 1/p0;

pQ = p0 * p_constant;

nQ = (pQ*rho)/(1-rho);

end
