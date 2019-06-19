//C-Exercise30
//Jurian Kahl
//Phanrattinon Nattawut

funcprot(0)
function [V0, CIl, CIr] = BS_EuCall_MC_IS (S0, r, sigma, K, T, mu, N, alpha)
    
    X = zeros(N,1);  
    Y = grand(N, 1, "nor", mu, 1);
    
    X = exp((-1)*Y*mu+0.5*mu^2).*max((S0*exp((r-0.5*sigma^2)*T+sigma*sqrt(T)*Y))-K,0);
    V0 = exp(-r*T)*mean(X);

    sigma = exp(-r*T) * stdev(X);
    Z_value=cdfnor('X',0,1,1-((1-alpha)/2),(1-alpha)/2);
    CIl=V0 - (Z_value * sigma / sqrt(N));
    CIr=V0 + (Z_value * sigma / sqrt(N));

endfunction


// Test data.
S0 = 100;
r = 0.05;
sigma = 0.3;
K = 220;
T = 1;
N = 10000;
alpha=0.95;
d =( log(K/S0) - (r-1/2*sigma^2)*T ) / ( sigma*sqrt(T) );

function y=g(x)
    y=max(x-K,0)
endfunction

//Display the result for mu=d. 
[V0, CIl,CIr] =  BS_EuCall_MC_IS (S0, r, sigma, K, T, d, N, alpha);
disp("The approximation of European Call option price with importance sampling is  " + string(V0) + ",  with " + "[" +string(CIl)+","+string(CIr) + "]" + " 95% confidence interval");
