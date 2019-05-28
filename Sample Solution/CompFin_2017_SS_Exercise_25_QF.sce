//Write function for computing the price of a European option in the BS-model using Monte-Carlo with control variate S(T)
function [V0, epsilon, beta] = EuOption_BS_MC_CV (S0, r, sigma, T, K, M, g)
    
    // Determine beta = Cov(V(T),S(T)) / Var(S(T)) by Monte Carlo simulation.
    X = grand(M, 1, 'nor', 0, 1);
    ST = S0*exp( (r-0.5*sigma^2)*T + sigma*sqrt(T)*X );
    VT = g(ST);
    Covar = mean( (ST-exp(r*T)*S0) .* (VT-mean(VT)) );
    beta = Covar / variance(ST);
    // (In the Black-Scholes model, we could compute Var(S(T)) also
    // analytically.)
    
    // Compute Monte Carlo estimator using the (discounted) stock as control variate.
    X = grand(M, 1, 'nor', 0, 1);
    ST_hat = S0*exp( (r-0.5*sigma^2)*T + sigma*sqrt(T)*X );
    Y = g(ST_hat)-beta*ST_hat;
    V0 = exp(-r*T)*mean(Y) + beta*S0;
    
    // Compute radius of 95% confidence interval (not part of the exercise).
    epsilon = 1.96 * sqrt(variance(Y)/M);
    
endfunction

// test parameters
S0 = 120;
r = 0.05;
sigma = 0.2;
T = 1;
K = 100;
M = 1000000;

//Payoff function for the European call
function y=g(x)
    y=max(x-K,0);
endfunction


// display result for test parameters
[V0, epsilon, beta] = EuOption_BS_MC_CV (S0, r, sigma, T, K, M, g);
disp("Price of European call by use of Monte-Carlo simulation with control variate: " + string(V0) + ", radius of 95% confidence interval: " + string(epsilon) + " Beta chosen in the estimation procedure:" + string(beta));

//Old function that computes the price of a European option in the BS-model using plain Monte-Carlo from C-Exercise 21
function [V0, epsilon] = EuOption_BS_MC (S0, r, sigma, T, K, M, g)
        
    // Generate an Mx1-vector of independent samples from
    // standard normally distributed random variables.
    X = grand(M, 1, 'nor', 0, 1);
    ST = S0*exp( (r-0.5*sigma^2)*T + sigma*sqrt(T)*X );
    
    // Compute Monte-Carlo estimator.
    Y = exp(-r*T)*g(ST);
    V0 = mean(Y);
    
    // Compute radius of 95% confidence interval 
    epsilon = 1.96 * sqrt(variance(Y)/M);
    
endfunction

[V0, epsilon] = EuOption_BS_MC (S0, r, sigma, T, K, M, g);
disp("Price of European call by use of plain Monte-Carlo simulation: " + string(V0) + ", radius of 95% confidence interval: " + string(epsilon));

//Closed formula for the price of a European Call in the BS-model
function V0=BS_EuCall_Closed(S0, r, sigma, T, K)
    d1=(log(S0/K)+r*T+sigma^2/2*T)/(sigma*sqrt(T));
    d2=d1-sigma*sqrt(T);
    V0=S0*cdfnor("PQ",d1,0,1)-K*exp(-r*T)*cdfnor("PQ",d2,0,1)
endfunction

V0=BS_EuCall_Closed(S0, r, sigma, T, K);
disp("Price of European call by use of the closed formula: " + string(V0));
