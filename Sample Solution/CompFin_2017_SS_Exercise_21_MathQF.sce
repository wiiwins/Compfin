function [V0, c1, c2] = EuOption_BS_MC (S0, r, sigma, T, M, g)
        
    // Generate an Mx1-vector of independent samples from
    // standard normally distributed random variables.
    X = grand(M, 1, 'nor', 0, 1);
    ST = S0*exp( (r-0.5*sigma^2)*T + sigma*sqrt(T)*X );
    
    // Compute Monte-Carlo estimator.
    Y = exp(-r*T)*g(ST);
    V0 = mean(Y);
    
    // Compute 95% confidence interval 
    epsilon = 1.96 * sqrt(variance(Y)/M);
    c1=V0-epsilon;
    c2=V0+epsilon;
    
endfunction

// test parameters
S0 = 95;
r = 0.05;
sigma = 0.2;
T = 1;
K = 100;
M = 100000;
function y=g(x)
    y=max(100-x,0)
endfunction


[V0, c1, c2] = EuOption_BS_MC (S0, r, sigma, T, M, g);
disp("Price of European put by use of plain Monte-Carlo simulation: " + string(V0) + ", 95% confidence interval: [" + string(c1)+","+string(c2)+"]");


//Closed formula for Put price (not asked for in the exercise)
function V0 = BS_EuPut(t, S_t, r, sigma, T, K)
    
    // A function in a function: Provide the cumulative distribution function
    // of the standard normal distribution as function Phi, using the internal
    // scilab function cdfnor.
    function p = Phi(x)
        p = cdfnor("PQ", x, zeros(x), ones(x));
    endfunction
    
    // Implement Black-Scholes formula (3.23) and below
    d_1 = ( log(S_t./K) + (r+1/2*sigma^2)*(T-t) ) ./ ( sigma*sqrt(T-t) );
    d_2 = d_1 - sigma*sqrt(T-t);
    V0 = K.*exp(-r*(T-t)).*Phi(-d_2)-S_t.*Phi(-d_1);
    
endfunction

V0 = BS_EuPut (0, S0, r, sigma, T, K);
disp("Exact price of European Put by use of the BS-Formula: " + string(V0));
