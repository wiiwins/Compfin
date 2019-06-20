//C-Exercise 31
//Jurian Kahl
//Phanrattinon Nattawut

function V0 = BS_EuOption_MC_CV (S0, r, sigma, T, K, M, f, g)
    
    X = grand(M, 1, 'nor', 0, 1);
    ST = S0*exp( (r-0.5*sigma^2)*T + sigma*sqrt(T)*X );
    Covar = mean( (f(ST)-mean(f(ST))) .* (g(ST)-mean(g(ST))) );
    beta = Covar / variance(f(ST));

    // Using max(ST-K,0) as control variate.
    X = grand(M, 1, 'nor', 0, 1);
    ST = S0*exp( (r-0.5*sigma^2)*T + sigma*sqrt(T)*X );
    dif = g(ST)-beta*f(ST);
    V0 = exp(-r*T)*mean(dif) + mean(exp(-r*T)*f(ST));
    
endfunction

// test parameters
S0 = 100;
r = 0.05;
sigma = 0.3;
T = 1;
K = 110;
M = 100000;

// Define more function for Calculation
function y=g(x)
    y=(max(x-K,0)).*x;
endfunction

function y=f(x)
    y=max(x-K,0);
endfunction

// display result for test parameters
V0 = BS_EuOption_MC_CV (S0, r, sigma, T, K, M, f, g);
disp("European call price by useing Monte-Carlo simulation with control variate: " + string(V0));

//plain Monte-Carlo from C-Exercise 27
function V0 = Eu_Option_BS_MC (S0, r, sigma, T, K, M, g)
        
    X = grand(M, 1, 'nor', 0, 1);
    ST = S0*exp( (r-0.5*sigma^2)*T + sigma*sqrt(T)*X );
    
    Y = exp(-r*T)*g(ST);
    V0 = mean(Y);
    
endfunction

V0 = Eu_Option_BS_MC (S0, r, sigma, T, K, M, g);
disp("European call price by useing plain Monte-Carlo simulation: " + string(V0));
