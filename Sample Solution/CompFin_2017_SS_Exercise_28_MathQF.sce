
function [V0, c1, c2] = Heston_EuCall_MC_Euler (S0, r, nu0, kappa, lambda, sigma_tilde, T, g, M, m)
    
    // Time step.
    delta_t = T/m;

    // Initialize the vector for stock prices S, X = log(S/S0) and nu.
    S = S0 * ones(M, 1);
    X = zeros(M, 1);
    nu = nu0*ones(M, 1);
    
    // Loop over t_1,...t_m.
    for k=1:m
        
        // Simulate increments of Brownian motion.
        delta_WQ = grand(M, 1, 'nor', 0, sqrt(delta_t));
        delta_W_tilde = grand(M, 1, 'nor', 0, sqrt(delta_t));
        
        // Euler step: Simulate the next values for X and nu.
        X = X + (r - 0.5*nu)*delta_t + sqrt(nu).*delta_WQ;
        nu = nu + (kappa - lambda*nu)*delta_t + sigma_tilde*sqrt(nu).*delta_W_tilde;
        
        // Correct nu if the simulation is negative.
        nu = max(nu, 0);
    end
    
    // Stock prices at maturity.
    S = S0*exp(X);    
       
    // Compute (discounted) payoff.
    VT =(exp(-r*T)*g(S));
        
    // Compute Monte Carlo estimate.
    V0 = mean(VT);
    
    epsilon = 1.96 * sqrt(variance(VT)/M);
    c1 = V0 - epsilon;
    c2 = V0 + epsilon;
               
endfunction

exec("./CompFin_2017_SS_QF_Exercise_17.sce");

// test parameters
S0 = 100;
r = 0.05;
nu0 = 0.2^2;
kappa = 0.5;
lambda = 2.5;
sigma_tilde = 1;
T = 1;
K = 100;
M = 10000;
m = 2500;

function y=g(x)
    y = max(x-K,0);
endfunction

// display result for test parameters
[V0,c1,c2] = Heston_EuCall_MC_Euler (S0, r, nu0, kappa, lambda, sigma_tilde, T, g, M, m);
disp("Price of European Call in Heston model using Euler method: " + string(V0) +". Lower bound of 95% confidence interval: " + string(c1)+". Upper bound of 95% confidence interval: "+ string(c2)+".");
