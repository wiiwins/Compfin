function [V0, epsilon] = EuOption_BS_MC_AV (S0, r, sigma, T, M, g)
    
    // The initial price of the option is given by E_Q (f(Z))
    // for a N(0,1)-distributed random variable Z and the following
    // function f.
    function y = f(x)
        y = exp(-r*T) * g(S0*exp((r-0.5*sigma^2)*T + sigma*sqrt(T)*x));
    endfunction
    
    // Generate an Mx1-vector of independent samples from
    // standard normally distributed random variables.
    X = grand(M, 1, 'nor', 0, 1);
    
    // Compute Monte-Carlo estimator using antithetic variables.
    Y = (f(X)+f(-X)) / 2;
    V0 = mean(Y);
    
    // Compute radius of 95% confidence interval.
    epsilon = 1.96 * sqrt(variance(Y)/M);
    
endfunction

function [V0, epsilon] = EuOption_BS_MC (S0, r, sigma, T, M, g)
    
    // The initial price of the european option is given by E_Q (f(Z))
    // for a N(0,1)-distributed random variable Z and the following
    // function f.
    function y = f(x)
        y = exp(-r*T) * g(S0*exp( (r-0.5*sigma^2)*T + sigma*sqrt(T)*x));
    endfunction
    
    // Generate an Mx1-vector of independent samples from
    // standard normally distributed random variables.
    X = grand(M, 1, 'nor', 0, 1);
    
    // Compute Monte-Carlo estimator.
    Y = f(X);
    V0 = mean(Y);
    
    // Compute radius of 95% confidence interval 
    epsilon = 1.96 * sqrt(variance(Y)/M);
    
endfunction

function [V0, phi] = BS_EuPut(t, S_t, r, sigma, T, K)
    
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


// test parameters
S0 = 100;
r = 0.05;
sigma = 0.2;
T = 1;
K = 100;
M = 100000;

function y=g(x)
    y=max(K-x,0)
endfunction

// display result for test parameters 
V0 = BS_EuPut (0, S0, r, sigma, T, K);
disp("Exact price of the European put calculated with the BS-Formula: " + string(V0) );

[V0, epsilon] = EuOption_BS_MC_AV (S0, r, sigma, T, M, g);
disp("Price of European call by use of Monte-Carlo simulation with antithetic variables for M simulations: " + string(V0) + ", radius of 95% confidence interval: " + string(epsilon)+ char(10)); //char(10) adds a newline at the the end of the output

[V0, epsilon] = EuOption_BS_MC (S0, r, sigma, T, 2*M, g);
disp("Price of European call by use of plain Monte-Carlo simulation for 2M simulations: " + string(V0) + ", radius of 95% confidence interval: " + string(epsilon));

