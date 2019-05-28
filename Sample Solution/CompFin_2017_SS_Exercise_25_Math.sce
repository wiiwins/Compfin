// Implement closed-form pricing formula for a geometric average option
// in the Black-Scholes model as given in T-Exercise 14.
function [V0, epsilon] = BS_GeomAvgOption_CF (S0, r, sigma, K, T, M)
    
    // Cumulative distribution function of the standard normal distribution.
    function p = Phi(x)
        p = cdfnor("PQ", x, zeros(x), ones(x));
    endfunction
           
    a = (r-0.5*sigma^2)*T/2;
    b = sqrt( (2*M+1)/(6*(M+1))) * sigma * sqrt(T);
    
    V0 = exp(-r*T) * ( S0 * exp(a+b^2/2) * Phi( (log(S0/K) + a)/b + b ) - K * Phi( (log(S0/K) + a)/b) );
           
endfunction

// Function for the control variate Monte Carlo simulation
function [V0, epsilon] = BS_AsianCall_MC_CV (S0, r, sigma, K, T, M, N)
    
    // Allocate memory of N Monte-Carlo samples of the option payoff.
    samplesPayoff = zeros(N,1);
    
    // Generate N Monte-Carlo samples of the stock price at the monitoring times
    // and compute the difference of the payoffs of the Asian option and the control variate.
    for k=1:N
        path=S0*ones(1,M+1);
        for l=2:M+1
            path(l)=path(l-1)*exp((r-sigma^2/2)*T/M+sigma*sqrt(T/M)*grand(1, 1, 'nor', 0, 1))
        end
        
        samplesPayoff(k) =  max( mean(path)  - K, 0 ) - max(  prod(path)^(1/(M+1)) - K, 0);        
    end
        
    // Compute Monte-Carlo estimator with control variate.
    V0 = exp(-r*T) *mean(samplesPayoff) + BS_GeomAvgOption_CF (S0, r, sigma, K, T, M);
    
    // Determine asymptotic 95%-level confidence interval based on the sample variance
    // of the Monte-Carlo samples with control variate.
    sigma_hat = exp(-r*T) * stdev(samplesPayoff);
    epsilon = 1.96 * sigma_hat / sqrt(N);
       
endfunction


// Function for Asian Call option with plain Monte Carlo
function [V0,epsilon] = BS_AsianCall_MC (S0, r, sigma, K, T, M, N)
    
    // Allocate memory of N Monte-Carlo samples of the option payoff.
    samplesPayoff = zeros(N,1);
    
    // Generate N Monte-Carlo samples of the stock price at the monitoring times
    // and compute the difference of the payoffs of the Asian option and the control variate.
    for k=1:N
        path=S0*ones(1,M+1);
        for l=2:M+1
            path(l)=path(l-1)*exp((r-sigma^2/2)*T/M+sigma*sqrt(T/M)*grand(1, 1, 'nor', 0, 1))
        end
        
        samplesPayoff(k) =  max( mean(path)  - K, 0 );        
    end
        
    // Compute Monte-Carlo estimator with control variate.
    V0 = exp(-r*T) *mean(samplesPayoff);
    
    // Determine asymptotic 95%-level confidence interval based on the sample variance
    // of the Monte-Carlo samples with control variate.
    sigma_hat = exp(-r*T) * stdev(samplesPayoff);
    epsilon = 1.96 * sigma_hat / sqrt(N);
    
endfunction

// Test data.
S0 = 100;
r = 0.05;
sigma = 0.2;
K = 100;
T = 1;
M = 50;
N = 10000;

// Call function and display the result.
[V0, epsilon] = BS_AsianCall_MC_CV (S0, r, sigma, K, T, M, N);
disp("The Monte-Carlo approximation with control variate to the price of the Asian call option is given by " + string(V0) + ";  radius of 95% confidence interval: " + string(epsilon) );

//Comparison with plain Monte Carlo
[V0_plain, epsilon_plain] = BS_AsianCall_MC (S0, r, sigma, K, T, M, N);
disp("The Monte-Carlo approximation with plain Monte Carlo to the price of the Asian call option is given by " + string(V0_plain) + ";  radius of 95% confidence interval: " + string(epsilon_plain) );

// We observe a significantly tighter confidence interval compared to the plain Monte-Carlo estimator.
