function [V0, CIl, CIr, epsilon] = BS_Call_MC_IS (S0, r, sigma, T, K, N, mu, alpha)
    
    // Allocate memory of N Monte-Carlo samples of the option payoff.
    samplesPayoff = zeros(N,1);  
      
    // Generate N random variables with N(mu,1)-distribution
    Y = grand(N, 1, "nor", mu, 1);
    
    
    // Compute Monte-Carlo estimator.
    samplesPayoff = exp((-1)*Y*mu+0.5*mu^2).*max((S0*exp((r-0.5*sigma^2)*T+sigma*sqrt(T)*Y))-K,0);
    V0 = exp(-r*T)*mean(samplesPayoff);

     // Determine asymptotic 95%-level confidence interval based on the sample variance
    // of the Monte-Carlo samples with control variate.
    sigma_hat = exp(-r*T) * stdev(samplesPayoff);
    alpha_level=cdfnor('X',0,1,(alpha+1)/2,(1-alpha)/2);
    epsilon = alpha_level * sigma_hat / sqrt(N);
    CIl=V0-epsilon;
    CIr=V0+epsilon;


endfunction


function [V0, epsilon] = EuOption_BS_MC (S0, r, sigma, T, K, M, g)
    
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

function C = BS_EuCall(t, S_t, r, sigma, T, K)
    
    // A function in a function: Provide the cumulative distribution function
    // of the standard normal distribution as function Phi, using the internal
    // scilab function cdfnor.
    function p = Phi(x)
        p = cdfnor("PQ", x, zeros(x), ones(x));
    endfunction
    
    // Implement Black-Scholes formula as given in C-Exercise 04
    d_1 = ( log(S_t./K) + (r+1/2*sigma^2)*(T-t) ) ./ ( sigma*sqrt(T-t) );
    d_2 = d_1 - sigma*sqrt(T-t);
    C = S_t.*Phi(d_1) - K.*exp(-r*(T-t)).*Phi(d_2);
    
endfunction

// Test data.
S0 = 100;
r = 0.05;
sigma = 0.2;
T = 1;
K = 200;
N = 1000;
alpha=0.95;
d =( log(K/S0) - (r-1/2*sigma^2)*T ) / ( sigma*sqrt(T) );
delta = abs(d);
mu = (1:0.001:d + delta)';

function y=g(x)
    y=max(x-K,0)
endfunction

// Call function for mu=d and display the result.
[V0, c1,c2,epsilon] =  BS_Call_MC_IS (S0, r, sigma, T, K, N, d, alpha);
[V0_plain, epsilon_plain] = EuOption_BS_MC(S0, r, sigma, T, K, N, g);
V0_BS = BS_EuCall(0, S0, r , sigma, T, K);
disp("The Monte-Carlo approximation with importance sampling to the price of the European Call option is given by " + string(V0) + ";   radius of 95% confidence interval: " + string(epsilon)+char(10) );
disp("The Monte-Carlo approximation without variance reduction is given by " + string(V0_plain) + ";   radius of 95% confidence interval: " + string(epsilon_plain)+char(10) );
disp("The real option price calculated with the BS-Formula is "+string(V0_BS));


// Comute a vector of approximations with mu from the interval [d-delta,d+delta]
M = length(mu);
V0_mu = length(M);
epsilon=length(M);
for i = 1:M
    [V0_mu(i),c1,c2,epsilon(i)] =  BS_Call_MC_IS (S0, r, sigma, T, K, N, mu(i));
end

// Compare with the Black Scholes formula (3.23) and the plain Monte Carlo 
// estimator (C-Exercise 21)

// Plot of standard estimator and IS estimators for mu from [d-delta,d+delta] 
clf();
plot(mu, [V0_mu V0_plain*ones(V0_mu) V0_BS*ones(V0_mu) V0_BS*ones(V0_mu)+epsilon,V0_BS*ones(V0_mu)-epsilon] );

title("Comparison of standard and importance sampling estimators for a European call option");
xlabel("$\mu$");
ylabel("option price");
legend("Importance sampling estimation", "Standard estimator", "Black-Scholes formula","CIr","CIl");
