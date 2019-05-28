function [V0, S] = BS_EuPut_FiDi_Explicit (r, sigma, a, b, m, nu_max, T, K)
    
     // Allocate memory for result vector.
    w = zeros(m + 1, 1);
    
    // Time and space discretization step sizes.
    delta_t = 0.5*sigma^2*T / nu_max;
    delta_x = (b-a) / m;
    
    // Check for stability of explicit scheme.
    lambda = delta_t / delta_x^2;
    if (lambda >= 1/2) then
        error("Finite difference scheme unstable.");
    end
    
    // Determine grid in time and space.
    t_tilde = (0:1:nu_max) * delta_t;
    x_tilde = a + (0:1:m) * delta_x;
    
    // Compute auxiliary variables.
    q = 2*r/sigma^2;
    alpha = 0.5 * (q-1);
    beta = 0.5 * (q+1);
            
    // Boundary condition for t_tilde = 0, corresponding to payoff
    // of the put option at maturity.
    w = max( exp(alpha*x_tilde)-exp(beta*x_tilde), 0 );
              
    
    
    // Boundary condition for x = b. Stays the same for all time points
    w(m+1) = 0;
    
    // Iterate through time layers.
    for nu=1:nu_max
        w(2:$-1) =  lambda*[w(1:$-2)] + (1-2*lambda)*w(2:$-1) + lambda*[w(3:$)]
        
        // Boundary condition for x = a.
        w(1) = exp(alpha*a + alpha^2*t_tilde(nu+1))-exp( beta*a + beta^2*t_tilde(nu+1) );
        
    end        
    
    // retransformation of heat equation
    V0 =( K .* w .* exp(-alpha*x_tilde - 0.5* sigma^2*T*(alpha^2 + q)));
    
    // vector of initial stock prices 
    S = K.*exp(x_tilde);
    
endfunction

// test parameters
r = 0.05;
sigma = 0.2;
a = -0.7;
b = 0.4;
m = 100; 
nu_max = 2000; 
T = 1; 
K = 100; 

// Application of function
[V0,S] = BS_EuPut_FiDi_Explicit (r, sigma, a, b, m, nu_max, T, K)

// Black-Scholes formula for the European Put from the lecture notes.
function [C, phi] = BS_EuPut(t, S_t, r, sigma, T, K)
    
    // cdf of standatd normal distribution
    function p = Phi(x)
        p = cdfnor("PQ", x, zeros(x), ones(x));
    endfunction
    
    // Implement Black-Scholes formula as given in C-Exercise 04
    d_1 = ( log(S_t./K) + (r+1/2*sigma^2)*(T-t) ) ./ ( sigma*sqrt(T-t) );
    d_2 = d_1 - sigma*sqrt(T-t);
    C = K.*exp(-r*(T-t)).*Phi(-d_2)-S_t.*Phi(-d_1) ;
    phi = Phi(d_1);
    
endfunction

// Put Black-Scholes function for test parameters.
C = BS_EuPut(0, S, r, sigma, T, K);

// Comparison of the explicit scheme and the Black-Scholes formula. 
plot (S, V0,'.');
plot(S, C, 'r');
legend('Explicit finite difference scheme', 'BS formula')
