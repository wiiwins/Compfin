//C-Exercise19
//Jurain Kahl
//Nattawut Phanrattinon

function V0 = BS_EuCall_Laplace (St, r, sigma, T, K, R)
    
    function y = f_tilde(z)
        y = K^(1-z) / (z*(z-1));
    endfunction
    
    // Characteristic function of BS model
    function v = chi (u)
        v = exp( %i*u*(log(S0)+r*T) - (%i*u+u^2)*sigma^2/2*T );
        endfunction
    
    // Integrand for the Laplace transform 
    function v = integrand (u)
        v = exp(-r*(T))/%pi * real(f_tilde(R+%i*u)*chi(u-%i*R) );
    endfunction
    
     // Perform integration. To approximate the indefinite integral from 0 to +inf, we have to cut
    // the domain of integration. We use an adaptive approach:
    
    // threshold for relative error 
    epsilon = 0.0001; 
    
    // integration bounds
    N = 1; 
    
    // compute option prices with initial and integration bounds
    VN = intg(0, N, integrand);
    V2N = intg(0, 2*N, integrand);
    
    while (abs(VN - V2N)/V2N >= epsilon)
        VN = V2N;
        N = 2*N;
        V2N = intg(0, 2*N, integrand);
    end
        
    V0 = V2N;
    
    
endfunction

// test parameters
S0 = 100;
r = 0.05;
sigma = 0.2;
T = 1;
K = 100;
R = 1.1;

BS_EuCall_Laplace (S0, r, sigma, T, K, R)
