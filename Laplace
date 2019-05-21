function V0 = Heston_EuCall_Laplace (St, r, gamma_t, kappa, lambda, sigma_tilde, T, t, K, R)
    
    // Laplace transform of the function f(x) = (e^x - K)^+,
    // cf. (4.6).
    function y = f_tilde(z)
        y = K^(1-z) / (z*(z-1));
    endfunction
    
    // Characteristic function of log(S(T)) in the Heston model,
    // cf. (4.8).
    function v = chi (u)
        d = sqrt(lambda^2+sigma_tilde^2*(%i*u+u^2));
        n = cosh(d*(T-t)/2) + lambda*sinh(d*(T-t)/2)/d;
        z1 = exp(lambda*(T-t)/2);
        z2 = (%i*u+u^2)*sinh(d*(T-t)/2)/d;
        v = exp( %i*u*(log(St)+r*(T-t))) * (z1/n)^(2*kappa/sigma_tilde^2) * exp(-gamma_t*z2/n);
    endfunction
    
    // Integrand for the Laplace transform method (cf. (4.5)).
    function v = integrand (u)
        v = exp(-r*(T-t))/%pi * real( f_tilde(R+%i*u)*chi(u-%i*R) );
    endfunction
    
    // upper bound for integration
    ub = 40;
    
    V0 = intg(0, ub, integrand);
    
endfunction

funcprot(0);

// test parameters
St = 100;
r = 0.05;
gamma_t = 0.2^2;
kappa = 0.5;
lambda = 2.5;
sigma_tilde = 1;
T = 1;
t=0;
K = 100;
R = 3;

// display result for test parameters
disp("Price of European call by use of Laplace transform approach: " + string(Heston_EuCall_Laplace (St, r, gamma_t, kappa, lambda, sigma_tilde, T, t, K, R)) );
