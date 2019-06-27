//C-Exercise 36
//Jurian Kahl
//Phanrattinon Nattawut

function V0 = Heston_EuOption_MC(S0, r, gamma, T, g)

    delta_t = T/m;

    //Set initial price and create vector for ST and gamma.
    ST = S0  * ones(M, 1);
    gamma = gamma0 * ones(M, 1);
    
    a=0;  
    b=0;
    for k=1:m
        delta_W = grand(M, 1, 'nor', 0, sqrt(delta_t));
        gamma = gamma + (kappa - lambda*gamma)*delta_t + sigma_tilde*sqrt(gamma).*delta_W;
        a=a+r-gamma/2;
        b=b+gamma;
    end
    
    // Compute Stock prices from given equation
    X = grand(M, 1, 'nor', 0, 1);
    ST=S0*exp(a*delta_t + sqrt(delta_t*b).*X);
       
    // Compute payoff.
    VT =(exp(-r*T)*g(ST));
        
    // Compute Monte Carlo estimate.
    V0 = mean(VT);
    

endfunction

// test parameters
K = 90;
p=1.2;
S0 = 100;
r = 0.05;
T = 1;

//from C-Exercise33
gamma0 = 0.3^2;
kappa = 0.3^2;
lambda = 2.5;
sigma_tilde = 0.2;
M = 10000;
m = 100;
R = p + 0.1;


function y=g(x)
    y = max(x^p-K,0);
endfunction

//From Exercise20.sce

function V0 = Heston_PCall_Laplace (S0, r, nu0, kappa, lambda, sigma_tilde, T, K, R, p)
    
    function y = f_tilde(z)
        y = -K^(p-z) / (p-z) - K^(1-z) / z;
    endfunction
    
    function v = chi (u)
        d = sqrt(lambda^2+sigma_tilde^2*(%i*u+u^2));
        n = cosh(d*T/2) + lambda*sinh(d*T/2)/d;
        z1 = exp(lambda*T/2);
        z2 = (%i*u+u^2)*sinh(d*T/2)/d;
        v = exp( %i*u*(log(S0)+r*T)) * (z1/n)^(2*kappa/sigma_tilde^2) * exp(-gamma0*z2/n);
    endfunction
    
    function v = integrand (u)
        v = exp(-r*T)/%pi * real( f_tilde(R+%i*u)*chi(u-%i*R) );
    endfunction
    
    V0 = intg(0, 100, integrand);
    
endfunction

// display result 
V0 = Heston_EuOption_MC(S0, r, gamma, T, g);
disp("Price of European Call in Heston model : " + string(V0));

V0= Heston_PCall_Laplace (S0, r, nu0, kappa, lambda, sigma_tilde, T, K, R, p);
disp("Price of Power Call in Heston model using Laplace transform approach: " + string(V0));
