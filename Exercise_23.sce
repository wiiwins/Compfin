//C-Exercise23
//Jurian Kahl
//Nattawut Phanrattinon

funcprot(0)
function V0 = Heston_PCall_Laplace (S0, r, nu0, kappa, lambda, sigma_tilde, T, K, R, p)
        
    //laplace transform for f(x)=(e^(xp) -K)^+ according to page 39
    function y = f_tilde(z)
        y = K^(1-z/p) / ((z-p)*z);
    endfunction
    
    //the characteristic function given by (4.8)
    function y = chi(u)
        d = sqrt(lambda^2 + sigma_tilde^2 * (%i*u + u^2));
        a = cosh(d*T/2) + lambda*sinh(d*T/2)/d;
        e1 = exp(lambda *T/2);
        e2 = (-nu0)*(%i* u + u^2)*sinh(d*T/2) / d;
        
        y = exp(%i*u*log(S0)+r*T) * (e1 / a)^(2*kappa / (sigma_tilde^2)) * exp(e2 / a);
    endfunction
    
    //the integrand for the laplace transform given by (4.5)
    function y = integrand(u)
        y = (exp(-r*T) / %pi) * real(f_tilde(R+%i*u) * chi(u-%i*R));
    endfunction
    
    //the fair price at time 0
    V0 = intg(0, 50, integrand);
endfunction

//test parameter
S0 = (50:150);
r = 0.05; 
nu0 = 0.3^2;
sigma = nu0;
kappa = 0.3^2;
lambda = 2.5;
sigma_tilde = 0.2;
T = 1;
K = 100;
p = 2;
R = p + 0.1;

//testing for Heston_PCall_Laplace
V0 = zeros(S0);
i=1
for n=S0
    V0(i) = Heston_PCall_Laplace (S0(i), r, nu0, kappa, lambda, sigma_tilde, T, K, R, p);
    i=i+1;
end

//the BS formula
function V0 = BS_Price_Int (S0, r, sigma, T, g)
      
    // Specify integrand for integration formula.
    function y = integrand (x)
        y = 1/sqrt(2*%pi) * g( S0*exp((r-0.5*sigma^2)*T + sigma*sqrt(T)*x) ) * exp(-r*T) * exp(-x^2/2);
    endfunction
    
    // Perform integration. To approximate the indefinite integral from -inf to +inf, we have to cut
    // the domain of integration. We use an adaptive approach:
    
    // threshold for relative error 
    epsilon = 0.0001; 
    
    // integration bounds
    N = 1; 
    
    // compute option prices with initial and doubled integration bounds
    VN = intg(-N, N, integrand);
    V2N = intg(-2*N, 2*N, integrand);
    
    while (abs(VN - V2N)/V2N >= epsilon)
        VN = V2N;
        N = 2*N;
        V2N = intg(-2*N, 2*N, integrand);
    end
        
    V0 = V2N;
    
endfunction

//using the bs formula with the parameters
V1 = zeros(S0);
i=1;
for j=S0
    function y = g(x)
        y = max(0,x^p-K);
    endfunction
    V1(i) = BS_Price_Int (S0(i), r, sigma, T, g);
    i = i+1;
end

plot(S0, V0)
plot(S0, V1, 'red')
title("Comparison of the Laplace transform and Black-Scholes integration formula")
xlabel("S0")
ylabel("value of option at time 0")
legend("Laplace", "BS integration", 4)
