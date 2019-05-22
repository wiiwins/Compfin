funcprot(0);
//computes the initial price of European call options with identicall maturity T and strikes K = (K1, ..., Kn)
function V0 = BS_EuCall_FFT (S0, r, sigma, T, K, R, N, M, kappa1)
    //model parameters delta and kappa (4.15)
    delta = M/N;
    kappa = kappa1 + ((1:N)-1) *((2*%pi) / M );
    
    //function f_tilde according to (4.11) for kappa = 0
    function y = f_tilde(z)
        y = 1/ (z*(z-1));       
    endfunction
    
    //characteristic function at time 0 of the Black-Scholes model given on page 41
    function y = chi(u)
        y = exp(%i*u*(log(S0)+r*T) -(%i*u+u^2)*((sigma^2)/2)*T);
    endfunction
    
    //function g given by (4.13)
    function y = g(u)
        y = f_tilde(R+%i*u) * chi(u-%i*R);
    endfunction
    
    //setting x according to (4.17)
    x = ones(1,N)
    for n=1:N
        x(n)= g((n-1/2)*delta) * delta * exp(-%i*(n-1)*delta * kappa1);
    end
    
    //the fourier transform of x
    x_hat = fft(x);
    
    //computing (4.18)
    V_kappa = exp(-r*T + (1-R)*kappa)*(1/%pi).*real(x_hat.*exp(-(%i/2)*delta*kappa));
    
    //computing the price of the option using linear interpolation
    V0 = interpln([kappa;V_kappa],log(K));
endfunction

//testing
S0 = 100;
r = 0.05;
sigma = 0.2;
T = 1;
K = (80:180);
R = 1.1;
N = 2^11;
M = 50;
kappa1 = log(80);

V1 = BS_EuCall_FFT (S0, r, sigma, T, K, R, N, M, kappa1)
