//Defining a function that computes prices of European calls in the B-S model using the FFT
function V0 = BS_EuCall_FFT (S0, r, sigma, T, K, R, N, M, kappa1)
    Delta=M/N;
    kappa=kappa1+((1:N)-1)*2*%pi/M;
    
    function y =g(u)
        function x = f_tilde_0(z)
            x = 1 / (z*(z-1));
        endfunction
        
        function v = chi (u)
            v = exp( %i*u*(log(S0)+r*T) - (%i*u+u^2)*sigma^2/2*T );
        endfunction       
        
        y=ones(1,length(u))
        for i=1:length(u)
            y(i)=f_tilde_0(R+%i*u(i))*chi(u(i)-%i*R)
        end
         
    endfunction
    
    x=g(((1:N)-1/2)*Delta)*Delta.*exp(-%i*((1:N)-1)*Delta*kappa1);
        
    x_hat=fft(x);
    
    V_kappa=1/%pi*exp(-r*T+(1-R)*kappa).*real(x_hat.*exp(-%i/2*Delta*kappa));
    
    V0=interpln([kappa; V_kappa],log(K));
   
endfunction

// test parameters
S0 = 100;
r = 0.05;
sigma = 0.2;
T = 1;
K = 80:1:180;
R = 1.1;
N = 2^11;
M = 50;
kappa1=log(80)

V1 = BS_EuCall_FFT(S0, r, sigma, T, K, R, N, M)


exec("CompFin_2019_SS_BS_Price_Int.sce")
V2=zeros(K);
j=1;
for i=K
    function y=g(x)
        y=max(0,x-i)
    endfunction
    V2(j)=BS_Price_Int (r, sigma, S0, T, g);
    j=j+1;
end


plot(V1)
plot(V2)
