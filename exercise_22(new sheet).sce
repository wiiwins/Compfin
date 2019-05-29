//Exercise22
//Jurain Kahl
//Phanrattinon Nattawut

//define Function Black-Sholes price

function V0 = EuCall_BS(S0, r, sigma, T, K)
    
    function p = Phi(x)
        p = cdfnor("PQ", x, zeros(x), ones(x));
    endfunction
    
    d1 = ( log(S0./K) + (r+0.5*sigma^2)*T ) ./ ( sigma*sqrt(T) );
    d2 = d1 - sigma*sqrt(T);
    V0 = S0.*Phi(d1) - K.*exp(-r*T).*Phi(d2);
    
endfunction

// Define Calibrate Function
function sigma = EuCall_BS_Calibrate (S0, r, T, K, V, sigma0)
    
    function e = g(s)
        e = V - EuCall_BS(S0, r, s, T, K);
    endfunction
    
    [xopt, sigma] = leastsq(g, sigma0);    
        
endfunction

// Test data & Call option data

S0 = 12658;
r = 0;
sigma0 = 0.3;
data = csvRead("./Dax_CallPrices_Eurex.csv",';',',');
K = data(2:$,1);
V = data(2:$,3);
T = data(2:$,2);

// Calibrate data
sigma = EuCall_BS_Calibrate (S0, r, T, K, V, sigma0);
disp("Calibration to test data yields sigma = " + string(sigma) );

// Plot 
scf(0); clf();
plot(K, V, 'b+');
plot(K, EuCall_BS(S0, r, sigma_calib, T, K), 'ro');
xlabel ("strike"); ylabel("market prices / model prices");
legend("market price", "model price");
title("market prices vs. model prices");
