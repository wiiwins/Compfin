function V0 = EuCall_BS_ClosedForm(S0, r, sigma, T, K)
    
    // Cumulative distribution function of the standard normal
    // distribution.
    function p = Phi(x)
        p = cdfnor("PQ", x, zeros(x), ones(x));
    endfunction
    
    // Implement Black-Scholes formula as given in (3.23),
    // where K and T can be vectors.
    d_1 = ( log(S0./K) + (r+1/2*sigma^2)*T ) ./ ( sigma*sqrt(T) );
    d_2 = d_1 - sigma*sqrt(T);
    V0 = S0.*Phi(d_1) - K.*exp(-r*T).*Phi(d_2);
    
endfunction

function sigma = EuCall_BS_Calibrate (S0, r, T, K, V, sigma0)
    
    function e = g(s)
        e = V - EuCall_BS_ClosedForm(S0, r, s, T, K);
    endfunction
    
    [gopt, sigma] = leastsq(g, sigma0);    
        
endfunction

// Test data, Call option data

S0 = 12658;
r = 0;
sigma0 = 0.3;
dax_options = csvRead("./Dax_CallPrices_Eurex.csv",';',',');
K = dax_options(2:$,1);
V = dax_options(2:$,3);
T = dax_options(2:$,2);

// Calibrate to test data.
sigma_calib = EuCall_BS_Calibrate (S0, r, T, K, V, sigma0);
disp("Calibration to test data yields sigma = " + string(sigma_calib) );

// Plot fit of calibrated sigma.
scf(0); clf();
plot(K, V, 'ro');
plot(K, EuCall_BS_ClosedForm(S0, r, sigma_calib, T, K), 'b+');
xlabel ("strike"); ylabel("market prices / model prices");
legend("market price", "model price");
title("market prices vs. model prices");

// Plot implied volatilities of market prices
n_strike=11;
n_T=5;
sigma_impl = zeros(n_strike,n_T);
for i=1:n_strike
    for j=1:n_T
        sigma_impl(i,j) = EuCall_BS_Calibrate (S0, r, T(n_strike*j), K(i), V(n_strike*(j-1)+i), sigma0);
    end
end

scf(1); clf();
surf(T(n_strike*(1:n_T)),K(1:n_strike),sigma_impl);
c=jetcolormap(50)
gcf().color_map = c;
gce().thickness=0;
label=("Implied volatility surface");
xtitle(label,"maturity T", "Strike K","implied volatility")
