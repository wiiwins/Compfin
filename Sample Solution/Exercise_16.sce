
//a)
exec("./CompFin_2019_SS_BS_Price_Int.sce")
function [Delta, vega, gamma] = BS_Greeks_num(r, sigma, S0, T, g, eps)
    
    //price function for initial parameters 
    V0 = BS_Price_Int(r, sigma, S0, T, g);
    V_Delta = BS_Price_Int(r, sigma, (1+eps)*S0, T, g);
    V_vega = BS_Price_Int(r, (1+eps)*sigma, S0, T, g);
    V_gamma = BS_Price_Int(r, sigma, (1+eps)*S0, T, g);
    V_gamma2 = BS_Price_Int(r, sigma, (1-eps)*S0, T, g)

    //Compute greeks num
    Delta = (V_Delta-V0)/(eps*S0)
    vega = (V_vega-V0)/(eps*sigma)
    gamma = (V_gamma-2*V0+V_gamma2)/((eps*S0)^2)
    
endfunction
    
//b)
r=0.03;
sigma=0.2;
T=1;
S0=60:1:140;
eps=0.001;

function y = g(x)
    y = max(x-100, 0)
endfunction


Delta = ones(S0);
vega = ones(S0);
gamma = ones(S0);

for i = 1:length(S0)
    [Delta(i), vega(i), gamma(i)] = BS_Greeks_num(r, sigma, S0(i), T, g, eps);
end


plot(S0, vega)
xlabel("S0")
ylabel("Delta")


