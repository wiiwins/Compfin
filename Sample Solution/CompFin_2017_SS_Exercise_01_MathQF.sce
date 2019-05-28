// C-Exercise 01
clear;
funcprot(0);

function Vn = capital (V0, r, n, c)
    
    // simple rate
    if (c == 0) then    
        Vn = V0 * (1+r)^n;
        
    // continuous rate
    elseif (c == 1) then    
        Vn = V0 * exp(r*n);
        
    // wrong value for parameter s   
    else
        error("Error: Argument for compound type must be 0 (continuous) or 1 (simple).");
    end
    
endfunction

// test parameters as in C-Exercise 01
V0 = 1000;
r = 0.05;
n = 10;
c = 0;

// call function with test parameters
Vn = capital(V0, r, n, c)

disp(Vn)
