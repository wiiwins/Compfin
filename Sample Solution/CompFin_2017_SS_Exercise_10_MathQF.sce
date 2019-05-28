// C-Exercise 10

// Ensure that the function CRR_AmPut from C-Exercise 06 is known
// to scilab.
exec("CompFin_2017_SS_QF_Exercise_06.sce")

function V0 = CRR_AmPut_Adapt (S0, r, sigma, T, K, M, epsilon)
   
    // Compute prices with M and 2M periods. 
    V_M = CRR_AmPut (S0, r, sigma, T, K, M);
    V_2M = CRR_AmPut (S0, r, sigma, T, K, 2*M);
    
    // Double number of periods as long as termination condition is not
    // satisfied.        
    while ( abs(V_M-V_2M) / V_M >= epsilon)
        M = 2*M;
        V_M = V_2M;
        V_2M = CRR_AmPut (S0, r, sigma, T, K, 2*M);
        disp("Price in Iteration:"+string(V_2M)+" Numeber of time Steps:"+string(M))
    end
    
    V0 = V_2M;
         
endfunction

// test parameters
S0 = 100;
r = 0.03;
sigma = 0.24
T = 3/4;
K = 95;
M = 5;
epsilon = 0.001;

funcprot(0);

// Display price for test parameters.
disp("The put price with adaptive step size control for the test parameters is given by: " + string(CRR_AmPut_Adapt (S0, r, sigma, T, K, M, epsilon)) );
