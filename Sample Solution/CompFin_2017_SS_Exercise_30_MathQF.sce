function [V0, epsilon, epsilon2] = UpOutPut_BS_MC_Richardson (S0, r, sigma, T, K, B, M, m)
    
    // Time step on the fine grid.
    delta_t = T/(2*m);
    
    S_fine = S0 * ones(M, 1);
    S_coarse = S0 * ones(M, 1);
    no_barrier_hit_fine = ones(M,1);
    no_barrier_hit_coarse = ones(M,1);
    
    // Loop over points on coarse grid.
    for k=1:m
        
        // Simulate two increments of Brownian motion on the fine grid.
        delta_W_1 = grand(M, 1, 'nor', 0, sqrt(delta_t));
        delta_W_2 = grand(M, 1, 'nor', 0, sqrt(delta_t));
        
        // 1st Euler step on fine grid.
        S_fine = S_fine + r*S_fine*delta_t + sigma*S_fine.*delta_W_1;
        no_barrier_hit_fine = no_barrier_hit_fine .* (S_fine<B);
        
        // 2nd Euler step on fine grid.
        S_fine = S_fine + r*S_fine*delta_t + sigma*S_fine.*delta_W_2;
        no_barrier_hit_fine = no_barrier_hit_fine .* (S_fine<B);
        
        // Euler step on coarse grid.
        S_coarse = S_coarse + r*S_coarse*2*delta_t + sigma*S_coarse.*(delta_W_1+delta_W_2);
        no_barrier_hit_coarse = no_barrier_hit_coarse .* (S_coarse<B);

    end
    
    // Compute (discounted) payoffs for paths on fine and coarse grids.
    VT_fine = no_barrier_hit_fine .*  (exp(-r*T)*max( K - S_fine , 0));
    VT_coarse = no_barrier_hit_coarse .*  (exp(-r*T)*max( K- S_coarse , 0));
        
    // Compute Monte Carlo estimate.
    VT = 2*VT_fine - VT_coarse;
    V0 = mean(VT);
    
    // Compute radius of 95% confidence interval (not part of the exercise).
    epsilon = 1.96 * sqrt(variance(VT)/M);
               
endfunction

exec('./CompFin_2017_SS_QF_Exercise_12.sce')

// test parameters
S0 = 100;
r = 0.05;
sigma = 0.2;
T = 1;
K = 100;
B = 110;
M = 10000;
m = 250;

// display result for test parameters
[V0, epsilon] = UpOutPut_BS_MC_Richardson (S0, r, sigma, T, K, B, M, m);
disp("Price of up-and-out call by use of Monte-Carlo simulation with Euler scheme with Richardson extrapolation: " + string(V0) + ", radius of 95% confidence interval: " + string(epsilon));
