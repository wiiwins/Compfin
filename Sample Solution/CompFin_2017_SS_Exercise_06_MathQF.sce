function [V_0,t] = CRR_AmPut (S_0, r, sigma, T, K, M)
    
    // Compute values of u, d and q according to Equations (2.4)--(2.7).
    delta_t = T/M;
    alpha = exp(r*delta_t);
    beta = 1/2 * ( 1/alpha + alpha*exp(sigma^2*delta_t) );
    u = beta + sqrt(beta^2-1);
    d = 1/u;
    q = ( exp(r*delta_t)-d ) / ( u-d );
    
    S=ones(M+1,M+1);
    
    // auxiliary matrices to build the stock price matrix
    for i=1:M+1
        for j=1:i
            //S(j,i) is the stock price at time t_{i-1) (1 <= i <= M+1) if
            // j-1 (1<=j<=i) upward jumps occured
            S(j,i) = S_0*u^(j-1)*d^(i-j)
        end    
    end

    // V(j,i) will in the end contain the price of the American put at 
    // time t_{i-1) (1 <= i <= M+1) if j-1 (1<=j<=i) upward jumps occured.
    // Initialize first with -1es.
    V = -ones(M+1, M+1);
    
    // The prices of the put at time t_M are given by the exercise function only.
    V(:,M+1) = max(K-S(:,M+1), 0);
        
    // Compute the put prices at times t_i (0<=i<=M-1) via the Snell envelope by
    // backward recursion (cf. Section 2.3).
    for k=M:-1:1
        V(1:k,k) = max( K-S(1:k,k), exp(-r*delta_t) * ( q*V(2:k+1,k+1) + (1-q)*V(1:k,k+1) ) );
    end
    
    // Return the price of the put at time t_0 = 0.
    V_0 = V(1,1);
        
endfunction

// test parameters
S_0 = 100;
r = 0.03;
sigma = 0.24;
T = 3/4;
K = 95;
M = 500;

// apply the function
tic()
V_0 = CRR_AmPut(S_0, r, sigma, T, K, M)
toc()
// display the price for the test parameters
disp("The price of the American put option for the test parameters is given by: " + string(V_0) );
