function V_0 = UpOutPut_BinMod (S_0, r, sigma, T, K, B, M)
    
    // Compute values of u, d and q according to Equations (2.4)--(2.7).
    delta_t = T/M;
    alpha = exp(r*delta_t);
    beta = 1/2 * ( 1/alpha + alpha*exp(sigma^2*delta_t) );
    u = beta + sqrt(beta^2-1);
    d = 1/u;
    q = ( exp(r*delta_t)-d ) / ( u-d );
    
    // auxiliary matrices to build the stock price matrix
    u_Matrix = repmat( (0:1:M)', 1, M+1 );
    d_Matrix = repmat( 0:1:M, M+1, 1 ) - u_Matrix;
    
    // S(j,i) is the stock price at time t_{i-1) (1 <= i <= M+1) if
    // j-1 (1<=j<=i) upward jumps occured
    S = S_0 * u.^u_Matrix .* d.^d_Matrix;
    
    // V(j,i) will in the end contain the price of the option at 
    // time t_{i-1) (1 <= i <= M+1) if j-1 (1<=j<=i) upward jumps occured, GIVEN
    // that the stock price has not been hit B at times t_0, ..., t_{i-1}.
    // Initialize first with -1es.
    V = -ones(M+1, M+1);
    
    // The prices of the put at time t_M are given by the exercise function only.
    // (Note that if the terminal price is above B, it is automatically also above K,
    // and hence nothing is payed out anyway.)
    V(:,M+1) = max(K-S(:,M+1), 0);
        
    // Compute the put prices at times t_i (0<=i<=M-1) via
    // backward recursion (cf. Section 2.3).
    for k=M:-1:1
        V(1:k,k) = exp(-r*delta_t) * ( q*V(2:k+1,k+1) + (1-q)*V(1:k,k+1) ) .* (S(1:k,k) < B) ;
    end
    
    // Return the price of the put at time t_0 = 0.
    V_0 = V(1,1);
        
endfunction

// test parameters
S_0 = 100;
r = 0.05;
sigma = 0.2;
T = 1;
K = 100;
B = 110;
M = 1000;

// display the price for the test parameters
disp("The price of the up-and-out put option for the test parameters is given by: " + string(UpOutPut_BinMod (S_0, r, sigma, T, K, B, M)) );
