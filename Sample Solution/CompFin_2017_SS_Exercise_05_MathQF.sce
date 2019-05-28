// C-Exercise 06, SS 2017

funcprot(0);

//a)Writing function that computes European put price
    function V0 = BinMod_EuPut_CF (S0, r, sigma, T, M, K)
        
        // Compute values of u, d and p as in the lecture notes.
        delta_t = T/M;
        alpha = exp(r*delta_t);
        beta = 1/2 * ( 1/alpha + alpha*exp(sigma^2*delta_t) );
        u = beta + sqrt(beta^2-1);
        d = 1/u;
        p = ( exp(r*delta_t)-d ) / ( u-d );
        
        // Compute auxiliary values a and \tilde{p} as defined in the exercise.
        a = ceil( ( log(K/S0) - M*log(d) ) / log(u/d) );
        p_tilde = p * u / exp(-r*delta_t);
        
        // Implement put option price formula from the exercise.
        V0 =  K * exp(-r*T) * (cdfbin("PQ", a-1, M, p, 1-p)) - S0 * (cdfbin("PQ", a-1, M, p_tilde, 1-p_tilde));
        
    endfunction

//b)
    // Test parameters.
    S0 = 100;
    r = 0.05;
    sigma = 0.2;
    T = 1;
    M = 10:500;
    K = 100;
    
    // Allocate memory for put option prices.
    V0 = zeros(M);
    
    // Compute call option price for different numbers of steps in the binomial model.
    for j=1:length(M)
        V0(j) = BinMod_EuPut_CF (S0, r, sigma, T, M(j), K);
    end
    
    // Plot results.
    scf(2); clf();
    
    // Plot prices in dependence on M.
    subplot(2,1,1);
    plot(M, V0, 'b-');
    plot(M, mean(V0($-10:$))*ones(M), 'r-');
    title("Put option prices in the binomial model"); xlabel("number of steps"); ylabel("put option price");
    legend("Price in the binomial model", "Approximate B-S price", 4);

//c)[Only for students from mathematical programmes]
    // Allocate memory for errors compared to the final price.
    Err = zeros(M);
    mean_V0 = mean(V0($-10:$));
    
    for j=1:length(M)
        Err(j) = abs(V0(j)-mean_V0);
    end
    // Plot errors in dependence on M as log-log-plot.
    subplot(2,1,2);
    plot2d("ll", M, Err);
    title("Approximate errors of put option prices in the binomial model compared to the B-S price"); xlabel("number of steps"); ylabel("error");
    
    // Estimate order of convergence.
    disp("The slope of the line in the log-log-plot is approximately " + string( reglin(log(M),log(Err)) ) + "." );

