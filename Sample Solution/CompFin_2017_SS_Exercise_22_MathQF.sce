// C-Exercise 22, SS 2017

// Density of the Beta distribution.
function y = Beta_density (x, alpha1, alpha2)
    y = 1/beta(alpha1,alpha2) * x^(alpha1-1) .* (1-x)^(alpha2-1) .* (x>=0) .* (x<=1);
endfunction

function X = Sample_Beta_AR (alpha1, alpha2, N)
    
    // Mode of the density of the Beta distribution.
    x_max = (alpha1-1) / (alpha1+alpha2-2);
    
    // Constant C for the acceptance/rejection method.
    C = Beta_density(x_max, alpha1, alpha2);

    // Generate one sample by the acceptance/rejection method.
    function x = SingleSample()
    
        success = %F
        while ~success
            U = rand(2,1);
            success = ( C*U(2) <= Beta_density(U(1), alpha1, alpha2) );
        end
        x = U(1);
        
    endfunction
    
    // Generate N samples.
    X = zeros(N,1);
    for n=1:N
        X(n) = SingleSample();
    end
    
endfunction

// Test parameters.
alpha1 = 2;
alpha2 = 3;
N = 2000;

// Generate N samples.
X = Sample_Beta_AR (alpha1, alpha2, N);

// Plot histogram of samples and density of the Beta distribution.
scf(0); clf();
histplot(50,X);
// Divide density by the number of bins such that they have the same scale.
x = 0:0.01:1; plot(x, 1/50*Beta_density(x,alpha1,alpha2));
legend("Acceptance/rejection method", "True density *1/50");
xlabel("x"); ylabel("Frequency / density");
title("Illustration of the acceptance/rejection method to sample from the Beta distribution");


