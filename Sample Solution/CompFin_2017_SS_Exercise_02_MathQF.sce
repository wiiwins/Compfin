// C-Exercise 02
clear;
funcprot(0);

function [y_hat, theta_hat] = cubic_regression (x, y)
    
    // set matrix X as in C-Exercise 02 (x.^k is the component-wise k-th power of the vector x, 
    // ones(x) is a vetor of 1es with the same length as x)
    X = [x.^3, x.^2, x, ones(x)];
    
    // set theta_hat as in C-Exercise 02 (X' is the transpose of X, inv(Z) is the inverse of
    // a matrix Z, * performs a matrix multiplication)
    //theta_hat = inv(X'*X) * X' * y;
    theta_hat = (X'*X)\(X' * y);
    
    // set y_hat as in C-Exercise 02
    y_hat = X * theta_hat;
    
endfunction

// test data as in C-Exercise 02
x = [0; 1; 2; 3; 4];
y = [1; 0; 3; 5; 8];

// call function with test data
[y_hat, theta_hat] = cubic_regression(x, y);

// define the regression function
function y=cubic_function(x, theta_hat)
    y=theta_hat(1)*x^3+theta_hat(2)*x^2+theta_hat(3)*x+theta_hat(4) 
endfunction

// plot x against y and y_hat
scf(0);                  // use figure number 0
clf();                   // clear this figure
plot(x, y, 'ro');        // plot x against y using red circles
x2=min(x):0.01:max(x);
plot(x2, cubic_function(x2), 'blue');    // plot the regression function as a blue line
// set title, labels for x- and y-axis, and define a legend for the plot
title("Data and regression function");
xlabel("x");
ylabel("y");
legend("data", "reg. function","in_upper_left");
