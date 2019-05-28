// C-Exercise 03
clear;
funcprot(0);

function lr = log_returns (data)
    
    // log(data) returns a vector with the natural logarithm applied to each component
    // of data. diff(x) for a vector x = (x_1,...,x_n) returns the vector (x_2-x_1,...,x_n-x_{n-1}).
    lr = diff(log(data));
    
endfunction

// import time series form file in current directory
ts_dax = csvRead('time_series_dax.csv', ";", ",","double" , [], [],[], 1);

// extract data from time series
s = ts_dax(:,2);
// as an alternative to s = s($:-1:1) one can also use s = flimdim(s,1)
s = s($:-1:1);



// apply function to imported time series
lr = log_returns(s);

// plot log-returns for DAX time series
scf(0);                // use figure number 0
clf();                 // clear this figure
plot(lr);     // plot log-returns (the x-axis is 1,...,length(log_returns))

// set title and labels for x- and y-axsis of the plot
title("log-returns of DAX in the period 04.01.2016 - 20.04.2017");
xlabel("trading day");
ylabel("log-return");

// compute mean and standard deviation (= root of the variance) of log-returns
ev = mean(lr);
std_dev = sqrt( variance(lr) );

// display annualized mean and standard deviation of log-returns
// ("+" concatenates two string variables, and string(x) converts a number variable x to a string
// variable)
disp("DAX log-returns: annualized mean = " + string(ev*250) + ", annualized standard deviation = " + string(std_dev*sqrt(250)) + ".");
