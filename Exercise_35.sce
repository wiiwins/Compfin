//C-Exercise 35
//Jurian Kahl
//Phanrattinon Nattawut


function [X_exact, X_Euler, X_Milshtein] = Sim_Paths_GeoBM(X0, mu, sigma, T, N)
    delta_t=T/N;
    delta_W=grand(N, 1, 'nor', 0, sqrt(delta_t));
    
    //set initial value for each method 
    X_exact=X0*ones(1);
    X_Euler=X0*ones(1);
    X_Milshtein=X0*ones(1);
    
    for i=1:N
        X_exact(i+1)= X_exact(i)*exp((mu-0.5*sigma^2)*delta_t + sigma*delta_W(i));
        
        X_Euler(i+1)=X_Euler(i) + X_Euler(i)*mu*delta_t + X_Euler(i)*sigma*delta_W(i);
        
        X_Milshtein(i+1)=X_Milshtein(i) + X_Milshtein(i)*mu*delta_t + X_Milshtein(i)*sigma*delta_W(i)+0.5*sigma^2*X_Milshtein(i)*((delta_W(i))^2-delta_t);
    end
     
endfunction

//Set Parameters
X0=100;
mu=0.1;
sigma=0.3;
T=1;
N(1)=10;
N(2)=100;
N(3)=1000;
N(4)=10000;

//Plot
clf()
for i=1:4
    [X_exact, X_Euler, X_Milshtein]=Sim_Paths_GeoBM(X0,mu,sigma,T,N(i));
    subplot(2,2,i);
    plot(X_exact,'red');
    plot(X_Euler, 'blue');
    plot(X_Milshtein, 'k--');
    xlabel('T/delta_t');
    ylabel('X(t)');
    title('N ='+string(N(i)));
    legend('Exact Solution','Euler method','Milshtein method');
end
