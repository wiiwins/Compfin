//C-Exercise 27 Computational Finance SS2017

function [X_exact, X_Euler, X_Milshtein]=Sim_Paths_GeoBM(X0,mu,sigma,T,N)
    Delta_t=T/N;
    Delta_W=grand(N, 1, 'nor', 0, sqrt(Delta_t));
    
    //Initialize vectors with starting value
    X_exact=X0*ones(N+1);
    X_Euler=X0*ones(N+1);
    X_Milshtein=X0*ones(N+1);
    
    //Recursive simulation according to the algorithms in Section 5.5 using identical Delta_W
    for i=1:N
        X_exact(i+1)=X_exact(i)*exp((mu-sigma^2/2)*Delta_t+sigma*Delta_W(i));
        X_Euler(i+1)=X_Euler(i)*(1+mu*Delta_t+sigma*Delta_W(i));
        X_Milshtein(i+1)=X_Milshtein(i)*(1+mu*Delta_t+sigma*Delta_W(i)+sigma^2/2*((Delta_W(i))^2-Delta_t));
    end
     
endfunction

//Parameters
X0=100;
mu=0.05;
sigma=0.2;
T=2;
N=10^(1:4);

scf(1)
clf(1)
for i=1:4
    [X_exact, X_Euler, X_Milshtein]=Sim_Paths_GeoBM(X0,mu,sigma,T,N(i));
    subplot(2,2,i)
    plot((0:N(i))*T/N(i), X_exact)
    plot((0:N(i))*T/N(i), X_Euler, 'red')
    plot((0:N(i))*T/N(i), X_Milshtein, 'green')
    xlabel('t')
    ylabel('X(t)')
    title('N='+string(N(i)))
    legend('Exact simulation','Euler approximation','Milshtein approximation')
end
