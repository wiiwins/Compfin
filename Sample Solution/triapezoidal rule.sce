function [r] =Trapezoidal(a,b,f,n,r) 
%TRAPEZOIDALNumericalintegrationfromatobusingtrapezoids
function [f] = fcnchk(f); 
    h=(b-a)/n; 
    s=0; 
    x=a;
for i=1:n-1 
    x=x+h; 
    s=s+feval(f,x); 
end
 s=0.5*(feval(f,a)+feval(f,b))+s; 
 r=h*s;
