function x=TAUrednoise(t,tau,eps)
    rng('shuffle');
    L=length(t);
    x=zeros(L,1);
    x(1)=randn(1);
    
    for j=2:L;
       x(j)=x(j-1)*exp(-(t(j)-t(j-1))/tau)+eps(j-1);
    end; 
end