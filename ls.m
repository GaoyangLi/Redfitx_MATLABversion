function sumsq=ls(a,t,x)
    dt=diff(t,1);
    x2=x(2:end);
    x1=x(1:end-1);
    sumsq=sum((x2-x1*sign(a).*abs(a).^dt).^2);
end