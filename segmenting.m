function [indices,theo,rr]=segmenting(t,nseg,ratio)
    T=t(end)-t(1);
    dt=T/(1/(1-ratio)+(nseg-1));
    Tseg=dt/(1-ratio);
    
    theo=zeros(2,nseg);
    theo(1,:)=[t(1):dt:(nseg-1)*dt+1];
    theo(2,:)=theo(1,:)+Tseg;
    
    indices=zeros(2,nseg);
    j=0;
    k=0;
    for i=1:nseg;
        m=k+1;
        k=lookup(m,i,1,t,theo);
        indices(1,i)=k;
        
        m=j+1;
        j=lookup(m,i,2,t,theo);
        indices(2,i)=j;  
    end;
    for i=2:nseg;
       rr(i-1)=1-(indices(1,i)-indices(1,i-1))/(indices(2,i-1)-indices(1,i-1));
    end;
end

function j=lookup(m,i,row,t,theo)
   for j=m:length(t);
       if t(j)>theo(row,i);
           break;
       end;
   end;
   if abs(t(j-1)-theo(row,i))<(t(j)-theo(row,i));
       j=j-1;
   end;
end