function [res,tau,distribution]=tauest(t,x,nsim)
    %%%% Written by Gaoyang Li, 23/Nov./2017;  
    
    %%%% For an unevenly sampled first-order auto-regressive time series: x(1)~N(0,1);
    %%%% x(i+1)=x(i)*exp(-(t(i+1)-t(i))/tau)+eps  (*),   eps~N(0,exp(-2*(t(i+1)-t(i))/tau)), 
    %%%%this program estimates the parameter tau in equation *.
    
    %%%%INPUTS: sampled time t, corresponding value x, numbers of Monte
    %%%%Carlo Simulation runs nsim (which will determine a distribution of
    %%%%the estimation of tau);
    %%%%OUTPUTS: results matrix res: 
    %%%%                     col. 1:  x(i)
    %%%%                     col. 2:  modelled x(i+1)'=x(i)*exp(-(t(i+1)-t(i))/tau)
    %%%%                     col. 3:  residue: x(i+1)-x(i+1)
    %%%%        scalar estimated tau;
    %%%%        vector distribution: Tau values computed in Monte Carlo Simulation.
    
    %%%%For more information and theoretical details, please refer to 
    %%%% "Mudelsee, M. (2002). TAUEST: a computer program for estimating persistence 
    %%%% in unevenly spaced weather/climate time series. Pergamon Press, Inc."
    %%%%  Also refer to redfit-x.f90 written by Kristin B. Olafsdottir, Michael Schulz 
    %%%%and Manfred Mudelsee for details of the fortran implementation.
    
    x=zscore(x);
    x=detrend(x);
    x1=x(1:end-1);
    x2=x(2:end);
    
    rho=rhoest(x);   %%%% rough estimate of rho as if the series were evenly sampled; 
    if rho<0;
        rho=0.05;
        disp('Warning: rho smaller than 0 !');
    elseif rho>=1;
        rho=0.95;
        disp('Warning: rho larger than 1 !');
    end;
    dt1=(t(end)-t(1))/(length(t)-1);
    scalet=-log(rho)/dt1;    %%%% scaling factor of time; see reference and redfit-x.f90 for detail
    tt=t*scalet;
    dt=diff(t,1);
    dtt=diff(tt,1);
    
    a=[0.001:0.001:0.999];
    for i=1:length(a)
        Sa(i)=sum((x2-x1.*a(i).^(dtt)).^2);
    end;
    if nsim>0;
        plot(a,Sa);
    end;
    
    a_ar1=exp(-1);
    tol=3e-8;
    [xmin1,fxmin1]=brent(-2,a_ar1,2,tol,tt,x);   
    [xmin2,fxmin2]=brent(a_ar1,0.5*(a_ar1+1.0),2,tol,tt,x);
    [xmin3,fxmin3]=brent(-2,0.5*(a_ar1-1.0),a_ar1,tol,tt,x);
    
    [tmp,I]=sort([fxmin1,fxmin2,fxmin3],'ascend');
    counter=0;
    MIN=[xmin1,xmin2,xmin3];
    for i=1:3;
        if MIN(I(i))>0 && MIN(I(i))<1 && MIN(I(i))~=a_ar1;
            xmin4=MIN(I(i));
            break;
        end;
    end;
    
    NumOfSolutions=length(xmin4);
    if NumOfSolutions==1;
        tau_prim=-1/(scalet*log(xmin4));
        if nsim>0
            disp(['The uncorrected tau is ',num2str(tau_prim)]);
        end;
        x2model=x1.*exp(-dt./tau_prim);
        residual=x2-x2model;
        res=[x1,x2model,residual];
        rho_corr=rho+(1+3*rho)/(length(t)-1);
        tau=-dt1/log(rho_corr);
        tau2=max([tau_prim,tau]);
        if nsim>0;
            distribution=MCtau(tau2,t,nsim);
        else
            distribution=0;
        end;
        tau_med=median(distribution);
        if nsim>0;
            disp(['Median tau from MC simulation: ',num2str(tau_med),'Corrected tau: ',num2str(tau)]);
        end;
    elseif NumOfSolutions==0;
        res=0;
        rho_corr=rho+(1+3*rho)/(length(t)-1);
        tau=-dt1/log(rho_corr);
        if nsim>0;
            disp(['No solution...',num2str(tau)]);
        end;
        distribution=0;
    end;
end

function rho=rhoest(x)
    sum1=sum(x(1:end-1).*x(2:end));
    sum2=sum(x(1:end-1).^2);
    rho=sum1/sum2;
end

function samples=MCtau(tau,t,nsim)
  L=length(t);
  samples=[];
  Eps=zeros(nsim,L-1);
  for i=1:nsim
    %if mod(i,100)==0;
       disp(['Simulations: ', num2str(i), '/', num2str(nsim)]);
    %end;
    for j=2:L;
        pd=makedist('Normal','mu',0,'sigma',sqrt(1-exp(-2*(t(j)-t(j-1))/tau)));
        Eps(:,j)=random(pd,nsim,1);
    end;
   
    x=TAUrednoise(t,tau,Eps(i,:));
    [tmp1,a,tmp]=tauest(t,x,0);
    samples=[samples,a];
  end;
end