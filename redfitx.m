function [Gxx,Gyy,Gxy,Cxy,Pimagx,Prealx,Pimagy,Prealy,Phi_xy,ciGXX,ciGYY,ciCXY,ciPhi]=redfitx(tx,x,ty,y,nseg,ratio,window,beta,nsim)
    %%%% Written by Gaoyang Li, 26/Nov./2017;  
    
    %%%% A matlab "translation" of redfit-x.f90 written by Kristin B. Olafsdottir, Michael Schulz 
    %%%% and Manfred Mudelsee. The code computes the cross-spectrum  Gxy, coherency Cxy, quadrature 
    %%%% spectrum Phi_xy between unevenly sampled time series x and y, along with confidence 
    %%%% intervals determined from Monte Carlo Simulation.
    
    %%%% NOTE: MCS for quadurature spectrum is highly time-consuming, so this part of the code has been 
    %%%% commented out and replaced by a theorectical estimation of the confidence interval for evenly-spaced 
    %%%% data.
    
    %%%% For details of theoretical background and formulas, refer to:
    %%%% "?lafsd¨®ttir, K. B., Schulz, M., & Mudelsee, M. (2016). Redfit-x: cross-spectral analysis of 
    %%%% unevenly spaced paleoclimate time series. Computers & Geosciences, 91, 11-18."
    %%%% and to:
    %%%% Scargle, 1989; The Astrophysical Journal 343, 874-887;
    
    %%%% INPUTS: 
    %%%%        Vectors tx, ty: time at which x and y values were
    %%%%        observed;
    
    %%%%        Vectors x, y: observations;
    
    %%%%        Scalar nseg: number of overlapping segments for spectrum
    %%%%        estimation (recommended: 8);
    
    %%%%        Scalar ratio: proportion of overlapping (recommended: 0.5 and 0.75)
    
    %%%%        Scalar window: window type, 0 for hanning, 1 for Kaiser;
    %%%%        Hanning is sufficient for most purposes; however, if extra
    %%%%        accuracy is needed for locating the peak (i.e. a narrower
    %%%%        main lobe, at the cost of less side-lobe attenuation), then 
    %%%%        Kaiser window with small beta (i.e. 4 or 5) is
    %%%%        recommended; if extra accuracy is needed for amplitudes 
    %%%%        (i.e. more suppresed side lobes, at the cost of wider main lobe),
    %%%%        then Kaiser window with larger beta (i.e. 7,8,9) 
    
    %%%%        Scalar beta: parameter for the Kaiser window;
    
    %%%%        Scalar nsim: number of MCS runs; recommended: 1000~2000
    
    %%%% OUTPUTS:
    %%%%        Vectors Gxx, Gyy, Gxy, Cxy, Phi_xy:
    %%%%                Scargle-Lomb spectrum of x, y, cross spectrum, 
    %%%%                coherency and quadruature spectrum between them.
    %%%%        Vectors Prealx, Pimagx, Prealy, Pimagy:
    %%%%                real and imaginary parts of x or y;
    %%%%        Matrix  ciGXX, ciGYY, ciCXY:
    %%%%                Confidence intervals determined from MCS;
    %%%%                col.1: 90%, col.2: 95%; col.3: 99%;
    %%%%                as GXX, GYY conform to chi-squared distribution,
    %%%%                CXY to F-distribution, so here one-sided test is
    %%%%                used (similar as in ANOVA).
    %%%%       Vector   ciPhi:
    %%%%                99% uncertainty of phase difference.
    
    tx=tx-tx(1);
    ty=ty-ty(1);
    [indices_x,~,rr_x]=segmenting(tx,nseg,ratio);
    [indices_y,~,rr_y]=segmenting(ty,nseg,ratio);
    
    nx=(indices_x(2,:)-indices_x(1,:))';
    ny=(indices_y(2,:)-indices_y(1,:))';
    dtx=tx(indices_x(2,:))-tx(indices_x(1,:));
    dty=ty(indices_y(2,:))-ty(indices_y(1,:));
    
    fl_x=max(1./(dtx));
    fl_y=max(1./(dty));
    fh_x=min(nx./(2*dtx));
    fh_y=min(ny./(2*dty));
    
    fl_xy=max([fl_x,fl_y]);
    fh_xy=min([fh_x,fh_y]);
    
    fx=[fl_x:fl_x/4:fh_x]';
    fy=[fl_y:fl_y/4:fh_y]';
    fxy=[fl_xy:fl_xy/4:fh_xy]';
    
    [Gxx,Gyy,Gxy,Cxy,Pimagx,Prealx,Pimagy,Prealy,Phi_xy,ciGXX,ciGYY,ciCXY,ciPhi]=main(tx,ty,fx,fy,fxy,indices_x,indices_y,x,y,nseg,ratio,window,beta,nsim);
end

function [Gxx2,Gyy2,Gxy,Cxy,Pimagx,Prealx,Pimagy,Prealy,Phi_xy,ciGXX,ciGYY,ciCXY,ciPhi]=main(tx,ty,fx,fy,fxy,indices_x,indices_y,x,y,nseg,ratio,window,beta,nsim)
    %[fx,fy,fxy,indices_x,indices_y]=redfitx(tx,ty,nseg,ratio);
    nfx=length(fx);
    nfy=length(fy);
    nfxy=length(fxy);
    fl_x=fx(1);
    fl_y=fy(1);
    fl_xy=fxy(1);
    nx=indices_x(2,:)-indices_x(1,:)+1;
    ny=indices_y(2,:)-indices_y(1,:)+1;
    
    Px=zeros(nfx,nseg);
    Py=zeros(nfy,nseg);
    Px2=zeros(nfxy,nseg);
    Py2=zeros(nfxy,nseg);
    Prealx=zeros(nfxy,nseg);
    Pimagx=zeros(nfxy,nseg);
    Prealy=zeros(nfxy,nseg);
    Pimagy=zeros(nfxy,nseg);
    cx=zeros(1,nseg);
    cy=zeros(1,nseg);
    
    for i=1:nseg
        Rx=indices_x(1,i):1:indices_x(2,i);
        Ry=indices_y(1,i):1:indices_y(2,i);
        ttx=tx(Rx);
        tty=ty(Ry);
        px=(ttx-ttx(1))/(ttx(end)-ttx(1));
        py=(tty-tty(1))/(tty(end)-tty(1));
        Wx=unevenWindow(px,window,beta);
        Wy=unevenWindow(py,window,beta);
        
        [~, Px(:,i), ~, ~, probx(:,i)]=SLspectrum(ttx,x(Rx).*Wx,fx);
        [~, Px2(:,i), Prealx(:,i),Pimagx(:,i),~]=SLspectrum(ttx,x(Rx).*Wx,fxy);
        Px(:,i)=Px(:,i)/nx(i);
        Px2(:,i)=Px2(:,i)/nx(i);
        Prealx(:,i)=Prealx(:,i)/sqrt(nx(i));
        Pimagx(:,i)=Pimagx(:,i)/sqrt(nx(i));
        
        [~, Py(:,i), ~, ~, ~]=SLspectrum(tty,y(Ry).*Wy,fy);
        [~, Py2(:,i), Prealy(:,i),Pimagy(:,i),~]=SLspectrum(tty,y(Ry).*Wy,fxy);
        Py(:,i)=Py(:,i)/ny(i);
        Py2(:,i)=Py2(:,i)/ny(i);
        Prealy(:,i)=Prealy(:,i)/sqrt(ny(i));
        Pimagy(:,i)=Pimagy(:,i)/sqrt(ny(i));
        
        [cx(i),cy(i)]=EffNseg(tx(Rx),ty(Ry),Wx,Wy,window,ratio,beta);
    end;
    
    neffx=nseg./(1+2*cx.^2-2*cx.^2/nseg);
    neffy=nseg./(1+2*cy.^2-2*cy.^2/nseg);
    neff=min([neffx;neffy]);
    
    Gxx=2/nseg/fl_x*sum(Px,2);
    Gyy=2/nseg/fl_y*sum(Py,2);
    
    Gxx2=2/nseg/fl_xy*sum(Px2,2);
    Gyy2=2/nseg/fl_xy*sum(Py2,2);
    
    Gxy=2/nseg/fl_xy*sum((Prealx+sqrt(-1)*Pimagx).*(Prealy-sqrt(-1)*Pimagy),2);
    Cxy=abs(Gxy).^2./(Gxx2.*Gyy2);
   
    Cxy=Cxy-(1-Cxy).^2/mean(neff);
    
    %%%%phase spectrum
    Phi_xy=atan2d(imag(Gxy),real(Gxy));
    
    %%%%MC simulation
    if nsim==0
        ciGXX=0; ciGYY=0; ciCXY=0; ciPhi=0;
    end;
    if nsim>0
        [tmp,taux,tmp1]=tauest(tx,x,0);
        [tmp,tauy,tmp1]=tauest(ty,y,0);
        gxx=zeros(nfx,nsim);
        gyy=zeros(nfy,nsim);
        cxy=zeros(nfxy,nsim);
        for j=2:length(tx)-1;
            pd=makedist('Normal','mu',0,'sigma',sqrt(1-exp(-2*(tx(j)-tx(j-1))/taux)));
            Eps_x(:,j)=random(pd,nsim,1);
        end;
        for j=2:length(ty)-1;
            pd=makedist('Normal','mu',0,'sigma',sqrt(1-exp(-2*(ty(j)-ty(j-1))/tauy)));
            Eps_y(:,j)=random(pd,nsim,1);
        end;
        for i=1:nsim
            if mod(i,100)==0;
                disp(i);
            end
            xx=TAUrednoise(tx,taux,Eps_x(i,:));
            yy=TAUrednoise(ty,tauy,Eps_y(i,:));
            [gxx(:,i),gyy(:,i),~,cxy(:,i),~,~,~,~,~]=main(tx,ty,fx,fy,fxy,indices_x,indices_y,xx,yy,nseg,ratio,window,beta,0);
        end;
        
        gxx=sort(gxx,2,'ascend');
        gyy=sort(gyy,2,'ascend');
        cxy=sort(cxy,2,'ascend');
        ind1=[fix(nsim*0.90),fix(nsim*0.95),fix(nsim*0.99)];
        ciGXX=gxx(:,ind1);
        ciGYY=gyy(:,ind1);
        ciCXY=cxy(:,ind1);
            
        %%%%%%%%%%%%%%%%%%%%%%%%%
        ciPhi=0;
        ind_sigcxy=find(Cxy>=ciCXY(:,3));
        G=sqrt(-2*(sqrt(1-Cxy(ind_sigcxy))./Cxy(ind_sigcxy))+2./Cxy(ind_sigcxy)-1);
        fxy_sig=fxy(ind_sigcxy);
        
        tmp11=randn(fix(nsim/2),length(tx));
        rng('shuffle');
        tmp12=randn(fix(nsim/2),length(tx));
 
        tmp21=randn(fix(nsim/2),length(ty));
        rng('shuffle');
        tmp22=randn(fix(nsim/2),length(ty));
        
        phi_xy1=zeros(nfxy,fix(nsim/2));
        phi_xy2=zeros(nfxy,fix(nsim/2));
        ciPhi=zeros(2,length(fxy));    
        
        ind2=[fix(fix(nsim/2)*0.025),fix(fix(nsim/2)*0.975)];
        for k=1:length(fxy_sig);
           n=ind_sigcxy(k);
           sigma=sqrt(1-Cxy(n))/sqrt(Cxy(n))/sqrt(2*mean(neff));
           ciPhi(1,n)=-2.58*sigma;
           ciPhi(2,n)=2.58*sigma;
%            f1=fx(n);  f2=fy(n); f3=fxy(n);
%            Eps_x=tmp11+G(k)*tmp12;
%            Eps_y=tmp12+G(k)*tmp11;
%            for i=1:fix(nsim/2);   %%% Tn set to tx;
%                if mod(i,100)==0;
%                   disp(i);
%                end;
%                xx=TAUrednoise(tx,taux,Eps_x(i,:));
%                yy=TAUrednoise(tx,tauy,Eps_y(i,:));    
%                [~,~,~,phi_xy1(:,i),~,~,~,~]=main(tx,tx,f1,f2,f3,indices_x,indices_x,xx,yy,nseg,ratio,window,beta,0);
%            end;
%            phi_xy1=phi_xy1-repmat(mean(phi_xy1,2),1,fix(nsim/2));
%            phi_xy1=sort(phi_xy1,2,'ascend');
%            ciPhi1=phi_xy1(ind_sigcxy(k),ind2);
%            
%            Eps_x=tmp21+G(k)*tmp22;
%            Eps_y=tmp22+G(k)*tmp21;
%            for i=1:fix(nsim/2);   %%% Tn set to ty;
%                if mod(i,100)==0;
%                   disp(i);
%                end
%                xx=TAUrednoise(ty,taux,Eps_x(i,:));
%                yy=TAUrednoise(ty,tauy,Eps_y(i,:));
%                [~,~,~,phi_xy2(:,i),~,~,~,~]=main(ty,ty,f1,f2,f3,indices_y,indices_y,xx,yy,nseg,ratio,window,beta,0);
%            end;
%            phi_xy2=phi_xy2-repmat(mean(phi_xy2,2),1,fix(nsim/2));
%            phi_xy2=sort(phi_xy2,2,'ascend');
%            ciPhi2=phi_xy2(ind_sigcxy(k),ind2);
%            
%            ciPhi(:,k)=(ciPhi1'+ciPhi2')/2;
        end;
    end;
end

function W=unevenWindow(p,flag,beta)
    if flag==0
        W=0.5*(1-cos(2*pi*p));
    elseif flag==1
        W=besseli(0,beta*sqrt(1-(2*p-1).^2))/besseli(0,beta);    
    end;
end

function [cx,cy]=EffNseg(tx,ty,Wx,Wy,window,ratio,beta)
    px=(tx-tx(1))/(tx(end)-tx(1));
    py=(ty-ty(1))/(ty(end)-ty(1));
    px=px(px<=ratio);
    py=py(py<=ratio);
    pxx=px+(1-ratio);
    pyy=py+(1-ratio);
    
    Wxx=unevenWindow(pxx,window,beta);
    Wyy=unevenWindow(pyy,window,beta);
    
    cx=sum(Wx(1:length(px)).*Wxx)./sum(Wx.^2);
    cy=sum(Wy(1:length(py)).*Wyy)./sum(Wy.^2);
end