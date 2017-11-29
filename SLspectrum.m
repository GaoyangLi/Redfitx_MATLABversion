function [f, P, Preal, Pimag, prob]=SLspectrum(t,h,F)
% Written by Gaoyang Li, 23/Nov./2017.

% Estimation of the Scargle-Lomb power spectrum for unevenly-sampled data
% INPUTS:
%      Vector: time t at which the observation is obtained;
%      Vector: observation h; must be of same length as t;
%      Vector: F, a vetor of frequencies at which P is evaluated;
% OUTPUTS:
%      Vector: f, frequencies at which power is evaluated;
%      Vector: P, power values;
%      Vector: Preal, Pimag: the real and imaginary parts of the SL Fourier
%      Transform;
%      prob: False Alarm Probability.

% Theoretical background: Scargle, 1989; The Astrophysical Journal 343, 874-887;

N=length(h);
T=max(t)-min(t);

L=length(t);
sumx=sum(h);
mu=sumx/L;
s2=var(h);

f=F;
if length(f)>1
    Ofac = f(1)/(f(2)-f(1));   %%% oversampling factor
end;

w = 2*pi*f;
nf=length(f);
tol1=1e-4;
tol2=1e-8;

csum=sum(cos(2*w*t.'),2);
ssum=sum(sin(2*w*t.'),2);
sumtc=sum(repmat(t',nf,1).*cos(2*w*t.'),2);
sumts=sum(repmat(t',nf,1).*sin(2*w*t.'),2);
wtau=zeros(nf,1);
for i=1:nf
    if sumtc(i)>tol1 || sumts(i)>tol1
        wtau(i)=atan2(ssum(i),csum(i))/2;
    else
        wtau(i)=atan2(-sumtc(i),sumts(i))/2;
    end;
end;

const=1/sqrt(2*s2);
cterm = cos(w*t.' - repmat(wtau,1,L));
sterm = sin(w*t.' - repmat(wtau,1,L));
sumr = sum(cterm*diag(h-mu),2);
sumi = sum(sterm*diag(h-mu),2);

scos2 = sum(cterm.^2,2);
ssin2 = sum(sterm.^2,2);

ftid=zeros(nf,1);
ftrd=zeros(nf,1);

ftrd=const*sumr./sqrt(scos2);
%ftid=const*sumi./sqrt(ssin2);
for i=1:nf                        %%%%%%% Taken from Scargle, 1989. Here
                                  %%%%%%% the algorithm will not be
                                  %%%%%%% implemented for omega=0.
    cross=sum(t'.*cterm(i,:).*sterm(i,:));
    if ssin2(i)<=tol1
        ftid(i)=const*sumx/sqrt(L);
        if (abs(cross)>tol2)
            ftid(i)=0.0;
        end;
    else
        ftid(i)=const*sumi(i)/sqrt(ssin2(i));
    end;
end;
work=complex(ftrd,ftid).*exp(complex(zeros(nf,1),wtau));
P=abs(work).^2;
Preal=real(work);
Pimag=imag(work);

if length(f)>1
    M=2*length(f)/Ofac;
    prob = M*exp(-P);
    inds = prob > 0.01;
    prob(inds) = 1-(1-exp(-P(inds))).^M;
else
    prob=0;
end;
end