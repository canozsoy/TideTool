
function [tsFiltered]=fun_lPass(record,hertz,cutOff)

%c
%hertz=10;
dt=1/hertz;
%record=xlsread('wave.xls');
%record=dlmread('irrTimeSeries_4Hz_withSurfBeat_mean0.dat');
[N,col]=size(record);

df=(1/N/dt);
X=fft(record);
ak=real(X);
bk=imag(X);
phs=atan2(bk,ak);
Xc=conj(X);
S=X.*Xc;
fk=0;
Sfk=0;
for i=2:N/2+1
fk(i)=(i-1)*df;
Sfk(i)=S(i)*(2/N^2)/(df);
% Sfk(i)=S(i)/(df);
end

Sfk(fk > cutOff)=0;

for i=1:N/2+1
cphs=complex(0,phs);
hnk=complex(1,0);
Xk(i,1)=hnk*sqrt(.5*Sfk(i)*df*N^2)*exp(cphs(i));
% Xk(i,1)=hnk*sqrt(Sfk(i)*df)*exp(cphs(i));
end
for i=2:N/2+1
Xk(N-i+2,1)=conj(Xk(i,1));
end
tsFiltered=ifft(Xk);