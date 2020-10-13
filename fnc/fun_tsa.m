% WATER SURFACE PROFILE TIME SERIES ANALYSIS

function [nwaves,chnwaves,Hm,Tm,Hs,Ts,H10,T10,H2perc,Hrms,amean,skewness,kurtosis]=fun_tsa(record,hertz)
% c
% data=dlmread('sil.txt');
% hertz=4;
% record=data(:,2);%xlsread('wave.xls');

tinc=1/hertz;
[ntimestep,nchannel]=size(record);
amean=mean(record);
% ETA RMS
sumsqr=0;
for i=1:ntimestep
    sumsqr=sumsqr+record(i,1)*record(i,1);
end
nrms=sqrt(sumsqr/ntimestep);

% SKEWNESS AND KURTOSIS
sumskew=0;
sumkurt=0;
for i=1:ntimestep
    sumskew=sumskew+(record(i,1)-amean)^3;
    sumkurt=sumkurt+(record(i,1)-amean)^4;
end
skewness=(1/nrms^3)*(1/ntimestep)*sumskew;
kurtosis=(1/nrms^4)*(1/ntimestep)*sumkurt;


% ARITHMETIC CORRECTION
n(:,1)=record(:,1)-amean;

% % LINEAR CORRECTION
% N0=ntimestep;
% N1=0;
% N2=0;
% Y0=0;
% Y1=0;
% for i=1:ntimestep
%     N1=N1+i;
%     N2=N2+i*i;
%     Y0=Y0+record(i,1);
%     Y1=Y1+i*record(i,1);
% end
% A0=(N2*Y0-N1*Y1)/(N0*N2-N1*N1);
% A1=(N0*Y1-N1*Y0)/(N0*N2-N1*N1);
% lmean(:,1)=A0+A1*record(:,1);
% n(:,1)=record(:,1)-lmean(i,1);

% NUMBER OF WAVES
nwaves=1;
for i=1:1:ntimestep-2
    if n(i,1)<0 && n(i+1,1)>=0
        nstart(nwaves,1)=i+1;
        tstart(nwaves,1)=(i-1)*tinc + (tinc)*abs(n(i,1))/(abs(n(i,1))+abs(n(i+1,1)));
        nwaves=nwaves+1;
    end
end

% PERIODS
for k=1:nwaves-2
    period(k,1)=tstart(k+1,1)-tstart(k,1);
end

% WAVE HEIGHTS
for k=1:nwaves-2
    nmaxo=0;
    nmino=0;
    for i=nstart(k,1):1:nstart(k+1,1)-1
        if n(i-1,1)<=n(i,1) && n(i,1)>n(i+1,1)
            A=0.5*(n(i-1,1)-2*n(i,1)+n(i+1,1));
            B=0.5*(n(i+1,1)-n(i-1,1));
            C=n(i,1);
            nmaxn=C-B*B/A;
            if nmaxn > nmaxo
                nmaxo=nmaxn;
            end
            nmax(k,1)=nmaxo;
        end

        if n(i-1,1)>=n(i,1) && n(i,1)<n(i+1,1)
            A=0.5*(n(i-1,1)-2*n(i,1)+n(i+1,1));
            B=0.5*(n(i+1,1)-n(i-1,1));
            C=n(i,1);
            nmin(k,1)=C-B*B/A;
            nminn=C-B*B/A;
            if nminn < nmino
                nmino=nminn;
            end
            nmin(k,1)=nmino;
        end
    end
    height(k,1)=abs(nmax(k,1))+abs(nmin(k,1));
end

% SORTING AND REST OF IT
chnwaves(:,1)=height(:,1);
chnwaves(:,2)=period(:,1);
chnwaves_unsorted=chnwaves;
chnwaves = sortrows(chnwaves, 1);
[row,col]=size(chnwaves);
chnwaves(1:row-nwaves(1)+2,:)=[];
[row,col]=size(chnwaves);
hm_param=mean(chnwaves);
Hm=hm_param(1);
Tm=hm_param(2);
Hmax=chnwaves(end,1);
Tmax=chnwaves(end,2);
chnsigwaves=chnwaves(round(row*2/3)+1:row,:);
hs_param=mean(chnsigwaves);
Hs=hs_param(1);
Ts=hs_param(2);
H10=0;
T10=0;
H100=0;
T100=0;
sumrms=0;
for k=1:row
sumrms=sumrms+chnwaves(k,1)^2;
end
Hrms=sqrt(sumrms/row);
nwaves=nwaves-2;
if nwaves>=10;
    chn10waves=chnwaves(round(row*9/10)+1:row,:);
    [r10,c10]=size(chn10waves);
    if r10==1
        h10_param=chn10waves;
    else
    h10_param=mean(chn10waves);
    end
    H10=h10_param(1);
    T10=h10_param(2);
end
if nwaves>=100;
    chn100waves=chnwaves(round(row*99/100)+1:row,:);
    [r100,c100]=size(chn100waves);
    if r100==1
        h100_param=chn100waves;
    else
    h100_param=mean(chn100waves);
    end
    H100=h100_param(1);
    T100=h100_param(2);
end
H2perc=0;
T200=0;
if nwaves>=200;
    chn200waves=chnwaves(round(row*49/50)+1:row,:);
    [r200,c200]=size(chn200waves);
    if r200==1
        H2perc_param=chn200waves;
    else
    H2perc_param=mean(chn200waves);
    end
    H2perc=chnwaves(round(row*49/50)+1,1);%H2perc_param(1);
    T200=chnwaves(round(row*49/50)+1,2);%h200_param(1);
end