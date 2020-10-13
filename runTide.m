%{
 Simple Tide and sea level variation analysis code by C.Ozsoy (2019),
 please fill inputs arguments!

 v01
    - if else condition at start is removed
    - peak conditions are altered due to previous error.
    - functions are placed under fnc folder

%}
clc
clear
close all
tic

%% Inputs

fileName='input.txt';                                                       % Enter filename
dt=60*60;                                                                   % Enter recording time increment (seconds)
FontName='Calibri';                                                         % Enter fontname for outputs
FontSize=25;                                                                % Enter fontsize for outputs
tend=48*60*60;                                                              % Enter plot upper limit (seconds)
lFilt=10*60*60;                                                             % Enter low pass filter for tide analysis (seconds)
hFilt=48*60*60;                                                             % Enter high pass filter for tide analysis (seconds)
lsFilt=8760/12*3600;                                                        % Enter seasonal low pass filter (seconds)
threshold=21;                                                               % Minimum data number in oneday (if a day has less data than the value that day will not be considered!)
switchLongTerm=1;                                                           % Long term sea level rise calculation switch 1:open 2:close

%% General Calculations

addpath('fnc');
tideData=fun_tideReader(fileName);
[row,~]=size(tideData);
NumtideData = zeros(row, 2);

for i=1:row
    temp=tideData(i,1)+" "+tideData(i,2);
    NumtideData(i,1)=datenum(temp,'dd.mm.yyyy HH:MM:SS');
    NumtideData(i,2)=str2double(tideData(i,3));
end

initial=ceil(NumtideData(1,1));
index=find(NumtideData(:,1)<initial);
NumtideData(index,:)=[];
[row,col]=size(NumtideData);

%% Spectral Analysis

fData=fft(NumtideData(:,2))/row;
a=real(fData);
b=imag(fData);
c=2.*sqrt(a.^2+b.^2);
f=(0:1:row-1)/(row-1)/dt;
T=1./f;
TFlip=fliplr(T);
Inend=knnsearch(TFlip',tend);

figure('units','normalized','outerposition',[0,0,1,1]);
plot(TFlip(2:Inend),c(2:Inend));
xtix=[TFlip(1),12*60*60,24*60*60,36*60*60,tend];
set(gca,'XTick',xtix,'XTickLabel',xtix./60./60);
set(gca,'FontName',FontName,'FontSize',FontSize);
grid on
xlabel('Period (Hours)');
ylabel('Constituent Amplitude (m)');

IndM2=knnsearch(f',1/12.42/60/60);
M2=c(IndM2);
IndS2=knnsearch(f',1/12/60/60);
S2=c(IndS2);
IndK1=knnsearch(f',1/23.93/60/60);
K1=c(IndK1);
IndO1=knnsearch(f',1/25.82/60/60);
O1=c(IndO1);
type=(K1+O1)/(M2+S2);

fprintf('##### Spectral Analysis #####\nPlease Check Constituents Manually!\n');

if type>=0 && type<0.25
    fprintf('Tide is Semiduirnal\n');
    option=1;
elseif type>=0.25 && type<1.5
    fprintf('Tide is mixed, mainly semi-duirnal\n');
    if type<0.5
        option=1;
    else
        option=2;
    end
elseif type>=1.5 && type<3.0
    fprintf('Tide is mixed, mainly diurnal\n');
    option=2;
elseif type>=3
    fprintf('Tide is diurnal\n');
    option=2;
else
    error('Tide type cannot be determined!');
end

if option==1
    spec_MLWS=-(M2+S2);
    spec_MLWN=-abs(M2-S2);
    spec_MHWN=abs(M2-S2);
    spec_MHWS=(M2+S2);
    fprintf('MLWS = %.2f\nMLWN = %.2f\nMHWN = %.2f\nMHWS = %.2f\n',spec_MLWS,...
        spec_MLWN,spec_MHWN,spec_MHWS);
elseif option==2
    spec_MLLW=-(M2+K1+O1);
    spec_MHLW=-abs(M2-(K1+O1));
    spec_MLHW=abs(M2-(K1+O1));
    spec_MHHW=(M2+K1+O1);
    fprintf('MLLW = %.2f\nMHLW = %.2f\nMLHW = %.2f\nMHHW = %.2f\n',spec_MLLW,...
        spec_MHLW,spec_MLHW,spec_MHHW);
end

fprintf('#############################\n');

%% Time-Series Analysis

[tsFiltered]=fun_lPass(NumtideData(:,2),1/dt,1/lFilt);
[tsFiltered]=fun_hPass(tsFiltered,1/dt,1/hFilt);
[tsFilteredS]=fun_lPass(NumtideData(:,2),1/dt,1/lsFilt);
tsFilteredMslCor=tsFiltered+mean(NumtideData(:,2));
st_date=NumtideData(1,1);
end_date=NumtideData(row,1);

k=1;
for j=st_date:(end_date-1)
    m=1;
    for i=1:row
        if j<=NumtideData(i,1) && NumtideData(i,1)<j+1
            sep_tide{k,1}(m,:)=NumtideData(i,:);
            m=m+1;
        end
    end
    k=k+1;
end

[row2,col2]=size(sep_tide);
k=1;

for i=1:row2
    if isempty(sep_tide{i,1}) || numel(sep_tide{i,1})<(2*threshold)
        continue
    else
        [pks,lcs]=findpeaks(sep_tide{i,1}(:,2));
        [pks2,cls2]=findpeaks(-sep_tide{i,1}(:,2));
        pks2=-pks2;
        if isempty(pks)
            out=sort(sep_tide{i,1}(:,2),'descend');
            pks(1,1)=out(1);
            pks(1,2)=out(2);
        end
        if isempty(pks2)
            out2=sort(sep_tide{i,1}(:,2));
            pks2(1,1)=out2(1);
            pks2(1,2)=out2(2);
        end
        if size(pks)==1
            out=sort(sep_tide{i,1}(:,2),'descend');
            if pks>=max(out)
                pks(2)=out(2);
            else
                pks(2)=max(out);
            end
        end
        if size(pks2)==1
            out2=sort(sep_tide{i,1}(:,2));
            if pks2<min(out2)
                pks2(2)=out2(2);
            else
                pks2(2)=min(out2);
            end
        end
        sortedValue = sort(pks);
        sortedValue2 = sort(pks);
        hhw(k)=sortedValue(end);
        lhw(k)=sortedValue(end - 1);
        llw(k)=sortedValue2(1);
        hlw(k)=sortedValue2(2);
        k=k+1;
    end
end

mhhw=mean(hhw);
mlhw=mean(lhw);
mllw=mean(llw);
mhlw=mean(hlw);
msl=mean(NumtideData(:,2));
txt=["MHHW","MSL","MLLW"];
[~,chnwaves,Hm,Tm,Hs,Ts,~,~,~,~,~,~,~]=fun_tsa(tsFiltered,1/dt);
st_date_num=str2double(datestr(st_date,'yyyy'));
end_date_num=str2double(datestr(end_date,'yyyy'));
years=st_date_num:end_date_num;
nYears=numel(years);
years(nYears+1)=end_date_num+1;

for i=1:nYears
    temp1=datenum(num2str(years(i)),'yyyy');
    temp2=datenum(num2str(years(i+1)),'yyyy');
    index=find(NumtideData(:,1)>=temp1 & NumtideData(:,1)<temp2);
    MinMax(i,1)=min(tsFilteredS(index));
    MinMax(i,2)=max(tsFilteredS(index));
end

MinMax(:,3)=MinMax(:,2)-MinMax(:,1);
SeasonalChange=max(MinMax(:,3));
if switchLongTerm==1
    XPlot=linspace(years(1),years(end),numel(NumtideData(:,1)));
    [As,Bs]=postreg(tsFilteredS,NumtideData(:,1),'hide');
    [AsPlot,BsPlot]=postreg(tsFilteredS,XPlot','hide');
end


%% Plots

figure('units','normalized','outerposition',[0,0.3,1,0.5]);
hold on
pl=plot(NumtideData(:,1),NumtideData(:,2));
pl.Color(4)=0.8;
line([NumtideData(1,1),NumtideData(row,1)],[mhhw mhhw],'LineWidth',2,'Color','r','LineStyle','--');
line([NumtideData(1,1),NumtideData(row,1)],[msl msl],'LineWidth',2,'Color',[0.85,0.325,0.098],'LineStyle','--');
line([NumtideData(1,1),NumtideData(row,1)],[mllw mllw],'LineWidth',2,'Color','k','LineStyle','--');
legend('Data',txt(1),txt(2),txt(3),'Position',[0.95  .61 0 0])
% text(NumtideData(1,1),mhhw+0.05,txt(1),'FontName',FontName,'FontSize',FontSize);
% text(NumtideData(1,1),mllw+0.05,txt(3),'FontName',FontName,'FontSize',FontSize);
% text(NumtideData(1,1),msl+0.05,txt(2),'FontName',FontName,'FontSize',FontSize);
datetick('x','yyyy');
xlabel('Record Date');
ylabel({'Sea Surface'; 'Elevation (m)'});
set(gca,'FontName',FontName,'FontSize',FontSize);

figure('units','normalized','outerposition',[0,0.3,1,0.5]);
hold on
plot(NumtideData(:,1),tsFiltered);
datetick('x','yyyy');
xlabel('Record Date');
ylabel({'Sea Surface'; 'Elevation (m)'});
set(gca,'FontName',FontName,'FontSize',FontSize);

figure('units','normalized','outerposition',[0,0.3,1,0.5]);
hold on
plot(NumtideData(:,1),tsFilteredS);
if switchLongTerm==1
    plot(NumtideData(:,1),As*NumtideData(:,1)+Bs,'Color','k');
    tp1=num2str(AsPlot,'y=%.6f');
    tp2=num2str(BsPlot,'%.6f');
    tsPlot=strcat(tp1,tp2);
    text(NumtideData(round(row*0.8),1),As*NumtideData(row,1)+Bs+0.02,tsPlot,'FontName',FontName,'FontSize',FontSize);
end
datetick('x','yyyy');
grid on
xlabel('Record Date');
ylabel({'Sea Surface'; 'Elevation (m)'});
set(gca,'FontName',FontName,'FontSize',FontSize);


fprintf('\n##### Time-Series Analysis #####\n');
fprintf('MLLW/MLWS = %.2f\n',mllw);
fprintf('MHLW/MLWN = %.2f\n',mhlw);
fprintf('MSL = %.2f\n',msl);
fprintf('MLHW/MHWN = %.2f\n',mlhw);
fprintf('MHHW/MHWS = %.2f\n',mhhw);
fprintf('Hm = %.2f\n',Hm);
fprintf('Tm = %.2f hour\n',Tm/60/60);
fprintf('Seasonal Variation = %.2f \n',SeasonalChange);

rmpath('fnc');
toc


