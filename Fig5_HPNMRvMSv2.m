% Get kpl values and create figure to show relationship between
% values derived from ICMS and from NMR
clear
savefig=0;

msfns={fullfile('data','MSPKSweepHth83.mat'),...
    fullfile('data','MSPKSweepPC3.mat'),...
    fullfile('data','MSPKSweepLnCap.mat')};
hpfns={fullfile('data','HPNMRFits_Hth83.mat'),...
    fullfile('data','HPNMRFits_PC3.mat'),...
    fullfile('data','HPNMRFits_LnCap.mat')};

hpkpls=[];
mskpls=[];
NRuns=zeros(1,3);
Names={};
for ii=1:3
    %Get kpl value for mass spec:
    load(msfns{ii},'residvec','kplvec','fitsvec','fdv');
    kpli=find(strcmp('kpl',fdv.knowns));
    kecpi=find(strcmp('kecp',fdv.fitvars));
    kecli=find(strcmp('kecl',fdv.fitvars));
    klpi=find(strcmp('klp',fdv.fitvars));
    rv=residvec;
    kecps=fitsvec(:,kecpi)';
    kecls=fitsvec(:,kecli)';
    klps=fitsvec(:,klpi)';
%    sos=sqrt((kecps.^2)+(kecls.^2)+(klps.^2));
%    sos=sqrt((kecps.^2)+(kecls.^2)+(klps.^2)+(kplvec.^2));
    sos=sqrt(((kecps/0.2).^2)+((kecls/0.2).^2)+((klps).^2)+((kplvec/1.5).^2));
    [kplll,tmparms]=estparms(residvec,kplvec,fitsvec,fdv,0);
    sosll=sqrt((tmparms(kecpi)/0.2)^2+(tmparms(kecli)/0.2)^2+(tmparms(klpi))^2+(kplll/1.5)^2);
    tmpmskpl=kplll;
    ll=find(kplvec>kplll,1); % look above lower limit...
    [minsos,minsosi]=min(sos(ll:end));
    minsosi=minsosi+ll-1;
    if minsos<sosll
        tmpmskpl=kplvec(minsosi);
    end

    figure(528);
    tiledlayout(5,1,'TileSpacing','Compact','Padding','Compact');
    nexttile
    plot(kplvec,rv,'-',kplll,min(rv)+0.01*(max(rv)-min(rv)),'r>');
    ylabel('Resid');
    title('MS Fits')    
    nexttile
    plot(kplvec,kecps);
    ylabel('kecp')
    nexttile
    plot(kplvec,kecls);
    ylabel('kecl')
    nexttile
    plot(kplvec,klps);
    ylabel('klp')
    nexttile
    plot(kplvec,sos,'-',kplvec(minsosi),sos(minsosi),'ro');
    ylabel('SOS')
    xlabel('kpl')


    %Get kpl values for multiple NMR runs:
    load(hpfns{ii},'Group','kplvec');
    Names{ii}=Group(1).Name;
    NRuns(ii)=length(Group(1).Data);
    tmphpkpl=zeros(1,NRuns(ii));
    for jj=1:NRuns(ii)
        kpli=find(strcmp('kpl',Group(1).Data(jj).fdv.knowns));
        kecpi=find(strcmp('kecp',Group(1).Data(jj).fdv.fitvars));
        kecli=find(strcmp('kecl',Group(1).Data(jj).fdv.fitvars));
        klpi=find(strcmp('klp',Group(1).Data(jj).fdv.fitvars));
        rv=Group(1).Data(jj).residvec;
        kecps=Group(1).Data(jj).fitsvec(:,kecpi)';
        kecls=Group(1).Data(jj).fitsvec(:,kecli)';
        klps=Group(1).Data(jj).fitsvec(:,klpi)';
%        sos=sqrt((kecps.^2)+(kecls.^2)+(klps.^2));
%        sos=sqrt((kecps.^2)+(kecls.^2)+(klps.^2)+(kplvec.^2));
        sos=sqrt(((kecps/0.2).^2)+((kecls/0.2).^2)+((klps).^2)+((kplvec/1.5).^2));
        [kplll,tmparms]=estparms(Group(1).Data(jj).residvec,kplvec,Group(1).Data(jj).fitsvec,Group(1).Data(jj).fdv,0);
        sosll=sqrt((tmparms(kecpi)/0.2)^2+(tmparms(kecli)/0.2)^2+(tmparms(klpi))^2+(kplll/1.5)^2);
        tmphpkpl(jj)=kplll; % lowest possible kpl value
        ll=find(kplvec>kplll,1); % look above lower limit...
        [minsos,minsosi]=min(sos(ll:end));
        minsosi=minsosi+ll-1;
        if minsos<sosll
            tmphpkpl(jj)=kplvec(minsosi);
        end

        figure(529);
        tiledlayout(5,1,'TileSpacing','Compact','Padding','Compact');
        nexttile
        plot(kplvec,rv,'-',kplll,min(rv)+0.01*(max(rv)-min(rv)),'r>');
        ylabel('Resid');
        title('NMR Fits')    
        nexttile
        plot(kplvec,kecps);
        ylabel('kecp')
        nexttile
        plot(kplvec,kecls);
        ylabel('kecl')
        nexttile
        plot(kplvec,klps);
        ylabel('klp')
        nexttile
        plot(kplvec,sos,'-',kplvec(minsosi),sos(minsosi),'ro');
        ylabel('SOS')
        xlabel('kpl')
        %pause

    end
    hpkpls=[hpkpls tmphpkpl];
    mskpls=[mskpls tmpmskpl*ones(1,NRuns(ii))];
    fprintf('%s: MS Kpl = %5.3f; HP Kpl = %5.3f \x00B1 %5.3f\n',...
        Group(1).Name,tmpmskpl, mean(tmphpkpl), std(tmphpkpl));
end

%What happens if we exclude LnCap run 1 as outlier, given ?
%mskpls=mskpls([1:8 10:12]);
%hpkpls=hpkpls([1:8 10:12]);
%NRuns(3)=NRuns(3)-1;

% Answer: LnCap run 1 does not significantly affect outcome either way.

f5h=figure(5);
plot(mskpls(1:NRuns(1)),hpkpls(1:NRuns(1)),'m+',...
    mskpls(NRuns(1)+[1:NRuns(2)]),hpkpls(NRuns(1)+[1:NRuns(2)]),'bx',...
    mskpls(NRuns(1)+NRuns(2)+[1:NRuns(3)]),hpkpls(NRuns(1)+NRuns(2)+[1:NRuns(3)]),'r*',...
    'LineWidth',1,'MarkerSize',8);
xlabel('Mass Spectrometry k_{PL} (s^{-1})','FontSize',14,'FontWeight','bold')
ylabel('HP NMR k_{PL} (s^{-1})','FontSize',14,'FontWeight','bold')
%axis tight
grid on

%Test for correlation:
[rho,pcor]=corrcoef(mskpls,hpkpls);
%Result: rho=0.667, p=0.0179

%Linear Regression
lr=fitlm(mskpls,hpkpls);
lrce=lr.Coefficients.Estimate;
lrmci=lr.coefCI;
%Result: y = 0.89765x - 0.07337
fprintf('Regression Slope = %5.3f \x00B1 %5.3f\n',lrce(2),lr.Coefficients.SE(2));
fprintf('Regression Slope Confidence Interval = [%5.3f, %5.3f]\n',lrmci(2,1),lrmci(2,2));

%slope = 0.89765 +/- 0.31735
R2=lr.Rsquared.Ordinary;
% R-Squared = 0.444; Adjusted R-Squared 0.389
p=lr.Coefficients.pValue(2);
fprintf('p = %5.3f\n',p)

hold on
plot(mskpls,lrce(1)+mskpls*lrce(2),'k-','linewidth',1);
hold off
Names{4}='Linear Regression';
f5lh=legend(Names,'Location','NorthWest','FontSize',10);
if lrce(1)>0
    tmps="+";
else
    tmps="-";
end
text(0.25,0.3,sprintf('y = %5.3fx %s %5.3f',lrce(2),tmps,abs(lrce(1))));
%text(0.25,0.1,sprintf('Pearson \\rho = %5.3f',rho(1,2)));
text(0.3,0.8,sprintf('R^2 = %5.3f',R2));
%fprintf('Pearson rho = %5.3f (P=%5.3f)\n',rho(1,2),pcor(1,2));

%ofn='TestNorms_SOS3.png';
tmp=datetime('today');
ofn=sprintf('%s.%d%02d%02d',mfilename,tmp.Year,tmp.Month,tmp.Day);
if savefig
    print(ofn,'-depsc');
    exportgraphics(gcf,[ofn '.png']);
end

%Copy figure to clipboard
copygraphics(gcf)

%% Results:
% 
% using kpl @ residual reduced to within 1% of min (from max):
% and min sos=sqrt(((kecps/0.2).^2)+((kecls/0.2).^2)+((klps).^2)+((kplvec/1.5).^2));
%
% Hth83: MS Kpl = 0.559; HP Kpl = 0.603 ± 0.140
% PC3: MS Kpl = 0.158; HP Kpl = 0.252 ± 0.024
% LnCap: MS Kpl = 0.890; HP Kpl = 0.934 ± 0.432
% Regression Slope = 0.932 ± 0.432
% Regression Slope Confidence Interval = [0.343, 1.521], delta=1.178
% p = 0.005
% R^2 = 0.554 
