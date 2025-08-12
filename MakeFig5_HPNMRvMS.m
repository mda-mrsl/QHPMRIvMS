% Get kpl values and create figure to show relationship between
% values derived from ICMS and from NMR
clear
resave=1;

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
    [tmpmskpl,~]=estparms(residvec,kplvec,fitsvec,fdv,0);

    %Get kpl values for multiple NMR runs:
    load(hpfns{ii},'Group','kplvec');
    Names{ii}=Group(1).Name;
    NRuns(ii)=length(Group(1).Data);
    tmphpkpl=zeros(1,NRuns(ii));
    for jj=1:NRuns(ii)
        [tmp,~]=estparms(Group(1).Data(jj).residvec,kplvec,Group(1).Data(jj).fitsvec,Group(1).Data(jj).fdv,0);
        tmphpkpl(jj)=tmp;
    end
    hpkpls=[hpkpls tmphpkpl];
    mskpls=[mskpls tmpmskpl*ones(1,NRuns(ii))];
    fprintf('%s: MS Kpl = %5.3f; HP Kpl = %5.3f \x00B1 %5.3f\n',...
        Group(1).Name,tmpmskpl, mean(tmphpkpl), std(tmphpkpl));
end

%What happens if we exclude LnCap run 1, kpl ~0.934 as outlier?
%mskpls=mskpls([1:8 10:12]);
%hpkpls=hpkpls([1:8 10:12]);
%NRuns(3)=NRuns(3)-1;
%[rho,p]=corrcoef(mskpls([1:8 10:12]),hpkpls([1:8 10:12]))
% % Result: rho=0.797, p=0.0033; 
% % Regression: y = 0.618x - 0.006
% % Regression slope: 0.618 +/- .0156
% % RSquared = 0.635
% Answer: LnCap run 1 does not significantly affect outcome either way.

figure(5)
plot(mskpls(1:NRuns(1)),hpkpls(1:NRuns(1)),'m+',...
    mskpls(NRuns(1)+[1:NRuns(2)]),hpkpls(NRuns(1)+[1:NRuns(2)]),'bx',...
    mskpls(NRuns(1)+NRuns(2)+[1:NRuns(3)]),hpkpls(NRuns(1)+NRuns(2)+[1:NRuns(3)]),'r*',...
    'LineWidth',1,'MarkerSize',8);
xlabel('Mass Spectrometry k_{PL} (s^{-1})','FontSize',18,'FontWeight','bold')
ylabel('HP NMR k_{PL} (s^{-1})','FontSize',18,'FontWeight','bold')
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
text(0.2,0.3,sprintf('y = %5.3fx %s %5.3f',lrce(2),tmps,abs(lrce(1))));
%text(0.25,0.1,sprintf('Pearson \\rho = %5.3f',rho(1,2)));
text(0.35,0.5,sprintf('R^2 = %5.3f',R2));
%fprintf('Pearson rho = %5.3f (P=%5.3f)\n',rho(1,2),pcor(1,2));

if resave||(~exist('Fig5-Corr.jpg','file'))
    print('Fig5-Corr','-depsc');
    exportgraphics(gcf,'Fig5-Corr.png');
end

%Copy figure to clipboard
copygraphics(gcf)


