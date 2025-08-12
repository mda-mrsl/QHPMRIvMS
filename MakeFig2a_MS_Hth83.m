clear

%Choose cancer cell line: Hth83 ATC cells, PC3 prostate cancer cells, LnCap
%prostate cancer cells 
tumor={'Hth83','PC3','LnCap'};
tflag=1;

%Load analysis file
tmpfn=fullfile('data',sprintf('MSPKSweep%s.mat',tumor{tflag}));
load(tmpfn);

%Find the best fit @ minimum residual:
[minresid minresi]=min(residvec);
%use kpl values [0,1]:
minresi=min([45,minresi]);
minresid=residvec(minresi);
%Estimate the first kpl value at which residual decreases 
%by less than one percent of full range
tmpi=1;
thresh=(max(residvec)-min(residvec))/100 + min(residvec);
while residvec(tmpi)>thresh
    tmpi=tmpi+1;
end
%assume residual is linear over this interval:
slope=(residvec(tmpi)-residvec(tmpi-1))/(kplvec(tmpi)-kplvec(tmpi-1));
kplest=kplvec(tmpi-1) + (thresh-residvec(tmpi-1))/slope;
%Same with rest of the fit parameters:
parmsest=zeros(1,size(fitsvec,2));
for qq=1:size(fitsvec,2)
    slope=(residvec(tmpi)-residvec(tmpi-1))/(fitsvec(tmpi,qq)-fitsvec(tmpi-1,qq));
    parmsest(qq)=fitsvec(tmpi-1,qq) + (thresh-residvec(tmpi-1))/slope;
end

%Set up plots
fdvi=fdv;
ntpi=101;
fdvi.taxis=linspace(fdv.taxis(1),fdv.taxis(end),101);
fdvi.ntp=ntpi;
fdvi.NFlips=ntpi;
fdvi.TR=(fdv.taxis(end)-fdv.taxis(1))/(ntpi-1)*ones(1,ntpi);
fdvi.VIFP=zeros(1,ntpi);
fdvi.VIFL=zeros(1,ntpi);
fdvi.FlipAngle=zeros(4,ntpi);

% get curves for best fit:
fdvi.knownvals(find(strcmp(fdv.knowns,'kpl')))=kplvec(minresi);
[EV,IV,vols,Mzevmin,Mzivmin] = P2LCMSv4(fitsvec(minresi,:),fdvi);

%Get interpolated curves at 1% threshold:
fdvi.knownvals(find(strcmp(fdv.knowns,'kpl')))=kplest;
[EV,IV,vols,Mzev,Mziv] = P2LCMSv4(parmsest,fdvi);
fdv.knownvals(end)=kplest;
residatest=P2LCMSv4Err(parmsest,fdv);
%Caclulate residual using these estimated parameters:
resnatest=sum(residatest.*residatest); 


endi = 45; % dont need to plot kpl>1 here...

figure(20+tflag)
set(gcf,'Name',sprintf('Mass Spec - %s',tumor{tflag}));
subplot(2,1,1)
f3a1h=plot(fdv.taxis,fdv.data(3,:),'gx',...
    fdv.taxis,fdv.data(4,:),'bx',...
    fdvi.taxis,Mzev(3,:),'g-',...
    fdvi.taxis,Mzev(4,:),'b-',...
    fdvi.taxis,Mzevmin(3,:),'g:',...
    fdvi.taxis,Mzevmin(4,:),'b:',...
    'linewidth',2,'MarkerSize',8);
f3a1h(5).LineWidth=1;
f3a1h(6).LineWidth=1;
if tflag==1
    ylabel('Label Count (arb)','FontSize',18,'FontWeight','bold')
end
xlabel('Time (s)','FontSize',18,'FontWeight','bold')
title(sprintf('%s',tumor{tflag}),'FontSize',22,'FontWeight','bold')
if tflag==2
    legend('IC M+3 Pyruvate Observed','IC M+3 Lactate Observed',...
        'IC M+3 Pyruvate Fit @ Thresh','IC M+3 Lactate Fit @ Tresh', ...
        'IC M+3 Pyruvate Best Fit','IC M+3 Lactate Best Fit', ...
        'Location','Best')
end
grid on
axis tight

subplot(2,1,2)

%kpl vs residual 
f3a2h=plot(kplvec(1:endi),residvec(1:endi),'linewidth',2);
grid on
if tflag==1
    ylabel('Residual (arb)','FontSize',18,'FontWeight','bold')
end
xlabel('k_P_L (s^{-1})','FontSize',18,'FontWeight','bold')
axis tight

%note kpl at which residual residual decreases by less than one percent
hold on
plot(kplest,resnatest,'rx','linewidth',3,'MarkerSize',16)
text(kplest,resnatest+max(residvec)/9,...
    sprintf('k_{PL} = %5.3f s^{-1}',kplest))
hold off

%Set figure size on scren:
set(gcf,'position',[0 60 600 600])

grid on
print(sprintf('Fig3-MS-%s',tumor{tflag}),'-depsc')
exportgraphics(gcf,sprintf('Fig3-MS-%s.png',tumor{tflag}))
xlim([0 1])

%Copy figure to clipboard
copygraphics(gcf)


