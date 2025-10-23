%% Make Fig 7
clear
close all

%To prepare fig 7 - import overlays, resize to 2" width @ 300dpi

hpdataprefn=fullfile('data','HPMRI_ATC_pre.mat');
pk1pcprefn=fullfile('data','HPMRIPK_pre_1PC.20250915.mat');
pk2pcprefn=fullfile('data','HPMRIPK_pre_2PC.20250915.mat');
pk3pcprefn=fullfile('data','HPMRIPK_pre_3PC.20250915.mat');

%Position of figures for proper formatting:
figpos=[1930 210 900 750];
figposh=[1921 200 560 420]; %for histograms
figpos=[1930 -200 900 750]; % at work
figposh=[1921 -700 560 420];
resave=1;

%% Figure 7b1 - baseline precursor/product (kpla):

load(hpdataprefn,'hptroif','hph','bgims','bgi','bgh');
load(pk1pcprefn,'kpla');

%Representative slice:
sno=3;
%Set zoom region:
zrx=[154-5 384+25];
zry=[113-25 315+10];
%Ignore negative values (noise) and those voxels less than 40% tumor:
kplathresh=0.001;
kplamask=((kpla>kplathresh).*(hptroif>0.4))>0; 
nbins=25;
maxpix=[8 10];

figure(721)
set(gcf,'Name','Baseline kpla overlay')
valnames=[];%'k''''_{PL}';
registeredoverlay2(valnames,kpla(:,:,sno).*kplamask(:,:,sno),hph,1,...
    bgims(:,:,bgi(sno)),bgh,zrx,zry,1);
%truesize(gcf,[450 450])
set(gcf,'position',figpos)
ofn='Fig7b1-Base-kpla-overlay';
if ~exist([ofn '.tif'],'file')||(resave==1)
    exportgraphics(gca,[ofn '.tif']);
    exportgraphics(gcf,[ofn '.png']);
end


figure(722)
set(gcf,'Name','Baseline kpla histogram')
textflag=0;
subplot(2,1,1);p1=0.88;p2=p1-0.18;p3=p2-0.18;fs=12;fw='normal';
histogram(kpla(kplamask),nbins,'BinLimits',[0,0.13],'FaceColor',[0 0 0.5]);
meankpla=mean(kpla(kplamask));
mediankpla=median(kpla(kplamask));
maxkpla=max(kpla(kplamask));
kplaaxis=axis % [-0.0063 0.1313 0 7.0]
kplaaxis(4)=8; % match scale on post
axis(kplaaxis);
if textflag
    text(kplaaxis(2)*0.525,kplaaxis(4)*p1,sprintf('Mean k''''_{PL} = %5.3f s^{-1}',meankpla),...
        'fontsize',fs,'fontweight',fw)
    text(kplaaxis(2)*0.525,kplaaxis(4)*p2,sprintf('Median k''''_{PL} = %5.3f s^{-1}',mediankpla),...
        'fontsize',fs,'fontweight',fw)
    text(kplaaxis(2)*0.525,kplaaxis(4)*p3,sprintf('Max k''''_{PL} = %5.3f s^{-1}',maxkpla),...
        'fontsize',fs,'fontweight',fw)
end
ylabel('Count','fontsize',14)
xlabel('k''''_{PL} (s^{-1})','fontsize',14)
set(gca,'fontsize',12)
set(gcf,'position',figposh)
ofn='FigS321v3-Base-kpla-histogram';
if ~exist([ofn '.eps'],'file')||(resave==1)
    print(ofn,'-depsc');
    exportgraphics(gcf,[ofn '.png']);
end

%% Figure 7a1 - baseline normalized lactate ratio
load(hpdataprefn,'spyr','pyri','slac','laci');

nLac=sum(slac(:,:,:,laci(1):laci(2)),4)./(sum(spyr(:,:,:,pyri(1):pyri(2)),4)+sum(slac(:,:,:,laci(1):laci(2)),4));

nlacmask=(kplamask.*(nLac>0))>0;

figure(711)
set(gcf,'Name','nlac overlay')
valnames=[];%'nLac';
registeredoverlay2(valnames,nLac(:,:,sno).*kplamask(:,:,sno),hph,1,...
    bgims(:,:,bgi(sno)),bgh,zrx,zry,1);
set(gcf,'position',figpos)
ofn='Fig7a1-Base-nlac-overlay';
if ~exist([ofn '.tif'],'file')||(resave==1)
    exportgraphics(gcf,[ofn '.tif']);
    exportgraphics(gcf,[ofn '.png']);
end

figure(712)
set(gcf,'Name','nLac histogram')
%textflag=0;
subplot(2,1,1);
p1=0.92;p2=p1-0.18;p3=p2-0.18;
f22h=histogram(nLac(nlacmask),nbins,'BinLimits',[0.0,1.0],'FaceColor',[0 0 0.5]);
meannlac=mean(nLac(nlacmask));
mediannlac=median(nLac(nlacmask));
maxnlac=max(nLac(nlacmask));
nlacaxis=axis % [0.0025 0.8275 0 4.0]
%nlacaxis(2)=1.0;
%axis(nlacaxis);
if textflag
    text(nlacaxis(2)*0.0,nlacaxis(4)*p1,sprintf('Mean nLac = %5.3f s^{-1}',meannlac),...
        'fontsize',fs,'fontweight',fw)
    text(nlacaxis(2)*0.525,nlacaxis(4)*p1,sprintf('Median nLac = %5.3f s^{-1}',mediannlac),...
        'fontsize',fs,'fontweight',fw)
    text(nlacaxis(2)*0.525,nlacaxis(4)*p2,sprintf('Max nLac = %5.3f s^{-1}',maxnlac),...
        'fontsize',fs,'fontweight',fw)
end
title('Baseline','fontsize',16)
ylabel('Count','fontsize',14)
xlabel('nLac (s^{-1})','fontsize',14)
set(gca,'fontsize',12)
set(gcf,'position',figposh)
ofn='FigS311v3-Base-nlac-histogram';
if ~exist([ofn '.eps'],'file')||(resave==1)
    print(ofn,'-depsc');
    exportgraphics(gcf,[ofn '.png']);
end

%% Fig 7c1 - baseline two-compartment model.

load(pk2pcprefn,'fdvb','kplbfits','kplbresids','kplvec');

%kpli=find(strcmp('kpl',fdvb.knowns));
vbi=find(strcmp('vb',fdvb.fitvars));
kvei=find(strcmp('kve',fdvb.fitvars));
klpi=find(strcmp('klp',fdvb.fitvars));
%kplb=kplbfits(:,:,:,kpli);
%vb=kplbfits(:,:,:,vbi);

kplb=zeros(size(kplamask));
kplbparms=zeros(16,16,8,4);
vb=zeros(size(kplb));

for iPE=1:size(kplbfits,1)
    for iRO=1:size(kplbfits,2)
        for iSL=1:size(kplbfits,3)
            if kplamask(iPE,iRO,iSL) %SNR,tumor
                %Find point where resid(kpl) is changing by less than 1%:
                residvec=squeeze(kplbresids(iPE,iRO,iSL,:))';
                tmpi=1;
                thresh=(max(residvec)-min(residvec))/100+min(residvec);
                while residvec(tmpi)>thresh
                    tmpi=tmpi+1;
                end
                %Assume residual is linear over this interval. Estimate kpl at the
                %crossing point:
                slope=(residvec(tmpi)-residvec(tmpi-1))/(kplvec(tmpi)-kplvec(tmpi-1));
                kplest=kplvec(tmpi-1) + (thresh-residvec(tmpi-1))/slope;
                kpll=kplest;
                %And find the rest of the parameters too:
                NParm=size(kplbfits,5);
                parmsest=zeros(1,NParm);
                fitsvec=squeeze(kplbfits(iPE,iRO,iSL,:,:));
                for qq=1:NParm
                    slope=(residvec(tmpi)-residvec(tmpi-1))/(fitsvec(tmpi,qq)-fitsvec(tmpi-1,qq));
                    parmsest(qq)=fitsvec(tmpi-1,qq) + (thresh-residvec(tmpi-1))/slope;
                end
                %Calculate SOS of nuisance parms for this parm set:
                sosll=sqrt((parmsest(kvei)/0.2)^2 + (parmsest(klpi))^2 + (kpll/1.5)^2);
                %Calculate SOS of nuisance parms for equivalent parm sets:
                kves=squeeze(kplbfits(iPE,iRO,iSL,:,kvei))';
                klps=squeeze(kplbfits(iPE,iRO,iSL,:,klpi))';
                sos=sqrt(((kves/0.2).^2)+(klps.^2)+((kplvec/1.5).^2));
                [minsos,minsosi]=min(sos(tmpi:end));
                if minsos<sosll
                    kplest=kplvec(minsosi+tmpi-1);
                    parmsest=squeeze(kplbfits(iPE,iRO,iSL,minsosi+tmpi-1,:))';
                else
                    minsos=sosll;
                end
                kplb(iPE,iRO,iSL)=kplest;
                vb(iPE,iRO,iSL)=parmsest(vbi);
                kplbparms(iPE,iRO,iSL,:)=(parmsest);
            end %mask
        end % Slice Loop
    end % RO Loop
end % PE Loop

kplbmask=(kplamask.*(vb<0.5))>0;

figure(731)
set(gcf,'Name','Baseline kplb overlay')
valnames=[];%'k''_{PL}';
registeredoverlay2(valnames,kplb(:,:,sno).*kplbmask(:,:,sno),hph,1,...
    bgims(:,:,bgi(sno)),bgh,zrx,zry,1);
set(gcf,'position',figpos)
ofn='Fig7c1v3-Base-kplb-overlay';
if ~exist([ofn '.tif'],'file')||(resave==1)
    exportgraphics(gcf,[ofn '.tif']);
    exportgraphics(gcf,[ofn '.png']);
end

figure(732)
set(gcf,'Name','Baseline kplb histogram')
%textflag=0;
subplot(2,1,1);
p1=0.88;p2=p1-0.20;p3=p2-0.20;
histogram(kplb(kplbmask),nbins,'BinLimits',[0,1.01],'FaceColor',[0 0 0.5]);
meankplb=mean(kplb(kplbmask));
mediankplb=median(kplb(kplbmask));
maxkplb=max(kplb(kplbmask));
kplbaxis=axis; % [-.006 0.346 0 4.0]
%make axes match between pre/post
kplbaxis(4)=10.0;
axis(kplbaxis);
if textflag
    text(kplbaxis(2)*0.52,kplbaxis(4)*p1,sprintf('Mean k''_{PL} = %5.3f s^{-1}',meankplb),...
        'fontsize',fs,'fontweight',fw)
    text(kplbaxis(2)*0.52,kplbaxis(4)*p2,sprintf('Median k''_{PL} = %5.3f s^{-1}',mediankplb),...
        'fontsize',fs,'fontweight',fw)
    text(kplbaxis(2)*0.52,kplbaxis(4)*p3,sprintf('Max k''_{PL} = %5.3f s^{-1}',maxkplb),...
        'fontsize',fs,'fontweight',fw)
else
    meankplb
    mediankplb
    maxkplb
end
ylabel('Count','fontsize',14)
xlabel('k''_{PL} (s^{-1})','fontsize',14)
set(gca,'fontsize',12)
set(gcf,'position',figposh)
ofn='FigS331v3-Base-kplb-histogram';
if ~exist([ofn '.eps'],'file')||(resave==1)
    print(ofn,'-depsc');
    exportgraphics(gcf,[ofn '.png']);
end


%% Fig 7c1 - baseline three-compartment model.

% load('230427_ModelC003Av4.mat');
% kplveclen=length(kplvec);

%load('S3PC3kplsweep5000.mat'); 
load(pk3pcprefn,'fdvc','kplcfits','kplcresids','kplvec');

kplcmask=kplbmask;
kplc=zeros(16,16,8);
kplcparms=zeros(16,16,8,7);

vbi=find(strcmp(fdvc.fitvars,'vb'));
kvei=find(strcmp(fdvc.fitvars,'kve'));
kecpi=find(strcmp(fdvc.fitvars,'kecp'));
kecli=find(strcmp(fdvc.fitvars,'kecl'));
klpi=find(strcmp(fdvc.fitvars,'klp'));
vefi=find(strcmp(fdvc.fitvars,'vef'));

%Extract kpl value from sweeps at each voxel:
for ii=1:size(kplcfits,1)
    for jj=1:size(kplcfits,2)
        for kk=1:size(kplcfits,3)
            if kplcmask(ii,jj,kk)
                %Find point where resid(kpl) is changing by less than 1%:
                maxr=max(kplcresids(ii,jj,kk,:));
                minr=min(kplcresids(ii,jj,kk,:));
                thresh=minr+(maxr-minr)/100;
                tmpi=1;
                while(kplcresids(ii,jj,kk,tmpi)>(thresh))&&(tmpi<length(kplvec))
                    tmpi=tmpi+1;
                end
                if tmpi==1
                    kplc(ii,jj,kk)=kplvec(tmpi);
                    parmsest=kplcfits(ii,jj,kk,tmpi);
                elseif tmpi==length(kplvec)
                    kplc(ii,jj,kk)=kplvec(tmpi);
                    parmsest=kplcfits(ii,jj,kk,tmpi);
                else
                    %Assume residual is linear over this interval
                    slope=(kplcresids(ii,jj,kk,tmpi)-kplcresids(ii,jj,kk,tmpi-1))/(kplvec(tmpi)-kplvec(tmpi-1));
                    %kplc(ii,jj,kk)=kplvec(tmpi-1) + (thresh-kplcresids(ii,jj,kk,tmpi-1))/slope;
                    kpll=kplvec(tmpi-1) + (thresh-kplcresids(ii,jj,kk,tmpi-1))/slope;
                    %Same estimate for individual parms:
                    NParm=size(kplcfits,5);
                    parmsest=zeros(1,NParm);
                    fitsvec=squeeze(kplcfits(ii,jj,kk,:,:));
                    for qq=1:NParm
                        slope=(kplcresids(ii,jj,kk,tmpi)-kplcresids(ii,jj,kk,tmpi-1))/(fitsvec(tmpi,qq)-fitsvec(tmpi-1,qq));
                        parmsest(qq)=fitsvec(tmpi-1,qq) + (thresh-kplcresids(ii,jj,kk,tmpi-1))/slope;
                    end
                    sosll=sqrt(((parmsest(kvei)/0.2).^2)+((parmsest(kecpi)/0.2).^2)+((parmsest(kecli)/0.2).^2)+...
                        ((parmsest(vefi)-0.5).^2)+((parmsest(klpi)/0.2).^2)+((kpll/1.5).^2));
                    kves=squeeze(kplcfits(ii,jj,kk,:,kvei))';
                    kecps=squeeze(kplcfits(ii,jj,kk,:,kecpi))';
                    kecls=squeeze(kplcfits(ii,jj,kk,:,kecli))';
                    klps=squeeze(kplcfits(ii,jj,kk,:,klpi))';
                    vefs=squeeze(kplcfits(ii,jj,kk,:,vefi))';
                    sos=sqrt(((kves/0.2).^2)+((kecps/0.2).^2)+((kecls/0.2).^2)+((vefs-0.5).^2)+((klps).^2)+((kplvec/1.5).^2));
                    [minsos,minsosi]=min(sos(tmpi:end));
                    if minsos<sosll
                        kplc(ii,jj,kk)=kplvec(minsosi+tmpi-1);
                        parmsest=squeeze(kplcfits(ii,jj,kk,minsosi+tmpi-1,:))';
                    else
                        kplc(ii,jj,kk)=kpll;
                        minsos=sosll;
                    end
                    kplcparms(ii,jj,kk,:)=parmsest;
                end %calculation of kplc
            end % mask
        end %slice
    end %rows
end %columns

figure(741)
set(gcf,'Name','Basline kplc overlay')
valnames=[];%'k_{PL}'
registeredoverlay2(valnames,kplc(:,:,sno).*kplcmask(:,:,sno),hph,1,...
    bgims(:,:,bgi(sno)),bgh,zrx,zry,1);
%truesize(gcf,[450 450])
%print('test','-djpeg')
set(gcf,'position',figpos)
ofn='Fig7d1v3-Base-kplc-overlay';
if ~exist([ofn '.tif'],'file')||(resave==1)
    exportgraphics(gcf,[ofn '.tif']);
    exportgraphics(gcf,[ofn '.png']);
end

figure(742)
set(gcf,'Name','Basline kplc histogram')
subplot(2,1,1)
p1=0.89;p2=p1-0.20;p3=p2-0.20;
histogram(kplc(kplcmask),nbins,'BinLimits',[0,2.0],'FaceColor',[0 0 0.5]);
meankplc=mean(kplc(kplcmask));
mediankplc=median(kplc(kplcmask));
maxkplc=max(kplc(kplcmask));
kplcaxis=axis; % [-.006 0.346 0 4.0]
kplcaxis(4)=11;
axis(kplcaxis)
if textflag
    text(kplcaxis(2)*0.45,kplcaxis(4)*p1,sprintf('Mean k_{PL} = %5.3f s^{-1}',meankplc),...
        'fontsize',fs,'fontweight',fw)
    text(kplcaxis(2)*0.45,kplcaxis(4)*p2,sprintf('Median k_{PL} = %5.3f s^{-1}',mediankplc),...
        'fontsize',fs,'fontweight',fw)
    text(kplcaxis(2)*0.45,kplcaxis(4)*p3,sprintf('Max k_{PL} = %5.3f s^{-1}',maxkplc),...
        'fontsize',fs,'fontweight',fw)
else
    meankplc
    mediankplc
    maxkplc
end
ylabel('Count','fontsize',14)
xlabel('k_{PL} (s^{-1})','fontsize',14)
set(gca,'fontsize',12)
set(gcf,'position',figposh)
ofn='FigS341v3-Base-kplc-histogram';
if ~exist([ofn '.eps'],'file')||(resave==1)
    print(ofn,'-depsc');
    exportgraphics(gcf,[ofn '.png']);
end

save('QHPMRIPre.mat','kpla','kplamask','nLac','nlacmask',...
    'fdvb','kplb','kplbparms','kplbmask',...
    'fdvc','kplc','kplcparms','kplcmask');
