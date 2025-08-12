%% Make Fig 7
clear
close all

%To prepare fig 7 - import overlays, resize to 2" width @ 300dpi

hpdatapostfn=fullfile('data','HPMRI_ATC_post.mat');
pk1pcpostfn=fullfile('data','HPMRIPK_post_1PC.mat');
pk2pcpostfn=fullfile('data','HPMRIPK_post_2PC.mat');
pk3pcpostfn=fullfile('data','HPMRIPK_post_3PC.mat');

%Position of figures for proper formatting:
figpos=[1930 210 900 750];
figposh=[1921 200 560 420]; %for histograms
resave=1;

%% Figure 7b2 - +8d precursor/product (kpla):

load(hpdatapostfn,'hptroif','hph','bgims','bgi','bgh');
load(pk1pcpostfn,'kpla');

%Representative slice:
sno=3;
%Set zoom region:
zrx=[154-5 384+25]-25;
zry=[113-25 315+10];
%Ignore negative values (noise) and those voxels less than 40% tumor:
kplathresh=0.001;
kplamask=((kpla>kplathresh).*(hptroif>0.4))>0; 
nbins=25;
maxpix=[6 9];

figure(723)
set(gcf,'Name','+8d kpla overlay')
valnames=[];%'k''''_{PL}';
registeredoverlay2(valnames,kpla(:,:,sno).*kplamask(:,:,sno),hph,1,...
    bgims(:,:,bgi(sno)),bgh,zrx,zry,1);
%truesize(gcf,[450 450])
set(gcf,'position',figpos)
ofn='Fig7b3-Post-kpla-overlay';
if ~exist([ofn '.tif'],'file')||(resave==1)
    exportgraphics(gca,[ofn '.tif']);
    exportgraphics(gcf,[ofn '.jpg']);
end


figure(724)
set(gcf,'Name','+8d kpla histogram')
textflag=0;
subplot(2,1,1);p1=0.77;p2=p1-0.18;;p3=p2-0.18;fs=12;fw='normal';
histogram(kpla(kplamask),nbins,'BinLimits',[0,0.13],'FaceColor',[0 0 0.5]);
meankpla=mean(kpla(kplamask));
mediankpla=median(kpla(kplamask));
maxkpla=max(kpla(kplamask));
kplaaxis=axis; 
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
ofn='FigS322-Post-kpla-histogram';
if ~exist([ofn '.eps'],'file')||(resave==1)
    print(ofn,'-depsc');
    exportgraphics(gcf,[ofn '.jpg']);
end

%% Figure 7a2 - +8d normalized lactate ratio
load(hpdatapostfn,'spyr','pyri','slac','laci');

nLac=sum(slac(:,:,:,laci(1):laci(2)),4)./(sum(spyr(:,:,:,pyri(1):pyri(2)),4)+sum(slac(:,:,:,laci(1):laci(2)),4));

nlacmask=(kplamask.*(nLac>0))>0;

figure(713)
set(gcf,'Name','+8d nLac overlay')
valnames=[];%'nLac';
registeredoverlay2(valnames,nLac(:,:,sno).*kplamask(:,:,sno),hph,1,...
    bgims(:,:,bgi(sno)),bgh,zrx,zry,1);
set(gcf,'position',figpos)
ofn='Fig7a3-Post-nLac-overlay';
if ~exist([ofn '.tif'],'file')||(resave==1)
    exportgraphics(gcf,[ofn '.tif']);
    exportgraphics(gcf,[ofn '.jpg']);
end

figure(714)
set(gcf,'Name','+8d nLac histogram')
%textflag=0;
subplot(2,1,1);
p1=0.92;p2=p1-0.18;p3=p2-0.18;
histogram(nLac(nlacmask),nbins,'BinLimits',[0.0,1.0],'FaceColor',[0 0 0.5]);
meannlac=mean(nLac(nlacmask));
mediannlac=median(nLac(nlacmask));
maxnlac=max(nLac(nlacmask));
nlacaxis=axis % [0.0025 0.8275 0 4.0]
nlacaxis(4)=5; % fix problem - y axis scales when x label set..?
if textflag
    text(nlacaxis(2)*0.55,nlacaxis(4)*p1,sprintf('Mean nLac = %5.3f s^{-1}',meannlac),...
        'fontsize',fs,'fontweight',fw)
    text(nlacaxis(2)*0.55,nlacaxis(4)*p2,sprintf('Median nLac = %5.3f s^{-1}',mediannlac),...
        'fontsize',fs,'fontweight',fw)
    text(nlacaxis(2)*0.55,nlacaxis(4)*p3,sprintf('Max nLac = %5.3f s^{-1}',maxnlac),...
        'fontsize',fs,'fontweight',fw)
end
title('+8 Days','fontsize',16)
ylabel('Count','fontsize',14)
xlabel('nLac (s^{-1})','fontsize',14)
set(gca,'fontsize',12)
set(gcf,'position',figposh)
axis(nlacaxis);
ofn='FigS312-Post-nLac-histogram';
if ~exist([ofn '.eps'],'file')||(resave==1)
    print(ofn,'-depsc');
    exportgraphics(gcf,[ofn '.jpg']);
end

%% Fig 7c2 - +8d two-compartment model.

load(pk2pcpostfn,'fdvb','kplbfits');

kpli=find(strcmp(fdvb.fitvars,'kpl'));
kplb=kplbfits(:,:,:,kpli);
vbi=find(strcmp(fdvb.fitvars,'vb'));
vb=kplbfits(:,:,:,vbi);

kplbmask=(kplamask.*(vb<0.5))>0;

figure(733)
set(gcf,'Name','+8d kplb overlay')
valnames=[];%'k''_{PL}';
registeredoverlay2(valnames,kplb(:,:,sno).*kplbmask(:,:,sno),hph,1,...
    bgims(:,:,bgi(sno)),bgh,zrx,zry,1);
set(gcf,'position',figpos)
ofn='Fig7c3-Post-kplb-overlay';
if ~exist([ofn '.tif'],'file')||(resave==1)
    exportgraphics(gcf,[ofn '.tif']);
    exportgraphics(gcf,[ofn '.jpg']);
end

figure(734)
set(gcf,'Name','+8d kplb histogram')
%textflag=0;
subplot(2,1,1);
p1=0.89;p2=p1-0.20;p3=p2-0.20;
histogram(kplb(kplbmask),nbins,'BinLimits',[0,1.01],'FaceColor',[0 0 0.5]);
meankplb=mean(kplb(kplbmask));
mediankplb=median(kplb(kplbmask));
maxkplb=max(kplb(kplbmask));
kplbaxis=axis; % [-.05 1.05 0 10.0]



if textflag
    text(kplbaxis(2)*0.525,kplbaxis(4)*p1,sprintf('Mean k''_{PL} = %5.3f s^{-1}',meankplb),...
        'fontsize',fs,'fontweight',fw)
    text(kplbaxis(2)*0.525,kplbaxis(4)*p2,sprintf('Median k''_{PL} = %5.3f s^{-1}',mediankplb),...
        'fontsize',fs,'fontweight',fw)
    text(kplbaxis(2)*0.525,kplbaxis(4)*p3,sprintf('Max k''_{PL} = %5.3f s^{-1}',maxkplb),...
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
ofn='FigS332-Post-kplb-histogram';
if ~exist([ofn '.eps'],'file')||(resave==1)
    print(ofn,'-depsc');
    exportgraphics(gcf,[ofn '.jpg']);
end


%% Fig 7c2 - +8d three-compartment model.

load(fullfile('data','HPMRIPK_post_3PC.mat'));

kplveclen=20;
kplvec=logspace(log10(0.1),log10(2),kplveclen);

kplcmask=kplbmask;
kplc=zeros(16,16,8);

vbi=find(strcmp(fdvc.fitvars,'vb'));
kvei=find(strcmp(fdvc.fitvars,'kve'));
kecpi=find(strcmp(fdvc.fitvars,'kecp'));

%Extract kpl value from sweeps at each voxel:
for ii=1:size(kplcfits,1)
    for jj=1:size(kplcfits,2)
        for kk=1:size(kplcfits,3)
            if kplcmask(ii,jj,kk)
                maxr=max(kplcresids(ii,jj,kk,:));
                minr=min(kplcresids(ii,jj,kk,:));
                thresh=minr+(maxr-minr)/100;
                tmpi=1;
                while(kplcresids(ii,jj,kk,tmpi)>(thresh))&&(tmpi<kplveclen)
                    tmpi=tmpi+1;
                end
                if tmpi==1
                    kplc(ii,jj,kk)=kplvec(1);
                elseif tmpi==kplveclen
                    kplc(ii,jj,kk)=kplvec(kplveclen);
                else
                    %Assume residual is linear over this interval
                    slope=(kplcresids(ii,jj,kk,tmpi)-kplcresids(ii,jj,kk,tmpi-1))/(kplvec(tmpi)-kplvec(tmpi-1));
                    kplc(ii,jj,kk)=kplvec(tmpi-1) + (thresh-kplcresids(ii,jj,kk,tmpi-1))/slope;
                end %calculation of kplc
            end % mask
        end %slice
    end %rows
end %columns

figure(743)
set(gcf,'Name','+8d kplc overlay')
valnames=[];%'k_{PL}';
registeredoverlay2(valnames,kplc(:,:,sno).*kplcmask(:,:,sno),hph,1,...
    bgims(:,:,bgi(sno)),bgh,zrx,zry,1);
set(gcf,'position',figpos)
ofn='Fig7d3-Post-kplc-overlay';
if ~exist([ofn '.tif'],'file')||(resave==1)
    exportgraphics(gcf,[ofn '.tif']);
    exportgraphics(gcf,[ofn '.jpg']);
end

figure(744)
set(gcf,'Name','+8d kplc histogram')
subplot(2,1,1)
p1=0.89;p2=p1-0.20;p3=p2-0.20;
histogram(kplc(kplcmask),nbins,'BinLimits',[0,2.0],'FaceColor',[0 0 0.5]);
meankplc=mean(kplc(kplcmask));
mediankplc=median(kplc(kplcmask));
maxkplc=max(kplc(kplcmask));
kplcaxis=axis; % [-.006 0.346 0 4.0]
if textflag
    text(kplcaxis(2)*0.5,kplcaxis(4)*p1,sprintf('Mean k_{PL} = %5.3f s^{-1}',meankplc),...
        'fontsize',fs,'fontweight',fw)
    text(kplcaxis(2)*0.5,kplcaxis(4)*p2,sprintf('Median k_{PL} = %5.3f s^{-1}',mediankplc),...
        'fontsize',fs,'fontweight',fw)
    text(kplcaxis(2)*0.5,kplcaxis(4)*p3,sprintf('Max k_{PL} = %5.3f s^{-1}',maxkplc),...
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
ofn='FigS342-Post-kplc-histogram';
if ~exist([ofn '.eps'],'file')||(resave==1)
    print(ofn,'-depsc');
    exportgraphics(gcf,[ofn '.jpg']);
end