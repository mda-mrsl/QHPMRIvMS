%% Compare HP MRI (baseline vs +8d) 2020-0482-003
clear

hpdataprefn=fullfile('data','HPMRI_ATC_pre.mat');
pk1pcprefn=fullfile('data','HPMRIPK_pre_1PC.mat');
pk2pcprefn=fullfile('data','HPMRIPK_pre_2PC.mat');
pk3pcprefn=fullfile('data','HPMRIPK_pre_3PC.mat');

hpdatapostfn=fullfile('data','HPMRI_ATC_post.mat');
pk1pcpostfn=fullfile('data','HPMRIPK_post_1PC.mat');
pk2pcpostfn=fullfile('data','HPMRIPK_post_2PC.mat');
pk3pcpostfn=fullfile('data','HPMRIPK_post_3PC.mat');

sno=3;
kplathresh=0.001;

%% nLac, kpla

%Baseline
load(hpdataprefn,'slac','spyr','hptroif','pyri','laci','maxpyrsnr');
load(pk1pcprefn,'kpla');
kpla1=kpla;
kplamask1=((kpla1>kplathresh).*(hptroif>0.4))>0;
nLac1=sum(slac(:,:,:,laci(1):laci(2)),4)./(sum(spyr(:,:,:,pyri(1):pyri(2)),4)+sum(slac(:,:,:,laci(1):laci(2)),4));
AUC1=sum(slac(:,:,:,laci(1):laci(2)),4)./(sum(spyr(:,:,:,pyri(1):pyri(2)),4));
nlacmask1=(kplamask1.*(nLac1>0))>0;

%+8d
load(hpdatapostfn,'slac','spyr','hptroif','pyri','laci','maxpyrsnr');
load(pk1pcpostfn,'kpla');
kpla2=kpla;
kplamask2=((kpla2>kplathresh).*(hptroif>0.4))>0;
nLac2=sum(slac(:,:,:,laci(1):laci(2)),4)./(sum(spyr(:,:,:,pyri(1):pyri(2)),4)+sum(slac(:,:,:,laci(1):laci(2)),4));
AUC2=sum(slac(:,:,:,laci(1):laci(2)),4)./(sum(spyr(:,:,:,pyri(1):pyri(2)),4));
nlacmask2=(kplamask2.*(nLac2>0))>0;
fprintf('nlac mask: %d vs %d\n', sum(sum(sum(nlacmask1))),sum(sum(sum(nlacmask2))));
fprintf('kpla mask: %d vs %d\n\n', sum(sum(sum(kplamask1))),sum(sum(sum(kplamask2))));
[~,paucs]=ttest2(AUC1(nlacmask1), AUC2(nlacmask2),'tail','right');
%paucs = .055
rspaucs=ranksum(AUC1(nlacmask1), AUC2(nlacmask2),'tail','right');
%rspaucs = 0.125
% % p > 0.05 for both tests.
fprintf('AUCr: ttest2 pval = %5.3f; ranksum pval = %5.3f\n\n',paucs,rspaucs)


nlac_meds=[median(nLac1(nlacmask1)) median(nLac2(nlacmask2))];
nlac_maxs=[max(nLac1(nlacmask1)) max(nLac2(nlacmask2))];
nlac_means=[mean(nLac1(nlacmask1)) mean(nLac2(nlacmask2))];
[~,pnlac]=ttest2(nLac1(nlacmask1),nLac2(nlacmask2),'tail','right');
%pnlac = 0.124
rspnlac=ranksum(nLac1(nlacmask1),nLac2(nlacmask2),'tail','right');
%rspnlac = 0.125
% % p > 0.05 for both tests.
fprintf('nLac: mean %5.3f vs %5.3f; med %5.3f vs %5.3f; max %5.3f vs %5.3f\n ',...
    nlac_means,nlac_meds,nlac_maxs);
fprintf('nLac: ttest2 pval = %5.3f; ranksum pval = %5.3f\n\n',pnlac,rspnlac)

kpla_means=[mean(kpla1(kplamask1)) mean(kpla2(kplamask2))];
kpla_meds=[median(kpla1(kplamask1)) median(kpla2(kplamask2))];
kpla_maxs=[max(kpla1(kplamask1)) max(kpla2(kplamask2))];
fprintf('kpla: mean %5.3f vs %5.3f; med %5.3f vs %5.3f; max %5.3f vs %5.3f\n ',...
    kpla_means,kpla_meds,kpla_maxs);
[~,pkpla]=ttest2(kpla1(kplamask1),kpla2(kplamask2),'tail','right');
%pkpla = 0.005
rspkpla=ranksum(kpla1(kplamask1),kpla2(kplamask2),'tail','right');
%rspkpla = 0.032
% % p < 0.05 for both tests.
fprintf('kpla: ttest2 pval = %5.3f; ranksum pval = %5.3f\n\n',pkpla,rspkpla)

%% two-compartment model

%Baseline:
load(pk2pcprefn,'fdvb','kplbfits');
kplb1=kplbfits(:,:,:,find(strcmp(fdvb.fitvars,'kpl')));
vbb1=kplbfits(:,:,:,find(strcmp(fdvb.fitvars,'vb')));
kveb1=kplbfits(:,:,:,find(strcmp(fdvb.fitvars,'kve')));
kplbmask1=(kplamask1.*(vbb1<0.5))>0; % down to N=31

%+8d:
load(pk2pcpostfn,'fdvb','kplbfits');
kplb2=kplbfits(:,:,:,find(strcmp(fdvb.fitvars,'kpl')));
vbb2=kplbfits(:,:,:,find(strcmp(fdvb.fitvars,'vb')));
kveb2=kplbfits(:,:,:,find(strcmp(fdvb.fitvars,'kve')));
kplbmask2=(kplamask2.*(vbb2<0.5))>0;
%kplbmask2=(kplamask2.*(vbb2<0.5).*(kplb2<0.99))>0;
fprintf('kplb mask: %d vs %d\n', sum(sum(sum(kplbmask1))),sum(sum(sum(kplbmask2))));

kplb_means=[mean(kplb1(kplbmask1)) mean(kplb2(kplbmask2))];
kplb_meds=[median(kplb1(kplbmask1)) median(kplb2(kplbmask2))];
kplb_maxs=[max(kplb1(kplbmask1)) max(kplb2(kplbmask2))];
fprintf('kplb: mean %5.3f vs %5.3f; med %5.3f vs %5.3f; max %5.3f vs %5.3f\n ',...
    kplb_means,kplb_meds,kplb_maxs);
[~,pkplb]=ttest2(kplb1(kplbmask1),kplb2(kplbmask2),'tail','right');
rspkplb=ranksum(kplb1(kplbmask1),kplb2(kplbmask2),'tail','right');
fprintf('kplb: ttest2 pval = %5.3f; ranksum pval = %5.3f\n\n',pkplb,rspkplb)


%% three-compartment model

%Baseline:
load(pk3pcprefn,'fdvc','kplcfits','kplcresids');

kplveclen=20;
kplvec=logspace(log10(0.1),log10(2),kplveclen);

kplcmask1=kplbmask1;
kplc1=zeros(size(kplb1));
vbc1=kplc1;
kvec1=kplc1;

for ii=1:size(kplcfits,1)
    for jj=1:size(kplcfits,2)
        for kk=1:size(kplcfits,3)
            if kplcmask1(ii,jj,kk)
                maxr=max(kplcresids(ii,jj,kk,:));
                minr=min(kplcresids(ii,jj,kk,:));
                thresh=minr+(maxr-minr)/100;
                tmpi=1;
                while(kplcresids(ii,jj,kk,tmpi)>thresh)&&(tmpi<kplveclen)
                    tmpi=tmpi+1;
                end
                if tmpi==1
                    kplc1(ii,jj,kk)=kplvec(1);
                elseif tmpi==kplveclen
                    kplc1(ii,jj,kk)=kplvec(kplveclen);
                else
                    %Assume residual is linear over this interval
                    slope=(kplcresids(ii,jj,kk,tmpi)-kplcresids(ii,jj,kk,tmpi-1))/(kplvec(tmpi)-kplvec(tmpi-1));
                    kplc1(ii,jj,kk)=kplvec(tmpi-1) + (thresh-kplcresids(ii,jj,kk,tmpi-1))/slope;
                end %calculation of kplc1
            end
        end
    end
end

%+8d
load(pk3pcpostfn,'fdvc','kplcfits','kplcresids');

kplveclen=20;
kplvec=logspace(log10(0.1),log10(2),kplveclen);

kplcmask2=kplbmask2;
kplc2=zeros(size(kplb2));

for ii=1:size(kplcfits,1)
    for jj=1:size(kplcfits,2)
        for kk=1:size(kplcfits,3)
            if kplcmask2(ii,jj,kk)
                maxr=max(kplcresids(ii,jj,kk,:));
                minr=min(kplcresids(ii,jj,kk,:));
                thresh=minr+(maxr-minr)/100;
                tmpi=1;
                while(kplcresids(ii,jj,kk,tmpi)>thresh)&&(tmpi<kplveclen)
                    tmpi=tmpi+1;
                end
                if tmpi==1
                    kplc2(ii,jj,kk)=kplvec(1);
                elseif tmpi==kplveclen
                    kplc2(ii,jj,kk)=kplvec(kplveclen);
                else
                    %Assume residual is linear over this interval
                    slope=(kplcresids(ii,jj,kk,tmpi)-kplcresids(ii,jj,kk,tmpi-1))/(kplvec(tmpi)-kplvec(tmpi-1));
                    kplc2(ii,jj,kk)=kplvec(tmpi-1) + (thresh-kplcresids(ii,jj,kk,tmpi-1))/slope;
                end %calculation of kplc1
            end
        end
    end
end

kplc_means=[mean(kplc1(kplcmask1)) mean(kplc2(kplcmask2))];
kplc_meds=[median(kplc1(kplcmask1)) median(kplc2(kplcmask2))];
kplc_maxs=[max(kplc1(kplcmask1)) max(kplc2(kplcmask2))];
fprintf('kplc: mean %5.3f vs %5.3f; med %5.3f vs %5.3f; max %5.3f vs %5.3f\n ',...
    kplc_means,kplc_meds,kplc_maxs);
% Left tail because metabolism is slightly increased in this model, in enhancing voxels.
[~,pkplc]=ttest2(kplc1(kplcmask1),kplc2(kplcmask2),'tail','left');
rspkplc=ranksum(kplc1(kplcmask1),kplc2(kplcmask2),'tail','left');
fprintf('kplc: ttest2 pval = %5.3f; ranksum pval = %5.3f\n\n',pkplc,rspkplc)
