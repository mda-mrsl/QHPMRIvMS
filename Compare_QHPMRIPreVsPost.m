%% Compare HP MRI (baseline vs +8d) 2020-0482-003
clear

infn1='QHPMRIPre.mat';
infn2='QHPMRIPost.mat';

%% nLac:

load(infn1,'nLac','nlacmask')
nLac1=nLac;
nlacmask1=nlacmask; 

load(infn2,'nLac','nlacmask')

nlac_voxels=[sum(nlacmask1,'all') sum(nlacmask,'all')]
% 32 vs 32
nlac_means=[mean(nLac1(nlacmask1)) mean(nLac(nlacmask))] 
% 0.3610 vs 0.3113
nlac_meds=[median(nLac1(nlacmask1)) median(nLac(nlacmask))]
% 0.4179 vs 0.3152
nlac_maxs=[max(nLac1(nlacmask1)) max(nLac(nlacmask))]
% 0.6844 vs 0.6601
[~,pnlac]=ttest2(nLac1(nlacmask1),nLac(nlacmask),'tail','right')
% pnlac = 0.1240
rspnlac=ranksum(nLac1(nlacmask1),nLac(nlacmask),'tail','right')
% rspnlac = 0.1255

%% AUC Ratio:
AUCr1=nLac1(nlacmask1)./(1-nLac1(nlacmask1));
AUCr2=nLac(nlacmask)./(1-nLac(nlacmask));

%voxels: 32 vs 32 voxels, same as nLac
AUCr_means=[mean(AUCr1) mean(AUCr2)]
% 0.7179 vs 0.5263
AUCr_meds=[median(AUCr1) median(AUCr2)]
% 0.7184 vs 0.4604
AUCr_maxs=[max(AUCr1) max(AUCr2)]
% 2.1687 vs 1.9420
[~,paucr]=ttest2(AUCr1,AUCr2,'tail','right')
% paucr = 0.0553
rspaucr=ranksum(AUCr1,AUCr2,'tail','right')
% rspaucr = 0.1255


%% kpla

load(infn1,'kpla','kplamask')
kpla1=kpla;
kplamask1=kplamask; 

load(infn2,'kpla','kplamask')
kpla_voxels=[sum(kplamask1,'all') sum(kplamask,'all')]
% 33 vs 32
kpla_means=[mean(kpla1(kplamask1)) mean(kpla(kplamask))]
% 0.0329 vs 0.0183
kpla_meds=[median(kpla1(kplamask1)) median(kpla(kplamask))]
% 0.0341 vs 0.0171
kpla_maxs=[max(kpla1(kplamask1)) max(kpla(kplamask))]
% 0.1150 vs 0.0452
[~,pkpla]=ttest2(kpla1(kplamask1),kpla(kplamask),'tail','right')
% pkpla = 0.0046 **
rspkpla=ranksum(kpla1(kplamask1),kpla(kplamask),'tail','right')
% rspkpla = 0.0316 *

%% two-compartment model
load(infn1,'fdvb','kplb','kplbparms','kplbmask')
kplb1=kplb;
vbb1=kplbparms(:,:,:,find(strcmp(fdvb.fitvars,'vb')));
kveb1=kplbparms(:,:,:,find(strcmp(fdvb.fitvars,'kve')));
kplbmask1=kplbmask;

load(infn2,'fdvb','kplb','kplbparms','kplbmask')
kplb2=kplb;
vbb2=kplbparms(:,:,:,find(strcmp(fdvb.fitvars,'vb')));
kveb2=kplbparms(:,:,:,find(strcmp(fdvb.fitvars,'kve')));
kplbmask2=kplbmask;

kplb_voxels=[sum(kplbmask1,'all') sum(kplbmask2,'all')]
% 31 vs 30
kplb_means=[mean(kplb1(kplbmask1)) mean(kplb2(kplbmask2))]
% 0.1401 vs 0.1228
kplb_meds=[median(kplb1(kplbmask1)) median(kplb2(kplbmask2))]
% 0.1269 vs 0.0905
kplb_maxs=[max(kplb1(kplbmask1)) max(kplb2(kplbmask2))]
% 0.6406 vs 0.6550
[~,pkplb]=ttest2(kplb1(kplbmask1),kplb2(kplbmask2),'tail','right')
% pkplb = 0.2995
rspkplb=ranksum(kplb1(kplbmask1),kplb2(kplbmask2),'tail','right')
% rspkplb = 0.1702

%% three-compartment model

load(infn1,'fdvc','kplc','kplcparms','kplcmask')
kplc1=kplc;
vbc1=kplcparms(:,:,:,find(strcmp(fdvc.fitvars,'vb')));
kvec1=kplcparms(:,:,:,find(strcmp(fdvc.fitvars,'kve')));
kplcmask1=kplcmask;

load(infn2,'fdvc','kplc','kplcparms','kplcmask')
kplc2=kplc;
vbc2=kplcparms(:,:,:,find(strcmp(fdvc.fitvars,'vb')));
kvec2=kplcparms(:,:,:,find(strcmp(fdvc.fitvars,'kve')));
kplcmask2=kplcmask;

kplc_voxels=[sum(kplcmask1,'all') sum(kplcmask2,'all')]
% 31 vs 30
kplc_means=[mean(kplc1(kplcmask1)) mean(kplc2(kplcmask2))]
% 1.2992 vs 0.9320
kplc_meds=[median(kplc1(kplcmask1)) median(kplc2(kplcmask2))]
% 1.3231 vs 0.5790
kplc_maxs=[max(kplc1(kplcmask1)) max(kplc2(kplcmask2))]
% 2.00 vs 2.00
[~,pkplc]=ttest2(kplc1(kplcmask1),kplc2(kplcmask2),'tail','right')
%  pkplc = 0.0212 *

rspkplc=ranksum(kplc1(kplcmask1),kplc2(kplcmask2),'tail','right')
% rspkplc = 0.0174 *

%% Perfusion terms from HP

vbb_means=[mean(vbb1(kplbmask1)) mean(vbb2(kplbmask2))]
% 0.111 vs 0.0861
[~,pvbb]=ttest2(vbb1(kplbmask1),vbb2(kplbmask2),'tail','right')
% pvbb=0.1790
rspvbb=ranksum(vbb1(kplbmask1),vbb2(kplbmask2),'tail','right')
% rspvbb = 0.2747

kveb_means=[mean(kveb1(kplbmask1)) mean(kveb2(kplbmask2))]
% 0.0206 vs 0.0101
[~,pkveb]=ttest2(kveb1(kplbmask1),kveb2(kplbmask2),'tail','right')
% p=0.0115 *
rspkveb=ranksum(kveb1(kplbmask1),kveb2(kplbmask2),'tail','right')
% p=0.0088 **

vbc_means=[mean(vbc1(kplcmask1)) mean(vbc2(kplcmask2))]
% 0.0443 vs 0.0287
[~,pvbc]=ttest2(vbc1(kplcmask1),vbc2(kplcmask2),'tail','right')
% p=0.3236
rspvbc=ranksum(vbc1(kplcmask1),vbc2(kplcmask2),'tail','right')
% p=0.5877

kvec_means=[mean(kvec1(kplcmask1)) mean(kvec2(kplcmask2))]
% 0.0632 vs 0.0263
[~,pkvec]=ttest2(kvec1(kplcmask1),kvec2(kplcmask2),'tail','right')
% p=0.0099 **
rspkvec=ranksum(kvec1(kplcmask1),kvec2(kplcmask2),'tail','right')
% p=0.0234 *

%% Perfusion terms from DCE:

dcesno=7;
load('..\PreDCE.mat');
ktrans1=ktrans_pre_high_res.data(:,:,dcesno).';
dcevb1=vp_pre_high_res.data(:,:,dcesno).';
dceve1=ve_pre_high_res.data(:,:,dcesno).';
load('..\PreDCEROI.mat','dceroi1');
dcemask1=(dceroi1.*(dcevb1<0.5).*(dceve1<1))>0;

load('..\PostDCE.mat');
ktrans2=ktrans_post_high_res.data(:,:,dcesno).';
dcevb2=vp_post_high_res.data(:,:,dcesno).';
dceve2=ve_post_high_res.data(:,:,dcesno).';
load('..\PostDCEROI.mat','dceroi2');
dcemask2=(dceroi2.*(dcevb2<0.5).*(dceve2<1))>0;

dcevb_means=[mean(dcevb1(dcemask1)) mean(dcevb2(dcemask2))]
% 0.0273 vs 0.0024
[~,p_dce_vb]=ttest2(dcevb1(dcemask1),dcevb2(dcemask2),'tail','right')
% p<<0.001
rsp_dce_vb=ranksum(dcevb1(dcemask1),dcevb2(dcemask2),'tail','right')
% p<<0.001

ktrans_means=[mean(ktrans1(dcemask1)) mean(ktrans2(dcemask2))]
% 0.7516 vs 0.1500
[~,p_ktrans]=ttest2(ktrans1(dcemask1),ktrans2(dcemask2),'tail','right')
% p<<0.001
rsp_ktrans=ranksum(ktrans1(dcemask1),ktrans2(dcemask2),'tail','right')
% p<<0.001