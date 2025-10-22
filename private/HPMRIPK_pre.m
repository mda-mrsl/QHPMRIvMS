%% Prepare for analysis

clear
close all

Name='HPMRI_ATC_Pre';
infn=fullfile('data',[Name '.mat']);
vifi=[9 11 3]; %coordinates of pyr vascular input function
pyrfa=20;
lacfa=30;
TR=3.0;

%No need to modify below here...
load(infn);
[NPE,NOP,NSL,NTP]=size(spyr);
dtt=datetime('today');
mfn=mfilename;
wf=pwd;
%create a copy of the analysis script - note handling of private folder
if ~strncmp(mfn,'Live',4) % running section by section - probably debugging
    copyfile(fullfile(wf,'private',sprintf('%s.m',mfn)),fullfile(wf,'private',sprintf('%s.%d%02d%02d.m',mfn,dtt.Year,dtt.Month,dtt.Day)));
else 
    fprintf('%s.m is running in LiveEditor, and was not copied.\n',mfn)
end

%% Precursor-product (aka one-compartment, 1PC, "Model A")
ofn1pc=fullfile('data','HPMRIPK_pre_1PC.mat');
ofn1pc=fullfile('data',sprintf('%s.%d%02d%02d.mat','HPMRIPK_pre_1PC',dtt.Year,dtt.Month,dtt.Day));


% Set "forward variables" that the PK script will need to know:
fdva.fitvars={'kpl'};
fdva.knowns={'T1Lac' 'klp' 'L0'};
fdva.knownvals=[55, 0, 0]; 
fdva.ntp=NTP;
fdva.NSeg=1;
fdva.NFlips=NTP;
fdva.TR=ones(1,NTP)*TR;
fdva.taxis=(0:NTP-1)*TR;
fdva.FlipAngle=repmat([pyrfa;lacfa],[1 NTP]);
fdva.VIFP=zeros(1,NTP);
fdva.Name=Name;
fdva.verbose=0;

%Run the fits:
if ~exist(ofn1pc,'file')
    kpla=zeros(NPE,NOP,NSL);
    for kk=1:NSL
        fprintf('1PC Slice %d\n',kk)
        for jj=1:NOP
            for ii=1:NPE
                if maxpyrsnr(ii,jj,kk)>5
                    pixeldata=[squeeze(spyr(ii,jj,kk,:))' ; squeeze(slac(ii,jj,kk,:))'];
                    fdva.data=pixeldata;
                    [kpl,resid] = lsqnonlin(@(x) P2L1Err(x,fdva),0.01,0,Inf,optimset('display','off'));
                    kpla(ii,jj,kk)=kpl;                
                end
            end
        end
    end
    save(ofn1pc,'fdva','kpla');
end

%% Two-Compartment Model (2PC, Model B)
ofn2pc=fullfile('data','HPMRIPK_pre_2PC.mat');
ofn2pc=fullfile('data',sprintf('%s.%d%02d%02d.mat','HPMRIPK_pre_2PC',dtt.Year,dtt.Month,dtt.Day));

fdvb=fdva;
fdvb.fitvars={'VIFScale' 'vb'    'kve'    'klp'};
%Bounds for optimization:
LB=           [0          0.001   0.001    1e-5];
UB=           [Inf        1.0     0.2      1.0]; 
%Bounds for initial guesses (ignore VIFScale which is set later):
LGB =         [0          0.2     0.2      1e-4];
UGB =         [0          0.2     0.2      0.05];
fdvb.knowns=  {'L0' 'P0' 'T1Pyr' 'T1Lac' 'Gam1' 'Gam2' 'tdel' 'kpl'};
fdvb.knownvals=[ 0    0    64      55      0      0      0     0.2];
%Note VIF is z-magnetization:
fdvb.VIFP=squeeze(spyr(vifi(1),vifi(2),vifi(3),:))'./sind(fdvb.FlipAngle(1,:));
%fdvb.VIFP=fdvb.VIFP/max(fdvb.VIFP);
fdvb.VIFL=zeros(1,NTP);
fdvb.UseVIF=1;

nparmsb=length(fdvb.fitvars);

nguessperkpl=5000;
kplvec=logspace(log10(0.01),log(1),30);

if ~exist(ofn2pc,'file')
    %Run the fits:
    kplbfits=zeros(NPE,NOP,NSL,length(kplvec),nparmsb);
    kplbresids=zeros(NPE,NOP,NSL,length(kplvec));
    maxpyrvif=max(spyr(vifi(1),vifi(2),vifi(3),:));
    for kk=1:NSL
        fprintf('2PC Slice %d\n',kk)
        for jj=1:NOP
            for ii=1:NPE
                if maxpyrsnr(ii,jj,kk)>5
                    pixeldata=[squeeze(spyr(ii,jj,kk,:))' ; squeeze(slac(ii,jj,kk,:))'];
                    fdvb.data=pixeldata;
                    maxpyrsig=max(pixeldata(1,:));
                    for kplvi=1:length(kplvec)
                        fdvb.knownvals(end)=kplvec(kplvi);
                        pfresids=Inf(1,nguessperkpl);
                        pffits=zeros(nguessperkpl,nparmsb);
                        parfor fiti=1:nguessperkpl
                            Guess=10.^(log10(LGB)+rand(1,length(LGB)).*(log10(UGB)-log10(LGB)));
                            Guess(1)=(1+rand(1))*maxpyrsig/Guess(2)/maxpyrvif;
                            try
                                [fits,resid] = lsqnonlin(@(x) P2L2Err(x,fdvb),Guess,LB,UB,optimset('display','off'));
                            catch materr
                                fprintf('Warning: %s\n',materr.identifier)
                                resid=Inf;
                            end    
                            pffits(fiti,:)=fits;
                            pfresids(fiti)=resid;
                        end %nguessperkpl
                        [bestresid,bestresidi]=min(pfresids);
                        kplbresids(ii,jj,kk,kplvi)=bestresid;
                        kplbfits(ii,jj,kk,kplvi,:)=pffits(bestresidi,:);
                    end %kplvec
                end % if maxpyrsnr
            end % ii
        end % jj
    end % kk
    save(ofn2pc,'fdvb','kplvec','kplbfits','kplbresids');
end


%% Three-Compartment Model (3PC, Model C) 

ofn3pc=fullfile('data','HPMRIPK_pre_3PC.mat');
ofn3pc=fullfile('data',sprintf('%s.%d%02d%02d.mat','HPMRIPK_pre_3PC',dtt.Year,dtt.Month,dtt.Day));

fdvc=fdvb;
fdvc.fitvars={'VIFScale' 'vb'    'kve'    'kecp'    'klp'    'kecl'   'vef'};
%Bounds for optimization:
LB=           [0          0.001   0.001    0.001     1e-5      1e-5     0.2];
UB=           [Inf        1.0     0.2      0.5       2         1        0.8]; 
%Bounds for initial guesses: (note VIFScale guess set later)
LGB=          [0.5        0.02    0.01     0.03      1e-3      1e-3     0.2];
UGB=          [2          0.2     0.2      0.3       0.05      0.7      0.8];
linscale =    [1          0       0        0         0         0        1];
fdvc.knowns=  {'Le0' 'Pe0' 'T1Pyr' 'T1Lac' 'Gam1' 'Gam2' 'tdel' 'Pi0' 'Li0' 'kpl'};
fdvc.knownvals=[ 0    0     64      55       0      0      0       0     0     0];
fdvc.FlipAngle(3:4,:)=fdvc.FlipAngle(1:2,:);
%Note VIF is assumed to be z-magnetization:
fdvc.VIFP=squeeze(spyr(vifi(1),vifi(2),vifi(3),:))'./sind(fdvc.FlipAngle(1,:));
%fdvc.VIFP=fdvc.VIFP/max(fdvc.VIFP);
fdvc.VIFL=zeros(1,NTP);
fdvc.UseVIF=1;

nparmsc=length(linscale);

kplveclen=30;
nguessperkpl=5000;
kplvec=logspace(log10(0.1),log10(2),kplveclen);
kplcfits=zeros(NPE,NOP,NSL,kplveclen,nparmsc);
kplcresids=zeros(NPE,NOP,NSL,kplveclen);

%Mask for voxels that are at least 40% tumor by cross-section
tmpmask=(hptroif>0.4).*(maxpyrsnr>5);

if ~exist(ofn3pc,'file')
    voxcount=0;
    maxpyrvif=max(spyr(vifi(1),vifi(2),vifi(3),:));
    for ii=1:NPE
        for jj=1:NOP
            for kk=1:NSL
                if tmpmask(ii,jj,kk)
                    voxcount=voxcount+1;
                    for kplvi=1:kplveclen
                        fprintf('%s - Voxel %d/%d, kpl %d/%d\n',datestr(now),voxcount,sum(sum(sum(tmpmask))),kplvi,kplveclen);
                        fdvc.knownvals(end)=kplvec(kplvi);
                        fdvc.data=[squeeze(spyr(ii,jj,kk,:))' ; squeeze(slac(ii,jj,kk,:))'];
                        maxpyrvox=max(squeeze(spyr(ii,jj,kk,:)));
                        pfresids=ones(1,nguessperkpl)*Inf;
                        pffits=zeros(nguessperkpl,nparmsc);
                        parfor fiti=1:nguessperkpl %parfor goes here
                            linguess=LGB+rand(1,nparmsc).*(UGB-LGB);
                            expguess=10.^(log10(LGB)+rand(1,nparmsc).*(log10(UGB)-log10(LGB)));
                            Guess=linguess.*linscale+expguess.*(1-linscale);
                            Guess(1)=(1+rand(1))*maxpyrvox/Guess(2)/maxpyrvif;
                            try
                                [fits,resid] = lsqnonlin(@(x) P2L3Err(x,fdvc),Guess,LB,UB,optimset('display','off'));
                            catch materr
                                fprintf('Warning: %s\n',materr.identifier)
                                fits=Guess;
                                resid=sum(P2L3Err(Guess,fdvc).^2);
                                if isnan(resid)
                                    resid=Inf;
                                end
                            end
                            %if resid<bestresids(fiti)
                            pffits(fiti,:)=fits;
                            pfresids(fiti)=resid;
                            %end
                        end % nguessperkpl
                        [bestresid, bestresidi]=min(pfresids);
                        kplcresids(ii,jj,kk,kplvi)=bestresid;
                        kplcfits(ii,jj,kk,kplvi,:)=pffits(bestresidi,:);                
                        fprintf('kpl = %8.3e --- bestresid = %7.5f\n\n',kplvec(kplvi),bestresid)
                    end % kplvec
                end % if voxel is in tumor
            end % slice loop (kk)
        end % jj
    end % ii
    save(ofn3pc,'fdvc','kplvec','kplcfits','kplcresids')
end
