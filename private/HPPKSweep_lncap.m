%%Set up the structure for analyzing HP NMR data.
clear

ofn = fullfile('data','HPNMRData_LnCap.mat'); % save spectra here
dfn = fullfile('data','HPNMRFits_LnCap.mat'); % save fits here

NG=1;
ND=4;

if exist(ofn,'file')
    load(ofn)
else
    Group=struct('Name','LnCap','Data',struct('Indir','tmp','NCells',0,'tstarti',0,'pyri',0,'laci',0));

    % Group 1: LnCap
    
    % NOTE Viability 29% viability.  Prob with heaters?
    Group(1).Data(1).Indir  =fullfile('data','20211005-Lai_PC3_lncap_hth83_cells','3');
    Group(1).Data(1).NCells =(3.05e6)*0.9; % count(cells/mL)*mL
    Group(1).Data(1).lastz=7;   %last point before signal;
    Group(1).Data(1).tstarti=10; % Stable by here.
    Group(1).Data(1).pyri=16345;
    Group(1).Data(1).laci=14689;

    % NOTE Viability 41%
    Group(1).Data(2).Indir  =fullfile('data','20211006-Lai_PC3_lncap_hth83_cells','4');
    Group(1).Data(2).NCells =(5.64e6)*0.9; 
    Group(1).Data(2).lastz=7;   %last point before signal;
    Group(1).Data(2).tstarti=11; % Stable by here.
    Group(1).Data(2).pyri=16359;
    Group(1).Data(2).laci=14700;
    
    % NOTE Viability 32%
    Group(1).Data(3).Indir  =fullfile('data','20211006-Lai_PC3_lncap_hth83_cells','7');
    Group(1).Data(3).NCells =(4.52e6)*0.9;
    Group(1).Data(3).lastz=7;   %last point before signal;
    Group(1).Data(3).tstarti=11; % Stable by here.
    Group(1).Data(3).pyri=16362;
    Group(1).Data(3).laci=14703;
        
    % NOTE Viability 69.7%
    Group(1).Data(4).Indir  =fullfile('data','20211007-Lai_PC3_lncap_hth83_cells','8');
    Group(1).Data(4).NCells =(8.7e6)*0.9;
    Group(1).Data(4).lastz=11;   %last point before signal;
    Group(1).Data(4).tstarti=13; % Stable by here.
    Group(1).Data(4).pyri=16354;
    Group(1).Data(4).laci=14696;
    
    for gn=1:NG
        for dn=1:ND
            Group(gn).Data(dn).Indir

            LB=3; 
            hdr=readBrukerHeader(fullfile(Group(gn).Data(dn).Indir,'acqus'));
            Group(gn).Data(dn).hdr=hdr;
            NOP=hdr.TD/2;
            fax=(0:NOP-1)*hdr.SW_h/NOP-hdr.SW_h/2;
            tax=(0:NOP-1)/hdr.SW_h;
            TR = 0.1612065+NOP/hdr.SW_h+0.03; %hdr.D(2)+NOP/hdr.SW_h+0.03;
            infil=fullfile(Group(gn).Data(dn).Indir,'ser');
            fh=fopen(infil);
            big=fread(fh,inf,'int32');
            fclose(fh);
            big=big(1:2:end)+1i*big(2:2:end);
            big=reshape(big,NOP,[])';
            NR=size(big,1);
            taxis=(0:NR-1)*TR;
            grpdly=round(hdr.GRPDLY)+2;
            big(:,1:end-grpdly+1)=big(:,grpdly:end);
            big=big.*repmat(exp(-LB*tax),NR,1);
            specs=fftshift(fft(big,[],2),2);

            fax=(0:NOP-1)*hdr.SW_h/NOP-hdr.SW_h/2;
            sumspec=sum(specs,1);

            figure(1)
            if Group(gn).Data(dn).pyri==1
                [pyrm, pyri]=max(abs(sumspec).*[zeros(1,15000) ones(1,32768-15000)]);
                plot(1:NOP,abs(sumspec),pyri,abs(sumspec(pyri)),'ro')
                title(sprintf('Pyr max at pyri=%d',pyri));
                Group(gn).Data(dn).pyri=pyri;
                pause;
            end
            if Group(gn).Data(dn).laci==1
                [lacm, laci]=max(abs(sumspec(1:15000)))
                plot(1:NOP,abs(sumspec),laci,abs(sumspec(laci)),'ro')
                title(sprintf('Lac max at laci=%d',laci));
                Group(gn).Data(dn).laci=laci;
                pause;
            end
            if Group(gn).Data(dn).tstarti==1
                [pyrm, tstarti]=max(abs(specs(:,Group(gn).Data(dn).pyri)))
                plot(1:NR,abs(specs(:,Group(gn).Data(dn).pyri)),tstarti,...
                    abs(specs(tstarti,Group(gn).Data(dn).pyri)),'ro');
                title(sprintf('Pyr max at tstarti=%d',tstarti));
                Group(gn).Data(dn).tstarti=tstarti;
                pause;
            end
            plot(1:NOP,abs(sumspec),[Group(gn).Data(dn).pyri Group(gn).Data(dn).laci],...
                abs(sumspec([Group(gn).Data(dn).pyri Group(gn).Data(dn).laci])),'ro');


            %Calculate pyr signal
            spyr=zeros(1,NR);
            tmp=specs*exp(-1i*angle(sumspec(Group(gn).Data(dn).pyri)));
            [~, ~, hhx]=hhfw(real(sumspec*exp(-1i*angle(sumspec(Group(gn).Data(dn).pyri)))),fax);
            for jj=1:NR
                hhfwa = hhfw(real(tmp(jj,:)),fax,hhx);
                spyr(jj)=hhfwa;
            end
            Group(gn).Data(dn).SPyr=spyr/max(spyr);

            %Calculate lac signal
            slac=zeros(1,NR);
            tmp=specs*exp(-1i*angle(sumspec(Group(gn).Data(dn).laci)));
            tmpss=sum(tmp,1);
            [hhfwa, hhfwd, hhx]=hhfw(real(tmpss(Group(gn).Data(dn).laci+(-100:100))),fax(Group(gn).Data(dn).laci+(-100:100)));
            hhx=hhx+Group(gn).Data(dn).laci-100-1;
            % pause
            for jj=1:NR
                hhfwa = hhfw(real(tmp(jj,:)),fax,hhx);
                slac(jj)=hhfwa;
            end
            Group(gn).Data(dn).SLac=slac/max(spyr);
            Group(gn).Data(dn).NLac=sum(slac)/(sum(slac)+sum(spyr));
            Group(gn).Data(dn).taxis=taxis;

            figure(2)
            plot(taxis,spyr,'g',taxis,slac,'b',taxis(Group(gn).Data(dn).tstarti),spyr(Group(gn).Data(dn).tstarti),'ro',...
                taxis(Group(gn).Data(dn).tstarti),slac(Group(gn).Data(dn).tstarti),'ro',...
                taxis(Group(gn).Data(dn).lastz),spyr(Group(gn).Data(dn).lastz),'rx','linewidth',2)
            axis tight;grid on
            ylabel('Relative Amplitude (arb)','fontsize',12)
            xlabel('Time (s)','fontsize',12)
            title(sprintf('Expt #%d: %s',dn,replace(Group(gn).Data(dn).Indir,'_','-')))
            %pause
        end % for dn
    end % for gn
    save(ofn,'Group')
end

    
%% Background variables for Pharmacokinetic Analysis:

celldum=mean([17.42 19.65 14.1 14.8]);
% 1um=1e-4cm
cellrcm=(1e-4)*celldum/2;
vhtcell = 4/3*pi*((cellrcm)^3); % expressed in mL

%According to PROSOL, 90-deg pulse with this probe is 15us @ 57.5W.
PLWRef=57.5; % Watts
%PRef=15; % us
PRef=33.75/2; %from our 1/15/2020 calibration tests
faff = 0.9285; % correction based on 9/17/21 recalibration


%% PK Analyses

if exist(dfn,'file')
    load(dfn)
    refit=questdlg('Found previous fits.  Run more fits?');
else
    refit='Yes';
    totruns=0;
end

if strcmp(refit,'Yes')
    oldruns=totruns;
    nruns=5000; % per kpl value - use 5000
    totruns=oldruns+nruns;
    
    nkpls=38; % close to same kpl resolution as above
    kplvec=logspace(-3,0,nkpls); % from 0.001 to 1
   
    for gn=1:NG
        for dn=1:ND
            fprintf('Starting Group %d Dataset %2d at %s\n',gn,dn,datestr(now))
            FA=(90*sqrt(Group(gn).Data(dn).hdr.PLW(2)/PLWRef)*Group(gn).Data(dn).hdr.P(2)/PRef);

            %Calculate the extracellular volume fraction from data:
            vc= (Group(gn).Data(dn).NCells*vhtcell)/1.181;
            vef=1-vc;

            %Unknowns to be fit:
            fdv.fitvars={'Pe0' 'kecp' 'kecl' 'T1Lac' 'klp'};
            LB         =[ 0     0      0      35      1e-10];
            UB         =[ Inf    Inf    Inf   55      Inf  ]; 

            LGB=[Group(gn).Data(dn).SPyr(Group(gn).Data(dn).tstarti)/vef/sind(FA) ... % Pe0
                LB(2:end)];
            UGB=[Group(gn).Data(dn).SPyr(Group(gn).Data(dn).tstarti)/vef/sind(FA) ...% Pe0
                1 1 UB(4:end-1) 1];

            %Describe acquisition scheme
            NR=length(Group(gn).Data(dn).taxis);
            fdv.ntp=NR-(Group(gn).Data(dn).lastz)+1; % Number of timepoints
            fdv.NSeg=1; % Segments per timepoint
            fdv.NFlips=(fdv.ntp)*(fdv.NSeg); % Total number of excitations

            %Describe temporal sampling scheme
            %TR=2; % defined above

            NOP=Group(gn).Data(dn).hdr.TD/2;
            TR = Group(gn).Data(dn).hdr.D(2)+NOP/Group(gn).Data(dn).hdr.SW_h+0.03;
            fdv.TR=ones(1,fdv.NFlips)*TR;
            fdv.taxis=cumsum(fdv.TR)-fdv.TR(1);

            %Describe excitation scheme
            fdv.FlipAngle=FA*faff*ones(4,fdv.NFlips);
 
            %Describe "vascular input function"
            %not real, but allows use of 3-compart scripts
            fdv.UseVIF=1;
            fdv.VIFP=zeros(1,fdv.NFlips);
            fdv.VIFL=zeros(1,fdv.NFlips);

            %Placeholder for data to be fit
            fdv.data=[Group(gn).Data(dn).SPyr(Group(gn).Data(dn).lastz:end); ...
                Group(gn).Data(dn).SLac(Group(gn).Data(dn).lastz:end)];
            fdv.Name='PK Analysis - LnCap Prostate';

            %Mask - only fit to data that has stabilized:
            tmp=zeros(2,fdv.NFlips);
            tmp(:,(Group(gn).Data(dn).tstarti-Group(gn).Data(dn).lastz+1):end)=1;
            fdv.Mask=tmp;
            
            %Other miscellaneous
            fdv.verbose=0;


            % Fit to 3-comparment model: vif =0 but nonzero IC for PyrEE, PyrIC, LacIC

            hppkopts=optimset('display','off');
            
            Group(gn).Data(dn).residvec=zeros(1,nkpls);
            Group(gn).Data(dn).fitsvec =zeros(nkpls,length(fdv.fitvars));
            
            %partmp=Group(1).Data(1);
            for kpli=1:nkpls
                fprintf('G%d D%d - %s\n',gn, dn, datestr(now));
                fprintf('kpli = %2d - kpl = %7.5f\n',kpli,kplvec(kpli));
            
                %Initialize 'nuisance parameters' - known or estimated eslewhere.
                % Note vef the extracellular fraction fraction:
                %      ve=(1-vb)*vef
                %      vc=1-vb-ve

                fdv.knowns=  {'Pi0','Li0','kve','vb', 'vef', 'Gam1','Gam2','tdel','Le0','VIFScale','T1Pyr','kpl'};
                fdv.knownvals=[0.0,  0.0,  0.0, 0.0, vef,   2.8,   4.5,   10,    0,    0,           64,     kplvec(kpli)];

                if oldruns>0
                    bestresid=Group(gn).Data(dn).resid;
                    bestfits =Group(gn).Data(dn).fits;
                    residhist=Inf*ones(1,totruns);
                    residhist(1:oldruns)=Group(gn).Data(dn).residhist;
                    fithist  =zeros(totruns,length(fdv.fitvars));
                    fithist(1:oldruns,:)=Group(gn).Data(dn).fithist;
                else
                    bestresid=Inf;
                    bestfits = zeros(1,length(fdv.fitvars));
                    residhist= Inf*ones(1,totruns); % totruns=nruns when oldruns=0
                    fithist =  zeros(totruns,length(fdv.fitvars));
                end

                if exist(dfn,'file')
                    Guess=Group(gn).Data(dn).fits(kpli);
                    [fits,resid] = lsqnonlin(@(x) P2LCv4Err(x,fdv),Guess,LB,UB,hppkopts);
                    bestfits=fits;
                    bestresid=resid;
                end
            
                parallelize=1;
                if parallelize
                    pfresids=zeros(1,totruns-oldruns);
                    pffits  =zeros(totruns-oldruns,length(fdv.fitvars));
                    parfor ii=(oldruns+1):(totruns)
                        Guess=LGB+(UGB-LGB).*rand(1,length(fdv.fitvars));
                        try
                            [fits,resid] = lsqnonlin(@(x) P2LCv4Err(x,fdv),Guess,LB,UB,hppkopts);
                        catch ME
                            fprintf('Warning! Exception Caught Line 310: %s.\n',ME.identifier)
                            resid=Inf;
                            fits = zeros(1,length(fdv.fitvars));
                        end
                        pffits(ii-oldruns,:)=fits;
                        pfresids(ii-oldruns)=resid;
                    end        
                    [bestresid,bri]=min(pfresids);
                    bestfits=pffits(bri,:);
                    fprintf('\n')
                else
                    for ii=(oldruns+1):(totruns)
                        Guess=LGB+(UGB-LGB).*rand(1,length(fdv.fitvars));
                        try
                            [fits,resid] = lsqnonlin(@(x) P2LCv4Err(x,fdv),Guess,LB,UB,hppkopts);
                        catch ME
                            fprintf('Warning! Exception Caught Line 327: %s.\n',ME.identifier)
                            resid=Inf;
                        end
                        residhist(ii)=resid;
                        fithist(ii,:)=fits;
                        if resid<bestresid
                            fprintf('%d.. ',ii)
                            bestresid=resid;
                            bestfits=fits;
                            fh=P2LCv4Plot(bestfits,fdv);shg
                        end
                    end        
%                     fprintf('\n')
                end

                for ii=1:length(fdv.fitvars)
                    fprintf('%s = %7.5e\n',fdv.fitvars{ii},bestfits(ii))
                end
                fprintf('bestresid = %7.5f\n\n',bestresid)
                Group(gn).Data(dn).fitsvec(kpli,:) = bestfits;
                Group(gn).Data(dn).residvec(kpli)= bestresid;
            end % kpli loop
            Group(gn).Data(dn).fdv  = fdv;
        end
    end
    save(dfn,'Group','totruns','kplvec')
end

fprintf('\nCalculations completed at %s\n',datestr(now));

