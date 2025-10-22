function [telapsed] = MSPK4v2(celsel,MassSpecData)

saveout=1;
verbose=0;

tstart=tic;

nruns=5000; % per kpl value
nexc=0;

% Sampling times:
taxis=MassSpecData.Hth83_20200305.Parental.texp; % [0 5 15 30 45 60]

%Select dataset to analyze: 1=Hth83; 2=PC3; 3=LnCap
%celsel=1;
switch celsel
    case 1 % Hth83
        name='PK Analysis of Mass Spec - Hth83 ATC'
        ofn='MSPKSweepHth83';

        % Mass spec values from cells:
        Pm0 = (MassSpecData.Hth83_20200305.Parental.Pm0 + MassSpecData.Hth83_20200305.Run2.Pm0 + MassSpecData.Hth83_20200305.Run3_20190207.Pm0)/3;
        Pm1 = (MassSpecData.Hth83_20200305.Parental.Pm1 + MassSpecData.Hth83_20200305.Run2.Pm1 + MassSpecData.Hth83_20200305.Run3_20190207.Pm1)/3;
        Pm2 = (MassSpecData.Hth83_20200305.Parental.Pm2 + MassSpecData.Hth83_20200305.Run2.Pm2 + MassSpecData.Hth83_20200305.Run3_20190207.Pm2)/3;
        Pm3 = (MassSpecData.Hth83_20200305.Parental.Pm3 + MassSpecData.Hth83_20200305.Run2.Pm3 + MassSpecData.Hth83_20200305.Run3_20190207.Pm3)/3;
        PmT = Pm0+Pm1+Pm2+Pm3;
        
        Lm0 = (MassSpecData.Hth83_20200305.Parental.Lm0 + MassSpecData.Hth83_20200305.Run2.Lm0 + MassSpecData.Hth83_20200305.Run3_20190207.Lm0)/3;
        Lm1 = (MassSpecData.Hth83_20200305.Parental.Lm1 + MassSpecData.Hth83_20200305.Run2.Lm1 + MassSpecData.Hth83_20200305.Run3_20190207.Lm1)/3;
        Lm2 = (MassSpecData.Hth83_20200305.Parental.Lm2 + MassSpecData.Hth83_20200305.Run2.Lm2 + MassSpecData.Hth83_20200305.Run3_20190207.Lm2)/3;
        Lm3 = (MassSpecData.Hth83_20200305.Parental.Lm3 + MassSpecData.Hth83_20200305.Run2.Lm3 + MassSpecData.Hth83_20200305.Run3_20190207.Lm3)/3;
        LmT = Lm0+Lm1+Lm2+Lm3;
        
        % Mass spec values from media @ t=60s
        % Note per 6/10/2020 email from Lin Tan: no lactate in fresh media at t=0. No unlabeled pyruvate.
        MLm0 = [0 0 0 0 0 458622132];
        MLm1 = [0 0 0 0 0 12471136];
        MLm2 = [0 0 0 0 0 151794];
        MLm3 = [0 0 0 0 0 20377337];
        MLmT = MLm0+MLm1+MLm2+MLm3;
        
        MPm0 = [0 0 0 0 0 18824185];
        MPm1 = [0 0 0 0 0 573813];
        MPm2 = [0 0 0 0 0 37574456];
        MPm3 = [0 0 0 0 0 1665402997];
        MPmT = MPm0+MPm1+MPm2+MPm3;
        
        %Volume of cells, assuming 25 million cells/dish with 22um diameter cells
        vcells=(25e6)*(4/3*pi*((11e-4)^3)); % in mL
        %per yun analysis/email (11/6/2020) there were approx 1.32e6
        %cells/plate:
    
    case 2 %PC3
        name='PK Analysis of Mass Spec - PC3 Prostate'
        ofn='MSPKSweepPC3';

        % Mass spec values from cells:
        Pm0 = MassSpecData.PC3_20200305.Average.Pm0;
        Pm1 = MassSpecData.PC3_20200305.Average.Pm1;
        Pm2 = MassSpecData.PC3_20200305.Average.Pm2;
        Pm3 = MassSpecData.PC3_20200305.Average.Pm3;
        PmT = Pm0+Pm1+Pm2+Pm3;
        
        Lm0 = MassSpecData.PC3_20200305.Average.Lm0;
        Lm1 = MassSpecData.PC3_20200305.Average.Lm1;
        Lm2 = MassSpecData.PC3_20200305.Average.Lm2;
        Lm3 = MassSpecData.PC3_20200305.Average.Lm3;
        LmT = Lm0+Lm1+Lm2+Lm3;
        
        % Mass spec values from media @ t=60s
        % Note per 6/10/2020 email from Lin Tan: no lactate in fresh media at t=0. No unlabeled pyruvate.
        MLm0 = [0 0 0 0 0 338983773];
        MLm1 = [0 0 0 0 0 8438420];
        MLm2 = [0 0 0 0 0 10392];
        MLm3 = [0 0 0 0 0 4037156];
        MLmT = MLm0+MLm1+MLm2+MLm3;
        
        MPm0 = [0 0 0 0 0 7584071];
        MPm1 = [0 0 0 0 0 297256];
        MPm2 = [0 0 0 0 0 23463675];
        MPm3 = [0 0 0 0 0 1159809726];
        MPmT = MPm0+MPm1+MPm2+MPm3;
        
        %Volume of PC3 cells estimated as 2785.9 fL/cell 
        % from: Oprea-Lager DE, van Kanten MP, van Moorselaar RJA, van den
        %   Eertwegh AJM, van de Ven PM, Bijnsdorp IV, Hoekstra OS, Geldof AA.
        %   [18F]Fluoromethylcholine as a chemotherapy response read-out in
        %   prostate cancer cellse.  Mol Imaging Biol 17(3):319-27, 2015.
        %Number of cells per plate estimated as 1.75 million (per yun email)
        
        vcells=(1.75e6)*(2785.9e-15); % in L
        vcells=vcells*1000; % in mL

    case 3 %LnCap
        name='PK Analysis of Mass Spec - LnCap Prostate'
        ofn='MSPKSweepLnCap';

        % Mass spec values from cells:
        Pm0 = MassSpecData.LnCap_20200305.Average.Pm0;
        Pm1 = MassSpecData.LnCap_20200305.Average.Pm1;
        Pm2 = MassSpecData.LnCap_20200305.Average.Pm2;
        Pm3 = MassSpecData.LnCap_20200305.Average.Pm3;
        PmT = Pm0+Pm1+Pm2+Pm3;
        
        Lm0 = MassSpecData.LnCap_20200305.Average.Lm0;
        Lm1 = MassSpecData.LnCap_20200305.Average.Lm1;
        Lm2 = MassSpecData.LnCap_20200305.Average.Lm2;
        Lm3 = MassSpecData.LnCap_20200305.Average.Lm3;
        LmT = Lm0+Lm1+Lm2+Lm3;
        
        % Mass spec values from media @ t=60s
        % Note per 6/10/2020 email from Lin Tan: no lactate in fresh media at t=0. No unlabeled pyruvate.
        MLm0 = [0 0 0 0 0 355920065];
        MLm1 = [0 0 0 0 0 9684838];
        MLm2 = [0 0 0 0 0 12071];
        MLm3 = [0 0 0 0 0 5875454];
        MLmT = MLm0+MLm1+MLm2+MLm3;
        
        MPm0 = [0 0 0 0 0 11174177];
        MPm1 = [0 0 0 0 0 299719];
        MPm2 = [0 0 0 0 0 14935140];
        MPm3 = [0 0 0 0 0 577019730];
        MPmT = MPm0+MPm1+MPm2+MPm3;
        
        %Volume of LnCap cells estimated as 3787.6 fL/cell 
        % from: Oprea-Lager DE, van Kanten MP, van Moorselaar RJA, van den
        %   Eertwegh AJM, van de Ven PM, Bijnsdorp IV, Hoekstra OS, Geldof AA.
        %   [18F]Fluoromethylcholine as a chemotherapy response read-out in
        %   prostate cancer cellse.  Mol Imaging Biol 17(3):319-27, 2015.
        %Number of cells per plate estimated as 2.22 million (per Yun email)
        
        vcells=(2.22e6)*(3787.6e-15); % in L
        vcells=vcells*1000; % in mL
end

mmask=[0 0 0 0 0 0];
cmask=[1 1 1 1 1 1];

vmedia=2; % mL

%% Set forward variables
eps=1e-100;
vef=vmedia/(vmedia+vcells);

fdv.fitvars={'klp'  'kecp' 'kecl' 'Pe0'};
LB=         [ 0      eps    eps    0];
UB=         [ Inf    1.0    1.0    Inf];
UG=         [ 0.2    0.1    0.1    5*MPm3(6)];

fdv.knowns=  {'kve','vb', 'vef', 'Gam1','Gam2','tdel','Le0','Pi0','Li0','kpx','klx','VIFScale','kpl'};
fdv.knownvals=[0.00, eps,  vef,   2.8,   4.5,   10,    0,    0,    0,    0,    0,    1,         0];

fdv.Name=name;
fdv.mask=[mmask;mmask;cmask;cmask];
fdv.ntp=length(Pm3);
fdv.NSeg=1;
fdv.NFlips=(fdv.ntp)*(fdv.NSeg); 
fdv.FlipAngle=zeros(4,fdv.NFlips);
fdv.TR=taxis(2:end)-taxis(1:end-1);
fdv.taxis=taxis;
fdv.UseVIF=1;
fdv.VIFP=zeros(1,fdv.NFlips);
fdv.VIFL=zeros(1,fdv.NFlips);
fdv.data=[MPm3; MLm3; Pm3; Lm3];
fdv.verbose=0;

%% Run the fits

hppkopts=optimset('display','off');

kplvec=logspace(log10(0.01),log10(1.5),30);
residvec=zeros(1,length(kplvec));
fitsvec =zeros(length(kplvec),length(fdv.fitvars));

for jj=1:length(kplvec)
    if verbose
        fprintf('jj = %d/%d:\n',jj,length(kplvec));
    end
    fdv.knownvals(end)=kplvec(jj);

    bestresid=Inf;
    bestfits = zeros(1,length(fdv.fitvars));

    for ii=1:nruns
        Guess=LB+(UG-LB).*rand(1,length(fdv.fitvars));
        try
            [fits,resid] = lsqnonlin(@(x) P2LCMSv4Err(x,fdv),Guess,LB,UB,hppkopts);
        catch ME
            fprintf('Warning! %s: %s\n',ofn,ME.identifier)
            resid=Inf;
            nexc=nexc+1;
        end
        if resid<bestresid
            if verbose
                fprintf('%04d: kpl = %6.4f; Resid = %6.4e\n',ii,fits(1),resid)
                fh=P2LCMSv4Plot(bestfits,fdv);shg
            end
            bestresid=resid;
            bestfits=fits;
        end
    end        
    residvec(jj)=bestresid;
    fitsvec(jj,:)=bestfits;
end

telapsed=toc(tstart);

if saveout
    save(fullfile('data',ofn),'fdv','bestfits','kplvec','residvec','fitsvec','LB','UB')
end