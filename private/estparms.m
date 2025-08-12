function [kplest parmsest]=estparms(residvec,kplvec,fitsvec,fdv,showme)

if ~exist('showme','var')
    showme=0;
end

maxres=max(residvec);
[minres minresi]=min(residvec);


thresh=minres+((maxres-minres)/100);

%Find point where resid(kpl) is changing by less than 1%:
tmpi=1;
thresh=(max(residvec)-min(residvec))/100+min(residvec);
while residvec(tmpi)>thresh
    tmpi=tmpi+1;
end

%Assume residual is linear over this interval. Estimate kpl at the
%crossing point:
slope=(residvec(tmpi)-residvec(tmpi-1))/(kplvec(tmpi)-kplvec(tmpi-1));
kplest=kplvec(tmpi-1) + (thresh-residvec(tmpi-1))/slope;

if showme
    figure(720)
    plot(kplvec,residvec,'-',kplest,thresh,'ro');
end

%And find the rest of the parameters too:
nparms=size(fitsvec,2);
parmsest=zeros(1,nparms);
for qq=1:nparms
    slope=(residvec(tmpi)-residvec(tmpi-1))/(fitsvec(tmpi,qq)-fitsvec(tmpi-1,qq));
    parmsest(qq)=fitsvec(tmpi-1,qq) + (thresh-residvec(tmpi-1))/slope;
    if showme
        figure(720+qq)
        plot(fitsvec(:,qq),residvec,'-',parmsest(qq),thresh,'ro');
        ylabel('residual');xlabel(fdv.fitvars{qq});
    end
end

if showme
    dresdkpl=(residvec(3:end)-residvec(1:end-2))./(kplvec(3:end)-kplvec(1:end-2));
    figure(799)
    plot(kplvec(2:end-1),dresdkpl);
    ylabel('dresn/dkpl');
    xlabel('kpl')
end
    
