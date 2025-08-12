function [resid] = P2LCv4Err(parms,fdv)

[EV,IV,vols] = P2LCv4(parms,fdv);
tot=(IV(1:2,:)*vols(1))+(EV(1:2,:)*vols(2))+(EV(3:4,:)*vols(3));
if isfield(fdv,'Mask')
    resid=(fdv.data-tot).*fdv.Mask;
else
    resid=(fdv.data-tot);
end
%Normalize to max of measured data. May need to weight...later
resid(1,:)=resid(1,:)/max(fdv.data(1,:));
resid(2,:)=resid(2,:)/max(fdv.data(2,:));

resid=reshape(resid,1,prod(size(fdv.data))); 