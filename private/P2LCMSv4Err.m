function [resid] = P2LCMSv4Err(parms,fdv)

[EVt,IVt,vols,EVz,IVz] = P2LCMSv4(parms,fdv);
resid=(fdv.data-EVz).*fdv.mask;
%Normalize to max of measured data. May need to weight...later
%resid(1,:)=resid(1,:)/max(fdv.data(1,:));
%resid(2,:)=resid(2,:)/max(fdv.data(2,:));

resid=reshape(resid,1,prod(size(fdv.data))); 