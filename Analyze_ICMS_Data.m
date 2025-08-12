
clear
load('data/MassSpecData.mat')

MediaFlag=0; % Do not include ICMS data from media

for CellLine=1:3
    MSPKSweep(CellLine,MassSpecData);
end











