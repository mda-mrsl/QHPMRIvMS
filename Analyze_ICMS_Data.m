
clear
load('data/MassSpecData.mat')

parfor CellLine=1:3
    MSPKSweep(CellLine,MassSpecData);
end












