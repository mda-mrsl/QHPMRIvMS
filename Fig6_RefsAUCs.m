%% Create anatomic reference images and AUC images
%Note overall figure will by 7.3in wide, 300dpi
%Baseline first:
clear
resave=1;

Name='HPMRI_ATC_Pre';
load(fullfile('data',[Name '.mat']));

%Select representative slice:
sno=3; 

%Position of figure for proper formatting:
figpos=[1940 210 900 750];

AUCp=squeeze(sum(spyr(:,:,sno,pyri(1):pyri(2)),4));
AUCl=squeeze(sum(slac(:,:,sno,laci(1):laci(2)),4));

%Set area to 'zoom' in on - pixel indices for background T2w images
zrx=[154-5 384+25];
zry=[113-25 315+10];
tmp=zeros(size(bgims,1),size(bgims,2));
tmp(zry(1):zry(2),zrx(1):zrx(2))=1;

f11h=figure(11);
set(gcf,'Name','Fig 6A: FOV Inset Ref 1 ')
imshow(bgims(:,:,bgi(sno))/max(max(squeeze(bgims(:,:,bgi(sno))))))
ylabel('Baseline','Fontsize',20,'Fontweight','Bold');
hold on
contour(troi(:,:,sno),[0.9 0.9],'r-','linewidth',2);
contour(tmp,[0.9 0.9],'y-','linewidth',2)
f11h.Position=figpos;
hold off
ofn='Fig6a_Pre_Ref.tif';
if (~exist(ofn,'file')||resave)
    exportgraphics(gca,ofn);
end

f12h=figure(12);
set(gcf,'Name','Fig 6B: Pyr AUC 1')
registeredoverlay2(sprintf('Pyr AUC, Slice %d',sno),AUCp,hph,4,...
    bgims(:,:,bgi(sno)),bgh,zrx,zry);
%truesize(gcf,[450 450])
f12h.Position=figpos;
ofn='Fig6b_Pre_PyrAUC.tif';
if (~exist(ofn,'file')||resave)
    exportgraphics(gca,ofn);
end

f13h=figure(13);
set(gcf,'Name','Fig 6C: Lac AUC 1')
registeredoverlay2(sprintf('Lac AUC, Slice %d',sno),AUCl,hph,4,...
    bgims(:,:,bgi(sno)),bgh,zrx,zry);
%truesize(gcf,[450 450])
f13h.Position=figpos;
ofn='Fig6c_Pre_LacAUC.tif';
if (~exist(ofn,'file')||resave)
    exportgraphics(gca,ofn);
end

%% Second timepoint:

Name='HPMRI_ATC_Post';
load(fullfile('data',[Name '.mat']));

AUCp=squeeze(sum(spyr(:,:,sno,pyri(1):pyri(2)),4));
AUCl=squeeze(sum(slac(:,:,sno,laci(1):laci(2)),4));

%Set area to zoom for background T2w images
zrx=[154-5 384+25]-25;
zry=[113-25 315+10];
tmp=zeros(size(bgims,1),size(bgims,2));
tmp(zry(1):zry(2),zrx(1):zrx(2))=1;

f21h=figure(21);
set(gcf,'Name','Fig 6C: FOV Inset Ref +8d')
imshow(bgims(:,:,bgi(sno))/max(max(squeeze(bgims(:,:,bgi(sno))))))
ylabel('+8 Days','Fontsize',20,'Fontweight','Bold');
hold on
contour(troi(:,:,sno),[0.9 0.9],'r-','linewidth',2);
contour(tmp,[0.9 0.9],'y-','linewidth',2)
hold off
f21h.Position=figpos;
ofn='Fig6d_Post_Ref.tif';
if (~exist(ofn,'file')||resave)
    exportgraphics(gca,ofn);
end

f22h=figure(22);
set(gcf,'Name','Fig 6D: Pyr AUC +8d')
registeredoverlay2(sprintf('Pyr AUC, Slice %d',sno),AUCp,hph,4,...
    bgims(:,:,bgi(sno)),bgh,zrx,zry);
%truesize(gcf,[450 450])
f22h.Position=figpos;
ofn='Fig6e_Post_PyrAUC.tif';
if (~exist(ofn,'file')||resave)
    exportgraphics(gca,ofn);
end

f23h=figure(23);
set(gcf,'Name','Fig 6E: Lac AUC +8d')
registeredoverlay2(sprintf('Lac AUC, Slice %d',sno),AUCl,hph,4,...
    bgims(:,:,bgi(sno)),bgh,zrx,zry);
%truesize(gcf,[450 450])
f23h.Position=figpos;
ofn='Fig6f_Post_LacAUC.tif';
if (~exist(ofn,'file')||resave)
    exportgraphics(gca,ofn);
end
