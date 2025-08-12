% This will be one panel for Fig 4, which will overall be 7.3" wide 
% (2 cols, 300dpi)
clear

load(fullfile('data','HPNMRFits_Hth83.mat')); 
NGrp=length(Group);
NSamp=length(Group(1).Data);
NParm=size(Group(1).Data(1).fitsvec,2);

for gn=1%:NGrp
    for dn= 4%1:NSamp
        fdv=Group(gn).Data(dn).fdv;
        fdv.Name=sprintf('%s Sample %d',Group(1).Name,dn);
        fitsvec=Group(gn).Data(dn).fitsvec;
        residvec=Group(gn).Data(dn).residvec;
        
        %Find the best fit @ minimum residual:
        [minresid minresi]=min(residvec);

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
        %And find the rest of the parameters too:
        parmsest=zeros(1,NParm);
        for qq=1:NParm
            slope=(residvec(tmpi)-residvec(tmpi-1))/(fitsvec(tmpi,qq)-fitsvec(tmpi-1,qq));
            parmsest(qq)=fitsvec(tmpi-1,qq) + (thresh-residvec(tmpi-1))/slope;
        end
        
        %Set up plots:
        f4ah=figure(410+dn); % Fig 4b
        set(gcf,'Name',fdv.Name)
        
        % get curves for best fit:
        fdv.knownvals(find(strcmp(fdv.knowns,'kpl')))=kplvec(minresi);
        [EV,IV,vols,Mzev,Mziv] = P2LCv4(fitsvec(minresi,:),fdv);
        totmin=(IV(1:2,:)*vols(1))+(EV(1:2,:)*vols(2))+(EV(3:4,:)*vols(3));
        
        % get curves at 1% threshold:
        fdv.knownvals(find(strcmp(fdv.knowns,'kpl')))=kplest;
        [EV,IV,vols,Mzev,Mziv] = P2LCv4(parmsest,fdv);
        tot=(IV(1:2,:)*vols(1))+(EV(1:2,:)*vols(2))+(EV(3:4,:)*vols(3));
        %calculate residual at estimated parameter set:
        residatest=P2LCv4Err(parmsest,fdv);
        resnatest=sum(residatest.^2);
        
        plotveci=zeros(1,length(fdv.taxis));
        plotveci(1:5)=1;
        plotveci(5:3:101)=1;
        %Plot!
        subplot(2,1,1)
        f4a1h=plot(fdv.taxis(plotveci==1),fdv.data(1,plotveci==1),'gx',...
            fdv.taxis(plotveci==1),fdv.data(2,plotveci==1),'bx',...
            fdv.taxis(1:101),tot(1,1:101),'g-',...
            fdv.taxis(1:101),tot(2,1:101),'b-',...
            fdv.taxis(1:101),totmin(1,1:101),'g:',...
            fdv.taxis(1:101),totmin(2,1:101),'b:',...
            Group(gn).Data(dn).taxis(Group(gn).Data(dn).tstarti+1)-Group(gn).Data(dn).taxis(Group(gn).Data(dn).lastz+1),...
                Group(gn).Data(dn).SPyr(Group(gn).Data(dn).tstarti),'ro',...
            'linewidth',2,'MarkerSize',8);
        f4a1h(5).LineWidth=1;
        f4a1h(6).LineWidth=1;
        grid on
        axis tight
        title(Group(1).Name(1:5),'Fontsize',22,'FontWeight','bold');
        ylabel('Signal Amplitude (arb)','FontSize',18,'FontWeight','bold'); % only for Fig4a
        xlabel('Time (s)','Fontsize',18,'FontWeight','bold');
        subplot(2,1,2)
        f4a2h=plot(kplvec,residvec,'-',kplest,resnatest,'rx','linewidth',2,'MarkerSize',12);
        ylabel('Residual (arb)','FontSize',18,'FontWeight','bold'); % Only Fig4a
        xlabel('k_{PL} (s^{-1})','FontSize',18,'FontWeight','bold')
        text(kplest,resnatest+max(residvec)/9,...
            sprintf('k_{PL} = %5.3f s^{-1}',kplest));
        grid on
        axis tight
        
        %Save to files
        set(gcf,'position',[0 60 600 600]);
        print(sprintf('Fig4-NMR-%s',Group.Name),'-depsc')
        exportgraphics(gcf,sprintf('Fig4-NMR-%s.png',Group(1).Name(1:5)))
        
        %Copy to clipboard:
        copygraphics(gcf)
    end
end

