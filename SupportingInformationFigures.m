%% Fig S1; plot sea level over time
figure
for RSLRA=[0.1592,0.1106,0.08023,0.03536,0.01881]
    load(['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Sill0_SSCReduction0_silldelay551_TLP0_TLPthick0_TLPfreq0_TLPdelay551_Erosion5.3e-09/meansealevel.mat']);
    plot(1:701,msl(1:701),'LineWidth',1.8)
    %ylim([-10 10])
    ylim([-1,5])
    xlim([0 825])
    leftcolors = [0 0 1; 0.1 0.20 0.9; 0.2 0.40 0.8; 0.3 0.60 0.7; 0.4 0.80 0.6];
    colororder(leftcolors)
    hold on
end
xline(300,'LineStyle',':');
xline(400,'LineStyle',':');
xline(550,'LineStyle',':');
hold on

drawbrace([0,4.3],[550,4.3],10,'Color','k')
drawbrace([551,4.3],[701,4.3],10,'Color','k')

text(80,4,'SLR = 1 mm yr^-^1','FontSize',8)
start=[0 3.5];
stop=[295 3.5];
arrow(start,stop,'width',0.3,'ends',[1 2],'baseangle',10,'tipangle',12,'length',10)
hold on
text(350,4,{'Ramp to','3 mm yr^-^1'},'FontSize',8,'HorizontalAlignment', 'center')
start=[305 3.5];
stop=[395 3.5];
arrow(start,stop,'width',0.3,'ends',[1 2],'baseangle',10,'tipangle',12,'length',10)
text(475,4,{'SLR =','3 mm yr^-^1'},'FontSize',8,'HorizontalAlignment', 'center')
start=[405 3.5];
stop=[550 3.5];
arrow(start,stop,'width',0.3,'ends',[1 2],'baseangle',10,'tipangle',12,'length',10)
text(630,4,{'Accelerating','SLR'},'FontSize',8,'HorizontalAlignment', 'center')
text(250,4.77,'Spinup','FontSize',8)
text(580,4.77,'Experiments','FontSize',8)
text(710,4,'High','FontSize',8,'Color',[0 0 1])
text(710,2.94,'Inter.-high','FontSize',8,'Color',[0.1 0.20 0.9])
text(710,2.27,'Intermediate','FontSize',8,'Color',[0.2 0.40 0.8])
text(710,1.3,'Inter.-low','FontSize',8,'Color',[0.3 0.60 0.7])
text(710,0.88,'Low','FontSize',8,'Color',[0.4 0.80 0.6])

xlabel('Time (yr)')
%ylabel('Rate of sea-level rise (mm yr^-^1)')
ylabel('Mean sea level (m)')

% l=legend('high','intermediate high','intermediate','intermediate low','low');
% l.Location='southwest';
% l.Color=[0.8 0.8 0.8];
% l.EdgeColor='none';
% l.BackgroundAlpha=0.3;

%% Fig S2: erosion rates through time for each erosion scenario 

RSLRA=0.08023;
color=[0.99608 0.78824 0.89412]; %start with color for low erosion scenario and will increase for every iteration of for-loop
for k=[1.65e-09,2.82e-09,5.3e-09,7.8e-09,9.95e-09]
    filename=['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Sill0_SSCReduction0_silldelay551_TLP0_TLPthick0_TLPfreq0_TLPdelay551_Erosion' num2str(k) '/'];
    myVars={'elev.mat','marsh edge.mat','meansealevel.mat','AnnualAccretion.mat','SLR.mat','resiliencymetrics.mat'};
    for i=1:6
        load([filename myVars{i}]);
    end
    erosion=zeros(701,1);
    accretiontoSLR=zeros(701,1);
    AvgAnnualAccretion=zeros(701,1);
    for yr=2:701
        erosion(yr)=Marsh_edge(yr)-Marsh_edge(yr-1);
        AvgAnnualAccretion(yr)=mean(AnnualAccretion(yr,Marsh_edge(yr):5501));
        accretiontoSLR(yr)=AvgAnnualAccretion(yr)/SLR(yr);
    end
    for yr=551:701
        if erosion(yr)>10
            erosion(yr)=NaN;
        end
    end
    marshage=find(resiliencymetrics(551:701,11)==0,1);
    if isempty(marshage)==1
        marshage=150;
    end
    year=551:marshage+550;
    if k==1.65e-09
        lowerosion=erosion(551:marshage+550);
    else if k==2.82e-09
        intermlowerosion=erosion(551:marshage+550);
    else if k==5.3e-09
        intermerosion=erosion(551:marshage+550);
    else if k==7.8e-09
        intermhigherosion=erosion(551:marshage+550);
    else if k==9.95e-09
        higherosion=erosion(551:marshage+550);
        color=[0.44706 0 0.22353];
    end
    end
    end
    end
    end
    figure
    x=(1:marshage-2)';
    scatter(x,erosion(552:marshage+549),10,'filled','MarkerFaceColor',color,'MarkerEdgeColor','none','MarkerFaceAlpha',0.8)
    hold on
    ylim([0 10])
    y=erosion(552:marshage+549);

    [xData, yData] = prepareCurveData(x,y);

    % Set up fittype and options.
    ft = fittype( 'smoothingspline' );
    opts = fitoptions( 'Method', 'SmoothingSpline' );
    opts.SmoothingParam = 0.01;

    % Fit model to data.
    [fitresult,q]= fit( xData, yData, ft, opts );

    % Plot fit with data.
    p=plot(fitresult,xData,yData);
    set(p,'LineWidth',3.5,'Color',color);
    hold on
    color=color+[-0.14 -0.2 -0.15];
end
ylabel('Erosion Rate (m yr^-^1)')
xlabel('Model Run Year')
ylim([0 8.5])
xlim([0 170])
text(115,1.4,'Erosion_{initial} = 0.5 m yr^-^1','FontSize',10,'Color',[0.99608 0.78824 0.89412])
text(112,2.2,'Erosion_{initial} = 1 m yr^-^1','FontSize',10,'Color',[0.85608 0.58824 0.74412])
text(107,4.2,'Erosion_{initial} = 2 m yr^-^1','FontSize',10,'Color',[0.71608 0.38824 0.69312])
text(100,5.8,'Erosion_{initial} = 3 m yr^-^1','FontSize',10,'Color',[0.67508 0.18824 0.54312])
text(93,7.2,'Erosion_{initial} = 4 m yr^-^1','FontSize',10,'Color',[0.44706 0 0.22353])
box on

%% Fig S3: TLP effect on elevation relative to MSL 
CO=10;RSLRA=0.08023;RSLR=3;Erosion=5.3e-09;
filename = ['Outputs_CO' num2str(CO) '_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR' num2str(RSLR) '_Sill0_SSCReduction0_silldelay551_TLP1_TLPthick0.05_TLPfreq15_TLPdelay551_Erosion' num2str(Erosion) '\'];
    myVars={'elev.mat','marsh edge.mat','meansealevel.mat','AnnualAccretion.mat','SLR.mat','resiliencymetrics.mat','VASEAmatrix.mat'};
    for i=1:7
        load([filename myVars{i}]);
    end
figure
p(1)=plot(0:150,msl(551:701),'Color','#4393c3','LineWidth',5);
hold on
for yr=551:701
    avgelev(yr)=mean(elev(yr,Marsh_edge(yr):5500));
end
p(2)=plot(0:150,avgelev(551:701),'Color','#91cf60','LineWidth',5);
hold on
% tidalrange(:,1)=msl(551:701)+0.5;
% tidalrange(:,2)=msl(551:701)-0.5;
p(3)=plot(0:150,msl(551:701)+0.5,'Color','k','LineWidth',2,'LineStyle','--');
hold on
p(4)=plot(0:150,msl(551:701)-0.5,'Color','k','LineWidth',2,'LineStyle','--');
hold on
scatter(0:15:136,2.96,80,'^','filled','MarkerFaceColor',[0.5,0.5,0.5],'MarkerEdgeColor','none');
s(1)=scatter(0,2.96,80,'^','filled','MarkerFaceColor',[0.5,0.5,0.5],'MarkerEdgeColor','none');
xlabel('Model Run Year')
ylabel('Elevation (m)')
legend([p(1) p(2) p(4) s(1)],'Mean sea level','Average marsh elevation','Tidal range','TLP Year')

%% Fig S4: sensitivity analysis of shoreline stabilization SSC reduction
figure
%graph VEC for each SSC reduction value for SS
subplot(1,2,1)
s=subaxis(1,2,1,'SpacingHoriz',0.04,'SpacingVert',0.04); %limits the padding between subplots
color=[0.553 0.859 0.212]; %=[0.553 0.859 0.212] start with color for low erosion scenario and will increase for every iteration of for-loop
CO=10; RSLRA=0.08023;RSLR=3;Erosion=5.3e-09;
for SSCReduction=[0,10,20,30,40]
    filename = ['Outputs_CO' num2str(CO) '_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR' num2str(RSLR) '_Sill1_SSCReduction' num2str(SSCReduction) '_silldelay551_TLP0_TLPthick0_TLPfreq0_TLPdelay551_Erosion' num2str(Erosion) '\'];
    myVars={'elev.mat','marsh edge.mat','meansealevel.mat','AnnualAccretion.mat','SLR.mat','resiliencymetrics.mat','VASEAmatrix.mat'};
    for i=1:7
        load([filename myVars{i}]);
    end
    marshage=find(resiliencymetrics(551:701,11)==0,1);
    if isempty(marshage)==1
        marshage=150;
    end
    Dmax = (0.237*1)-0.092+0.5;
    yrDmax=zeros(701,1);
    yrDmax=msl+0.5-Dmax;
    marshvolumeaboveDmax=zeros(701,1);
    for yr=1:701
        pond=find(elev(yr, Marsh_edge(yr)+1:5500) < yrDmax(yr),1)-1;
    if isempty(pond)==1
        pond=5500-Marsh_edge(yr)-1;
    end
    diff=elev(yr,Marsh_edge(yr):pond+Marsh_edge(yr))-yrDmax(yr);
    marshvolumeaboveDmax(yr)=sum(diff);
    clear diff
    end
    p1=plot(0:150,marshvolumeaboveDmax(551:701),'LineWidth',5,'Color',color);%#66a61e
    hold on
    color=color+[-0.1 -0.2 -0.05];
end
ax=gca;
ax.FontSize=16;
xlabel('Model Run Year')
ylabel('Volumetric elevation capital (m^3)')
text(3,285,'a  Breakwater','FontSize',16,'FontWeight','bold')
ylim([-10 300])
l=legend('SSC Reduction = 0%','SSC Reduction = 10%','SSC Reduction = 20%','SSC Reduction = 30%','SSC Reduction = 40%');
l.EdgeColor='none'; %grey legend outline
l.Color=[0.94 0.94 0.94 0.94];
%--------------------------------------------------------------------------------------------------------------------
%graph VEC for each SSC reduction value for SS + TLP
subplot(1,2,2)
s=subaxis(1,2,2,'SpacingHoriz',0.04,'SpacingVert',0.04); %limits the padding between subplots
color=[0.647 0.635 0.808]; %start with color for low erosion scenario and will increase for every iteration of for-loop
CO=10; RSLRA=0.08023;RSLR=3;Erosion=5.3e-09;
for SSCReduction=[0,10,20,30,40]
    filename = ['Outputs_CO' num2str(CO) '_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR' num2str(RSLR) '_Sill1_SSCReduction' num2str(SSCReduction) '_silldelay551_TLP1_TLPthick0.05_TLPfreq15_TLPdelay551_Erosion' num2str(Erosion) '\'];
    myVars={'elev.mat','marsh edge.mat','meansealevel.mat','AnnualAccretion.mat','SLR.mat','resiliencymetrics.mat','VASEAmatrix.mat'};
    for i=1:7
        load([filename myVars{i}]);
    end
    marshage=find(resiliencymetrics(551:701,11)==0,1);
    if isempty(marshage)==1
        marshage=150;
    end
    Dmax = (0.237*1)-0.092+0.5;
    yrDmax=zeros(701,1);
    yrDmax=msl+0.5-Dmax;
    marshvolumeaboveDmax=zeros(701,1);
    for yr=1:701
        pond=find(elev(yr, Marsh_edge(yr)+1:5500) < yrDmax(yr),1)-1;
    if isempty(pond)==1
        pond=5500-Marsh_edge(yr)-1;
    end
    diff=elev(yr,Marsh_edge(yr):pond+Marsh_edge(yr))-yrDmax(yr);
    marshvolumeaboveDmax(yr)=sum(diff);
    clear diff
    end
    p1=plot(0:150,marshvolumeaboveDmax(551:701),'LineWidth',5,'Color',color);%#66a61e
    hold on
    %[0.647 0.635 0.808]
    color=color+[-0.14 -0.15 -0.15];
end
s1=scatter(0:15:150,297,80,'^','filled','MarkerFaceColor',[0.5,0.5,0.5],'MarkerEdgeColor','none');
ax=gca;
ax.FontSize=16;
xlabel('Model Run Year')
text(3,285,'b  TLP + Breakwater','FontSize',16,'FontWeight','bold')
ylim([-10 300])
l=legend('SSC Reduction = 0%','SSC Reduction = 10%','SSC Reduction = 20%','SSC Reduction = 30%','SSC Reduction = 40%','TLP Year');
l.EdgeColor='none'; %grey legend outline
l.Color=[0.94 0.94 0.94 0.94];
%% Fig S5: sensitivity analysis of SSC on vertical resilience
% fig a: sensitivity graph with SSC on x-axis, time to interior drown on y-axis, with different restoration strategies as 4 diff lines
figure
subplot(1,2,1)
% set up for-loop to pull matrices of various SSC files for intermediate RSLRA for no restoration
RSLRA=0.1592;
k=5.3e-09;
%sill
color='#66a61e';
for CO=[10,20,30,40,50]
    filename=['Outputs_CO' num2str(CO) '_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Sill1_SSCReduction0_silldelay551_TLP0_TLPthick0_TLPfreq0_TLPdelay551_Erosion' num2str(k) '/'];
    myVars={'elev.mat','marsh edge.mat','meansealevel.mat','AnnualAccretion.mat','SLR.mat','resiliencymetrics.mat'};
    for i=1:6
        load([filename myVars{i}]);
    end

    % define onset of interior drowning
    marshage=find(resiliencymetrics(551:701,11)==0,1)+551;
    if isempty(marshage)==1
        marshage=701;
    end
    for yr=552:701
        widthslope(yr)=(resiliencymetrics(yr,7)-resiliencymetrics(yr-1,7))/1; %change in width each year
    end
    threshold=find(widthslope(552:701)<widthslope(551)-5,1)-1+551;

    if isempty(threshold)==1
        threshold=marshage;
        color='k';
    end

    thresholdsmatrix(1,CO/10)=threshold;
    scatter(CO,threshold-551,110,'MarkerFaceColor',color,'MarkerEdgeColor','none')
    hold on
    color='#66a61e';
end
p1=plot(10:10:50,thresholdsmatrix(1:5)-551,'LineWidth',3,'Color','#66a61e');
hold on

%TLP
color='#7570b3';
for CO=[10,20,30,40,50]
    filename=['Outputs_CO' num2str(CO) '_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Sill0_SSCReduction0_silldelay551_TLP1_TLPthick0.05_TLPfreq15_TLPdelay551_Erosion' num2str(k) '/'];
    myVars={'elev.mat','marsh edge.mat','meansealevel.mat','AnnualAccretion.mat','SLR.mat','resiliencymetrics.mat'};
    for i=1:6
        load([filename myVars{i}]);
    end

    % define onset of interior drowning
    marshage=find(resiliencymetrics(551:701,11)==0,1)+551;
    if isempty(marshage)==1
        marshage=701;
    end
    for yr=552:701
        widthslope(yr)=(resiliencymetrics(yr,7)-resiliencymetrics(yr-1,7))/1; %change in width each year
    end
    threshold=find(widthslope(552:701)<widthslope(551)-5,1)-1+551;

    if isempty(threshold)==1
        threshold=marshage;
        color='k';
    end

    thresholdsmatrix(1,CO/10)=threshold;
    scatter(CO,threshold-551,110,'MarkerFaceColor',color,'MarkerEdgeColor','none')
    hold on
    color='#7570b3';

end
p2=plot(10:10:50,thresholdsmatrix(1:5)-551,'LineWidth',3,'Color','#7570b3');
hold on

%TLP and sill
color='#e6ab02';
for CO=[10,20,30,40,50]
    filename=['Outputs_CO' num2str(CO) '_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Sill1_SSCReduction0_silldelay551_TLP1_TLPthick0.05_TLPfreq15_TLPdelay551_Erosion' num2str(k) '/'];
    myVars={'elev.mat','marsh edge.mat','meansealevel.mat','AnnualAccretion.mat','SLR.mat','resiliencymetrics.mat'};
    for i=1:6
        load([filename myVars{i}]);
    end

    % define onset of interior drowning
    marshage=find(resiliencymetrics(551:701,11)==0,1)+551;
    if isempty(marshage)==1
        marshage=701;
    end
    for yr=552:701
        widthslope(yr)=(resiliencymetrics(yr,7)-resiliencymetrics(yr-1,7))/1; %change in width each year
    end
    threshold=find(widthslope(552:701)<widthslope(551)-5,1)-1+551;

    if isempty(threshold)==1
        threshold=marshage;
        color='k';
    end

    thresholdsmatrix(1,CO/10)=threshold;
    scatter(CO,threshold-551,110,'MarkerFaceColor',color,'MarkerEdgeColor','none')
    hold on
    color='#e6ab02';

end
p3=plot(10:10:50,thresholdsmatrix(1:5)-551,'LineWidth',3,'Color','#e6ab02');
hold on

% no restoration
color='#e7298a';
for CO=[10,20,30,40,50]
    filename=['Outputs_CO' num2str(CO) '_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Sill0_SSCReduction0_silldelay551_TLP0_TLPthick0_TLPfreq0_TLPdelay551_Erosion' num2str(k) '/'];
    myVars={'elev.mat','marsh edge.mat','meansealevel.mat','AnnualAccretion.mat','SLR.mat','resiliencymetrics.mat'};
    for i=1:6
        load([filename myVars{i}]);
    end

    % define onset of interior drowning
    marshage=find(resiliencymetrics(551:701,11)==0,1)+551;
    if isempty(marshage)==1
        marshage=701;
    end
    for yr=552:701
        widthslope(yr)=(resiliencymetrics(yr,7)-resiliencymetrics(yr-1,7))/1; %change in width each year
    end
    threshold=find(widthslope(552:701)<widthslope(551)-5,1)-1+551;

    if isempty(threshold)==1
        threshold=marshage;
        color='k';
    end

    thresholdsmatrix(1,CO/10)=threshold;
    scatter(CO,threshold-551,110,'MarkerFaceColor',color,'MarkerEdgeColor','none')
    hold on
    color='#e7298a';
end
p4=plot(10:10:50,thresholdsmatrix(1:5)-551,'LineWidth',3,'Color','#e7298a','LineStyle',':');
hold on

xlabel('External suspended sediment concentration (mg L^-^1)')
ylabel('Year at which interior drowning initiates (yr)')
xticks([10 20 30 40 50])
text(11,118,'a','FontWeight','bold')
l=legend([p4,p1,p2,p3],'No restoration','Breakwater','TLP','TLP and breakwater');
l.Location='southeast';
l.EdgeColor='none';
l.Color=[0.97 0.97 0.97 0.97];

% -----------------------------------------------------------------------------------------------------------------------------------------------------------
% fig b: sensitivity graph with SSC on x-axis, time to collapse on y-axis, with different restoration strategies as 4 diff lines
subplot(1,2,2)


%sill
color='#66a61e';
for CO=[10,20,30,40,50]
    filename=['Outputs_CO' num2str(CO) '_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Sill1_SSCReduction0_silldelay551_TLP0_TLPthick0_TLPfreq0_TLPdelay551_Erosion' num2str(k) '/'];
    myVars={'elev.mat','marsh edge.mat','meansealevel.mat','AnnualAccretion.mat','SLR.mat','resiliencymetrics.mat'};
    for i=1:6
        load([filename myVars{i}]);
    end

    % define onset of interior drowning
    marshage=find(resiliencymetrics(551:701,11)==0,1)+551;
    if isempty(marshage)==1
        marshage=701;
        color='k';
    end

    marshagematrix(1,CO/10)=marshage;
    scatter(CO,marshage-551,110,'MarkerFaceColor',color,'MarkerEdgeColor','none')
    hold on
    color='#66a61e';
end
p1=plot(10:10:50,marshagematrix(1:5)-551,'LineWidth',3,'Color','#66a61e');
hold on

%TLP
color='#7570b3';
for CO=[10,20,30,40,50]
    filename=['Outputs_CO' num2str(CO) '_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Sill0_SSCReduction0_silldelay551_TLP1_TLPthick0.05_TLPfreq15_TLPdelay551_Erosion' num2str(k) '/'];
    myVars={'elev.mat','marsh edge.mat','meansealevel.mat','AnnualAccretion.mat','SLR.mat','resiliencymetrics.mat'};
    for i=1:6
        load([filename myVars{i}]);
    end

    % define onset of interior drowning
    marshage=find(resiliencymetrics(551:701,11)==0,1)+551;
    if isempty(marshage)==1
        marshage=701;
        color='k';
    end

    marshagematrix(1,CO/10)=marshage;
    scatter(CO,marshage-551,110,'MarkerFaceColor',color,'MarkerEdgeColor','none')
    hold on
    color='#7570b3';
end
p2=plot(10:10:50,marshagematrix(1:5)-551,'LineWidth',3,'Color','#7570b3');
hold on

%TLP and sill
color='#e6ab02';
for CO=[20,30,40,50]
    filename=['Outputs_CO' num2str(CO) '_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Sill1_SSCReduction0_silldelay551_TLP1_TLPthick0.05_TLPfreq15_TLPdelay551_Erosion' num2str(k) '/'];
    myVars={'elev.mat','marsh edge.mat','meansealevel.mat','AnnualAccretion.mat','SLR.mat','resiliencymetrics.mat'};
    for i=1:6
        load([filename myVars{i}]);
    end

    % define onset of interior drowning
    marshage=find(resiliencymetrics(551:701,11)==0,1)+551;
    if isempty(marshage)==1
        marshage=701;
        color='k';
    end

    marshagematrix(1,CO/10)=marshage;
    scatter(CO,marshage-551,110,'MarkerFaceColor',color,'MarkerEdgeColor','none')
    hold on
    color='#e6ab02';

end
p3=plot(10:10:50,marshagematrix(1:5)-551,'LineWidth',3,'Color','#e6ab02');
hold on

%no restoration
color='#e7298a';
for CO=[10,20,30,40,50]
    filename=['Outputs_CO' num2str(CO) '_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Sill0_SSCReduction0_silldelay551_TLP0_TLPthick0_TLPfreq0_TLPdelay551_Erosion' num2str(k) '/'];
    myVars={'elev.mat','marsh edge.mat','meansealevel.mat','AnnualAccretion.mat','SLR.mat','resiliencymetrics.mat'};
    for i=1:6
        load([filename myVars{i}]);
    end

    % define onset of interior drowning
    marshage=find(resiliencymetrics(551:701,11)==0,1)+551;
    if isempty(marshage)==1
        marshage=701;
        color='k';
    end
    marshagematrix(1,CO/10)=marshage;
    scatter(CO,marshage-551,110,'MarkerFaceColor',color,'MarkerEdgeColor','none')
    hold on
    color='#e7298a';
end
p4=plot(10:10:50,marshagematrix(1:5)-551,'LineWidth',3,'Color','#e7298a','LineStyle',':');
hold on

xlabel('External suspended sediment concentration (mg L^-^1)')
ylabel('Year at which marsh collapses (yr)')
l2=legend([p4,p1,p2,p3],'No restoration','Breakwater','TLP','TLP and breakwater');
l2.Location='southeast';
l2.EdgeColor='none';
l2.Color=[0.97 0.97 0.97 0.97];
xticks([10 20 30 40 50])
text(11,148,'b','FontWeight','bold')

%% Fig S6: sensitivity analysis of SLRA on vertical resilience
CO=10;
scatterRSLRA=[0.01881,0.03536,0.08023,0.1106,0.1592];
k=5.3e-09;
figure
subplot(1,2,1)
%sill
color='#66a61e';
holder=1;
for RSLRA=[0.01881,0.03536,0.08023,0.1106,0.1592]
    filename=['Outputs_CO' num2str(CO) '_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Sill1_SSCReduction0_silldelay551_TLP0_TLPthick0_TLPfreq0_TLPdelay551_Erosion' num2str(k) '/'];
    myVars={'elev.mat','marsh edge.mat','meansealevel.mat','AnnualAccretion.mat','SLR.mat','resiliencymetrics.mat'};
    for i=1:6
        load([filename myVars{i}]);
    end

    % define onset of interior drowning
    marshage=find(resiliencymetrics(551:701,11)==0,1)+551;
    if isempty(marshage)==1
        marshage=701;
    end
    for yr=552:701
        widthslope(yr)=(resiliencymetrics(yr,7)-resiliencymetrics(yr-1,7))/1; %change in width each year
    end
    threshold=find(widthslope(552:701)<widthslope(551)-5,1)-1+551;

    if isempty(threshold)==1
        threshold=marshage;
        color='k';
    end

    thresholdsmatrix(1,holder)=threshold;
    scatter(RSLRA,threshold-551,110,'MarkerFaceColor',color,'MarkerEdgeColor','none')
    hold on
    color='#66a61e';
    holder=holder+1;
    clear elev marsh edge meansealevel AnnualAccretion SLR.mat resiliencymetrics
end
p1=plot(scatterRSLRA,thresholdsmatrix(1:5)-551,'LineWidth',3,'Color','#66a61e');
hold on

%TLP
color='#7570b3';
holder=1;
for RSLRA=[0.01881,0.03536,0.08023,0.1106,0.1592]
    filename=['Outputs_CO' num2str(CO) '_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Sill0_SSCReduction0_silldelay551_TLP1_TLPthick0.05_TLPfreq15_TLPdelay551_Erosion' num2str(k) '/'];
    myVars={'elev.mat','marsh edge.mat','meansealevel.mat','AnnualAccretion.mat','SLR.mat','resiliencymetrics.mat'};
    for i=1:6
        load([filename myVars{i}]);
    end

    % define onset of interior drowning
    marshage=find(resiliencymetrics(551:701,11)==0,1)+551;
    if isempty(marshage)==1
        marshage=701;
    end
    for yr=552:701
        widthslope(yr)=(resiliencymetrics(yr,7)-resiliencymetrics(yr-1,7))/1; %change in width each year
    end
    threshold=find(widthslope(552:701)<widthslope(551)-5,1)-1+551;

    if isempty(threshold)==1
        threshold=marshage;
        color='k';
    end

    thresholdsmatrix(1,holder)=threshold;
    scatter(RSLRA,threshold-551,110,'MarkerFaceColor',color,'MarkerEdgeColor','none')
    hold on
    color='#7570b3';
    holder=holder+1;
    clear elev marsh edge meansealevel AnnualAccretion SLR.mat resiliencymetrics
end
p2=plot(scatterRSLRA,thresholdsmatrix(1:5)-551,'LineWidth',3,'Color','#7570b3');
hold on

%TLP and sill
color='#e6ab02';
holder=1;
for RSLRA=[0.01881,0.03536,0.08023,0.1106,0.1592]
    filename=['Outputs_CO' num2str(CO) '_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Sill1_SSCReduction0_silldelay551_TLP1_TLPthick0.05_TLPfreq15_TLPdelay551_Erosion' num2str(k) '/'];
    myVars={'elev.mat','marsh edge.mat','meansealevel.mat','AnnualAccretion.mat','SLR.mat','resiliencymetrics.mat'};
    for i=1:6
        load([filename myVars{i}]);
    end

    % define onset of interior drowning
    marshage=find(resiliencymetrics(551:701,11)==0,1)+551;
    if isempty(marshage)==1
        marshage=701;
    end
    for yr=552:701
        widthslope(yr)=(resiliencymetrics(yr,7)-resiliencymetrics(yr-1,7))/1; %change in width each year
    end
    threshold=find(widthslope(552:701)<widthslope(551)-5,1)-1+551;

    if isempty(threshold)==1
        threshold=marshage;
        color='k';
    end

    thresholdsmatrix(1,holder)=threshold;
    scatter(RSLRA,threshold-551,110,'MarkerFaceColor',color,'MarkerEdgeColor','none')
    hold on
    color='#e6ab02';
    holder=holder+1;
    clear elev marsh edge meansealevel AnnualAccretion SLR.mat resiliencymetrics
end
p3=plot(scatterRSLRA,thresholdsmatrix(1:5)-551,'LineWidth',3,'Color','#e6ab02');
hold on

% no restoration
color='#e7298a';
holder=1;
for RSLRA=[0.01881,0.03536,0.08023,0.1106,0.1592]
    filename=['Outputs_CO' num2str(CO) '_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Sill0_SSCReduction0_silldelay551_TLP0_TLPthick0_TLPfreq0_TLPdelay551_Erosion' num2str(k) '/'];
    myVars={'elev.mat','marsh edge.mat','meansealevel.mat','AnnualAccretion.mat','SLR.mat','resiliencymetrics.mat'};
    for i=1:6
        load([filename myVars{i}]);
    end

    % define onset of interior drowning
    marshage=find(resiliencymetrics(551:701,11)==0,1)+551;
    if isempty(marshage)==1
        marshage=701;
    end
    for yr=552:701
        widthslope(yr)=(resiliencymetrics(yr,7)-resiliencymetrics(yr-1,7))/1; %change in width each year
    end
    threshold=find(widthslope(552:701)<widthslope(551)-5,1)-1+551;

    if isempty(threshold)==1
        threshold=marshage;
        color='k';
    end

    thresholdsmatrix(1,holder)=threshold;
    scatter(RSLRA,threshold-551,110,'MarkerFaceColor',color,'MarkerEdgeColor','none')
    hold on
    color='#e7298a';
    holder=holder+1;
    clear elev marsh edge meansealevel AnnualAccretion SLR.mat resiliencymetrics
end
p4=plot(scatterRSLRA,thresholdsmatrix(1:5)-551,'LineWidth',3,'Color','#e7298a','LineStyle',':');
hold on

xlabel('Sea-level rise acceleration rate (mm yr^-^2)')
ylabel('Year at which interior drowning initiates (yr)')
xticks([0 0.04 0.08 0.12 0.16])
text(0.005,148,'a','FontWeight','bold')
l=legend([p4,p1,p2,p3],'No restoration','Breakwater','TLP','TLP and breakwater');
l.Location='southwest';
l.EdgeColor='none'; %grey legend outline
l.Color=[0.97 0.97 0.97 0.97];
xlim([0 0.18])
ylim([40 150])

% -----------------------------------------------------------------------------------------------------------------------------------------------------------
% fig b: sensitivity graph with SLRA on x-axis, time to collapse on y-axis, with different restoration strategies as 4 diff lines

subplot(1,2,2)


%sill
color='#66a61e';
holder=1;
for RSLRA=[0.01881,0.03536,0.08023,0.1106,0.1592]
    filename=['Outputs_CO' num2str(CO) '_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Sill1_SSCReduction0_silldelay551_TLP0_TLPthick0_TLPfreq0_TLPdelay551_Erosion' num2str(k) '/'];
    myVars={'elev.mat','marsh edge.mat','meansealevel.mat','AnnualAccretion.mat','SLR.mat','resiliencymetrics.mat'};
    for i=1:6
        load([filename myVars{i}]);
    end

    % define onset of interior drowning
    marshage=find(resiliencymetrics(551:701,11)==0,1)+551;
    if isempty(marshage)==1
        marshage=701;
        color='k';
    end

    marshagematrix(1,holder)=marshage;
    scatter(RSLRA,marshage-551,110,'MarkerFaceColor',color,'MarkerEdgeColor','none')
    hold on
    color='#66a61e';
    holder=holder+1;
    clear elev marsh edge meansealevel AnnualAccretion SLR.mat resiliencymetrics
end
p1=plot(scatterRSLRA,marshagematrix(1:5)-551,'LineWidth',3,'Color','#66a61e');
hold on

%TLP
color='#7570b3';
holder=1;
for RSLRA=[0.01881,0.03536,0.08023,0.1106,0.1592]
    filename=['Outputs_CO' num2str(CO) '_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Sill0_SSCReduction0_silldelay551_TLP1_TLPthick0.05_TLPfreq15_TLPdelay551_Erosion' num2str(k) '/'];
    myVars={'elev.mat','marsh edge.mat','meansealevel.mat','AnnualAccretion.mat','SLR.mat','resiliencymetrics.mat'};
    for i=1:6
        load([filename myVars{i}]);
    end

    % define onset of interior drowning
    marshage=find(resiliencymetrics(551:701,11)==0,1)+551;
    if isempty(marshage)==1
        marshage=701;
        color='k';
    end

    marshagematrix(1,holder)=marshage;
    scatter(RSLRA,marshage-551,110,'MarkerFaceColor',color,'MarkerEdgeColor','none')
    hold on
    color='#7570b3';
    holder=holder+1;
    clear elev marsh edge meansealevel AnnualAccretion SLR.mat resiliencymetrics
end
p2=plot(scatterRSLRA,marshagematrix(1:5)-551,'LineWidth',3,'Color','#7570b3');
hold on

%TLP and sill
color='#e6ab02';
holder=1;
for RSLRA=[0.01881,0.03536,0.08023,0.1106,0.1592]
    filename=['Outputs_CO' num2str(CO) '_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Sill1_SSCReduction0_silldelay551_TLP1_TLPthick0.05_TLPfreq15_TLPdelay551_Erosion' num2str(k) '/'];
    myVars={'elev.mat','marsh edge.mat','meansealevel.mat','AnnualAccretion.mat','SLR.mat','resiliencymetrics.mat'};
    for i=1:6
        load([filename myVars{i}]);
    end

    % define onset of interior drowning
    marshage=find(resiliencymetrics(551:701,11)==0,1)+551;
    if isempty(marshage)==1
        marshage=701;
        color='k';
    end

    marshagematrix(1,holder)=marshage;
    scatter(RSLRA,marshage-551,110,'MarkerFaceColor',color,'MarkerEdgeColor','none')
    hold on
    color='#e6ab02';
    holder=holder+1;
    clear elev marsh edge meansealevel AnnualAccretion SLR.mat resiliencymetrics
end
p3=plot(scatterRSLRA,marshagematrix(1:5)-551,'LineWidth',3,'Color','#e6ab02');
hold on

%no restoration
color='#e7298a';
holder=1;
for RSLRA=[0.01881,0.03536,0.08023,0.1106,0.1592]
    filename=['Outputs_CO' num2str(CO) '_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Sill0_SSCReduction0_silldelay551_TLP0_TLPthick0_TLPfreq0_TLPdelay551_Erosion' num2str(k) '/'];
    myVars={'elev.mat','marsh edge.mat','meansealevel.mat','AnnualAccretion.mat','SLR.mat','resiliencymetrics.mat'};
    for i=1:6
        load([filename myVars{i}]);
    end

    % define onset of interior drowning
    marshage=find(resiliencymetrics(551:701,11)==0,1)+551;
    if isempty(marshage)==1
        marshage=701;
        color='k';
    end
    marshagematrix(1,holder)=marshage;
    scatter(RSLRA,marshage-551,110,'MarkerFaceColor',color,'MarkerEdgeColor','none')
    hold on
    color='#e7298a';
    holder=holder+1;
    clear elev marsh edge meansealevel AnnualAccretion SLR.mat resiliencymetrics
end
p4=plot(scatterRSLRA,marshagematrix(1:5)-551,'LineWidth',3,'Color','#e7298a','LineStyle',':');
hold on
xlabel('Sea-level rise acceleration rate (mm yr^-^2)')
ylabel('Year at which marsh collapses (yr)')
l2=legend([p4,p1,p2,p3],'No restoration','Breakwater','TLP','TLP and breakwater');
l2.Location='southwest';
l2.EdgeColor='none'; %grey legend outline
l2.Color=[0.97 0.97 0.97 0.97];
xticks([0 0.04 0.08 0.12 0.16])
xlim([0 0.18])
text(0.005,148,'b','FontWeight','bold')

