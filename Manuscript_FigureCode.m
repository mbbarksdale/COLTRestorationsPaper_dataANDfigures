%% fig 1; OM depth profile from spinup
percentorganic_spinup=zeros(550,5891);
percentorganic_spinup=(organic_dep_autoch+organic_dep_alloch)./(organic_dep_autoch+organic_dep_alloch+mineral_dep)*100;

figure
set(gca,'color','#3F4F4D','XColor','w','YColor','w')
c=colorbar;
colormap copper
c.Color="w";
clim([0, 45]);
xlim([1 500])
ylim([-0.2 1.4])
xlabel('Distance from marsh edge (m)')

for yr=2:550
    X = [1:500; 2:501; 2:501; 1:500];
    Y = [elev(yr-1, 5001:5500); elev(yr,5001:5500); elev(yr,5002:5501); elev(yr-1,5002:5501)];
    C = repmat(percentorganic_spinup(yr,5001:5500),4,1);
    patch(X,Y,C,'EdgeColor','interp');
    hold on
end

%% Fig 2; validation of OM-depth profiles

table=readtable("GoodwinEdge.csv");
GillenGoodwinEdge=table2array(table(:,3:4));
table=readtable("GoodwinInterior.csv");
GillenGoodwinInterior=table2array(table(:,3:4));
table=readtable("CatlettEdge.csv");
GillenCatlettEdge=table2array(table(:,3:4));
table=readtable("CatlettInterior.csv");
GillenCatlettInterior=table2array(table(:,3:4));

figure
scatter(GillenGoodwinInterior(:,2),GillenGoodwinInterior(:,1),'filled','LineWidth',3,'MarkerFaceColor','#276419')
hold on
scatter(GillenCatlettInterior(:,2),GillenCatlettInterior(:,1),'filled','LineWidth',3,'MarkerFaceColor','#4d9221')
hold on
scatter(GillenGoodwinEdge(:,2),GillenGoodwinEdge(:,1),'filled','LineWidth',3,'MarkerFaceColor','#7fbc41')
hold on
scatter(GillenCatlettEdge(:,2),GillenCatlettEdge(:,1),'filled','LineWidth',3,'MarkerFaceColor','#b8e186')
hold on

table=readtable("Holmquist_2018_depth_series_data.csv");
coredata=table2array(table(:,3:6));
avgcore=zeros(10,1);
coredata(:,5)=mean(coredata(:,1:2),2);
coredataMidAtlantic=coredata([1:472,10220:10249,12795:12967,13432:13671,14776:15374],:); %just takes data from cores from VA, DE, or MD, and those that are ID'ed as salt marsh
run=1;
for i=1:10:99
    increment=find(coredataMidAtlantic(:,5)>i & coredataMidAtlantic(:,5)<i+9); 
    avgcore(run,1)=mean(coredataMidAtlantic(increment,4),'omitmissing');
    %scatter(avgcore(run,1),i+4,'filled','MarkerFaceColor','k')
    hold on
    %str = sprintf('%.f',length(increment)); %uncomment these 3 lines of  IF CODE IF WANT TO KNOW HOW MANY datapoints (n=) each point averaged
    %totalstr = strjoin(string(str));
    %text(avgcore(run,1)+0.01,i+4,totalstr,'FontSize',8)
    %hold on
    run=run+1;
end
scatter(avgcore(1:10,1)*100,5:10:95,40,'MarkerEdgeColor','k')

load("MarshStrat_transectspinup_initialRSLR1_rampedfinalRSLR3_CO10.mat")

fractionorganic_spinup=zeros(550,length(mineral_dep));
fractionorganic_spinup=(organic_dep_autoch+organic_dep_alloch)./(organic_dep_autoch+organic_dep_alloch+mineral_dep);

for yr=2:550 %changed from 2:550
    run=1;
    for profile=[5001,5250,5500] %[5001,5250,5500]
        elev_OM(yr-1,run)=fractionorganic_spinup(yr,profile);
        elev_OM(yr-1,run+1)=elev(549,profile)*100-(mean(elev(yr-1:yr,profile))*100); %=elev(550,profile)*100-(mean(elev(yr-1:yr,profile))*100);
        run=run+2;
    end
end

%different matrices for 3 different profiles. average elevation b/w year and yr-1 in 1st column and OM in second  
p1=plot(elev_OM(1:549,1)*100,elev_OM(1:549,2),'LineWidth',3,'Color','#8e0152'); %plot(elev_OM(:,1),elev_OM(:,2),'LineWidth',3)
hold on
p2=plot(elev_OM(1:549,3)*100,elev_OM(1:549,4),'LineWidth',3,'Color','#de77ae');
hold on
plot(elev_OM(1:549,5)*100,elev_OM(1:549,6),'LineWidth',3,'Color','#f1b6da')
hold on
xlabel('Organic Matter (Fraction)')
ylabel('Depth (cm)')
hold on

ylim([0 160])
l=legend('Observed Marsh Interior (Goodwin Island, VA)','Observed Marsh Interior (Catlett Island, VA)','Observed Marsh Edge (Goodwin Island, VA)','Observed Marsh Edge (Catlett Island, VA)','Mid-Atlantic Salt Marsh Average','Spinup Marsh Edge','Spinup Middle Marsh','Spinup Marsh-Upland Boundary','Location','southeast');
l.EdgeColor='none'; %grey legend outline
l.Color=[0.94 0.94 0.94 0.94];
set(gca, 'YDir','reverse')
%% Fig 3; SLR vs accretion ratio
figure
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]); %makes really large figure box (takes up most of computer screen)
subplot(2,2,[1,3])
RSLRA=0.08023;

%create first panel of SLR vs accretion excess 
s=subaxis(2,2,[1,3],'SpacingHoriz',0.1,'SpacingVert',0.1); %make subpanels closer together
color=[0.99608 0.78824 0.89412]; %start with color for low erosion scenario and will increase for every iteration of for-loop
for k=[1.65e-09,2.82e-09,5.3e-09,7.8e-09,9.95e-09] %creates for-loop to cycle through 5 erosion scenarios (changes Be, or Ke, the erosion coefficient)
    filename=['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Breakwater0_SSCReduction0_breakwaterdelay551_TLP0_TLPthick0_TLPfreq0_TLPdelay551_Erosion' num2str(k) '/']; %calls the folder
    myVars={'elev.mat','marsh edge.mat','meansealevel.mat','AnnualAccretion.mat','SLR.mat','resiliencymetrics.mat'}; %establishes the files to pull from folder above
    for i=1:6 %creates a for-loop to call each of the files
        load([filename myVars{i}]); %calls each file as part of myVars
    end
    %create empty matrices
    erosion=zeros(701,1);
    accretiontoSLR=zeros(701,1);
    AvgAnnualAccretion=zeros(701,1);
    marshage=find(resiliencymetrics(551:701,7)==0,1)-1+551; %find marsh age where resiliencymetrics 7th column (marsh width) equals zero, then add 1 (for marsh lifespan)
    if isempty(marshage)==1 %if the marsh width never collapses to 0...
        marshage=701; %sets marsh age to the oldest it can be
    end
    for yr=2:marshage %start for-loop that runs year by year
        erosion(yr)=Marsh_edge(yr)-Marsh_edge(yr-1); %erosion for each year is the marsh edge of that year minus the marsh edge position from the year before
        AvgAnnualAccretion(yr)=mean(AnnualAccretion(yr,Marsh_edge(yr):5501));%find the avg annual accretion for each year
        accretiontoSLR(yr)=AvgAnnualAccretion(yr)/SLR(yr); %rate of accretion to sea-level rise
    end
    if k==9.95e-09
        marshage=marshage-1;
        color=[0.44706 0 0.22353];
    end
    p=plot(SLR(551:marshage)*1000,accretiontoSLR(551:marshage),'LineWidth',7,'Color',color);
    hold on
    color=color+[-0.14 -0.2 -0.15];
end
ylim([-0.8 1.3])
xlim([3 12.5])
ax=gca;
ax.FontSize=22;
% ax.XColor=[0.5 0.5 0.5];
% ax.YColor=[0.5 0.5 0.5];
yline(1,'--',{'Accretion surplus'},'LineWidth',2,'FontSize',22);
dim=[0.336 0.685 .1 .1]; %sets dimensions to write accretion deficit
a=annotation('textbox',dim,'String','Accretion deficit','EdgeColor','none','FontSize',22);
dim2=[0.1 0.806 .1 .1]; %sets dimensions to write "a" in upper left of plot
b=annotation('textbox',dim2,'String','a','EdgeColor','none','FontSize',22,'FontWeight','bold');
l=legend('Erosion_{initial} = 0.5 m yr^-^1','Erosion_{initial} = 1 m yr^-^1','Erosion_{initial} = 2 m yr^-^1','Erosion_{initial} = 3 m yr^-^1','Erosion_{initial} = 4 m yr^-^1');
l.Location='southwest';
l.EdgeColor='none'; %grey legend outline
l.Color=[0.97 0.97 0.97 0.97];
l.FontSize=22;
ylabel('Accretion excess (dimensionless)')
xlabel('Rate of sea-level rise (mm yr^-^1)')
yticks([-1 -0.5 0 0.5 1]) %sets the ytick labels

%-----------------------------------------------------------------------------------------------------------------------------
% figure 3b graph erosion on x vs SLR threshold rate on y

subplot(2,2,2)
s=subaxis(2,2,2,'SpacingHoriz',0.1,'SpacingVert',0.1); %limits the padding between subplots
color=[0.99608 0.78824 0.89412];
edgecolor='none';
acceleration=1; %will cycle through acceleration
RSLRA=0.08023;
erosion=1; %start with erosion=1 and will increase for every iteration of for-loop, used to designate correct matrix row to place values in each iteration of loop
for k=[1.65e-09,2.82e-09,5.3e-09,7.8e-09,9.95e-09,9.900E-10,3.812E-09,6.300E-09,9.090E-09] %,1.884E-09, 2.118E-09, 2.352E-09, 2.586E-09, 3.812E-09, 4.804E-09, 5.800E-09, 6.300E-09, 6.800E-09, 7.300E-09, 8.230E-09, 8.660E-09, 9.090E-09, 9.520E-09, 7.300E-09, 8.230E-09, 8.660E-09, 9.090E-09, 9.520E-09, 3.300E-10, 6.600E-10, 9.900E-10, 1.320E-09, 3.316E-09, 4.308E-09]
    filename=['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Breakwater0_SSCReduction0_breakwaterdelay551_TLP0_TLPthick0_TLPfreq0_TLPdelay551_Erosion' num2str(k) '/'];
    myVars={'elev.mat','marsh edge.mat','meansealevel.mat','AnnualAccretion.mat','SLR.mat','resiliencymetrics.mat','fetch.mat','bgb.mat'};
    for i=1:8
        load([filename myVars{i}]);
    end
    if k==9.95e-09
        color=[0.44706 0 0.22353]; %designates special color for high erosion rate
    end

    marshage=find(resiliencymetrics(551:701,7)==0,1)-1; %find marsh age where resiliencymetrics 7th column (marsh width) equals zero, then add 1 (for marsh lifespan)
    if isempty(marshage)==1
        marshage=150; %here marshage is 150 (as opposed to how calculated for Fig 2a, but doesn't matter because just using these marsh age values to correlate to SLR at that year)
    end
    drowningfound=0;
    for yr=552:701 %start at year 552 because year 550 to 551 typically has large slope break as marsh suddenly begins to erode
        widthslope(yr)=(resiliencymetrics(yr,7)-resiliencymetrics(yr-1,7))/1; %calculates the slope of the marsh width curve through time (for each year)
        marshedgeerosion(yr)=Marsh_edge(yr)-Marsh_edge(yr-1);
        drowningstart=find(bgb(yr,Marsh_edge(yr)+1:5500)==0);
        if isempty(drowningstart)==0 && drowningfound==0
            drowningstartmatrix(1,erosion)=yr;
            drowningfound=1;
        end
    end

    threshold=find(widthslope(552:701)<widthslope(551)-5,1)-1; %find where a change of 10+ occurs in widthslope from one year to the next
    meanerosion=mean(marshedgeerosion(552:marshage+549)); %calculate the average erosion rate for the entire scenario
    if erosion==1 && isempty(threshold)==1 %if lowest erosion rate (first loop through) and there is no slope break...
        threshold=marshage; %...then assume that threshold doesn't occur until marsh age (this isn't the case for any intermediate SLR scenarios)
    else if erosion>1 && isempty(threshold)==1 && SLR(marshage+551)>thresholdsmatrix(1,erosion-1) %else for erosion greater than lowest, if threshold is empty and SLR of final marsh age is greater than SLR of threshold for low erosion scenario...
            threshold=marshage; %...then threshold equals marshage again (this does not occur for any interm. SLRA scenarios shown right now in Fig 2, but is a failsafe to continue to show that SLR threshold rate increases for increasing erosion
    else if isempty(threshold)==1 %otherwise, if threshold is empty (but does not meet the conditions listed above...
            continue %...then continue to next iteration of loop
    end
    end
    end

    if isempty(drowningstartmatrix(1,erosion))==1
        drowningstartmatrix(1,erosion)=701;
        color=[0 0 0];
    end

    if erosion>5
        color=[1 1 1];
        edgecolor=[0.5 0.5 0.5];
    end

    meanerosionmatrix(1,erosion)=meanerosion; %fill in meanerosionmatrix with avg erosion for that erosion scnenario
    thresholdsmatrix(1,erosion)=SLR(drowningstartmatrix(1,erosion))*1000; %fill in thresholdsmatrix with SLR rate at which threshold occurs
    %scatter(meanerosion,SLR(threshold+551)*1000,100,'MarkerFaceColor','k','MarkerEdgeColor','none');
    scatter(meanerosion,thresholdsmatrix(erosion),130,'MarkerFaceColor',color,'MarkerEdgeColor',edgecolor);
    color=color+[-0.14 -0.2 -0.15]; %change color matrix for next iteration of loop
    erosion = erosion + 1; %add 1 to erosion for next iteration of loop
    hold on
    clear AnnualAccretion elev msl resiliencymetrics SLR threshold widthslope yr meanerosion meanerosion1 fetcherosion marshedgeerosion drowningfound
end

[xData, yData] = prepareCurveData(meanerosionmatrix, thresholdsmatrix);

% Set up fittype and options.
% ft = fittype( 'exp2' );
% opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% opts.Display = 'Off';
%
% % Fit model to data.
% [fitresult, gof] = fit(xData, yData, ft, opts )

% Plot fit with data.
% h=plot(fitresult);
% set(h,'LineWidth',3,'Color',[.5 .5 .5]);
% uistack(h,'bottom') %brings the square above the line for the legend
% l=legend('example');
% set(l,'visible','off')

% beta0=[fitresult.a,fitresult.b,fitresult.c,fitresult.d];
% modelFun=@(b,x) b(1) * exp(b(2) * x) + b(3)*exp(b(4)*x);
% mdl=fitnlm(xData,yData,modelFun,beta0)


coefficients = polyfit(meanerosionmatrix,thresholdsmatrix, 2);
% Create a new x axis with exactly 1000 points (or whatever you want).
xFit = linspace(min(meanerosionmatrix), max(meanerosionmatrix), 1000);
% Get the estimated yFit value for each of those 1000 new x locations.
yFit = polyval(coefficients , xFit);
X=[meanerosionmatrix',meanerosionmatrix'.^2];
mdl = fitlm(X, thresholdsmatrix)
p=plot(xFit, yFit, 'color', [0.7 0.7 0.7], 'LineWidth', 7); % Plot fitted line.
uistack(p,'bottom') %pushes line below markers

ylim([8.5 9.3])
hold on; % Set hold on so the next plot does not blow away the one we just drew.
coefficients=table2array(mdl.Coefficients);
str = sprintf('SLR_{drowning.start} = %.2f(ero.)^2 %.2f(ero.) + %.2f', round(coefficients(3,1), 3), round(coefficients(2,1), 3), round(coefficients(1,1), 3));
%str= {'SLR_{drowning.start} = ' round(coefficients(3,1),2) '(erosion)^2' round(coefficients(2,1)) '(erosion) +' round(coefficients(1,1))};
totalstr = strjoin(string(str));
%text(0.2,8.55,totalstr,'FontSize',18,'Color','#48CAE4');
text(0.3,9.22,totalstr,'FontSize',22,'Color',[.5 .5 .5]);
dim2=[0.55 0.808 .1 .1];
b=annotation('textbox',dim2,'String','b','EdgeColor','none','FontSize',22,'FontWeight','bold');
ylabel({'SLR rate at which interior';'drowning initiates (mm yr^-^1)'})

ax=gca;
ax.FontSize=22;
box on
xlim([0 5])
hold on %-------------------

%add breakwaters dotted line
%     filename=['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Breakwater1_SSCReduction20_breakwaterdelay551_TLP0_TLPthick0_TLPfreq0_TLPdelay551_Erosion5.3e-09/'];
%     myVars={'elev.mat','marsh edge.mat','meansealevel.mat','AnnualAccretion.mat','SLR.mat','resiliencymetrics.mat','VASEAmatrix.mat','fetch.mat'};
%     for i=1:8
%         load([filename myVars{i}]);
%     end
%
%     marshage=find(resiliencymetrics(551:701,7)==0,1)-1; %find marsh age where resiliencymetrics 7th column (marsh width) equals zero, then add 1 (for marsh lifespan)
%     if isempty(marshage)==1
%         marshage=150; %here marshage is 150 (as opposed to how calculated for Fig 2a, but doesn't matter because just using these marsh age values to correlate to SLR at that year)
%     end
%     for yr=552:701 %start at year 552 because year 550 to 551 typically has large slope break as marsh suddenly begins to erode
%         widthslope(yr)=(resiliencymetrics(yr,7)-resiliencymetrics(yr-1,7))/1; %calculates the slope of the marsh width curve through time (for each year)
%         marshedgeerosion(yr)=Marsh_edge(yr)-Marsh_edge(yr-1);
%     end
%     threshold=find(widthslope(552:701)<widthslope(551)-5,1)-1; %find where a change of 10+ occurs in widthslope from one year to the next
%     if erosion==1 && isempty(threshold)==1 %if lowest erosion rate (first loop through) and there is no slope break...
%         threshold=marshage; %...then assume that threshold doesn't occur until marsh age (this isn't the case for any intermediate SLR scenarios)
%     else if erosion>1 && isempty(threshold)==1 && SLR(marshage+551)>thresholdsmatrix(1,erosion-1) %else for erosion greater than lowest, if threshold is empty and SLR of final marsh age is greater than SLR of threshold for low erosion scenario...
%         threshold=marshage; %...then threshold equals marshage again (this does not occur for any interm. SLRA scenarios shown right now in Fig 2, but is a failsafe to continue to show that SLR threshold rate increases for increasing erosion
%     % else if isempty(threshold)==1 %otherwise, if threshold is empty (but does not meet the conditions listed above...
%     %     continue %...then continue to next iteration of loop
%     % end
%     end
%     end
% scatter(0,SLR(threshold+551)*1000,100,'MarkerFaceColor','#66a61e','MarkerEdgeColor','none')
%-----------------------------------------------------------------------------------------------------------------------------
%figure 3c relating erosion to time of drowning and survival
subplot(2,2,4)
s=subaxis(2,2,4,'SpacingHoriz',0.1,'SpacingVert',0.1); %limits the padding between subplots
% yline(threshold+551,'LineWidth',4,'Color','#66a61e','LineStyle',':')
% hold on
% yline(marshage+551,'LineWidth',4,'Color','#66a61e','LineStyle',':')
% hold on
%creates legend for the two lines
p1=plot([0.2 0.5],[90 90],'Color','#006CA5','LineWidth',7);
hold on
p2=plot([0.2 0.5],[81.7 81.7],'Color','#48CAE4','LineWidth',7);
hold on
s1=scatter(0.35,90,145,'Marker','square','MarkerEdgeColor',[0.71 0.38824 0.59412],'MarkerFaceColor','w','LineWidth',2.5);
uistack(s1,'top') %brings the square above the line for the legend
hold on
s2=scatter(0.35,81.7,135,'Marker','^','MarkerEdgeColor',[0.71 0.38824 0.59412],'MarkerFaceColor','w','LineWidth',2.5);
hold on
text(0.6,90,'Years until total collapse','EdgeColor','none','FontSize',22);
text(0.6,81.7,'Years until initiation of interior drowning','EdgeColor','none','FontSize',22);
hold on
rectangle('Position',[0.13 77.5 4.3 15.5],'FaceColor',[0.8 0.8 0.8],'FaceAlpha',0.2,'EdgeColor','none')

edgecolor='none';
color=[0.99608 0.78824 0.89412]; %for no restoration color (erosion changes)
acceleration=1;
RSLRA=0.08023;
erosion=1;
for k=[1.65e-09,2.82e-09,5.3e-09,7.8e-09,9.95e-09,9.900E-10,3.812E-09,6.300E-09,9.090E-09]
    filename=['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Breakwater0_SSCReduction0_breakwaterdelay551_TLP0_TLPthick0_TLPfreq0_TLPdelay551_Erosion' num2str(k) '/'];
    myVars={'elev.mat','marsh edge.mat','meansealevel.mat','AnnualAccretion.mat','SLR.mat','resiliencymetrics.mat','fetch.mat','bgb.mat'};
    for i=1:8
        load([filename myVars{i}]);
    end
    if k==9.95e-09
        color=[0.44706 0 0.22353];
    end

    marshage=find(resiliencymetrics(551:701,11)==0,1);
    if isempty(marshage)==1
        marshage=150;
    end

    drowningfound=0;

    for yr=552:701
        widthslope(yr)=(resiliencymetrics(yr,7)-resiliencymetrics(yr-1,7))/1; %change in width each year
        marshedgeerosion(yr)=Marsh_edge(yr)-Marsh_edge(yr-1);
        drowningstart=find(bgb(yr,Marsh_edge(yr)+1:5500)==0)-551;
        if isempty(drowningstart)==0 && drowningfound==0
            drowningstartmatrix(1,erosion)=yr-551;
            drowningfound=1;
        end
    end

    threshold=find(widthslope(552:701)<widthslope(551)-5,1)-1;
    meanerosion=mean(marshedgeerosion(552:marshage+551-2));

    if erosion==1 && isempty(threshold)==1
        threshold=marshage;
    else if erosion>1 && isempty(threshold)==1 && marshage>thresholdsmatrix(1,erosion-1)
            threshold=marshage;
    else if isempty(threshold)==1
            continue
    end
    end
    end

    if isempty(drowningstartmatrix(1,erosion))==1
        drowningstartmatrix(1,erosion)=150;
        color=[0 0 0];
    end

    if erosion>5
        color=[1 1 1];
        edgecolor=[0.5 0.5 0.5];
    end
    thresholdsmatrix(1,erosion)=drowningstartmatrix(1,erosion); %fill in thresholdsmatrix with SLR rate at which threshold occurs
    meanerosionmatrix(1,erosion)=meanerosion;
    %thresholdsmatrix(1,erosion)=threshold;
    scatter(meanerosion,thresholdsmatrix(erosion),130,'MarkerFaceColor',color,'MarkerEdgeColor',edgecolor,'Marker','^')
    color=color+[-0.14 -0.2 -0.15]; %for TLP / erosion scenarios
    erosion = erosion + 1;
    hold on
    clear AnnualAccretion elev msl resiliencymetrics SLR threshold widthslope yr meanerosion meanerosion1 fetcherosion marshedgeerosion bgb
end

% Get coefficients of a line fit through the data.
if length(meanerosionmatrix) == length(thresholdsmatrix)
    meanerosionmatrix1=meanerosionmatrix;
else
    meanerosionmatrix1=meanerosionmatrix(1:length(thresholdsmatrix));
end
coefficients = polyfit(meanerosionmatrix1,thresholdsmatrix, 2);
% Create a new x axis with exactly 1000 points (or whatever you want).
xFit = linspace(min(meanerosionmatrix1), max(meanerosionmatrix1), 1000);
% Get the estimated yFit value for each of those 1000 new x locations.
yFit = polyval(coefficients , xFit);
X=[meanerosionmatrix',meanerosionmatrix'.^2];
mdl2 = fitlm(X, thresholdsmatrix)

% Plot everything.
hold on; % Set hold on so the next plot does not blow away the one we just drew.
p=plot(xFit, yFit, 'color', '#48CAE4', 'LineWidth', 7); % Plot fitted line.

coefficients=table2array(mdl2.Coefficients);
str = sprintf('Year_{drowning.start} = %.2f(ero.)^2 %.2f(ero.) + %.2f', round(coefficients(3,1), 3), round(coefficients(2,1), 3), round(coefficients(1,1), 3));
totalstr = strjoin(string(str));
text(0.3,62,totalstr,'FontSize',22,'Color','#20bad9'); %#20bad9 %was #48CAE4
dim2=[0.55 0.355 .1 .1];
b=annotation('textbox',dim2,'String','c','EdgeColor','none','FontSize',22,'FontWeight','bold');
xlabel('Average rate of erosion (m yr^-^1)')
ylabel('Model Run Year')
ax=gca;
ax.FontSize=22;
uistack(p,'bottom') %moves the linear regression line to the bottom, so that markers are above it

%add breakwater marker in green
%     filename=['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Breakwater1_SSCReduction20_breakwaterdelay551_TLP0_TLPthick0_TLPfreq0_TLPdelay551_Erosion5.3e-09/'];
%     myVars={'elev.mat','marsh edge.mat','meansealevel.mat','AnnualAccretion.mat','SLR.mat','resiliencymetrics.mat','fetch.mat'};
%     for i=1:7
%         load([filename myVars{i}]);
%     end
%
%     marshage=find(resiliencymetrics(551:701,7)==0,1)-1; %find marsh age where VASEAmatrix third column (marsh width) equals zero, then add 1 (for marsh lifespan)
%     if isempty(marshage)==1
%         marshage=150; %here marshage is 150 (as opposed to how calculated for Fig 2a, but doesn't matter because just using these marsh age values to correlate to SLR at that year)
%     end
%     for yr=552:701 %start at year 552 because year 550 to 551 typically has large slope break as marsh suddenly begins to erode
%         widthslope(yr)=(resiliencymetrics(yr,7)-resiliencymetrics(yr-1,7))/1; %calculates the slope of the marsh width curve through time (for each year)
%         marshedgeerosion(yr)=Marsh_edge(yr)-Marsh_edge(yr-1);
%     end
%     threshold=find(widthslope(552:701)<widthslope(551)-5,1)-1+551;
%     if erosion==1 && isempty(threshold)==1 %if lowest erosion rate (first loop through) and there is no slope break...
%         threshold=marshage; %...then assume that threshold doesn't occur until marsh age (this isn't the case for any intermediate SLR scenarios)
%     else if erosion>1 && isempty(threshold)==1 && SLR(marshage+551)>thresholdsmatrix(1,erosion-1) %else for erosion greater than lowest, if threshold is empty and SLR of final marsh age is greater than SLR of threshold for low erosion scenario...
%         threshold=marshage; %...then threshold equals marshage again (this does not occur for any interm. SLRA scenarios shown right now in Fig 2, but is a failsafe to continue to show that SLR threshold rate increases for increasing erosion
%     % else if isempty(threshold)==1 %otherwise, if threshold is empty (but does not meet the conditions listed above...
%     %     continue %...then continue to next iteration of loop
%     % end
%     end
%     end
% scatter(0,threshold,100,'MarkerFaceColor','#66a61e','MarkerEdgeColor','none')

hold on; %make second line indicating time of survival
clear all

%make second scatterplot/linear fit for marshage
edgecolor='none';
color=[0.99608 0.78824 0.89412]; %for no restoration color (erosion changes)
acceleration=1;
RSLRA=0.08023;
erosion=1;
for k=[1.65e-09,2.82e-09,5.3e-09,7.8e-09,9.95e-09,9.900E-10,3.812E-09,6.300E-09,9.090E-09]
    filename=['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Breakwater0_SSCReduction0_breakwaterdelay551_TLP0_TLPthick0_TLPfreq0_TLPdelay551_Erosion' num2str(k) '/'];
    myVars={'elev.mat','marsh edge.mat','meansealevel.mat','AnnualAccretion.mat','SLR.mat','resiliencymetrics.mat','fetch.mat','bgb.mat'};
    for i=1:8
        load([filename myVars{i}]);
    end
    if k==9.95e-09
        color=[0.44706 0 0.22353];
    end

    marshage=find(resiliencymetrics(551:701,11)==0,1);
    if isempty(marshage)==1
        marshage=150;
    end
    for yr=552:701
        marshedgeerosion(yr)=Marsh_edge(yr)-Marsh_edge(yr-1);
    end

    if erosion>5
        color=[1 1 1];
        edgecolor=[0.5 0.5 0.5];
    end
    meanerosion=mean(marshedgeerosion(552:marshage+551-2));
    meanerosionmatrix(1,erosion)=meanerosion; %inserts the avg erosion into the new meanerosion matrix which holds avg erosion for each erosion scenario even as for-loop turns to next iteration
    marshagematrix(1,erosion)=marshage;
    scatter(meanerosion,marshage,130,'MarkerFaceColor',color,'MarkerEdgeColor',edgecolor,'Marker','square')
    color=color+[-0.14 -0.2 -0.15];
    erosion = erosion + 1;
    hold on
    clear AnnualAccretion elev msl resiliencymetrics SLR threshold widthslope yr meanerosion meanerosion1 fetcherosion marshedgeerosion bgb
end
% Get coefficients of a line fit through the data.
coefficients = polyfit(meanerosionmatrix,marshagematrix, 2);
% Create a new x axis with exactly 1000 points (or whatever you want).
xFit = linspace(min(meanerosionmatrix), max(meanerosionmatrix), 1000);
% Get the estimated yFit value for each of those 1000 new x locations.
yFit = polyval(coefficients , xFit);
X=[meanerosionmatrix',meanerosionmatrix'.^2];
mdl3 = fitlm(X, marshagematrix)
% Plot everything.
hold on; % Set hold on so the next plot does not blow away the one we just drew.
p1=plot(xFit, yFit, 'color','#006CA5', 'LineWidth', 7); % Plot fitted line.
hold on
uistack(p1,'bottom') %moves the linear regression line to the bottom, so that markers are above it
coefficients=table2array(mdl3.Coefficients);
str = sprintf('Year_{collapse} = %.2f(ero.)^2 %.2f(ero.) + %.2f', round(coefficients(3,1), 3), round(coefficients(2,1), 3), round(coefficients(1,1), 3));
totalstr = strjoin(string(str));
text(0.3,123,totalstr,'FontSize',22,'Color','#006CA5');
ax=gca;
ax.FontSize=22;
% ax.XColor=[0.5 0.5 0.5];
% ax.YColor=[0.5 0.5 0.5];
ylim([55 130])
xlim([0 5])
hold on

%add breakwater marker in green
%     filename=['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Breakwater1_SSCReduction20_breakwaterdelay551_TLP0_TLPthick0_TLPfreq0_TLPdelay551_Erosion5.3e-09/'];
%     myVars={'elev.mat','marsh edge.mat','meansealevel.mat','AnnualAccretion.mat','SLR.mat','resiliencymetrics.mat','fetch.mat'};
%     for i=1:7
%         load([filename myVars{i}]);
%     end
%
%     marshage=find(resiliencymetrics(551:701,3)==0,1)-1; %find marsh age where resiliencymetrics 7th column (marsh width) equals zero, then add 1 (for marsh lifespan)
%     if isempty(marshage)==1
%         marshage=150; %here marshage is 150 (as opposed to how calculated for Fig 2a, but doesn't matter because just using these marsh age values to correlate to SLR at that year)
%     end
%     for yr=552:701 %start at year 552 because year 550 to 551 typically has large slope break as marsh suddenly begins to erode
%         widthslope(yr)=(resiliencymetrics(yr,7)-resiliencymetrics(yr-1,7))/1; %calculates the slope of the marsh width curve through time (for each year)
%         marshedgeerosion(yr)=Marsh_edge(yr)-Marsh_edge(yr-1);
%     end
%     threshold=find(widthslope(552:701)<widthslope(551)-5,1)-1+551;
%     if erosion==1 && isempty(threshold)==1 %if lowest erosion rate (first loop through) and there is no slope break...
%         threshold=marshage; %...then assume that threshold doesn't occur until marsh age (this isn't the case for any intermediate SLR scenarios)
%     else if erosion>1 && isempty(threshold)==1 && SLR(marshage+551)>thresholdsmatrix(1,erosion-1) %else for erosion greater than lowest, if threshold is empty and SLR of final marsh age is greater than SLR of threshold for low erosion scenario...
%         threshold=marshage; %...then threshold equals marshage again (this does not occur for any interm. SLRA scenarios shown right now in Fig 2, but is a failsafe to continue to show that SLR threshold rate increases for increasing erosion
%     % else if isempty(threshold)==1 %otherwise, if threshold is empty (but does not meet the conditions listed above...
%     %     continue %...then continue to next iteration of loop
%     % end
%     end
%     end
% scatter(0,marshage+551,100,'MarkerFaceColor','#66a61e','MarkerEdgeColor','none')

%% Fig 5; marsh width vs VEC
%5a
k=5.3e-09;
RSLRA=0.08023;
figure
%set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.6, 0.96],'InnerPosition',[0.09,0.6,0.4,0.4],'Position',[0.09,0.6,0.4,0.4]); %makes really large figure box (takes up most of computer screen)
subplot(1,2,1)
s=subaxis(1,2,1,'SpacingHoriz',0.1,'SpacingVert',0.04); %limits the padding between subplots
%graph shoreline stabilization first
load(['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Breakwater1_SSCReduction0_breakwaterdelay551_TLP0_TLPthick0_TLPfreq0_TLPdelay551_Erosion' num2str(k) '/resiliencymetrics.mat']);
plot(0:150,resiliencymetrics(551:701,7),'lineWidth',5,'Color','#66a61e'); %#66a61e
hold on
%graph TLP next
load(['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Breakwater0_SSCReduction0_breakwaterdelay551_TLP1_TLPthick0.05_TLPfreq15_TLPdelay551_Erosion' num2str(k) '/resiliencymetrics.mat']);
plot(0:150,resiliencymetrics(551:701,7),'lineWidth',5,'Color','#7570b3'); %#7570b3
hold on
%graph TLP+SS next
load(['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Breakwater1_SSCReduction0_breakwaterdelay551_TLP1_TLPthick0.05_TLPfreq15_TLPdelay551_Erosion' num2str(k) '/resiliencymetrics.mat']);
plot(0:150,resiliencymetrics(551:701,7),'lineWidth',5,'Color','#e6ab02'); %#e6ab02
hold on
load(['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Breakwater0_SSCReduction0_breakwaterdelay551_TLP0_TLPthick0_TLPfreq0_TLPdelay551_Erosion' num2str(k) '/resiliencymetrics.mat']);
plot(0:150,resiliencymetrics(551:701,7),'lineWidth',5,'Color','#e7298a','LineStyle',':'); %#e7298a
hold on
ylim([0 525])
hold on
s1=scatter(0:15:150,523,80,'^','filled','MarkerFaceColor',[0.5,0.5,0.5],'MarkerEdgeColor','none');
% l2=legend([p4(1),p1(1),p2(1),p3(1),s1(1)],'No restoration','SS','TLP','TLP and SS','TLP Year');
% l2.Location='southwest';
% l2.EdgeColor='none';
% l2.Color=[0.97 0.97 0.97 0.97];
% l2.FontSize=12;
ax=gca;
ax.FontSize=16;
%xlabel('Model Run Year')
ylabel('Width of vegetated marsh (m)')
xlabel('Model Run Year')

text(3,515,'a','FontSize',16,'FontWeight','bold')

yticks([0 100 200 300 400 500]) %sets the ytick labels

% ----------------------------------------------------------------------------------------------------------------------------
% Fig5b w volumetric elevation capital for diff restoration strategies

% graph breakwater first
subplot(1,2,2)
s=subaxis(1,2,2,'SpacingHoriz',0.1,'SpacingVert',0.04); %limits the padding between subplots

ylim([0 525])
for i=0:15:150
    plot([i i],[525 450],'LineStyle','--')
end
hold on
filename=['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Breakwater1_SSCReduction0_breakwaterdelay551_TLP0_TLPthick0_TLPfreq0_TLPdelay551_Erosion' num2str(k) '/'];
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
p1=plot(0:150,marshvolumeaboveDmax(551:701),'LineWidth',5,'Color','#66a61e');%#66a61e
hold on

%graph TLP next
filename=['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Breakwater0_SSCReduction0_breakwaterdelay551_TLP1_TLPthick0.05_TLPfreq15_TLPdelay551_Erosion' num2str(k) '/'];
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
p2=plot(0:150,marshvolumeaboveDmax(551:701),'LineWidth',5,'Color','#7570b3'); %#7570b3
hold on

%graph TLP and breakwater next
filename=['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Breakwater1_SSCReduction0_breakwaterdelay551_TLP1_TLPthick0.05_TLPfreq15_TLPdelay551_Erosion' num2str(k) '/'];
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
p3=plot(0:150,marshvolumeaboveDmax(551:701),'LineWidth',5,'Color','#e6ab02'); %#e6ab02
hold on

%finally graph no restoration last as dashed line on top
filename=['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Breakwater0_SSCReduction0_breakwaterdelay551_TLP0_TLPthick0_TLPfreq0_TLPdelay551_Erosion' num2str(k) '/'];
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
yr=0:150;
p4=plot(yr,marshvolumeaboveDmax(551:701),'LineWidth',5,'Color','#e7298a','LineStyle',':');
hold on
%end
xlim([0 150])
ylim([0 275])
ylabel('Volumetric elevation capital (m^3)')
%xlabel('Model run year')
hold on
s1=scatter(0:15:150,273,80,'^','filled','MarkerFaceColor',[0.5,0.5,0.5],'MarkerEdgeColor','none');
l=legend([p4(1),p1(1),p2(1),p3(1),s1(1)],'No restoration','Bw','TLP','TLP and Bw','TLP Year');
l.Location='southwest';
l.Orientation='horizontal';
l.EdgeColor='k';
%l.Color=[0.97 0.97 0.97 0.97];
l.FontSize=16;
ax=gca;
ax.FontSize=16;
xlabel('Model Run Year')
text(3,268,'b','FontSize',16,'FontWeight','bold')


%% Fig 6 of accretion rates/components for all 4 restoration strategies
% *keep org and min separate until the very end
%figure out percentage of Fb_org and Fb_min made up of by eroded material (saved in fluxes matrix)
k=5.3e-09;
RSLRA=0.08023;
figure
%legend('Marsh accretion','External & mudflat sediment','Recycled eroded sediment')

filename=['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Breakwater0_SSCReduction0_breakwaterdelay551_TLP0_TLPthick0_TLPfreq0_TLPdelay551_Erosion' num2str(k) '/'];
myVars={'fluxes.mat','marsh edge.mat','meansealevel.mat','AnnualAccretion.mat','SLR.mat','VASEAmatrix.mat','organic deposition.mat','mineral deposition.mat','savedcompaction.mat','bgb.mat'}; % org. dep and min. dep in m; savedcompaction in m2.
for i=1:10
    load([filename myVars{i}]);
end

erodedminpercent=(fluxes(1,:)./fluxes(7,:)).'; %[kg/yr] divided by [kg/yr] = percent
erodedorgpercent=(fluxes(2,:)./fluxes(8,:)).'; %[kg/yr] divided by [kg/yr] = percent

for yr=551:701
    if fluxes(7,yr)<0 %if more sediment exported from bay than imported...
        erodedminpercent(yr)=erodedminpercent(yr-1); %then assume same % eroded material from year before...
        if isnan(erodedminpercent(yr)==1)
            erodedminpercent(yr)=0;
        end
    end
    if erodedminpercent(yr)>1
        erodedminpercent(yr)=1;
    end
    if fluxes(8,yr)<0  %if more sediment exported from bay than imported...
        erodedorgpercent(yr)=erodedorgpercent(yr-1); %then assume same % eroded material from year before...
        if isnan(erodedorgpercent(yr)==1)
            erodedorgpercent(yr)=0;
        end
    end
    if erodedorgpercent(yr)>1
        erodedorgpercent(yr)=1;
    end
end

%multiply the percentage of Fe_org and Fe_min that makes up organic & mineral deposition

erodedmindeposition=zeros(701,6152);
erodedorgdeposition=zeros(701,6152);
for yr=551:701
    erodedmindeposition(yr,:)=erodedminpercent(yr)*mineral_dep(yr,:); %[g/m]
    erodedorgdeposition(yr,:)=erodedorgpercent(yr)*organic_dep_alloch(yr,:); %[g/m]
end

%translate [g/m] into [m] by dividing by bulk density of org and mineral matter
rhos=2000*1000;%bulk density of mineral matter converted from kg/m3 into [g/m3]
rhoo=85.0*1000;%bulk density of organic matter converted from kg/m3 into [g/m3]
volume_erodedmindeposition=erodedmindeposition/rhos; %[m2]
volume_erodedorgdeposition=erodedorgdeposition/rhoo; %[m2]

%average across marsh
avgvolume_erodedmindeposition=zeros(701,1);
avgvolume_erodedorgdeposition=zeros(701,1);

for yr=551:701
    avgvolume_erodedmindeposition(yr)=mean(volume_erodedmindeposition(yr,Marsh_edge(yr):5501)); %[m2]
    avgvolume_erodedorgdeposition(yr)=mean(volume_erodedorgdeposition(yr,Marsh_edge(yr):5501)); %[m2]
end

%save decomposition of OM every year as average volume per meter

avgdecomp=zeros(701,1);
for yr=551:701
    avgdecomp(yr)=mean(savedcompaction(yr,Marsh_edge(yr):5501)); %[m2]
end

%save average belowground biomass deposition every year and redefine in [m]

avgbgbdeposition=zeros(701,1);
for yr=551:701
    avgbgbdeposition(yr)=mean(bgb(yr,Marsh_edge(yr):5501)); %[g/m]
end

volume_bgbdeposition=avgbgbdeposition/rhoo; %[m2]

%save the external/mudflat mineral & OM contribution and redefine in [m]
%calculate it as the remaining percentage after taking out erosion component
externalandmudflatmin=zeros(701,6152);
externalandmudflatorg=zeros(701,6152);
externalandmudflatmin2=zeros(701,6152);
for yr=551:701
    externalandmudflatmin(yr,:)=(1-erodedminpercent(yr))*mineral_dep(yr,:); %[g/m]
    externalandmudflatmin2(yr,:)=mineral_dep(yr,:)-erodedmindeposition(yr,:);
    externalandmudflatorg(yr,:)=(1-erodedorgpercent(yr))*organic_dep_alloch(yr,:); %[g/m]
end

%translate [g/m] into [m] by dividing by bulk density of org and mineral matter
rhos=2000*1000;%bulk density of mineral matter [g/m3]
rhoo=85.0*1000;%bulk density of organic matter [g/m3]
volume_externalandmudflatmin=externalandmudflatmin/rhos; %[m2]
volume_externalandmudflatorg=externalandmudflatorg/rhoo; %[m2]

%average across marsh
avgvolume_externalandmudflatmin=zeros(701,1);
avgvolume_externalandmudflatorg=zeros(701,1);

for yr=551:701
    avgvolume_externalandmudflatmin(yr)=mean(volume_externalandmudflatmin(yr,Marsh_edge(yr):5501)); %[m2]
    avgvolume_externalandmudflatorg(yr)=mean(volume_externalandmudflatorg(yr,Marsh_edge(yr):5501)); %[m2]
end

%add all accretion components and subtract by decomposition for organic
totalminaccretion=avgvolume_externalandmudflatmin+avgvolume_erodedmindeposition; %[m2]
totalorgaccretion=avgvolume_externalandmudflatorg+avgvolume_erodedorgdeposition+volume_bgbdeposition; %[m2]
totaleroded=avgvolume_erodedmindeposition+avgvolume_erodedorgdeposition; %[m2]
totalexternal=avgvolume_externalandmudflatmin+avgvolume_externalandmudflatorg; %[m2]

totalaccretion=totalminaccretion+totalorgaccretion-avgdecomp; %[m2]

%check how closely match the AnnualAccretion matrix
AvgAnnualAccretion=zeros(701,1);
for yr=551:701
    AvgAnnualAccretion(yr)=mean(AnnualAccretion(yr,Marsh_edge(yr):5501));
    accretiontoSLR(yr)=AvgAnnualAccretion(yr)/SLR(yr);
end
comparison=totalaccretion./AvgAnnualAccretion;

%make figure
accretioncomponents=zeros(701,4);
accretioncomponents(:,1)=totalexternal*1000;
accretioncomponents(:,2)=totaleroded*1000;
accretioncomponents(:,3)=(volume_bgbdeposition-avgdecomp)*1000;
accretioncomponents(:,4)=volume_bgbdeposition*1000;
accretioncomponents(:,5)=-avgdecomp*1000;

marshage=find(VASEAmatrix(551:701,3)==0,1)-2;
if isempty(marshage)==1
    marshage=150;
end

switchyr=find(accretioncomponents(551:701,3)<=0,2);
switchyr=switchyr(find(switchyr>30));
if isempty(switchyr)==1
    switchyr=marshage;
end

subplot(2,2,1)
s=subaxis(2,2,1,'SpacingHoriz',0.04,'SpacingVert',0.04); %limits the padding between subplots

f=area(switchyr-3:switchyr,accretioncomponents(switchyr+548:switchyr+551,3));
f(1).FaceColor='#ece2f0';
f(1).EdgeColor='none';
hold on
a=area(0:switchyr-2,accretioncomponents(551:switchyr+549,1:3));
ylabel({'Average accretion rate','(mm yr^-^1)'})
%xlabel('Model Run Year')

a(1).FaceColor='#838996';
a(1).EdgeColor='none';
a(2).FaceColor='#4c516d';
a(2).EdgeColor='none';
a(3).FaceColor='#ece2f0';
a(3).EdgeColor='none';
hold on
e=area(switchyr-2:marshage,accretioncomponents(switchyr+549:marshage+551,1:2));
e(1).FaceColor='#838996';
e(1).EdgeColor='none';
e(2).FaceColor='#4c516d';
e(2).EdgeColor='none';
hold on
d=area(switchyr:marshage,accretioncomponents(switchyr+551:marshage+551,3));
d(1).FaceColor='#ece2f0';
d(1).EdgeColor='none';

hold on
plot(0:marshage,totalaccretion(551:marshage+551)*1000,'Color','#e7298a','LineWidth',4)
hold on
plot(0:marshage,SLR(551:marshage+551)*1000,'Color','k','LineWidth',4,'LineStyle',':')
hold on
p=plot(0:150,zeros(size(0:150)),'Color','k','LineWidth',0.2);
uistack(p,'top')
r=rectangle('Position',[marshage -1 150-marshage 2],'FaceColor','w','EdgeColor','none');
uistack(r,'top');
xlim([0 150])
ylim([-15 20])
xline(marshage,'-r','Marsh collapse','FontSize',16,'LineWidth',3,'LabelVerticalAlignment','middle')

ax=gca;
ax.FontSize=16;

text(3,19,'a','FontSize',16,'FontWeight','bold')
text(8,18,'No restoration','FontSize',16,'Color','#e7298a')
text(75,10,'SLR','FontSize',16,'Rotation',18)
legend('Net belowground biomass','External & mudflat sediment','Recycled eroded sediment','Location','northoutside','Orientation','Horizontal','FontSize',16,'EdgeColor','w')

text(rand,rand,'\fontsize{40}\color[rgb]{0.906,0.161,0.541}     —\fontsize{40}\color[rgb]{0.400,0.651,0.118}—\fontsize{40}\color[rgb]{0.459,0.439,0.702}—\fontsize{40}\color[rgb]{0.902,0.671,0.008}—\fontsize{20}\color[rgb]{0,0,0}^{ Total accretion}\color[rgb]{0.5,0.5,0.5}^{    ▲}\fontsize{20}\color[rgb]{0,0,0}^{Year in which marsh receives TLP                    }','BackgroundColor','w','EdgeColor','k','LineWidth',0.5)
plottools

% ----------------------------------------------------------------------------------------------------------------------------
% fig 6b; create accretion components graph
% *keep org and min separate until the very end
%figure out percentage of Fb_org and Fb_min made up of by eroded material (saved in fluxes matrix)

filename=['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Breakwater1_SSCReduction0_breakwaterdelay551_TLP0_TLPthick0_TLPfreq0_TLPdelay551_Erosion' num2str(k) '/'];
myVars={'fluxes.mat','marsh edge.mat','meansealevel.mat','AnnualAccretion.mat','SLR.mat','VASEAmatrix.mat','organic deposition.mat','mineral deposition.mat','savedcompaction.mat','bgb.mat'}; % org. dep and min. dep in m; savedcompaction in m2.
for i=1:10
    load([filename myVars{i}]);
end

erodedminpercent=(fluxes(1,:)./fluxes(7,:)).'; %[kg/yr] divided by [kg/yr] = percent
erodedorgpercent=(fluxes(2,:)./fluxes(8,:)).'; %[kg/yr] divided by [kg/yr] = percent

for yr=551:701
    if fluxes(7,yr)<0 %if more sediment exported from bay than imported...
        erodedminpercent(yr)=erodedminpercent(yr-1); %then assume same % eroded material from year before...
        if isnan(erodedminpercent(yr)==1)
            erodedminpercent(yr)=0;
        end
    end
    if erodedminpercent(yr)>1
        erodedminpercent(yr)=1;
    end
    if fluxes(8,yr)<0  %if more sediment exported from bay than imported...
        erodedorgpercent(yr)=erodedorgpercent(yr-1); %then assume same % eroded material from year before...
        if isnan(erodedorgpercent(yr)==1)
            erodedorgpercent(yr)=0;
        end
    end
    if erodedorgpercent(yr)>1
        erodedorgpercent(yr)=1;
    end
end

%multiply the percentage of Fe_org and Fe_min that makes up organic & mineral deposition

erodedmindeposition=zeros(701,6152);
erodedorgdeposition=zeros(701,6152);
for yr=551:701
    erodedmindeposition(yr,:)=erodedminpercent(yr)*mineral_dep(yr,:); %[g/m]
    erodedorgdeposition(yr,:)=erodedorgpercent(yr)*organic_dep_alloch(yr,:); %[g/m]
end

%translate [g/m] into [m] by dividing by bulk density of org and mineral matter
rhos=2000*1000;%bulk density of mineral matter converted from kg/m3 into [g/m3]
rhoo=85.0*1000;%bulk density of organic matter converted from kg/m3 into [g/m3]
volume_erodedmindeposition=erodedmindeposition/rhos; %[m2]
volume_erodedorgdeposition=erodedorgdeposition/rhoo; %[m2]

%average across marsh
avgvolume_erodedmindeposition=zeros(701,1);
avgvolume_erodedorgdeposition=zeros(701,1);

for yr=551:701
    avgvolume_erodedmindeposition(yr)=mean(volume_erodedmindeposition(yr,Marsh_edge(yr):5501)); %[m2]
    avgvolume_erodedorgdeposition(yr)=mean(volume_erodedorgdeposition(yr,Marsh_edge(yr):5501)); %[m2]
end

%save decomposition of OM every year as average volume per meter

avgdecomp=zeros(701,1);
for yr=551:701
    avgdecomp(yr)=mean(savedcompaction(yr,Marsh_edge(yr):5501)); %[m2]
end

%save average belowground biomass deposition every year and redefine in [m]

avgbgbdeposition=zeros(701,1);
for yr=551:701
    avgbgbdeposition(yr)=mean(bgb(yr,Marsh_edge(yr):5501)); %[g/m]
end

volume_bgbdeposition=avgbgbdeposition/rhoo; %[m2]

%save the external/mudflat mineral & OM contribution and redefine in [m]
%calculate it as the remaining percentage after taking out erosion component
externalandmudflatmin=zeros(701,6152);
externalandmudflatorg=zeros(701,6152);
externalandmudflatmin2=zeros(701,6152);
for yr=551:701
    externalandmudflatmin(yr,:)=(1-erodedminpercent(yr))*mineral_dep(yr,:); %[g/m]
    externalandmudflatmin2(yr,:)=mineral_dep(yr,:)-erodedmindeposition(yr,:);
    externalandmudflatorg(yr,:)=(1-erodedorgpercent(yr))*organic_dep_alloch(yr,:); %[g/m]
end

%translate [g/m] into [m] by dividing by bulk density of org and mineral matter
rhos=2000*1000;%bulk density of mineral matter [g/m3]
rhoo=85.0*1000;%bulk density of organic matter [g/m3]
volume_externalandmudflatmin=externalandmudflatmin/rhos; %[m2]
volume_externalandmudflatorg=externalandmudflatorg/rhoo; %[m2]

%average across marsh
avgvolume_externalandmudflatmin=zeros(701,1);
avgvolume_externalandmudflatorg=zeros(701,1);

for yr=551:701
    avgvolume_externalandmudflatmin(yr)=mean(volume_externalandmudflatmin(yr,Marsh_edge(yr):5501)); %[m2]
    avgvolume_externalandmudflatorg(yr)=mean(volume_externalandmudflatorg(yr,Marsh_edge(yr):5501)); %[m2]
end

%add all accretion components and subtract by decomposition for organic
totalminaccretion=avgvolume_externalandmudflatmin+avgvolume_erodedmindeposition; %[m2]
totalorgaccretion=avgvolume_externalandmudflatorg+avgvolume_erodedorgdeposition+volume_bgbdeposition; %[m2]
totaleroded=avgvolume_erodedmindeposition+avgvolume_erodedorgdeposition; %[m2]
totalexternal=avgvolume_externalandmudflatmin+avgvolume_externalandmudflatorg; %[m2]

totalaccretion=totalminaccretion+totalorgaccretion-avgdecomp; %[m2]

%check how closely match the AnnualAccretion matrix
AvgAnnualAccretion=zeros(701,1);
for yr=551:701
    AvgAnnualAccretion(yr)=mean(AnnualAccretion(yr,Marsh_edge(yr):5501));
    accretiontoSLR(yr)=AvgAnnualAccretion(yr)/SLR(yr);
end
comparison=totalaccretion./AvgAnnualAccretion;

%make figure
accretioncomponents=zeros(701,4);
accretioncomponents(:,1)=totalexternal*1000;
accretioncomponents(:,2)=totaleroded*1000;
accretioncomponents(:,3)=(volume_bgbdeposition-avgdecomp)*1000;
accretioncomponents(:,4)=volume_bgbdeposition*1000;
accretioncomponents(:,5)=-avgdecomp*1000;

marshage=find(VASEAmatrix(551:701,3)==0,1)-2;
if isempty(marshage)==1
    marshage=150;
end

switchyr=find(accretioncomponents(551:701,3)<=0,2);
switchyr=switchyr(find(switchyr>30));
if isempty(switchyr)==1
    switchyr=marshage;
end

subplot(2,2,2)
s=subaxis(2,2,2,'SpacingHoriz',0.04,'SpacingVert',0.04); %limits the padding between subplots

f=area(switchyr-3:switchyr,accretioncomponents(switchyr+548:switchyr+551,3));
f(1).FaceColor='#ece2f0';
f(1).EdgeColor='none';
hold on
a=area(0:switchyr-2,accretioncomponents(551:switchyr+549,1:3));

a(1).FaceColor='#838996';
a(1).EdgeColor='none';
a(2).FaceColor='#4c516d';
a(2).EdgeColor='none';
a(3).FaceColor='#ece2f0';
a(3).EdgeColor='none';
hold on
e=area(switchyr-2:marshage,accretioncomponents(switchyr+549:marshage+551,1:2));
e(1).FaceColor='#838996';
e(1).EdgeColor='none';
e(2).FaceColor='#4c516d';
e(2).EdgeColor='none';
hold on
d=area(switchyr:marshage,accretioncomponents(switchyr+551:marshage+551,3));
d(1).FaceColor='#ece2f0';
d(1).EdgeColor='none';

hold on
plot(0:marshage,totalaccretion(551:marshage+551)*1000,'Color','#66a61e','LineWidth',4)
hold on
plot(0:marshage,SLR(551:marshage+551)*1000,'Color','k','LineWidth',4,'LineStyle',':')
hold on
p=plot(0:150,zeros(size(0:150)),'Color','k','LineWidth',0.2);
uistack(p,'top')
r=rectangle('Position',[marshage -1 150-marshage 2],'FaceColor','w','EdgeColor','none');
uistack(r,'top');
xlim([0 150])
ylim([-15 20])

ax=gca;
ax.FontSize=16;

text(3,19,'b','FontSize',16,'FontWeight','bold')
text(8,18,'Breakwater','FontSize',16,'Color','#66a61e')
text(75,9.9,'SLR','FontSize',16,'Rotation',18)

xline(marshage,'-r','Marsh collapse','FontSize',16,'LineWidth',3,'LabelVerticalAlignment','middle')

%-------------------------------------------------------------------------
% fig 6c
subplot(2,2,3)
s1=subaxis(2,2,3,'SpacingHoriz',0.04,'SpacingVert',0.04); %limits the padding between subplots
k=5.3e-09;
RSLRA=0.08023;
filename=['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Breakwater0_SSCReduction0_breakwaterdelay551_TLP1_TLPthick0.05_TLPfreq15_TLPdelay551_Erosion' num2str(k) '/'];
myVars={'fluxes.mat','marsh edge.mat','meansealevel.mat','AnnualAccretion.mat','SLR.mat','VASEAmatrix.mat','organic deposition.mat','mineral deposition.mat','savedcompaction.mat','bgb.mat'}; % org. dep and min. dep in m; savedcompaction in m2.
for i=1:10
    load([filename myVars{i}]);
end

erodedminpercent=(fluxes(1,:)./fluxes(7,:)).'; %[kg/yr] divided by [kg/yr] = percent
erodedorgpercent=(fluxes(2,:)./fluxes(8,:)).'; %[kg/yr] divided by [kg/yr] = percent

for yr=551:701
    if fluxes(7,yr)<0 %if more sediment exported from bay than imported...
        erodedminpercent(yr)=erodedminpercent(yr-1); %then assume same % eroded material from year before...
        if isnan(erodedminpercent(yr)==1)
            erodedminpercent(yr)=0;
        end
    end
    if erodedminpercent(yr)>1
        erodedminpercent(yr)=1;
    end
    if fluxes(8,yr)<0  %if more sediment exported from bay than imported...
        erodedorgpercent(yr)=erodedorgpercent(yr-1); %then assume same % eroded material from year before...
        if isnan(erodedorgpercent(yr)==1)
            erodedorgpercent(yr)=0;
        end
    end
    if erodedorgpercent(yr)>1
        erodedorgpercent(yr)=1;
    end
end

%multiply the percentage of Fe_org and Fe_min that makes up organic & mineral deposition

erodedmindeposition=zeros(701,6152);
erodedorgdeposition=zeros(701,6152);
for yr=551:701
    erodedmindeposition(yr,:)=erodedminpercent(yr)*mineral_dep(yr,:); %[g/m]
    erodedorgdeposition(yr,:)=erodedorgpercent(yr)*organic_dep_alloch(yr,:); %[g/m]
end


%translate [g/m] into [m] by dividing by bulk density of org and mineral matter
rhos=2000*1000;%bulk density of mineral matter converted from kg/m3 into [g/m3]
rhoo=85.0*1000;%bulk density of organic matter converted from kg/m3 into [g/m3]
volume_erodedmindeposition=erodedmindeposition/rhos; %[m2]
volume_erodedorgdeposition=erodedorgdeposition/rhoo; %[m2]

%average across marsh
avgvolume_erodedmindeposition=zeros(701,1);
avgvolume_erodedorgdeposition=zeros(701,1);

for yr=551:701
    avgvolume_erodedmindeposition(yr)=mean(volume_erodedmindeposition(yr,Marsh_edge(yr):5501)); %[m2]
    avgvolume_erodedorgdeposition(yr)=mean(volume_erodedorgdeposition(yr,Marsh_edge(yr):5501)); %[m2]
end

avgvolume_erodedmindeposition(551:15:701)=NaN;
[avgvolume_erodedmindeposition,TF]=fillmissing(avgvolume_erodedmindeposition,'linear');
avgvolume_erodedorgdeposition(551:15:701)=NaN;
[avgvolume_erodedorgdeposition,TF]=fillmissing(avgvolume_erodedorgdeposition,'linear');

%save decomposition of OM every year as average volume per meter

avgdecomp=zeros(701,1);
for yr=551:701
    avgdecomp(yr)=mean(savedcompaction(yr,Marsh_edge(yr):5501)); %[m2]
end

%save average belowground biomass deposition every year and redefine in [m]

avgbgbdeposition=zeros(701,1);
for yr=551:701
    avgbgbdeposition(yr)=mean(bgb(yr,Marsh_edge(yr):5501)); %[g/m]
end

volume_bgbdeposition=avgbgbdeposition/rhoo; %[m2]

%save the external/mudflat mineral & OM contribution and redefine in [m]
%calculate it as the remaining percentage after taking out erosion component
externalandmudflatmin=zeros(701,6152);
externalandmudflatorg=zeros(701,6152);
externalandmudflatmin2=zeros(701,6152);
for yr=551:701
    externalandmudflatmin(yr,:)=(1-erodedminpercent(yr))*mineral_dep(yr,:); %[g/m]
    externalandmudflatmin2(yr,:)=mineral_dep(yr,:)-erodedmindeposition(yr,:);
    externalandmudflatorg(yr,:)=(1-erodedorgpercent(yr))*organic_dep_alloch(yr,:); %[g/m]
end

%translate [g/m] into [m] by dividing by bulk density of org and mineral matter
rhos=2000*1000;%bulk density of mineral matter [g/m3]
rhoo=85.0*1000;%bulk density of organic matter [g/m3]
volume_externalandmudflatmin=externalandmudflatmin/rhos; %[m2]
volume_externalandmudflatorg=externalandmudflatorg/rhoo; %[m2]

%average across marsh
avgvolume_externalandmudflatmin=zeros(701,1);
avgvolume_externalandmudflatorg=zeros(701,1);

for yr=551:701
    avgvolume_externalandmudflatmin(yr)=mean(volume_externalandmudflatmin(yr,Marsh_edge(yr):5501)); %[m2]
    avgvolume_externalandmudflatorg(yr)=mean(volume_externalandmudflatorg(yr,Marsh_edge(yr):5501)); %[m2]
end

avgvolume_externalandmudflatmin(551:15:701)=NaN;
[avgvolume_externalandmudflatmin,TF]=fillmissing(avgvolume_externalandmudflatmin,'linear');
avgvolume_externalandmudflatorg(551:15:701)=NaN;
[avgvolume_externalandmudflatorg,TF]=fillmissing(avgvolume_externalandmudflatorg,'linear');

%add all accretion components and subtract by decomposition for organic
totalminaccretion=avgvolume_externalandmudflatmin+avgvolume_erodedmindeposition; %[m2]
totalorgaccretion=avgvolume_externalandmudflatorg+avgvolume_erodedorgdeposition+volume_bgbdeposition; %[m2]
totaleroded=avgvolume_erodedmindeposition+avgvolume_erodedorgdeposition; %[m2]
totalexternal=avgvolume_externalandmudflatmin+avgvolume_externalandmudflatorg; %[m2]

totalaccretion=totalminaccretion+totalorgaccretion-avgdecomp; %[m2]

%check how closely match the AnnualAccretion matrix
AvgAnnualAccretion=zeros(701,1);
for yr=551:701
    AvgAnnualAccretion(yr)=mean(AnnualAccretion(yr,Marsh_edge(yr):5501));
    accretiontoSLR(yr)=AvgAnnualAccretion(yr)/SLR(yr);
end
comparison=totalaccretion./AvgAnnualAccretion;

%make figure
accretioncomponents=zeros(701,4);
accretioncomponents(:,1)=totalexternal*1000;
accretioncomponents(:,2)=totaleroded*1000;
accretioncomponents(:,3)=(volume_bgbdeposition-avgdecomp)*1000;
accretioncomponents(:,4)=volume_bgbdeposition*1000;
accretioncomponents(:,5)=-avgdecomp*1000;

marshage=find(VASEAmatrix(551:701,3)==0,1)-2; %subtract 2 b/c of way matlab counts using find function
if isempty(marshage)==1
    marshage=150;
end

switchyr=find(accretioncomponents(551:701,3)<=0,2);
switchyr=switchyr(find(switchyr>30));
if isempty(switchyr)==1
    a=area(0:marshage,accretioncomponents(551:marshage+551,1:3));
    a(1).FaceColor='#838996'; %#838996
    a(1).EdgeColor='none';
    a(2).FaceColor='#4c516d'; %#4c516d
    a(2).EdgeColor='none';
    a(3).FaceColor='#ece2f0'; %#ece2f0
    a(3).EdgeColor='none';
else
    %switchyr=marshage;
    f=area(switchyr-3:switchyr,accretioncomponents(switchyr+548:switchyr+551,3));
    f(1).FaceColor='#ece2f0'; %#ece2f0
    f(1).EdgeColor='none';
    hold on
    a=area(0:switchyr-2,accretioncomponents(551:switchyr+549,1:3));

    a(1).FaceColor='#838996'; %#838996
    a(1).EdgeColor='none';
    a(2).FaceColor='#4c516d'; %#4c516d
    a(2).EdgeColor='none';
    a(3).FaceColor='#ece2f0'; %#ece2f0
    a(3).EdgeColor='none';
    hold on
    e=area(switchyr-2:marshage,accretioncomponents(switchyr+549:marshage+551,1:2));
    e(1).FaceColor='#838996'; %#838996
    e(1).EdgeColor='none';
    e(2).FaceColor='#4c516d'; %#4c516d
    e(2).EdgeColor='none';
    hold on
    d=area(switchyr:marshage,accretioncomponents(switchyr+551:marshage+551,3));
    d(1).FaceColor='#ece2f0'; %#ece2f0
    d(1).EdgeColor='none';
end

ylabel({'Average accretion rate','(mm yr^-^1)'})
xlabel('Model Run Year')

hold on
plot(0:marshage,totalaccretion(551:marshage+551)*1000,'Color','#7570b3','LineWidth',4) %#7570b3
hold on
plot(0:marshage,SLR(551:marshage+551)*1000,'Color','k','LineWidth',4,'LineStyle',':')
hold on

scatter(0:15:marshage,19,80,'^','filled','MarkerFaceColor',[0.5,0.5,0.5],'MarkerEdgeColor','none')
%yline(0)
p=plot(0:150,zeros(size(0:150)),'Color','k','LineWidth',0.2);
uistack(p,'top')
r=rectangle('Position',[marshage -1 150-marshage 2],'FaceColor','w','EdgeColor','none');
uistack(r,'top');
xlim([0 150])
ylim([-15 20])
ax=gca;
ax.FontSize=16;

text(8,18,'Thin-layer Placement','FontSize',16,'Color','#7570b3')
text(110,13,'SLR','FontSize',16,'Rotation',18)

xline(marshage,'-r','Marsh collapse','FontSize',16,'LineWidth',3,'LabelVerticalAlignment','middle')
text(3,19,'c','FontSize',16,'FontWeight','bold')

% ----------------------------------------------------------------------------------------------------------------------------
% fig 6d
subplot(2,2,4)
s2=subaxis(2,2,4,'SpacingHoriz',0.04,'SpacingVert',0.04); %limits the padding between subplots
filename=['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Breakwater1_SSCReduction0_breakwaterdelay551_TLP1_TLPthick0.05_TLPfreq15_TLPdelay551_Erosion' num2str(k) '/'];
myVars={'fluxes.mat','marsh edge.mat','meansealevel.mat','AnnualAccretion.mat','SLR.mat','VASEAmatrix.mat','organic deposition.mat','mineral deposition.mat','savedcompaction.mat','bgb.mat'}; % org. dep and min. dep in m; savedcompaction in m2.
for i=1:10
    load([filename myVars{i}]);
end

erodedminpercent=(fluxes(1,:)./fluxes(7,:)).'; %[kg/yr] divided by [kg/yr] = percent
erodedorgpercent=(fluxes(2,:)./fluxes(8,:)).'; %[kg/yr] divided by [kg/yr] = percent

for yr=551:701
    if fluxes(7,yr)<0 %if more sediment exported from bay than imported...
        erodedminpercent(yr)=erodedminpercent(yr-1); %then assume same % eroded material from year before...
        if isnan(erodedminpercent(yr)==1)
            erodedminpercent(yr)=0;
        end
    end
    if erodedminpercent(yr)>1
        erodedminpercent(yr)=1;
    end
    if fluxes(8,yr)<0  %if more sediment exported from bay than imported...
        erodedorgpercent(yr)=erodedorgpercent(yr-1); %then assume same % eroded material from year before...
        if isnan(erodedorgpercent(yr)==1)
            erodedorgpercent(yr)=0;
        end
    end
    if erodedorgpercent(yr)>1
        erodedorgpercent(yr)=1;
    end
end

%multiply the percentage of Fe_org and Fe_min that makes up organic & mineral deposition

erodedmindeposition=zeros(701,6152);
erodedorgdeposition=zeros(701,6152);
for yr=551:701
    erodedmindeposition(yr,:)=erodedminpercent(yr)*mineral_dep(yr,:); %[g/m]
    erodedorgdeposition(yr,:)=erodedorgpercent(yr)*organic_dep_alloch(yr,:); %[g/m]
end

%translate [g/m] into [m] by dividing by bulk density of org and mineral matter
rhos=2000*1000;%bulk density of mineral matter converted from kg/m3 into [g/m3]
rhoo=85.0*1000;%bulk density of organic matter converted from kg/m3 into [g/m3]
volume_erodedmindeposition=erodedmindeposition/rhos; %[m2]
volume_erodedorgdeposition=erodedorgdeposition/rhoo; %[m2]

%average across marsh
avgvolume_erodedmindeposition=zeros(701,1);
avgvolume_erodedorgdeposition=zeros(701,1);

for yr=551:701
    avgvolume_erodedmindeposition(yr)=mean(volume_erodedmindeposition(yr,Marsh_edge(yr):5501)); %[m2]
    avgvolume_erodedorgdeposition(yr)=mean(volume_erodedorgdeposition(yr,Marsh_edge(yr):5501)); %[m2]
end

avgvolume_erodedmindeposition(551:15:701)=NaN;
[avgvolume_erodedmindeposition,TF]=fillmissing(avgvolume_erodedmindeposition,'linear');
avgvolume_erodedorgdeposition(551:15:701)=NaN;
[avgvolume_erodedorgdeposition,TF]=fillmissing(avgvolume_erodedorgdeposition,'linear');

%save decomposition of OM every year as average volume per meter

avgdecomp=zeros(701,1);
for yr=551:701
    avgdecomp(yr)=mean(savedcompaction(yr,Marsh_edge(yr):5501)); %[m2]
end

%save average belowground biomass deposition every year and redefine in [m]

avgbgbdeposition=zeros(701,1);
for yr=551:701
    avgbgbdeposition(yr)=mean(bgb(yr,Marsh_edge(yr):5501)); %[g/m]
end

volume_bgbdeposition=avgbgbdeposition/rhoo; %[m2]

%save the external/mudflat mineral & OM contribution and redefine in [m]
%calculate it as the remaining percentage after taking out erosion component
externalandmudflatmin=zeros(701,6152);
externalandmudflatorg=zeros(701,6152);
externalandmudflatmin2=zeros(701,6152);
for yr=551:701
    externalandmudflatmin(yr,:)=(1-erodedminpercent(yr))*mineral_dep(yr,:); %[g/m]
    externalandmudflatmin2(yr,:)=mineral_dep(yr,:)-erodedmindeposition(yr,:);
    externalandmudflatorg(yr,:)=(1-erodedorgpercent(yr))*organic_dep_alloch(yr,:); %[g/m]
end

%translate [g/m] into [m] by dividing by bulk density of org and mineral matter
rhos=2000*1000;%bulk density of mineral matter [g/m3]
rhoo=85.0*1000;%bulk density of organic matter [g/m3]
volume_externalandmudflatmin=externalandmudflatmin/rhos; %[m2]
volume_externalandmudflatorg=externalandmudflatorg/rhoo; %[m2]

%average across marsh
avgvolume_externalandmudflatmin=zeros(701,1);
avgvolume_externalandmudflatorg=zeros(701,1);

for yr=551:701
    avgvolume_externalandmudflatmin(yr)=mean(volume_externalandmudflatmin(yr,Marsh_edge(yr):5501)); %[m2]
    avgvolume_externalandmudflatorg(yr)=mean(volume_externalandmudflatorg(yr,Marsh_edge(yr):5501)); %[m2]
end

avgvolume_externalandmudflatmin(551:15:701)=NaN;
[avgvolume_externalandmudflatmin,TF]=fillmissing(avgvolume_externalandmudflatmin,'linear');
avgvolume_externalandmudflatorg(551:15:701)=NaN;
[avgvolume_externalandmudflatorg,TF]=fillmissing(avgvolume_externalandmudflatorg,'linear');

%add all accretion components and subtract by decomposition for organic
totalminaccretion=avgvolume_externalandmudflatmin+avgvolume_erodedmindeposition; %[m2]
totalorgaccretion=avgvolume_externalandmudflatorg+avgvolume_erodedorgdeposition+volume_bgbdeposition; %[m2]
totaleroded=avgvolume_erodedmindeposition+avgvolume_erodedorgdeposition; %[m2]
totalexternal=avgvolume_externalandmudflatmin+avgvolume_externalandmudflatorg; %[m2]

totalaccretion=totalminaccretion+totalorgaccretion-avgdecomp; %[m2]

%check how closely match the AnnualAccretion matrix
AvgAnnualAccretion=zeros(701,1);
for yr=551:701
    AvgAnnualAccretion(yr)=mean(AnnualAccretion(yr,Marsh_edge(yr):5501));
    accretiontoSLR(yr)=AvgAnnualAccretion(yr)/SLR(yr);
end
comparison=totalaccretion./AvgAnnualAccretion;

%make figure
accretioncomponents=zeros(701,4);
accretioncomponents(:,1)=totalexternal*1000;
accretioncomponents(:,2)=totaleroded*1000;
accretioncomponents(:,3)=(volume_bgbdeposition-avgdecomp)*1000;
accretioncomponents(:,4)=volume_bgbdeposition*1000;
accretioncomponents(:,5)=-avgdecomp*1000;


marshage=find(VASEAmatrix(551:701,3)==0,1)-2;
if isempty(marshage)==1
    marshage=150;
end

switchyr=find(accretioncomponents(551:701,3)<=0,2);
switchyr=switchyr(find(switchyr>30));
if isempty(switchyr)==1
    switchyr=marshage;
end

f=area(switchyr-3:switchyr,accretioncomponents(switchyr+548:switchyr+551,3));
f(1).FaceColor='#ece2f0';
f(1).EdgeColor='none';
hold on
a=area(0:switchyr-2,accretioncomponents(551:switchyr+549,1:3));
xlabel('Model Run Year')

a(1).FaceColor='#838996';
a(1).EdgeColor='none';
a(2).FaceColor='#4c516d';
a(2).EdgeColor='none';
a(3).FaceColor='#ece2f0';
a(3).EdgeColor='none';
hold on
e=area(switchyr-2:marshage,accretioncomponents(switchyr+549:marshage+551,1:2));
e(1).FaceColor='#838996';
e(1).EdgeColor='none';
e(2).FaceColor='#4c516d';
e(2).EdgeColor='none';
hold on
d=area(switchyr:marshage,accretioncomponents(switchyr+551:marshage+551,3));
d(1).FaceColor='#ece2f0';
d(1).EdgeColor='none';

hold on
plot(0:marshage,totalaccretion(551:marshage+551)*1000,'Color','#e6ab02','LineWidth',4)
hold on
plot(0:marshage,SLR(551:marshage+551)*1000,'Color','k','LineWidth',4,'LineStyle',':')
hold on

scatter(0:15:marshage,19,80,'^','filled','MarkerFaceColor',[0.5,0.5,0.5],'MarkerEdgeColor','none')

p=plot(0:150,zeros(size(0:150)),'Color','k','LineWidth',0.2);
uistack(p,'top')
r=rectangle('Position',[marshage -1 150-marshage 2],'FaceColor','w','EdgeColor','none');
uistack(r,'top');
xlim([0 150])
ylim([-15 20])
ax=gca;
ax.FontSize=16;
text(8,18,'TLP and Breakwater','FontSize',16,'Color','#e6ab02')
text(100,12,'SLR','FontSize',16,'Rotation',18)
xline(marshage,'-r','Marsh collapse','FontSize',16,'LineWidth',3,'LabelVerticalAlignment','middle')
text(3,19,'d','FontSize',16,'FontWeight','bold')

%% Fig 7: comparing contribution of erosion to avg accretion rates for TLP vs no-restoration
figure
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.6, 0.96],'InnerPosition',[0.09,0.6,0.4,0.4],'Position',[0.09,0.6,0.4,0.4]); %makes really large figure box (takes up most of computer screen)
k=7.8e-09;
RSLRA=0.08023;
% create erosion component for no restoration
%figure out percentage of Fb_org and Fb_min made up of by eroded material (saved in fluxes matrix)
filename=['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Sill0_SSCReduction0_silldelay551_TLP0_TLPthick0_TLPfreq0_TLPdelay551_Erosion' num2str(k) '/'];
myVars={'fluxes.mat','marsh edge.mat','meansealevel.mat','AnnualAccretion.mat','SLR.mat','VASEAmatrix.mat','organic deposition.mat','mineral deposition.mat','savedcompaction.mat','bgb.mat','suspended sediment.mat'}; % org. dep and min. dep in m; savedcompaction in m2.
for i=1:11
    load([filename myVars{i}]);
end

erodedminpercent=(fluxes(1,:)./fluxes(7,:)).'; %[kg/yr] divided by [kg/yr] = percent
erodedorgpercent=(fluxes(2,:)./fluxes(8,:)).'; %[kg/yr] divided by [kg/yr] = percent

MHW=zeros(701,1);
yravgelevbelowMHW=zeros(701,1);

for yr=551:701
    MHW(yr)=msl(yr)+0.5;
    yravgelevbelowMHW(yr)=MHW(yr)-VASEAmatrix(yr,2);
    if fluxes(7,yr)<0 %if more sediment exported from bay than imported...
        erodedminpercent(yr)=erodedminpercent(yr-1); %then assume same % eroded material from year before...
        if isnan(erodedminpercent(yr)==1)
            erodedminpercent(yr)=0;
        end
    end
    if erodedminpercent(yr)>1
        erodedminpercent(yr)=1;
    end
    if fluxes(8,yr)<0  %if more sediment exported from bay than imported...
        erodedorgpercent(yr)=erodedorgpercent(yr-1); %then assume same % eroded material from year before...
        if isnan(erodedorgpercent(yr)==1)
            erodedorgpercent(yr)=0;
        end
    end
    if erodedorgpercent(yr)>1
        erodedorgpercent(yr)=1;
    end
end

%multiply the percentage of Fe_org and Fe_min that makes up organic & mineral deposition

erodedmindeposition=zeros(701,width(bgb));
erodedorgdeposition=zeros(701,width(bgb));
for yr=551:701
    erodedmindeposition(yr,:)=erodedminpercent(yr)*mineral_dep(yr,:); %[g/m]
    erodedorgdeposition(yr,:)=erodedorgpercent(yr)*organic_dep_alloch(yr,:); %[g/m]
end

%translate [g/m] into [m] by dividing by bulk density of org and mineral matter
rhos=2000*1000;%bulk density of mineral matter converted from kg/m3 into [g/m3]
rhoo=85.0*1000;%bulk density of organic matter converted from kg/m3 into [g/m3]
volume_erodedmindeposition=erodedmindeposition/rhos; %[m2]
volume_erodedorgdeposition=erodedorgdeposition/rhoo; %[m2]

%average across marsh
avgvolume_erodedmindeposition=zeros(701,1);
avgvolume_erodedorgdeposition=zeros(701,1);

for yr=551:701
    avgvolume_erodedmindeposition(yr)=mean(volume_erodedmindeposition(yr,Marsh_edge(yr):5501)); %[m2]
    avgvolume_erodedorgdeposition(yr)=mean(volume_erodedorgdeposition(yr,Marsh_edge(yr):5501)); %[m2]
end

totalerodedflux=(fluxes(1,:)+fluxes(2,:))/1000;
TEFzero=find(totalerodedflux==0);
totalerodedflux(TEFzero)=NaN;

marshage=find(VASEAmatrix(551:701,3)==0,1)-1;
if isempty(marshage)==1
    marshage=150;
end

% x=0:marshage;
% y=totalerodedflux(551:marshage+551);
% [xData, yData] = prepareCurveData(x,y);
%
% % Set up fittype and options.
% ft = fittype( 'smoothingspline' );
% opts = fitoptions( 'Method', 'SmoothingSpline' );
% opts.SmoothingParam = 0.45;
%
% % Fit model to data.
% [fitresult,q]= fit( xData, yData, ft, opts );
%
% % Plot fit with data.
% yyaxis right
% p=plot(fitresult,xData,yData);
% set(p,'LineWidth',8,'Color',[0.91 0.16 0.54 0.2]);
% hold on

subplot(1,2,1)
s1=subaxis(1,2,1,'SpacingHoriz',0.04,'SpacingVert',0.1); %limits the padding between subplots
%add all accretion components and subtract by decomposition for organic
totaleroded=avgvolume_erodedmindeposition+avgvolume_erodedorgdeposition; %[m2]
%yyaxis left
plot(0:marshage,totaleroded(551:marshage+551)*1000,'Color','#e7298a','LineWidth',5)
hold on

subplot(1,2,2)
s2=subaxis(1,2,2,'SpacingHoriz',0.04,'SpacingVert',0.1); %limits the padding between subplots
scatter(yravgelevbelowMHW(551:701),totaleroded(551:701)*1000,'MarkerFaceColor','#e7298a','MarkerEdgeColor','none')
hold on

% ----------------------------------------------------------------------------------------------------------------------------
% make TLP line; create erosion component for no restoration
%figure out percentage of Fb_org and Fb_min made up of by eroded material (saved in fluxes matrix)

filename=['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Sill0_SSCReduction0_silldelay551_TLP1_TLPthick0.05_TLPfreq15_TLPdelay551_Erosion' num2str(k) '/'];
myVars={'fluxes.mat','marsh edge.mat','meansealevel.mat','AnnualAccretion.mat','SLR.mat','VASEAmatrix.mat','organic deposition.mat','mineral deposition.mat','savedcompaction.mat','bgb.mat','suspended sediment.mat'}; % org. dep and min. dep in m; savedcompaction in m2.
for i=1:11
    load([filename myVars{i}]);
end

erodedminpercent=(fluxes(1,:)./fluxes(7,:)).'; %[kg/yr] divided by [kg/yr] = percent
erodedorgpercent=(fluxes(2,:)./fluxes(8,:)).'; %[kg/yr] divided by [kg/yr] = percent

MHW=zeros(701,1);
yravgelevbelowMHW=zeros(701,1);

for yr=551:701
    MHW(yr)=msl(yr)+0.5;
    yravgelevbelowMHW(yr)=MHW(yr)-VASEAmatrix(yr,2);
    if fluxes(7,yr)<0 %if more sediment exported from bay than imported...
        erodedminpercent(yr)=erodedminpercent(yr-1); %then assume same % eroded material from year before...
        if isnan(erodedminpercent(yr)==1)
            erodedminpercent(yr)=0;
        end
    end
    if erodedminpercent(yr)>1
        erodedminpercent(yr)=1;
    end
    if fluxes(8,yr)<0  %if more sediment exported from bay than imported...
        erodedorgpercent(yr)=erodedorgpercent(yr-1); %then assume same % eroded material from year before...
        if isnan(erodedorgpercent(yr)==1)
            erodedorgpercent(yr)=0;
        end
    end
    if erodedorgpercent(yr)>1
        erodedorgpercent(yr)=1;
    end
end

%multiply the percentage of Fe_org and Fe_min that makes up organic & mineral deposition

erodedmindeposition=zeros(701,width(bgb));
erodedorgdeposition=zeros(701,width(bgb));
for yr=551:701
    erodedmindeposition(yr,:)=erodedminpercent(yr)*mineral_dep(yr,:); %[g/m]
    erodedorgdeposition(yr,:)=erodedorgpercent(yr)*organic_dep_alloch(yr,:); %[g/m]
end


%translate [g/m] into [m] by dividing by bulk density of org and mineral matter
rhos=2000*1000;%bulk density of mineral matter converted from kg/m3 into [g/m3]
rhoo=85.0*1000;%bulk density of organic matter converted from kg/m3 into [g/m3]
volume_erodedmindeposition=erodedmindeposition/rhos; %[m2]
volume_erodedorgdeposition=erodedorgdeposition/rhoo; %[m2]

%average across marsh
avgvolume_erodedmindeposition=zeros(701,1);
avgvolume_erodedorgdeposition=zeros(701,1);

for yr=551:701
    avgvolume_erodedmindeposition(yr)=mean(volume_erodedmindeposition(yr,Marsh_edge(yr):5501)); %[m2]
    avgvolume_erodedorgdeposition(yr)=mean(volume_erodedorgdeposition(yr,Marsh_edge(yr):5501)); %[m2]
end

avgvolume_erodedmindeposition(551:15:701)=NaN;
[avgvolume_erodedmindeposition,TF]=fillmissing(avgvolume_erodedmindeposition,'linear');
avgvolume_erodedorgdeposition(551:15:701)=NaN;
[avgvolume_erodedorgdeposition,TF]=fillmissing(avgvolume_erodedorgdeposition,'linear');

totalerodedflux=(fluxes(1,:)+fluxes(2,:))/1000;
TEFzero=find(totalerodedflux==0);
totalerodedflux(TEFzero)=NaN;

marshage=find(VASEAmatrix(551:701,3)==0,1)-1;
if isempty(marshage)==1
    marshage=150;
end

%add all accretion components and subtract by decomposition for organic
totaleroded=avgvolume_erodedmindeposition+avgvolume_erodedorgdeposition; %[m2]
%yyaxis left
axes(s1);
plot(0:marshage,totaleroded(551:marshage+551)*1000,'Color','#7570b3','LineWidth',5,'LineStyle','-')
ylabel({'Average annual accretion of','eroded sediment (mm yr^-^1)'},'HorizontalAlignment', 'center')
hold on
scatter(0:15:136,6.96,160,'^','filled','MarkerFaceColor',[0.5,0.5,0.5],'MarkerEdgeColor','none');


xlabel('Model run year')
ax=gca;
ax.FontSize=20;
hold on
l=legend('No restoration','TLP','TLP Year');
l.Location='west';
l.Color=[0.8 0.8 0.8];
l.EdgeColor='none';
l.BackgroundAlpha=0.3;
l.FontSize=20;
xlim([0 150])
text(3,6.8,'a','FontWeight','bold','FontSize',20)

axes(s2);
scatter(yravgelevbelowMHW(551:701),totaleroded(551:701)*1000,'MarkerFaceColor','#7570b3','MarkerEdgeColor','none')
xlim([0 1.2])
xticks([0 0.2 0.4 0.6 0.8 1 1.2]) %sets the ytick labels
box on
text(0.02,6.8,'b','FontWeight','bold','FontSize',20)

%ylabel('Average annual accretion of eroded sediment (mm yr^-^1 m^-^1)')
xlabel('Depth below MHW (m)')
l=legend('No restoration','TLP');
l.Location='east';
l.Color=[0.8 0.8 0.8];
l.EdgeColor='none';
l.BackgroundAlpha=0.3;
l.FontSize=20;
ax=gca;
ax.FontSize=20;
%% fig 8; marsh width and VEC graphs for different erosion rates
% fig 8a
figure
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.6, 0.96],'InnerPosition',[0.09,0.6,0.4,0.4],'Position',[0.09,0.6,0.4,0.4]); %makes really large figure box (takes up most of computer screen)
s=subaxis(2,2,1,'SpacingHoriz',0.04,'SpacingVert',0.04); %limits the padding between subplots
%subplot(2,2,1)
RSLRA=0.08023;
color=[0.99608 0.78824 0.89412]; %no-restoration color
erosion=1;
for k=[1.65e-09,2.82e-09,5.3e-09,7.8e-09,9.95e-09]
    filename=['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Breakwater0_SSCReduction0_breakwaterdelay551_TLP0_TLPthick0_TLPfreq0_TLPdelay551_Erosion' num2str(k) '/'];
    myVars={'elev.mat','marsh edge.mat','meansealevel.mat','AnnualAccretion.mat','SLR.mat','resiliencymetrics.mat','VASEAmatrix.mat'};
    for i=1:7
        load([filename myVars{i}]);
    end
    if k==9.95e-09
        color=[0.44706 0 0.22353]; %for no-restoration
    end
    plot(0:150,VASEAmatrix(551:701,3),'LineWidth',4,'Color',color)
    hold on
    marshage=find(resiliencymetrics(551:701,11)==0,1)+551;
    if isempty(marshage)==1
        marshage=701;
    end
    for yr=552:701
        widthslope(yr)=(VASEAmatrix(yr,3)-VASEAmatrix(yr-1,3))/1; %change in width each year
        marshedgeerosion(yr)=Marsh_edge(yr)-Marsh_edge(yr-1);
    end
    threshold=find(widthslope(552:701)<widthslope(551)-5,1)-1;
    meanerosion=mean(marshedgeerosion(552:marshage-2));

    if erosion==1 && isempty(threshold)==1
        threshold=marshage;
    else if erosion>1 && isempty(threshold)==1 && marshage>thresholdsmatrix(1,erosion-1)
            threshold=marshage;
    else if isempty(threshold)==1
            continue
    end
    end
    end
    scatter(threshold,VASEAmatrix(threshold+551,3),'MarkerFaceColor',color,'MarkerEdgeColor','k')
    hold on
    plot([threshold,120],[VASEAmatrix(threshold+551,3),400],'LineStyle',':','Color',[0.5 0.5 0.5],'LineWidth',1.5)
    color=color+[-0.14 -0.2 -0.15]; %for no restoration color
    erosion=erosion+1;
end
filename=['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Breakwater1_SSCReduction0_breakwaterdelay551_TLP0_TLPthick0_TLPfreq0_TLPdelay551_Erosion' num2str(k) '/'];
myVars={'elev.mat','marsh edge.mat','meansealevel.mat','AnnualAccretion.mat','SLR.mat','resiliencymetrics.mat','VASEAmatrix.mat'};
for i=1:7
    load([filename myVars{i}]);
end
plot(0:150,VASEAmatrix(551:701,3),'Color','#66a61e','LineWidth',4,'LineStyle','--') %#e6ab02 for TLP+breakwater and #66a61e for breakwater
marshage=find(resiliencymetrics(551:701,11)==0,1)+551;
if isempty(marshage)==1
    marshage=701;
end
for yr=552:701
    widthslope(yr)=(VASEAmatrix(yr,3)-VASEAmatrix(yr-1,3))/1; %change in width each year
    marshedgeerosion(yr)=Marsh_edge(yr)-Marsh_edge(yr-1);
end
threshold=find(widthslope(552:701)<widthslope(551)-5,1)-1;
meanerosion=mean(marshedgeerosion(552:marshage-2));

if erosion==1 && isempty(threshold)==1
    threshold=marshage;
else if erosion>1 && isempty(threshold)==1 && marshage>thresholdsmatrix(1,erosion-1)
        threshold=marshage;
end
end
scatter(threshold,VASEAmatrix(threshold+551,3),'MarkerFaceColor','#66a61e','MarkerEdgeColor','k','Marker','^')
hold on
plot([threshold,120],[VASEAmatrix(threshold+551,3),400],'LineStyle',':','Color',[0.5 0.5 0.5],'LineWidth',1.5)
% l=legend('Erosion_{initial} = 0.5 m yr^-^1','Erosion_{initial} = 1 m yr^-^1','Erosion_{initial} = 2 m yr^-^1','Erosion_{initial} = 3 m yr^-^1','Erosion_{initial} = 4 m yr^-^1','Breakwater');
% l.Location='southwest';
% l.EdgeColor='none'; %grey legend outline
% l.Color=[0.97 0.97 0.97 0.97];
% l.FontSize=10;
ylim([0 525])
ylabel('Vegetated marsh width (m)')
%xlabel('Model run year')

ax=gca;
ax.FontSize=16;

text(3,515,'a','FontSize',16,'FontWeight','bold')
text(120,405,'Thresholds','FontSize',12,'Color',[0.5 0.5 0.5])
yticks([0 100 200 300 400 500])

% fig 8b------------------------------------------------------------------------------------------------
subplot(2,2,2)
s=subaxis(2,2,2,'SpacingHoriz',0.04,'SpacingVert',0.04); %limits the padding between subplots

hold on
RSLRA=0.08023;
color=[0.88039 0.796078 0.878431]; %for TLP color (erosion changes)
for k=[1.65e-09,2.82e-09,5.3e-09,7.8e-09,9.95e-09]
    filename=['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Breakwater0_SSCReduction0_breakwaterdelay551_TLP1_TLPthick0.05_TLPfreq15_TLPdelay551_Erosion' num2str(k) '/'];
    myVars={'elev.mat','marsh edge.mat','meansealevel.mat','AnnualAccretion.mat','SLR.mat','resiliencymetrics.mat','VASEAmatrix.mat'};
    for i=1:7
        load([filename myVars{i}]);
    end
    if k==9.95e-09
        color=[0.256862 0.1607843 0.286274]; %for TLP
    else if k==7.8e-09 %for TLP scenarios
            color=[0.3941176 0.30196078 0.53333333];
    else if k==5.3e-09 %for TLP scenarios
            color=[0.6333333 0.5411764 0.74117647];
    end
    end
    end
    plot(0:150,VASEAmatrix(551:701,3),'LineWidth',4,'Color',color)
    hold on
    marshage=find(resiliencymetrics(551:701,11)==0,1)+551;
    if isempty(marshage)==1
        marshage=701;
    end
    for yr=552:701
        widthslope(yr)=(VASEAmatrix(yr,3)-VASEAmatrix(yr-1,3))/1; %change in width each year
        marshedgeerosion(yr)=Marsh_edge(yr)-Marsh_edge(yr-1);
    end
    threshold=find(widthslope(552:701)<widthslope(551)-5,1)-1;
    meanerosion=mean(marshedgeerosion(552:marshage-2));

    if erosion==1 && isempty(threshold)==1
        threshold=marshage;
    end
    if k==7.8e-09
        scatter(0:15:marshage-551,520,80,'^','filled','MarkerFaceColor',[0.5,0.5,0.5],'MarkerEdgeColor','none')
    end
    scatter(threshold,VASEAmatrix(threshold+551,3),'MarkerFaceColor',color,'MarkerEdgeColor','k')
    hold on
    plot([threshold,120],[VASEAmatrix(threshold+551,3),400],'LineStyle',':','Color',[0.5 0.5 0.5],'LineWidth',1.5)
    color=color+[-0.14 -0.2 -0.15]; %for no restoration color
    erosion=erosion+1;
end
filename=['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Breakwater1_SSCReduction0_breakwaterdelay551_TLP1_TLPthick0.05_TLPfreq15_TLPdelay551_Erosion' num2str(k) '/'];
myVars={'elev.mat','marsh edge.mat','meansealevel.mat','AnnualAccretion.mat','SLR.mat','resiliencymetrics.mat','VASEAmatrix.mat'};
for i=1:7
    load([filename myVars{i}]);
end
plot(0:150,VASEAmatrix(551:701,3),'Color','#e6ab02','LineWidth',4,'LineStyle','--') %#e6ab02 for TLP+breakwater and #66a61e for breakwater
marshage=find(resiliencymetrics(551:701,11)==0,1)+551;
if isempty(marshage)==1
    marshage=701;
end
for yr=552:701
    widthslope(yr)=(VASEAmatrix(yr,3)-VASEAmatrix(yr-1,3))/1; %change in width each year
    marshedgeerosion(yr)=Marsh_edge(yr)-Marsh_edge(yr-1);
end
threshold=find(widthslope(552:701)<widthslope(551)-5,1)-1;
meanerosion=mean(marshedgeerosion(552:marshage-2));

if erosion==1 && isempty(threshold)==1
    threshold=marshage;
else if erosion>1 && isempty(threshold)==1 && marshage>thresholdsmatrix(1,erosion-1)
        threshold=marshage;
end
end
scatter(threshold,VASEAmatrix(threshold+551,3),'MarkerFaceColor','#e6ab02','MarkerEdgeColor','k','Marker','^')
hold on
plot([threshold,120],[VASEAmatrix(threshold+551,3),400],'LineStyle',':','Color',[0.5 0.5 0.5],'LineWidth',1.5)
% l=legend('Erosion_{initial} = 0.5 m yr^-^1','Erosion_{initial} = 1 m yr^-^1','Erosion_{initial} = 2 m yr^-^1','Erosion_{initial} = 3 m yr^-^1','Erosion_{initial} = 4 m yr^-^1','Breakwater');
% l.Location='southwest';
% l.EdgeColor='none'; %grey legend outline
% l.Color=[0.97 0.97 0.97 0.97];
% l.FontSize=10;
ylim([0 525])
text(120,405,'Thresholds','FontSize',12,'Color',[0.5 0.5 0.5])
%ylabel('Vegetated marsh width (m)')
%xlabel('Model run year')
ax=gca;
ax.FontSize=16;
text(3,513,'b','FontSize',16,'FontWeight','bold')
yticks([0 100 200 300 400 500]) %sets the ytick labels

%Fig 8c---------------------------------------------------------------------
subplot(2,2,3)
s=subaxis(2,2,3,'SpacingHoriz',0.04,'SpacingVert',0.04); %limits the padding between subplots
RSLRA=0.08023;
color=[0.99608 0.78824 0.89412]; %no-restoration color
for k=[1.65e-09,2.82e-09,5.3e-09,7.8e-09,9.95e-09]
    filename=['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Breakwater0_SSCReduction0_breakwaterdelay551_TLP0_TLPthick0_TLPfreq0_TLPdelay551_Erosion' num2str(k) '/'];
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
    if k==9.95e-09
        color=[0.44706 0 0.22353]; %for no-restoration
    end
    yr=0:150;
    plot(yr,marshvolumeaboveDmax(551:701),'LineWidth',5,'Color',color);
    hold on
    marshage=find(resiliencymetrics(551:701,11)==0,1)+551;
    if isempty(marshage)==1
        marshage=701;
    end
    for yr=552:701
        widthslope(yr)=(VASEAmatrix(yr,3)-VASEAmatrix(yr-1,3))/1; %change in width each year
        marshedgeerosion(yr)=Marsh_edge(yr)-Marsh_edge(yr-1);
    end
    threshold=find(widthslope(552:701)<widthslope(551)-5,1)-1;
    meanerosion=mean(marshedgeerosion(552:marshage-2));

    if erosion==1 && isempty(threshold)==1
        threshold=marshage;
    end

    scatter(threshold,marshvolumeaboveDmax(threshold+551),'MarkerFaceColor',color,'MarkerEdgeColor','k')
    hold on
    plot([threshold,110],[marshvolumeaboveDmax(threshold+551),80],'LineStyle',':','Color',[0.5 0.5 0.5],'LineWidth',1.5)
    hold on
    color=color+[-0.14 -0.2 -0.15]; %for no restoration color
    erosion=erosion+1;
end
filename=['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Breakwater1_SSCReduction0_breakwaterdelay551_TLP0_TLPthick0_TLPfreq0_TLPdelay551_Erosion' num2str(k) '/'];
myVars={'elev.mat','marsh edge.mat','meansealevel.mat','AnnualAccretion.mat','SLR.mat','resiliencymetrics.mat','VASEAmatrix.mat'};
for i=1:7
    load([filename myVars{i}]);
end

marshage=find(resiliencymetrics(551:701,11)==0,1)+551;
if isempty(marshage)==1
    marshage=701;
end
for yr=552:701
    widthslope(yr)=(VASEAmatrix(yr,3)-VASEAmatrix(yr-1,3))/1; %change in width each year
    marshedgeerosion(yr)=Marsh_edge(yr)-Marsh_edge(yr-1);
end
threshold=find(widthslope(552:701)<widthslope(551)-5,1)-1;
meanerosion=mean(marshedgeerosion(552:marshage-2));
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
if erosion==1 && isempty(threshold)==1
    threshold=marshage;
else if erosion>1 && isempty(threshold)==1 && marshage>thresholdsmatrix(1,erosion-1)
        threshold=marshage;
end
end
plot(0:150,marshvolumeaboveDmax(551:701),'Color','#66a61e','LineWidth',4,'LineStyle','--') %#e6ab02 for TLP+breakwater and #66a61e for breakwater

scatter(threshold,marshvolumeaboveDmax(threshold+551),'MarkerFaceColor','#66a61e','MarkerEdgeColor','k','Marker','^')
hold on
plot([threshold,110],[marshvolumeaboveDmax(threshold+551),80],'LineStyle',':','Color',[0.5 0.5 0.5],'LineWidth',1.5)
ylim([0 275])

ylabel('Volumetric elevation capital (m^3)')
xlabel('Model run year')


text(3,268,'c','FontSize',16,'FontWeight','bold')
text(111,81,'Inflection points','FontSize',12,'Color',[0.5 0.5 0.5])
ax=gca;
ax.FontSize=16;

%Fig 8d---------------------------------------------------------------------

subplot(2,2,4)
s=subaxis(2,2,4,'SpacingHoriz',0.04,'SpacingVert',0.04); %limits the padding between subplots
RSLRA=0.08023;
color=[0.88039 0.796078 0.878431]; %for TLP color (erosion changes)
for k=[1.65e-09,2.82e-09,5.3e-09,7.8e-09,9.95e-09]
    filename=['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Breakwater0_SSCReduction0_breakwaterdelay551_TLP1_TLPthick0.05_TLPfreq15_TLPdelay551_Erosion' num2str(k) '/'];
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
    if k==9.95e-09
        color=[0.256862 0.1607843 0.286274]; %for TLP
    else if k==7.8e-09 %for TLP scenarios
            color=[0.3941176 0.30196078 0.53333333];
    else if k==5.3e-09 %for TLP scenarios
            color=[0.6333333 0.5411764 0.74117647];
    end
    end
    end
    yr=0:150;
    plot(yr,marshvolumeaboveDmax(551:701),'LineWidth',5,'Color',color);
    hold on
    color=color+[-0.14 -0.2 -0.15]; %for no restoration color
    erosion=erosion+1;
    if k==7.8e-09
        scatter(0:15:marshage,272,80,'^','filled','MarkerFaceColor',[0.5,0.5,0.5],'MarkerEdgeColor','none')
    end
end
filename=['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Breakwater1_SSCReduction0_breakwaterdelay551_TLP1_TLPthick0.05_TLPfreq15_TLPdelay551_Erosion' num2str(k) '/'];
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
plot(0:150,marshvolumeaboveDmax(551:701),'Color','#e6ab02','LineWidth',4,'LineStyle','--') %#e6ab02 for TLP+breakwater and #66a61e for breakwater
ylim([0 275])
%ylabel('Volumetric elevation capital (m^3)')
xlabel('Model run year')
ax=gca;
ax.FontSize=16;
text(3,268,'d','FontSize',16,'FontWeight','bold')
%% fig 9 SLRA rates vs VEC of different restorations
figure
RSLRA=0.01881;
k=5.3e-09;

% graph sill first
subaxis(2,3,1,'SpacingHoriz',0.04,'SpacingVert',0.1); %limits the padding between subplots
hold on
filename=['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Sill1_SSCReduction0_silldelay551_TLP0_TLPthick0_TLPfreq0_TLPdelay551_Erosion' num2str(k) '/'];
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
p1=plot(SLR(551:701)*1000,marshvolumeaboveDmax(551:701),'LineWidth',5,'Color','#66a61e');%#66a61e
hold on

%graph TLP next
filename=['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Sill0_SSCReduction0_silldelay551_TLP1_TLPthick0.05_TLPfreq15_TLPdelay551_Erosion' num2str(k) '/'];
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
p2=plot(SLR(551:701)*1000,marshvolumeaboveDmax(551:701),'LineWidth',5,'Color','#7570b3'); %#7570b3
hold on

%graph TLP and sill next
filename=['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Sill1_SSCReduction0_silldelay551_TLP1_TLPthick0.05_TLPfreq15_TLPdelay551_Erosion' num2str(k) '/'];
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
p3=plot(SLR(551:701)*1000,marshvolumeaboveDmax(551:701),'LineWidth',5,'Color','#e6ab02'); %#e6ab02
hold on

%finally graph no restoration last as dashed line on top
filename=['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Sill0_SSCReduction0_silldelay551_TLP0_TLPthick0_TLPfreq0_TLPdelay551_Erosion' num2str(k) '/'];
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
p4=plot(SLR(551:701)*1000,marshvolumeaboveDmax(551:701),'LineWidth',5,'Color','#e7298a','LineStyle',':');
hold on
ylabel('Volumetric elevation capital (m^3)')
hold on
increment=(SLR(701)*1000-3)/10;
s1=scatter(3:increment:SLR(701)*1000,283,80,'^','filled','MarkerFaceColor',[0.5,0.5,0.5],'MarkerEdgeColor','none');
hold on
l=legend([p4(1),p1(1),p2(1),p3(1),s1(1)],'No restoration','Shoreline Stabilization','TLP','TLP and Shoreline Stabilization','TLP Year');
l.Location='southwest';
l.Orientation='horizontal';
l.EdgeColor='k';
l.Color=[0.97 0.97 0.97 0.97];
l.FontSize=16;
text((SLR(701)*1000-3)/25+3,278,'a','FontSize',16,'FontWeight','bold')
text((SLR(701)*1000-3)/25+2.9,18,'Low SLRA','FontSize',16,'FontWeight','bold')
xlim([3 SLR(701)*1000])
ylim([0 285])
ax=gca;
ax.FontSize=16;
xlabel('Rate of sea-level rise (mm yr^-^1)');

%-----------------------------------------------------
%figure 9B
RSLRA=0.03536;
k=5.3e-09;

% graph sill first
ax1=subaxis(2,3,2,'SpacingHoriz',0.04,'SpacingVert',0.1); %limits the padding between subplots
hold on
filename=['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Sill1_SSCReduction0_silldelay551_TLP0_TLPthick0_TLPfreq0_TLPdelay551_Erosion' num2str(k) '/'];
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
plot(SLR(551:701)*1000,marshvolumeaboveDmax(551:701),'LineWidth',5,'Color','#66a61e');%#66a61e
hold on

%graph TLP next
filename=['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Sill0_SSCReduction0_silldelay551_TLP1_TLPthick0.05_TLPfreq15_TLPdelay551_Erosion' num2str(k) '/'];
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
plot(SLR(551:701)*1000,marshvolumeaboveDmax(551:701),'LineWidth',5,'Color','#7570b3'); %#7570b3
hold on

%graph TLP and sill next
filename=['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Sill1_SSCReduction0_silldelay551_TLP1_TLPthick0.05_TLPfreq15_TLPdelay551_Erosion' num2str(k) '/'];
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
plot(SLR(551:701)*1000,marshvolumeaboveDmax(551:701),'LineWidth',5,'Color','#e6ab02'); %#e6ab02
hold on

%finally graph no restoration last as dashed line on top
filename=['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Sill0_SSCReduction0_silldelay551_TLP0_TLPthick0_TLPfreq0_TLPdelay551_Erosion' num2str(k) '/'];
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
plot(SLR(551:701)*1000,marshvolumeaboveDmax(551:701),'LineWidth',5,'Color','#e7298a','LineStyle',':');
hold on
%end
hold on
increment=(SLR(701)*1000-3)/10;
scatter(3:increment:SLR(701)*1000,283,80,'^','filled','MarkerFaceColor',[0.5,0.5,0.5],'MarkerEdgeColor','none');
ax1.FontSize=16;
text((SLR(701)*1000-3)/25+3,278,'b','FontSize',16,'FontWeight','bold')
text((SLR(701)*1000-3)/25+2.9,18,'Interm.-Low SLRA','FontSize',16,'FontWeight','bold')
hold on
xlim([3 SLR(701)*1000])
ylim([0 285])
ax=gca;
ax.FontSize=16;
xlabel('Rate of sea-level rise (mm yr^-^1)')

%----------------------------------------------------
%Fig 9c
RSLRA=0.08023;
ax1=subaxis(2,3,3,'SpacingHoriz',0.04,'SpacingVert',0.1); %limits the padding between subplots
ylim([0 525])
for i=0:15:150
    plot([i i],[525 450],'LineStyle','--')
end
hold on
filename=['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Sill1_SSCReduction0_silldelay551_TLP0_TLPthick0_TLPfreq0_TLPdelay551_Erosion' num2str(k) '/'];
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
plot(SLR(551:701)*1000,marshvolumeaboveDmax(551:701),'LineWidth',5,'Color','#66a61e');%#66a61e
hold on

%graph TLP next
filename=['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Sill0_SSCReduction0_silldelay551_TLP1_TLPthick0.05_TLPfreq15_TLPdelay551_Erosion' num2str(k) '/'];
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
plot(SLR(551:701)*1000,marshvolumeaboveDmax(551:701),'LineWidth',5,'Color','#7570b3'); %#7570b3
hold on
increment=(SLR(701)*1000-3)/10;
scatter(3:increment:SLR(marshage+550)*1000,283,80,'^','filled','MarkerFaceColor',[0.5,0.5,0.5],'MarkerEdgeColor','none');

%graph TLP and sill next
filename=['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Sill1_SSCReduction0_silldelay551_TLP1_TLPthick0.05_TLPfreq15_TLPdelay551_Erosion' num2str(k) '/'];
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
plot(SLR(551:701)*1000,marshvolumeaboveDmax(551:701),'LineWidth',5,'Color','#e6ab02'); %#e6ab02
hold on

%finally graph no restoration last as dashed line on top
filename=['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Sill0_SSCReduction0_silldelay551_TLP0_TLPthick0_TLPfreq0_TLPdelay551_Erosion' num2str(k) '/'];
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
plot(SLR(551:701)*1000,marshvolumeaboveDmax(551:701),'LineWidth',5,'Color','#e7298a','LineStyle',':');
hold on
ax1.FontSize=16;
text((SLR(701)*1000-3)/25+3,278,'c','FontSize',16,'FontWeight','bold')
text((SLR(701)*1000-3)/25+2.8,18,'Interm. SLRA','FontSize',16,'FontWeight','bold')
xlim([3 SLR(701)*1000])
ylim([0 285])
ax=gca;
ax.FontSize=16;
xlabel('Rate of sea-level rise (mm yr^-^1)')
%--------------------------------------------------------------------------------
%figure 9d
RSLRA=0.1106;
k=5.3e-09;

% graph sill first
ax1=subaxis(2,3,4,'SpacingHoriz',0.04,'SpacingVert',0.1); %limits the padding between subplots
hold on
filename=['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Sill1_SSCReduction0_silldelay551_TLP0_TLPthick0_TLPfreq0_TLPdelay551_Erosion' num2str(k) '/'];
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
plot(SLR(551:701)*1000,marshvolumeaboveDmax(551:701),'LineWidth',5,'Color','#66a61e');%#66a61e
hold on

%graph TLP next
filename=['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Sill0_SSCReduction0_silldelay551_TLP1_TLPthick0.05_TLPfreq15_TLPdelay551_Erosion' num2str(k) '/'];
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
plot(SLR(551:701)*1000,marshvolumeaboveDmax(551:701),'LineWidth',5,'Color','#7570b3'); %#7570b3
hold on
increment=(SLR(701)*1000-3)/10;
scatter(3:increment:SLR(marshage+550)*1000,283,80,'^','filled','MarkerFaceColor',[0.5,0.5,0.5],'MarkerEdgeColor','none');hold on

%graph TLP and sill next
filename=['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Sill1_SSCReduction0_silldelay551_TLP1_TLPthick0.05_TLPfreq15_TLPdelay551_Erosion' num2str(k) '/'];
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
plot(SLR(551:701)*1000,marshvolumeaboveDmax(551:701),'LineWidth',5,'Color','#e6ab02'); %#e6ab02
hold on

%finally graph no restoration last as dashed line on top
filename=['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Sill0_SSCReduction0_silldelay551_TLP0_TLPthick0_TLPfreq0_TLPdelay551_Erosion' num2str(k) '/'];
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
plot(SLR(551:701)*1000,marshvolumeaboveDmax(551:701),'LineWidth',5,'Color','#e7298a','LineStyle',':');
hold on
ylabel('Volumetric elevation capital (m^3)')
ax1.FontSize=16;
text((SLR(701)*1000-3)/25+3,278,'d','FontSize',16,'FontWeight','bold')
text((SLR(701)*1000-3)/25+2.5,18,'Interm.-High SLRA','FontSize',16,'FontWeight','bold')
xlim([3 SLR(701)*1000])
ylim([0 285])
ax=gca;
ax.FontSize=16;
xlabel('Rate of sea-level rise (mm yr^-^1)')
%-----------------------------------------------------
%Fig 9e
RSLRA=0.1592;
ax1=subaxis(2,3,5,'SpacingHoriz',0.04,'SpacingVert',0.1); %limits the padding between subplots
ylim([0 525])
hold on
filename=['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Sill1_SSCReduction0_silldelay551_TLP0_TLPthick0_TLPfreq0_TLPdelay551_Erosion' num2str(k) '/'];
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
plot(SLR(551:701)*1000,marshvolumeaboveDmax(551:701),'LineWidth',5,'Color','#66a61e');%#66a61e
hold on

%graph TLP next
filename=['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Sill0_SSCReduction0_silldelay551_TLP1_TLPthick0.05_TLPfreq15_TLPdelay551_Erosion' num2str(k) '/'];
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
plot(SLR(551:701)*1000,marshvolumeaboveDmax(551:701),'LineWidth',5,'Color','#7570b3'); %#7570b3
hold on

%graph TLP and sill next
filename=['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Sill1_SSCReduction0_silldelay551_TLP1_TLPthick0.05_TLPfreq15_TLPdelay551_Erosion' num2str(k) '/'];
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
plot(SLR(551:701)*1000,marshvolumeaboveDmax(551:701),'LineWidth',5,'Color','#e6ab02'); %#e6ab02
hold on
increment=(SLR(701)*1000-3)/10;
scatter(3:increment:SLR(marshage+550)*1000,283,80,'^','filled','MarkerFaceColor',[0.5,0.5,0.5],'MarkerEdgeColor','none');
hold on

%finally graph no restoration last as dashed line on top
filename=['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Sill0_SSCReduction0_silldelay551_TLP0_TLPthick0_TLPfreq0_TLPdelay551_Erosion' num2str(k) '/'];
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
plot(SLR(551:701)*1000,marshvolumeaboveDmax(551:701),'LineWidth',5,'Color','#e7298a','LineStyle',':');
hold on
ax1.FontSize=16;
text((SLR(701)*1000-3)/25+3,278,'e','FontSize',16,'FontWeight','bold')
text((SLR(701)*1000-3)/25+2.5,18,'High SLRA','FontSize',16,'FontWeight','bold')
xlim([3 SLR(701)*1000])
ylim([0 285])
ax=gca;
ax.FontSize=16;
xlabel('Rate of sea-level rise (mm yr^-^1)');
inspect
