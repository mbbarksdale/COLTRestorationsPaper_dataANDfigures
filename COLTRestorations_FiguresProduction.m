%% fig 1a OM profile from spinup; percentorganicforspinup
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


%% Fig 1b; plot sea level over time
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

text(100,4,'SLR = 1 mm yr^-^1','FontSize',8)
text(350,4,{'Ramp to','3 mm yr^-^1'},'FontSize',8,'HorizontalAlignment', 'center')
text(475,4,{'SLR =','3 mm yr^-^1'},'FontSize',8,'HorizontalAlignment', 'center')
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

%% Fig 2a; SLR vs accretion ratio
figure
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]); %makes really large figure box (takes up most of computer screen)
subplot(2,2,[1,3])
s=subaxis(2,2,[1,3],'SpacingHoriz',0.1,'SpacingVert',0.1); %make subpanels closer together
color=[0.99608 0.78824 0.89412]; %start with color for low erosion scenario and will increase for every iteration of for-loop
for k=[1.65e-09,2.82e-09,5.3e-09,7.8e-09,9.95e-09] %creates for-loop to cycle through 5 erosion scenarios (changes Be, or Ke, the erosion coefficient)
    filename=['Outputs_CO10_RSLRA0.1106_spinuprampedRSLR1_runstartRSLR3_Sill0_SSCReduction0_silldelay551_TLP0_TLPthick0_TLPfreq0_TLPdelay551_Erosion' num2str(k) '/']; %calls the folder
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
    p=plot(SLR(551:marshage)*1000,accretiontoSLR(551:marshage),'LineWidth',6,'Color',color);
    hold on
    color=color+[-0.14 -0.2 -0.15];
end
ylim([-1 1.2])
xlim([3 14])
ax=gca;
ax.FontSize=22;
% ax.XColor=[0.5 0.5 0.5];
% ax.YColor=[0.5 0.5 0.5];
yline(1,'--',{'Accretion surplus'},'LineWidth',2,'FontSize',18);
dim=[0.347 0.73 .1 .1]; %sets dimensions to write accretion deficit
a=annotation('textbox',dim,'String','Accretion deficit','EdgeColor','none','FontSize',18);
dim2=[0.1 0.806 .1 .1]; %sets dimensions to write "a" in upper left of plot
b=annotation('textbox',dim2,'String','a','EdgeColor','none','FontSize',18,'FontWeight','bold');
l=legend('Erosion_{initial} = 0.5 m yr^-^1','Erosion_{initial} = 1 m yr^-^1','Erosion_{initial} = 2 m yr^-^1','Erosion_{initial} = 3 m yr^-^1','Erosion_{initial} = 4 m yr^-^1');
l.Location='southwest';
l.EdgeColor='none'; %grey legend outline
l.Color=[0.97 0.97 0.97 0.97];
l.FontSize=18;
ylabel('Accretion excess (dimensionless)')
xlabel('Rate of sea-level rise (mm yr^-^1)')
yticks([-1 -0.5 0 0.5 1]) %sets the ytick labels 

%-----------------------------------------------------------------------------------------------------------------------------
% figure 2b graph erosion on x vs SLR threshold rate on y 

subplot(2,2,2)
s=subaxis(2,2,2,'SpacingHoriz',0.1,'SpacingVert',0.1); %limits the padding between subplots
color=[0.99608 0.78824 0.89412];
edgecolor='none';
acceleration=1; %will cycle through acceleration 
RSLRA=0.08023;
    erosion=1; %start with erosion=1 and will increase for every iteration of for-loop, used to designate correct matrix row to place values in each iteration of loop
for k=[1.65e-09,2.82e-09,5.3e-09,7.8e-09,9.95e-09,9.900E-10,3.812E-09,6.300E-09,9.090E-09] %,1.884E-09, 2.118E-09, 2.352E-09, 2.586E-09, 3.812E-09, 4.804E-09, 5.800E-09, 6.300E-09, 6.800E-09, 7.300E-09, 8.230E-09, 8.660E-09, 9.090E-09, 9.520E-09, 7.300E-09, 8.230E-09, 8.660E-09, 9.090E-09, 9.520E-09, 3.300E-10, 6.600E-10, 9.900E-10, 1.320E-09, 3.316E-09, 4.308E-09]  
    filename=['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Sill0_SSCReduction0_silldelay551_TLP0_TLPthick0_TLPfreq0_TLPdelay551_Erosion' num2str(k) '/'];
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
    scatter(meanerosion,thresholdsmatrix(erosion),100,'MarkerFaceColor',color,'MarkerEdgeColor',edgecolor);
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
p=plot(xFit, yFit, 'color', [0.7 0.7 0.7], 'LineWidth', 4); % Plot fitted line.
uistack(p,'bottom') %pushes line below markers

ylim([8.5 9.3])
hold on; % Set hold on so the next plot does not blow away the one we just drew.
coefficients=table2array(mdl.Coefficients);
str = sprintf('SLR_{drowning.start} = %.2f(erosion)^2 %.2f(erosion) + %.2f', round(coefficients(3,1), 3), round(coefficients(2,1), 3), round(coefficients(1,1), 3));
%str= {'SLR_{drowning.start} = ' round(coefficients(3,1),2) '(erosion)^2' round(coefficients(2,1)) '(erosion) +' round(coefficients(1,1))};
totalstr = strjoin(string(str));
%text(0.2,8.55,totalstr,'FontSize',18,'Color','#48CAE4');
text(0.3,9.25,totalstr,'FontSize',18,'Color',[.6 .6 .6]);
dim2=[0.55 0.808 .1 .1];
b=annotation('textbox',dim2,'String','b','EdgeColor','none','FontSize',18,'FontWeight','bold');
ylabel({'Rate of SLR at which interior';'drowning initiates (mm yr^-^1)'})

ax=gca;
ax.FontSize=22;
box on
xlim([0 5])
hold on %-------------------

%add sills dotted line
%     filename=['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Sill1_SSCReduction20_silldelay551_TLP0_TLPthick0_TLPfreq0_TLPdelay551_Erosion5.3e-09/'];
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

%figure 2c relating erosion to time of drowning and survival

subplot(2,2,4)
s=subaxis(2,2,4,'SpacingHoriz',0.1,'SpacingVert',0.1); %limits the padding between subplots
% yline(threshold+551,'LineWidth',4,'Color','#66a61e','LineStyle',':')
% hold on
% yline(marshage+551,'LineWidth',4,'Color','#66a61e','LineStyle',':')
% hold on
%creates legend for the two lines
p1=plot([0.2 0.5],[642.5 642.5],'Color','#006CA5','LineWidth',4);
hold on
p2=plot([0.2 0.5],[637.2 637.2],'Color','#48CAE4','LineWidth',4);
hold on
s1=scatter(0.35,642.5,125,'Marker','square','MarkerEdgeColor',[0.71 0.38824 0.59412],'MarkerFaceColor','none','LineWidth',2.5);
uistack(s1,'top') %brings the square above the line for the legend
hold on
s2=scatter(0.35,637.2,115,'Marker','^','MarkerEdgeColor',[0.71 0.38824 0.59412],'MarkerFaceColor','none','LineWidth',2.5);
hold on
text(0.6,642.5,'Age at total collapse','EdgeColor','none','FontSize',18);
text(0.6,637.2,'Age at initiation of interior drowning','EdgeColor','none','FontSize',18);
hold on
rectangle('Position',[0.13 633 3.3 13],'FaceColor',[0.8 0.8 0.8],'FaceAlpha',0.2,'EdgeColor','none')

edgecolor='none';
color=[0.99608 0.78824 0.89412]; %for no restoration color (erosion changes)
acceleration=1;
RSLRA=0.08023;
    erosion=1;
for k=[1.65e-09,2.82e-09,5.3e-09,7.8e-09,9.95e-09,9.900E-10,3.812E-09,6.300E-09,9.090E-09] 
    filename=['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Sill0_SSCReduction0_silldelay551_TLP0_TLPthick0_TLPfreq0_TLPdelay551_Erosion' num2str(k) '/'];
    myVars={'elev.mat','marsh edge.mat','meansealevel.mat','AnnualAccretion.mat','SLR.mat','resiliencymetrics.mat','fetch.mat','bgb.mat'};
    for i=1:8
        load([filename myVars{i}]);
    end
    if k==9.95e-09
        color=[0.44706 0 0.22353];
    end

    marshage=find(resiliencymetrics(551:701,11)==0,1)+551;
    if isempty(marshage)==1
        marshage=701;
    end

    drowningfound=0;

    for yr=552:701
        widthslope(yr)=(resiliencymetrics(yr,7)-resiliencymetrics(yr-1,7))/1; %change in width each year
        marshedgeerosion(yr)=Marsh_edge(yr)-Marsh_edge(yr-1);
        drowningstart=find(bgb(yr,Marsh_edge(yr)+1:5500)==0);
        if isempty(drowningstart)==0 && drowningfound==0
            drowningstartmatrix(1,erosion)=yr;
            drowningfound=1;
        end
    end

    threshold=find(widthslope(552:701)<widthslope(551)-5,1)-1+551;
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

    if isempty(drowningstartmatrix(1,erosion))==1
        drowningstartmatrix(1,erosion)=701;
        color=[0 0 0];
    end

    if erosion>5
        color=[1 1 1];
        edgecolor=[0.5 0.5 0.5];
    end
    thresholdsmatrix(1,erosion)=drowningstartmatrix(1,erosion); %fill in thresholdsmatrix with SLR rate at which threshold occurs
    meanerosionmatrix(1,erosion)=meanerosion;
    %thresholdsmatrix(1,erosion)=threshold;
    scatter(meanerosion,thresholdsmatrix(erosion),110,'MarkerFaceColor',color,'MarkerEdgeColor',edgecolor,'Marker','^')
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
p=plot(xFit, yFit, 'color', '#48CAE4', 'LineWidth', 4); % Plot fitted line.

coefficients=table2array(mdl2.Coefficients);
str = sprintf('Age_{drowning.start} = %.2f(erosion)^2 %.2f(erosion) + %.2f', round(coefficients(3,1), 3), round(coefficients(2,1), 3), round(coefficients(1,1), 3));
totalstr = strjoin(string(str));
text(0.3,617,totalstr,'FontSize',18,'Color','#48CAE4');
dim2=[0.55 0.355 .1 .1];
b=annotation('textbox',dim2,'String','c','EdgeColor','none','FontSize',18,'FontWeight','bold');
xlabel('Average rate of erosion (m yr^-^1)')
ylabel('Marsh age (yr)')
ax=gca;
ax.FontSize=18;
uistack(p,'bottom') %moves the linear regression line to the bottom, so that markers are above it

%add sill marker in green
%     filename=['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Sill1_SSCReduction20_silldelay551_TLP0_TLPthick0_TLPfreq0_TLPdelay551_Erosion5.3e-09/'];
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
    filename=['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Sill0_SSCReduction0_silldelay551_TLP0_TLPthick0_TLPfreq0_TLPdelay551_Erosion' num2str(k) '/'];
    myVars={'elev.mat','marsh edge.mat','meansealevel.mat','AnnualAccretion.mat','SLR.mat','resiliencymetrics.mat','fetch.mat','bgb.mat'};
    for i=1:8
        load([filename myVars{i}]);
    end
    if k==9.95e-09
         color=[0.44706 0 0.22353];
    end
 
    marshage=find(resiliencymetrics(551:701,11)==0,1)+551;
    if isempty(marshage)==1
        marshage=701;
    end
    for yr=552:701
        marshedgeerosion(yr)=Marsh_edge(yr)-Marsh_edge(yr-1);
    end

    if erosion>5
        color=[1 1 1];
        edgecolor=[0.5 0.5 0.5];
    end
    meanerosion=mean(marshedgeerosion(552:marshage-2));       
    meanerosionmatrix(1,erosion)=meanerosion; %inserts the avg erosion into the new meanerosion matrix which holds avg erosion for each erosion scenario even as for-loop turns to next iteration
    marshagematrix(1,erosion)=marshage;
    scatter(meanerosion,marshage,125,'MarkerFaceColor',color,'MarkerEdgeColor',edgecolor,'Marker','square')
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
p1=plot(xFit, yFit, 'color','#006CA5', 'LineWidth', 4); % Plot fitted line.
hold on
uistack(p1,'bottom') %moves the linear regression line to the bottom, so that markers are above it
coefficients=table2array(mdl3.Coefficients);
str = sprintf('Age_{marsh.collapse} = %.2f(erosion)^2 %.2f(erosion) + %.2f', round(coefficients(3,1), 3), round(coefficients(2,1), 3), round(coefficients(1,1), 3));
totalstr = strjoin(string(str));
text(0.3,675,totalstr,'FontSize',18,'Color','#006CA5');
ax=gca;
ax.FontSize=22;
% ax.XColor=[0.5 0.5 0.5];
% ax.YColor=[0.5 0.5 0.5];
ylim([610 680])
xlim([0 5])
hold on

%add sill marker in green
%     filename=['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Sill1_SSCReduction20_silldelay551_TLP0_TLPthick0_TLPfreq0_TLPdelay551_Erosion5.3e-09/'];
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



%% Fig3a; marsh width graph
k=5.3e-09;
RSLRA=0.08023;
figure
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.6, 0.96],'InnerPosition',[0.09,0.6,0.4,0.4],'Position',[0.09,0.6,0.4,0.4]); %makes really large figure box (takes up most of computer screen)
subplot(2,2,1)
s=subaxis(2,2,1,'SpacingHoriz',0.1,'SpacingVert',0.04); %limits the padding between subplots
filename=['CO10_RSLRA' num2str(RSLRA) '_Erosion' num2str(k) '/VASEAmatrix_'];
myVars={'control','sill','TLP','TLPandsill'};
for i=1:4
    load([filename myVars{i}]);
end
p1=plot(0:150,VASEAmatrix_sill(551:701,3),'lineWidth',5,'Color','#66a61e'); %#66a61e
hold on
p2=plot(0:150,VASEAmatrix_TLP(551:701,3),'lineWidth',5,'Color','#7570b3'); %#7570b3
hold on
p3=plot(0:150,VASEAmatrix_TLPandsill(551:701,3),'lineWidth',5,'Color','#e6ab02'); %#e6ab02
hold on
p4=plot(0:150,VASEAmatrix_control(551:701,3),'lineWidth',5,'Color','#e7298a','LineStyle',':'); %#e7298a
hold on
ylim([0 525])
for i=0:15:150
    plot([i i],[525 455],'LineStyle','--','Color','k')
end
l2=legend([p4(1),p1(1),p2(1),p3(1)],'No restoration','Sill','TLP','TLP and sill');
l2.Location='southwest';
l2.EdgeColor='none';
l2.Color=[0.97 0.97 0.97 0.97];
l2.FontSize=12;
ax=gca;
ax.FontSize=16;
%xlabel('Model Run Year')
ylabel('Width of vegetated marsh (m)')

text(3,515,'a','FontSize',16,'FontWeight','bold')
text(130,467,'TLP','FontSize',16,'Rotation',90)
text(115,467,'TLP','FontSize',16,'Rotation',90)


yticks([0 100 200 300 400 500]) %sets the ytick labels 

% ----------------------------------------------------------------------------------------------------------------------------
% Fig3b w volumetric elevation capital for diff restoration strategies 

% graph sill first 
subplot(2,2,2)
s=subaxis(2,2,2,'SpacingHoriz',0.1,'SpacingVert',0.04); %limits the padding between subplots

    ylim([0 525])
    for i=0:15:150
        plot([i i],[525 450],'LineStyle','--')
    end
hold on

%RSLRA=0.1592;
    filename=['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Sill1_SSCReduction20_silldelay551_TLP0_TLPthick0_TLPfreq0_TLPdelay551_Erosion' num2str(k) '/'];
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
    pond=find(elev(yr, Marsh_edge(yr)+1:5501) < yrDmax(yr),1);
    if isempty(pond)==1
        pond=5501-Marsh_edge(yr)+1;
    end
    diff=elev(yr,Marsh_edge(yr):pond+Marsh_edge(yr))-yrDmax(yr);
    marshvolumeaboveDmax(yr)=sum(diff);
    clear diff
end
    p1=plot(0:150,marshvolumeaboveDmax(551:701),'LineWidth',5,'Color','#66a61e');%#66a61e
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
    pond=find(elev(yr, Marsh_edge(yr)+1:5501) < yrDmax(yr),1);
    if isempty(pond)==1
        pond=5501-Marsh_edge(yr)+1;
    end
    diff=elev(yr,Marsh_edge(yr):pond+Marsh_edge(yr))-yrDmax(yr);
    marshvolumeaboveDmax(yr)=sum(diff);
    clear diff
end
    p2=plot(0:150,marshvolumeaboveDmax(551:701),'LineWidth',5,'Color','#7570b3'); %#7570b3
   hold on

%graph TLP and sill next
    filename=['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Sill1_SSCReduction20_silldelay551_TLP1_TLPthick0.05_TLPfreq15_TLPdelay551_Erosion' num2str(k) '/'];
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
    pond=find(elev(yr, Marsh_edge(yr)+1:5501) < yrDmax(yr),1);
    if isempty(pond)==1
        pond=5501-Marsh_edge(yr)+1;
    end
    diff=elev(yr,Marsh_edge(yr):pond+Marsh_edge(yr))-yrDmax(yr);
    marshvolumeaboveDmax(yr)=sum(diff);
    clear diff
end
    p3=plot(0:150,marshvolumeaboveDmax(551:701),'LineWidth',5,'Color','#e6ab02'); %#e6ab02
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
    pond=find(elev(yr, Marsh_edge(yr)+1:5501) < yrDmax(yr),1);
    if isempty(pond)==1
        pond=5501-Marsh_edge(yr)+1;
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
for i=0:15:150
    plot([i i],[275 240],'LineStyle','--','Color','k')
    hold on
end
l=legend([p4(1),p1(1),p2(1),p3(1)],'No restoration','Sill','TLP','TLP and sill');
l.Location='southwest';
l.EdgeColor='none';
l.Color=[0.97 0.97 0.97 0.97];
l.FontSize=12;
ax=gca;
ax.FontSize=16;

% text(3,268,'b','FontSize',16,'FontWeight','bold')
text(130,245,'TLP','FontSize',16,'Rotation',90)
text(115,245,'TLP','FontSize',16,'Rotation',90)
% ----------------------------------------------------------------------------------------------------------------------------
% fig 3c; create accretion components graph 
% *keep org and min separate until the very end
%figure out percentage of Fb_org and Fb_min made up of by eroded material (saved in fluxes matrix)
filename=['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Sill0_SSCReduction0_silldelay551_TLP0_TLPthick0_TLPfreq0_TLPdelay551_Erosion' num2str(k) '/'];
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

marshage=find(VASEAmatrix(551:701,3)==0,1)-1;
if isempty(marshage)==1
    marshage=150;
end

switchyr=find(accretioncomponents(551:701,3)<=0,2);
switchyr=switchyr(find(switchyr>30));
if isempty(switchyr)==1
    switchyr=marshage;
end

subplot(2,2,3)
s=subaxis(2,2,3,'SpacingHoriz',0.1,'SpacingVert',0.04); %limits the padding between subplots
f=area(switchyr-3:switchyr,accretioncomponents(switchyr+548:switchyr+551,3));
f(1).FaceColor='#ece2f0';
f(1).EdgeColor='none';
hold on
a=area(0:switchyr-2,accretioncomponents(551:switchyr+549,1:3));
ylabel({'Average annual accretion rate'; '(mm yr^-^1 m^-^1)'})
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
plot(0:marshage,totalaccretion(551:marshage+551)*1000,'Color','#e7298a','LineWidth',5)
hold on
plot(0:marshage,SLR(551:marshage+551)*1000,'Color','k','LineWidth',4,'LineStyle',':')
hold on
yline(0)
xlim([0 150])
ylim([-8 12])

ax=gca;
ax.FontSize=16;

text(3,11.4,'c','FontSize',16,'FontWeight','bold')
text(62,7.3,'total accretion','FontSize',16,'Color','#e7298a','Rotation',-24)
text(75,10,'SLR','FontSize',16,'Rotation',26)
text(9,11.4,'No restoration','FontSize',16,'Color','#e7298a')

text(4,-4,'Accretionary components','FontSize',12)
plot([3,74],[-4.5 -4.5],'Color','k')

text(10,-5,'Net belowground biomass','FontSize',12)
rectangle('Position',[3 -5.4 6 0.8],'FaceColor','#ece2f0','EdgeColor','none')
text(10,-6,'Recycled eroded sediment','FontSize',12)
rectangle('Position',[3 -6.4 6 0.8],'FaceColor','#4c516d','EdgeColor','none')
text(10,-7,'External & mudflat sediment','FontSize',12)
rectangle('Position',[3 -7.4 6 0.8],'FaceColor','#838996','EdgeColor','none')

yticks([-8 -4 0 4 8 12]) %sets the ytick labels 

% ----------------------------------------------------------------------------------------------------------------------------
% fig 3d; create accretion components graph 
% *keep org and min separate until the very end
%figure out percentage of Fb_org and Fb_min made up of by eroded material (saved in fluxes matrix)

filename=['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Sill1_SSCReduction20_silldelay551_TLP0_TLPthick0_TLPfreq0_TLPdelay551_Erosion' num2str(k) '/'];
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

marshage=find(VASEAmatrix(551:701,3)==0,1)-1;
if isempty(marshage)==1
    marshage=150;
end

switchyr=find(accretioncomponents(551:701,3)<=0,2);
switchyr=switchyr(find(switchyr>30));
if isempty(switchyr)==1
    switchyr=marshage;
end

subplot(2,2,4)
s=subaxis(2,2,4,'SpacingHoriz',0.1,'SpacingVert',0.04); %limits the padding between subplots

f=area(switchyr-3:switchyr,accretioncomponents(switchyr+548:switchyr+551,3));
f(1).FaceColor='#ece2f0';
f(1).EdgeColor='none';
hold on
a=area(0:switchyr-2,accretioncomponents(551:switchyr+549,1:3));
ylabel({'Average annual accretion rate'; '(mm yr^-^1 m^-^1)'})
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
plot(0:marshage,totalaccretion(551:marshage+551)*1000,'Color','#66a61e','LineWidth',5)
hold on
plot(0:marshage,SLR(551:marshage+551)*1000,'Color','k','LineWidth',4,'LineStyle',':')
hold on
yline(0)
xlim([0 150])
ylim([-8 12])

ax=gca;
ax.FontSize=16;

text(3,11.4,'d','FontSize',16,'FontWeight','bold')
text(48,6.2,'total accretion','FontSize',16,'Color','#66a61e','Rotation',-18)
text(75,9.9,'SLR','FontSize',16,'Rotation',26)
text(9,10.5,{'Sill = no erosional','component'},'FontSize',16,'Color','#66a61e')

text(4,-5,'Accretionary components','FontSize',12)
plot([3,74],[-5.5 -5.5],'Color','k')

text(10,-6,'Net belowground biomass','FontSize',12)
rectangle('Position',[3 -6.4 6 0.8],'FaceColor','#ece2f0','EdgeColor','none')
text(10,-7,'External & mudflat sediment','FontSize',12)
rectangle('Position',[3 -7.4 6 0.8],'FaceColor','#838996','EdgeColor','none')

yticks([-8 -4 0 4 8 12]) %sets the ytick labels 



%% fig 4; marsh width and VEC graphs for different erosion rates
% fig 4a
figure
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.6, 0.96],'InnerPosition',[0.09,0.6,0.4,0.4],'Position',[0.09,0.6,0.4,0.4]); %makes really large figure box (takes up most of computer screen)
s=subaxis(2,2,1,'SpacingHoriz',0.04,'SpacingVert',0.04); %limits the padding between subplots
%subplot(2,2,1)
RSLRA=0.08023;
color=[0.99608 0.78824 0.89412]; %no-restoration color
erosion=1;
for k=[1.65e-09,2.82e-09,5.3e-09,7.8e-09,9.95e-09]
    filename=['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Sill0_SSCReduction0_silldelay551_TLP0_TLPthick0_TLPfreq0_TLPdelay551_Erosion' num2str(k) '/'];
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
filename=['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Sill1_SSCReduction20_silldelay551_TLP0_TLPthick0_TLPfreq0_TLPdelay551_Erosion' num2str(k) '/'];
 myVars={'elev.mat','marsh edge.mat','meansealevel.mat','AnnualAccretion.mat','SLR.mat','resiliencymetrics.mat','VASEAmatrix.mat'};
    for i=1:7
        load([filename myVars{i}]);
    end
plot(0:150,VASEAmatrix(551:701,3),'Color','#66a61e','LineWidth',4,'LineStyle','--') %#e6ab02 for TLP+sill and #66a61e for sill
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
% l=legend('Erosion_{initial} = 0.5 m yr^-^1','Erosion_{initial} = 1 m yr^-^1','Erosion_{initial} = 2 m yr^-^1','Erosion_{initial} = 3 m yr^-^1','Erosion_{initial} = 4 m yr^-^1','Sill');
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

% fig 4b
% ------------------------------------------------------------------------------------------------
subplot(2,2,2)
s=subaxis(2,2,2,'SpacingHoriz',0.04,'SpacingVert',0.04); %limits the padding between subplots
for i=0:15:150
    plot([i i],[525 470],'LineStyle','--','Color','k')
    hold on
end

hold on
RSLRA=0.08023;
color=[0.88039 0.796078 0.878431]; %for TLP color (erosion changes)
for k=[1.65e-09,2.82e-09,5.3e-09,7.8e-09,9.95e-09]
        filename=['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Sill0_SSCReduction0_silldelay551_TLP1_TLPthick0.05_TLPfreq15_TLPdelay551_Erosion' num2str(k) '/'];
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

    scatter(threshold,VASEAmatrix(threshold+551,3),'MarkerFaceColor',color,'MarkerEdgeColor','k')
    hold on
    plot([threshold,120],[VASEAmatrix(threshold+551,3),400],'LineStyle',':','Color',[0.5 0.5 0.5],'LineWidth',1.5)
    color=color+[-0.14 -0.2 -0.15]; %for no restoration color
    erosion=erosion+1;
end
filename=['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Sill1_SSCReduction20_silldelay551_TLP1_TLPthick0.05_TLPfreq15_TLPdelay551_Erosion' num2str(k) '/'];
myVars={'elev.mat','marsh edge.mat','meansealevel.mat','AnnualAccretion.mat','SLR.mat','resiliencymetrics.mat','VASEAmatrix.mat'};
    for i=1:7
        load([filename myVars{i}]);
    end
plot(0:150,VASEAmatrix(551:701,3),'Color','#e6ab02','LineWidth',4,'LineStyle','--') %#e6ab02 for TLP+sill and #66a61e for sill
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
% l=legend('Erosion_{initial} = 0.5 m yr^-^1','Erosion_{initial} = 1 m yr^-^1','Erosion_{initial} = 2 m yr^-^1','Erosion_{initial} = 3 m yr^-^1','Erosion_{initial} = 4 m yr^-^1','Sill');
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
text(131,467,'TLP','FontSize',16,'Rotation',90)
text(116,467,'TLP','FontSize',16,'Rotation',90)
text(146,467,'TLP','FontSize',16,'Rotation',90)
%text(120,430,'TLP','FontSize',16)

% X=[130 130 130];
% Y=[432 432 432];
% U=[-15,0,15];
% V=[33 33 33];
% quiver(X,Y,U,V,'k')

yticks([0 100 200 300 400 500]) %sets the ytick labels 

%Fig 4c---------------------------------------------------------------------
subplot(2,2,3)
s=subaxis(2,2,3,'SpacingHoriz',0.04,'SpacingVert',0.04); %limits the padding between subplots
RSLRA=0.08023;
color=[0.99608 0.78824 0.89412]; %no-restoration color
for k=[1.65e-09,2.82e-09,5.3e-09,7.8e-09,9.95e-09]   
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
    pond=find(elev(yr, Marsh_edge(yr)+1:5501) < yrDmax(yr),1);
    if isempty(pond)==1
        pond=5501-Marsh_edge(yr)+1;
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
filename=['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Sill1_SSCReduction20_silldelay551_TLP0_TLPthick0_TLPfreq0_TLPdelay551_Erosion' num2str(k) '/'];
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
    pond=find(elev(yr, Marsh_edge(yr)+1:5501) < yrDmax(yr),1);
    if isempty(pond)==1
        pond=5501-Marsh_edge(yr)+1;
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
    plot(0:150,marshvolumeaboveDmax(551:701),'Color','#66a61e','LineWidth',4,'LineStyle','--') %#e6ab02 for TLP+sill and #66a61e for sill

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

%Fig 4d---------------------------------------------------------------------

subplot(2,2,4)
s=subaxis(2,2,4,'SpacingHoriz',0.04,'SpacingVert',0.04); %limits the padding between subplots
RSLRA=0.08023;
color=[0.88039 0.796078 0.878431]; %for TLP color (erosion changes)
for k=[1.65e-09,2.82e-09,5.3e-09,7.8e-09,9.95e-09]   
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
    pond=find(elev(yr, Marsh_edge(yr)+1:5501) < yrDmax(yr),1);
    if isempty(pond)==1
        pond=5501-Marsh_edge(yr)+1;
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
end
filename=['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Sill1_SSCReduction20_silldelay551_TLP1_TLPthick0.05_TLPfreq15_TLPdelay551_Erosion' num2str(k) '/'];
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
    pond=find(elev(yr, Marsh_edge(yr)+1:5501) < yrDmax(yr),1);
    if isempty(pond)==1
        pond=5501-Marsh_edge(yr)+1;
    end
    diff=elev(yr,Marsh_edge(yr):pond+Marsh_edge(yr))-yrDmax(yr);
    marshvolumeaboveDmax(yr)=sum(diff);
    clear diff
end 
plot(0:150,marshvolumeaboveDmax(551:701),'Color','#e6ab02','LineWidth',4,'LineStyle','--') %#e6ab02 for TLP+sill and #66a61e for sill
for i=0:15:150
    plot([i i],[275 242],'LineStyle','--','Color','k')
end
ylim([0 275])
%ylabel('Volumetric elevation capital (m^3)')
xlabel('Model run year')
ax=gca;
ax.FontSize=16;

text(3,268,'d','FontSize',16,'FontWeight','bold')
text(131,245,'TLP','FontSize',16,'Rotation',90)
text(116,245,'TLP','FontSize',16,'Rotation',90)
text(146,245,'TLP','FontSize',16,'Rotation',90)


%% SUPFIG 1; plot sea level over time
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


%% SUP Fig. 3 with accretionary components for TLP, TLP+sill
% ----------------------------------------------------------------------------------------------------------------------------
% fig Xa; create accretion components graph 
% *keep org and min separate until the very end
%figure out percentage of Fb_org and Fb_min made up of by eroded material (saved in fluxes matrix)
figure
s1=subaxis(1,2,1,'SpacingHoriz',0.04,'SpacingVert',0.1); %limits the padding between subplots
k=7.8e-09;
RSLRA=0.08023;
filename=['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Sill0_SSCReduction0_silldelay551_TLP1_TLPthick0.05_TLPfreq15_TLPdelay551_Erosion' num2str(k) '/'];
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

marshage=find(VASEAmatrix(551:701,3)==0,1)-1;
if isempty(marshage)==1
    marshage=150;
end

switchyr=find(accretioncomponents(551:701,3)<=0,2);
switchyr=switchyr(find(switchyr>30));
if isempty(switchyr)==1
    switchyr=marshage;
end

f=area(switchyr-3:switchyr,accretioncomponents(switchyr+548:switchyr+551,3));
f(1).FaceColor='none'; %#ece2f0
f(1).EdgeColor='none';
hold on
a=area(0:switchyr-2,accretioncomponents(551:switchyr+549,1:3));
ylabel('Average annual accretion rate (mm yr^-^1 m^-^1)')
xlabel('Model Run Year')

a(1).FaceColor='none'; %#838996
a(1).EdgeColor='none';
a(2).FaceColor='none'; %#4c516d
a(2).EdgeColor='none';
a(3).FaceColor='none'; %#ece2f0
a(3).EdgeColor='none';
hold on
e=area(switchyr-2:marshage,accretioncomponents(switchyr+549:marshage+551,1:2));
e(1).FaceColor='none'; %#838996
e(1).EdgeColor='none';
e(2).FaceColor='none'; %#4c516d
e(2).EdgeColor='none';
hold on
d=area(switchyr:marshage,accretioncomponents(switchyr+551:marshage+551,3));
d(1).FaceColor='none'; %#ece2f0
d(1).EdgeColor='none';

hold on
plot(0:marshage,totalaccretion(551:marshage+551)*1000,'Color','none','LineWidth',4) %#7570b3
hold on
plot(0:marshage,SLR(551:marshage+551)*1000,'Color','k','LineWidth',3,'LineStyle',':')
hold on
for i=0:15:150
    plot([i i],[22 19],'LineStyle','--','Color','k')
    hold on
end
yline(0)
xlim([0 150])
ylim([-15 22])
ax=gca;
ax.FontSize=16;

%text(3,21,'a','FontSize',16,'FontWeight','bold')
text(40,10,'total accretion','FontSize',16,'Color','none','Rotation',30) %#7570b3
%text(110,13,'SLR','FontSize',16,'Rotation',23)
%text(9,21,'TLP','FontSize',16,'Color','#7570b3')

text(5,-5.7,'Accretionary components','FontSize',16)
plot([4,85],[-6.6 -6.6],'Color','k')
text(11,19,'TLP','FontSize',16,'Rotation',90)
text(26,19,'TLP','FontSize',16,'Rotation',90)
text(41,19,'TLP','FontSize',16,'Rotation',90)

text(13,-8,'Net belowground biomass','FontSize',16)
rectangle('Position',[3 -9 8 2],'FaceColor','#ece2f0','EdgeColor','none')
text(13,-10.5,'Recycled eroded sediment','FontSize',16)
rectangle('Position',[3 -11.5 8 2],'FaceColor','#4c516d','EdgeColor','none')
text(13,-13,'External & mudflat sediment','FontSize',16)
rectangle('Position',[3 -14 8 2],'FaceColor','#838996','EdgeColor','none')
title(['\fontsize{16} a {\color[rgb]{0.46 0.44 0.70}TLP}']);
ax = gca;
ax.TitleHorizontalAlignment = 'left';

hold on
yyaxis right
plot(0:marshage,msl(551:marshage+551),'Color','#ADD8E6','LineWidth',3)
ylim([-1.55 2.45])
%%
% ----------------------------------------------------------------------------------------------------------------------------
% supplementary fig Xb; create accretion components graph 

% *keep org and min separate until the very end
%figure out percentage of Fb_org and Fb_min made up of by eroded material (saved in fluxes matrix)
s2=subaxis(1,2,2,'SpacingHoriz',0.04,'SpacingVert',0.1); %limits the padding between subplots
filename=['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Sill1_SSCReduction20_silldelay551_TLP1_TLPthick0.05_TLPfreq15_TLPdelay551_Erosion' num2str(k) '/'];
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


marshage=find(VASEAmatrix(551:701,3)==0,1)-1;
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
%ylabel('Accretion Rate (mm yr^-^1)')
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
plot(0:marshage,SLR(551:marshage+551)*1000,'Color','k','LineWidth',3,'LineStyle',':')
hold on
for i=0:15:150
    plot([i i],[22 19],'LineStyle','--','Color','k')
    hold on
end
yline(0)
xlim([0 150])
ylim([-15 22])

ax=gca;
ax.FontSize=16;

%text(3,21,'b','FontSize',16,'FontWeight','bold')
text(40,9,'total accretion','FontSize',16,'Color','#e6ab02','Rotation',18)
text(100,12,'SLR','FontSize',16,'Rotation',23)
%text(9,21,'TLP and sill','FontSize',16,'Color','#e6ab02','BackgroundColor','white')
text(11,19,'TLP','FontSize',16,'Rotation',90)
text(26,19,'TLP','FontSize',16,'Rotation',90)
text(41,19,'TLP','FontSize',16,'Rotation',90)


text(5,-5.7,'Accretionary components','FontSize',16)
plot([4,85],[-6.6 -6.6],'Color','k')

text(13,-8,'Net belowground biomass','FontSize',16)
rectangle('Position',[3 -9 8 2],'FaceColor','#ece2f0','EdgeColor','none')
text(13,-10.5,'External & mudflat sediment','FontSize',16)
rectangle('Position',[3 -11.5 8 2],'FaceColor','#838996','EdgeColor','none')
title(['\fontsize{16} b {\color[rgb]{0.9 0.67 0.01}TLP and sill}']);
ax = gca;
ax.TitleHorizontalAlignment = 'left';

%% Supplementary Fig. 4
% comparing contribution of erosion to avg accretion rates for TLP vs no-restoration
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
plot(0:marshage,totaleroded(551:marshage+551)*1000,'Color','#e7298a','LineWidth',4)
hold on
text(116,6.4,'TLP','FontSize',20,'Rotation',90)
text(131,6.4,'TLP','FontSize',20,'Rotation',90)
text(101,6.4,'TLP','FontSize',20,'Rotation',90)


subplot(1,2,2)
s2=subaxis(1,2,2,'SpacingHoriz',0.04,'SpacingVert',0.1); %limits the padding between subplots
scatter(yravgelevbelowMHW(551:701),totaleroded(551:701)*1000,'MarkerFaceColor','#e7298a','MarkerEdgeColor','none')
hold on

% yyaxis right
% scatter(0:marshage,totalerodedflux(551:marshage+551),10,'MarkerFaceColor','#e7298a','MarkerEdgeColor','none')
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
plot(0:marshage,totaleroded(551:marshage+551)*1000,'Color','#7570b3','LineWidth',4,'LineStyle','-')
ylabel({'Average annual accretion of','eroded sediment (mm yr^-^1 m^-^1)'},'HorizontalAlignment', 'center')
hold on
for i=0:15:150
    plot([i i],[7 6.4],'LineStyle','--','Color','k')
    hold on
end

xlabel('Model run year')
ax=gca;
ax.FontSize=20;
hold on
l=legend('No restoration','TLP');
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

%% SUPPLEMENTARY FIG S5
% fig a: sensitivity graph with SSC on x-axis, time to interior drown on y-axis, with different restoration strategies as 4 diff lines
figure
subplot(1,2,1)
% set up for-loop to pull matrices of various SSC files for intermediate RSLRA for no restoration
RSLRA=0.1592;
k=5.3e-09;

%sill
color='#66a61e';
for CO=[10,20,30,40,50]
    filename=['Outputs_CO' num2str(CO) '_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Sill1_SSCReduction20_silldelay551_TLP0_TLPthick0_TLPfreq0_TLPdelay551_Erosion' num2str(k) '/'];
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
    scatter(CO,threshold,110,'MarkerFaceColor',color,'MarkerEdgeColor','none')
    hold on
    color='#66a61e';
end
p1=plot(10:10:50,thresholdsmatrix(1:5),'LineWidth',3,'Color','#66a61e');
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
    scatter(CO,threshold,110,'MarkerFaceColor',color,'MarkerEdgeColor','none')
    hold on
    color='#7570b3';

end
p2=plot(10:10:50,thresholdsmatrix(1:5),'LineWidth',3,'Color','#7570b3');
hold on

%TLP and sill
color='#e6ab02';
for CO=[10,20,30,40,50]
    filename=['Outputs_CO' num2str(CO) '_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Sill1_SSCReduction20_silldelay551_TLP1_TLPthick0.05_TLPfreq15_TLPdelay551_Erosion' num2str(k) '/'];
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
    scatter(CO,threshold,110,'MarkerFaceColor',color,'MarkerEdgeColor','none')
    hold on
    color='#e6ab02';

end
p3=plot(10:10:50,thresholdsmatrix(1:5),'LineWidth',3,'Color','#e6ab02');
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
    scatter(CO,threshold,110,'MarkerFaceColor',color,'MarkerEdgeColor','none')
    hold on
    color='#e7298a';
end
p4=plot(10:10:50,thresholdsmatrix(1:5),'LineWidth',3,'Color','#e7298a','LineStyle',':');
hold on

xlabel('External suspended sediment concentration (mg L^-^1)')
ylabel('Age at which interior drowning initiates (yr)')
xticks([10 20 30 40 50])
text(11,668,'a','FontWeight','bold')
l=legend([p4,p1,p2,p3],'No restoration','Sill','TLP','TLP and sill');
l.Location='southeast';
l.EdgeColor='none';
l.Color=[0.97 0.97 0.97 0.97];

% -----------------------------------------------------------------------------------------------------------------------------------------------------------
% fig b: sensitivity graph with SSC on x-axis, time to collapse on y-axis, with different restoration strategies as 4 diff lines
subplot(1,2,2)


%sill
color='#66a61e';
for CO=[10,20,30,40,50]
    filename=['Outputs_CO' num2str(CO) '_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Sill1_SSCReduction20_silldelay551_TLP0_TLPthick0_TLPfreq0_TLPdelay551_Erosion' num2str(k) '/'];
    myVars={'elev.mat','marsh edge.mat','meansealevel.mat','AnnualAccretion.mat','SLR.mat','resiliencymetrics.mat'};
    for i=1:6
        load([filename myVars{i}]);
    end

    % define onset of interior drowning 
    marshage=find(resiliencymetrics(551:701,11)==0,1)+551;
    if isempty(marshage)==1
        marshage=701;
    end
   
    marshagematrix(1,CO/10)=marshage;
    scatter(CO,marshage,110,'MarkerFaceColor',color,'MarkerEdgeColor','none')
    hold on
    color='#66a61e';
end
p1=plot(10:10:50,marshagematrix(1:5),'LineWidth',3,'Color','#66a61e');
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
        
    marshagematrix(1,CO/10)=marshage;
    scatter(CO,marshage,110,'MarkerFaceColor',color,'MarkerEdgeColor','none')
    hold on
    color='#7570b3';
end
p2=plot(10:10:50,marshagematrix(1:5),'LineWidth',3,'Color','#7570b3')
hold on

%TLP and sill
color='#e6ab02';
for CO=[20,30,40,50]
    filename=['Outputs_CO' num2str(CO) '_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Sill1_SSCReduction20_silldelay551_TLP1_TLPthick0.05_TLPfreq15_TLPdelay551_Erosion' num2str(k) '/'];
    myVars={'elev.mat','marsh edge.mat','meansealevel.mat','AnnualAccretion.mat','SLR.mat','resiliencymetrics.mat'};
    for i=1:6
        load([filename myVars{i}]);
    end

    % define onset of interior drowning 
    marshage=find(resiliencymetrics(551:701,11)==0,1)+551;
    if isempty(marshage)==1
        marshage=701;
    end
        
    marshagematrix(1,CO/10)=marshage;
    scatter(CO,marshage,110,'MarkerFaceColor',color,'MarkerEdgeColor','none')
    hold on
    color='#e6ab02';

end
p3=plot(10:10:50,marshagematrix(1:5),'LineWidth',3,'Color','#e6ab02');
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
    end
    marshagematrix(1,CO/10)=marshage;
    scatter(CO,marshage,110,'MarkerFaceColor',color,'MarkerEdgeColor','none')
    hold on
    color='#e7298a';
end
p4=plot(10:10:50,marshagematrix(1:5),'LineWidth',3,'Color','#e7298a','LineStyle',':');
hold on

xlabel('External suspended sediment concentration (mg L^-^1)')
ylabel('Age at which marsh collapses (yr)')
l2=legend([p4,p1,p2,p3],'No restoration','Sill','TLP','TLP and sill');
l2.Location='southeast';
l2.EdgeColor='none';
l2.Color=[0.97 0.97 0.97 0.97];
xticks([10 20 30 40 50])
text(11,708,'b','FontWeight','bold')

%% SUPPLEMENTARY FIG S6
% fig a: sensitivity graph with SLRA on x-axis, time to interior drown on y-axis, with different restoration strategies as 4 diff lines
figure
subplot(1,2,1)
% set up for-loop to pull matrices of various SSC files for intermediate RSLRA for no restoration

k=7.8e-09;
CO=10;
scatterRSLRA=[0.01881,0.03536,0.08023,0.1106,0.1592];

%sill
color='#66a61e';
holder=1;
for RSLRA=[0.01881,0.03536,0.08023,0.1106,0.1592]
    filename=['Outputs_CO' num2str(CO) '_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Sill1_SSCReduction20_silldelay551_TLP0_TLPthick0_TLPfreq0_TLPdelay551_Erosion' num2str(k) '/'];
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
    scatter(RSLRA,threshold,110,'MarkerFaceColor',color,'MarkerEdgeColor','none')
    hold on
    color='#66a61e';
    holder=holder+1;
    clear elev marsh edge meansealevel AnnualAccretion SLR.mat resiliencymetrics
end
p1=plot(scatterRSLRA,thresholdsmatrix(1:5),'LineWidth',3,'Color','#66a61e');
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
    scatter(RSLRA,threshold,110,'MarkerFaceColor',color,'MarkerEdgeColor','none')
    hold on
    color='#7570b3';
    holder=holder+1;
    clear elev marsh edge meansealevel AnnualAccretion SLR.mat resiliencymetrics
end
p2=plot(scatterRSLRA,thresholdsmatrix(1:5),'LineWidth',3,'Color','#7570b3');
hold on

%TLP and sill
color='#e6ab02';
holder=1;
for RSLRA=[0.01881,0.03536,0.08023,0.1106,0.1592]
    filename=['Outputs_CO' num2str(CO) '_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Sill1_SSCReduction20_silldelay551_TLP1_TLPthick0.05_TLPfreq15_TLPdelay551_Erosion' num2str(k) '/'];
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
    scatter(RSLRA,threshold,110,'MarkerFaceColor',color,'MarkerEdgeColor','none')
    hold on
    color='#e6ab02';
    holder=holder+1;
    clear elev marsh edge meansealevel AnnualAccretion SLR.mat resiliencymetrics
end
p3=plot(scatterRSLRA,thresholdsmatrix(1:5),'LineWidth',3,'Color','#e6ab02');
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
    scatter(RSLRA,threshold,110,'MarkerFaceColor',color,'MarkerEdgeColor','none')
    hold on
    color='#e7298a';
    holder=holder+1;
    clear elev marsh edge meansealevel AnnualAccretion SLR.mat resiliencymetrics
end
p4=plot(scatterRSLRA,thresholdsmatrix(1:5),'LineWidth',3,'Color','#e7298a','LineStyle',':');
hold on

xlabel('Sea-level rise acceleration rate (mm yr^-^2)')
ylabel('Age at which interior drowning initiates (yr)')
xticks([0 0.04 0.08 0.12 0.16])
text(0.005,718,'a','FontWeight','bold')
l=legend([p4,p1,p2,p3],'No restoration','Sill','TLP','TLP and sill');
l.Location='southwest';
l.EdgeColor='none'; %grey legend outline
l.Color=[0.97 0.97 0.97 0.97];
xlim([0 0.18])


% -----------------------------------------------------------------------------------------------------------------------------------------------------------
% fig b: sensitivity graph with SLRA on x-axis, time to collapse on y-axis, with different restoration strategies as 4 diff lines
subplot(1,2,2)


%sill
color='#66a61e';
holder=1;
for RSLRA=[0.01881,0.03536,0.08023,0.1106,0.1592]
    filename=['Outputs_CO' num2str(CO) '_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Sill1_SSCReduction20_silldelay551_TLP0_TLPthick0_TLPfreq0_TLPdelay551_Erosion' num2str(k) '/'];
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
    scatter(RSLRA,marshage,110,'MarkerFaceColor',color,'MarkerEdgeColor','none')
    hold on
    color='#66a61e';
    holder=holder+1;
    clear elev marsh edge meansealevel AnnualAccretion SLR.mat resiliencymetrics
end
p1=plot(scatterRSLRA,marshagematrix(1:5),'LineWidth',3,'Color','#66a61e');
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
    scatter(RSLRA,marshage,110,'MarkerFaceColor',color,'MarkerEdgeColor','none')
    hold on
    color='#7570b3';
    holder=holder+1;
    clear elev marsh edge meansealevel AnnualAccretion SLR.mat resiliencymetrics
end
p2=plot(scatterRSLRA,marshagematrix(1:5),'LineWidth',3,'Color','#7570b3')
hold on

%TLP and sill
color='#e6ab02';
holder=1;
for RSLRA=[0.01881,0.03536,0.08023,0.1106,0.1592]
    filename=['Outputs_CO' num2str(CO) '_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Sill1_SSCReduction20_silldelay551_TLP1_TLPthick0.05_TLPfreq15_TLPdelay551_Erosion' num2str(k) '/'];
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
    scatter(RSLRA,marshage,110,'MarkerFaceColor',color,'MarkerEdgeColor','none')
    hold on
    color='#e6ab02';
    holder=holder+1;
    clear elev marsh edge meansealevel AnnualAccretion SLR.mat resiliencymetrics
end
p3=plot(scatterRSLRA,marshagematrix(1:5),'LineWidth',3,'Color','#e6ab02');
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
    scatter(RSLRA,marshage,110,'MarkerFaceColor',color,'MarkerEdgeColor','none')
    hold on
    color='#e7298a';
    holder=holder+1;
    clear elev marsh edge meansealevel AnnualAccretion SLR.mat resiliencymetrics
end
p4=plot(scatterRSLRA,marshagematrix(1:5),'LineWidth',3,'Color','#e7298a','LineStyle',':');
hold on

xlabel('Sea-level rise acceleration rate (mm yr^-^2)')
ylabel('Age at which marsh collapses (yr)')
l2=legend([p4,p1,p2,p3],'No restoration','Sill','TLP','TLP and sill');
l2.Location='southwest';
l2.EdgeColor='none'; %grey legend outline
l2.Color=[0.97 0.97 0.97 0.97];
xticks([0 0.04 0.08 0.12 0.16])
xlim([0 0.18])
text(0.005,708,'b','FontWeight','bold')

