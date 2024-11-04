function [resiliencymetrics,AnnualAccretion]=MARindex(amp,msl,SLR,startyear,endyear,Marsh_edge,elevation,realmarshend,B,bgb,Dmax)

%%function uses Raposa et al. (2016) to 

%COMMENT OUT WHEN RUNNING INSIDE FULL TRANSECT.M
% load('elev.mat')
% load('meansealevel.mat')
% load('real marsh end.mat') 
% load('SLR.mat') 
% load('marsh edge.mat') 
% load('bgb.mat') 
% startyear=551; %
% endyear=651; %
% amp=0.5; %
% B=length(elev); %
% elevation=elev; %
% tr=1;%
% amp=tr/2;%
% Dmax=(0.237*tr)-0.092+amp;%
%bgb(1:startyear-1,5001:5500)=bgb_25;%
%SLR(1:300)=0.001; %1 mm yr-1 SLR for first 399 years of spinup %
%SLR(400:startyear)=SLR(startyear); %add SLR rate from input value %

mhw=(msl+amp); % add tidal amplitude to msl matrix (with yearly values) to get mean high water level for every year

%make empty matrices for start
resiliencymetrics=zeros(endyear,11); 
totalmarshwidth=zeros(1,endyear);
vegmarshwidth=zeros(1,endyear);
PerMarBelBMaxDepth=zeros(1,endyear); %1st MAR (marsh resilience) metric
ElevAvg=zeros(1,endyear); %feed into 5th, 6th, 7th and 8th MAR metric
ElevChange=zeros(1,endyear); %feed into 5th, 6th, 7th and 8th MAR
AvgElevAboveMSL=zeros(1,endyear); %8th MAR metric
AnnualAccretion=zeros(endyear,B);

ElevAvg(1)=mean(elevation(bgb(1,Marsh_edge(1):5500)>0)); %average marsh elev for this year; define marsh as marsh edge to 5500
for i=2:endyear %start for-loop
    AnnualAccretion(i,:)=elevation(i,:)-elevation(i-1,:); %[m yr-1] finds elevation change from previous year for every cell in domain
    ElevAvg(i)=mean(elevation(i,bgb(i,Marsh_edge(i):5500)>0)); %[m] average marsh elev for this year; define marsh as marsh edge to 5500
    tempmarshcells=find(bgb(i,Marsh_edge(i):5500)>0); %creates vector of all vegetated marsh cells
    ElevChange(i)=mean(AnnualAccretion(i,tempmarshcells+Marsh_edge(i)-1)); %[m yr-1] calculates average annual elevation change across the vegetated marsh
end

for yr=startyear:endyear %starts for-loop for model run years
    tempmarshcells=find(bgb(yr,Marsh_edge(yr):5500)>0); %creates vector of all vegetated marsh cells
    tempmarshcells=tempmarshcells+Marsh_edge(yr)-1; %aligns vector with marsh edge position  
    totalmarshwidth(yr)=5500-Marsh_edge(yr); %finds total marsh width between current edge and initial marsh-forest boundary
    vegmarshwidth(yr)=length(tempmarshcells); %[m] finds total vegetated marsh width
    DepthBelMHW=zeros(1,totalmarshwidth(yr));  %creates vector for holding depth below mean high water values
    DepthBelMHW=mhw(yr)-elevation(yr,Marsh_edge(yr):5500); % [m] for each cell in marsh, gives the depth below mhw
    PerMarBelBMaxDepth=length(find(DepthBelMHW>=(Dmax/2)))/(totalmarshwidth(yr)+1)*100;%[%] calculates the percentage of marsh elevations below BMax (maximum biomass) depth below mean high water
    ThirdsMarshElev=(max(elevation(yr,Marsh_edge(yr):5500))-min(elevation(yr,Marsh_edge(yr):5500))/3)+min(elevation(yr,Marsh_edge(yr):5500)); %divides plant elevation distribution into thirds
    PerMarInLowestThird=length(find(elevation(yr,Marsh_edge(yr):5500)<=ThirdsMarshElev))/(totalmarshwidth(yr)+1)*100; %[%] calculates the percentage of marsh cells that have depths in thelowest third of depth distribution
    yrMarshBGB=bgb(yr,Marsh_edge(yr):5500); %findssthe belowground biomass of every marsh cell for that year
    ThirdsMarshProd=((max(yrMarshBGB)-min(yrMarshBGB))/3)+min(yrMarshBGB); %divides the distribution into thirds
    PerMarInLowestThirdProd=length(find(yrMarshBGB<=ThirdsMarshProd))/vegmarshwidth(yr)*100; %[%] calculates the percentage of marsh cells that have belowground biomass in the lowest third productivity range
    
    resiliencymetrics(yr,1)=PerMarBelBMaxDepth; % first resiliency metric
    resiliencymetrics(yr,2)=PerMarInLowestThird; % second resiliency metric
    resiliencymetrics(yr,3)=PerMarInLowestThirdProd; % third resiliency metric
  
    %Scores the first plant MAR metric
    if PerMarBelBMaxDepth<=20
        PlantScore1=5;
    elseif PerMarBelBMaxDepth>20 && PerMarBelBMaxDepth<=40
        PlantScore1=4;
    elseif PerMarBelBMaxDepth>40 && PerMarBelBMaxDepth<=60
        PlantScore1=3;
    elseif PerMarBelBMaxDepth>60 && PerMarBelBMaxDepth<=80
        PlantScore1=2;
    else 
        PlantScore1=1;
    end
    
    %Scores the second plant MAR metric
    if PerMarInLowestThird<=20
        PlantScore2=5;
    elseif PerMarInLowestThird>20 && PerMarInLowestThird<=40
        PlantScore2=4;
    elseif PerMarInLowestThird>40 && PerMarInLowestThird<=60
        PlantScore2=3;
    elseif PerMarInLowestThird>60 && PerMarInLowestThird<=80
        PlantScore2=2;
    else 
        PlantScore2=1;
    end

        %Scores the third plant MAR metric
    if PerMarInLowestThirdProd<=20
        PlantScore3=5;
    elseif PerMarInLowestThirdProd>20 && PerMarInLowestThirdProd<=40
        PlantScore3=4;
    elseif PerMarInLowestThirdProd>40 && PerMarInLowestThirdProd<=60
        PlantScore3=3;
    elseif PerMarInLowestThirdProd>60 && PerMarInLowestThirdProd<=80
        PlantScore3=2;
    else 
        PlantScore3=1;
    end
    
    %Averages scores and puts in new resiliencymetrix value
    PlantScore=(PlantScore1+PlantScore2+PlantScore3)/3;
    
    resiliencymetrics(yr,4)=mean(ElevChange(yr-30:yr))/mean(SLR(yr-30:yr)); %calculates long-term (30-year) accretion to SLR ratio
    resiliencymetrics(yr,5)=mean(ElevChange(yr-5:yr))/mean(SLR(yr-5:yr)); %calculates short-term (5-year) accretion to SLR ratio
    PerMarWithARdeficit=length(find(AnnualAccretion(yr,Marsh_edge(yr):5500)<SLR(yr)))/(5500-Marsh_edge(yr)+1)*100; %claculates percentage of marsh with accretion rate deficit
    resiliencymetrics(yr,6)=PerMarWithARdeficit; %6th MAR metric 

    %Scores 1st vertical resilience MAR metric
    if resiliencymetrics(yr,4)>=2
        VertScore3=5;
    elseif resiliencymetrics(yr,4)>1.5 && resiliencymetrics(yr,4)<2
        VertScore3=4;
    elseif resiliencymetrics(yr,4)>1 && resiliencymetrics(yr,4)<=1.5
        VertScore3=3;
    elseif resiliencymetrics(yr,4)>0.5 && resiliencymetrics(yr,4)<=1
        VertScore3=2;
    else 
        VertScore3=1;
    end
    
    %Scores 2nd vertical resilience MAR metric
    if resiliencymetrics(yr,5)>=2
        VertScore4=5;
    elseif resiliencymetrics(yr,5)>1.5 && resiliencymetrics(yr,5)<2
        VertScore4=4;
    elseif resiliencymetrics(yr,5)>1 && resiliencymetrics(yr,5)<=1.5
        VertScore4=3;
    elseif resiliencymetrics(yr,5)>0.5 && resiliencymetrics(yr,5)<=1
        VertScore4=2;
    else 
        VertScore4=1;
    end
    
    
    %Scores 3rd vertical resilience MAR metric
    if PerMarWithARdeficit<=20
        VertScore5=5;
    elseif PerMarWithARdeficit>20 && PerMarWithARdeficit<=40
        VertScore5=4;
    elseif PerMarWithARdeficit>40 && PerMarWithARdeficit<=60
        VertScore5=3;
    elseif PerMarWithARdeficit>60 && PerMarWithARdeficit<=80
        VertScore5=2;
    else 
        VertScore5=1;
    end

    VertScore=(VertScore3+VertScore4+VertScore5)/3; %averages the three vertical scores 

    resiliencymetrics(yr,7)=vegmarshwidth(yr); %7th MAR metric

    %Scores single lateral resilience MAR metric (marsh width)
    if vegmarshwidth(yr)>=500
        LatScore=5;
    elseif (vegmarshwidth(yr)>=400) && (vegmarshwidth(yr)<500)
        LatScore=4;
    elseif (vegmarshwidth(yr)>=300) && (vegmarshwidth(yr)<400)
        LatScore=3;
    elseif (vegmarshwidth(yr)>=200) && (vegmarshwidth(yr)<300)
        LatScore=2;
    else 
        LatScore=1;
    end
    
    resiliencymetrics(yr,8)=PlantScore;
    resiliencymetrics(yr,9)=VertScore;
    resiliencymetrics(yr,10)=LatScore;

    resiliencymetrics(yr,11)=PlantScore+VertScore+LatScore; %cumulative of all (15 is maximum)
    if vegmarshwidth(yr)==0 % if no vegetated marsh, designates remaining resilicney score metrics
        resiliencymetrics(yr:endyear,1:3)=100;
        resiliencymetrics(yr:endyear,6)=100;
        resiliencymetrics(yr:endyear,8:11)=0; 
        break %ends for-loop
    end
   
end
save('resiliencymetrics.mat','resiliencymetrics')
save('AnnualAccretion.mat','AnnualAccretion')