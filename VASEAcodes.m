%%streamline code for saving for VASEA
Nuptake=zeros(150,500);
totalNuptake=zeros(150,1);
yearlyCaccumulation=zeros(150,1);
ElevAvg=zeros(150,1);
for CO=10
%for runstartRSLR=3
    for RSLRA=[0.1592,0.1106,0.08023,0.03536,0.01881]
    for RSLR=3
    for Erosion=[1.65e-09,2.82e-09,5.3e-09,7.8e-09,9.95e-09]
            %for Erosion=2.984e-09
    filename = ['Outputs_CO' num2str(CO) '_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR' num2str(RSLR) '_Sill0_SSCReduction0_silldelay551_TLP0_TLPthick0_TLPfreq0_TLPdelay551_Erosion' num2str(Erosion) '\'];
    %filename = ['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Sill0_SSCReduction0_silldelay551_TLP0_TLPthick0_TLPfreq0_TLPdelay551_Erosion' num2str(Erosion) '\'];
    %filename = ['Outputs_CO10_RSLRA0.221_spinupinitialRSLR1_runstartRSLR3_Sill0_SSCReduction0_TLP0_TLPthick0_TLPfreq0_Erosion2.988E-09\'];
    myVars={'bgb','agb','organic deposition','real marsh end','marsh edge','elev','meansealevel','resiliencymetrics','totalmarshOMto5500','SLR'}; 
    for i=1:10
        load([filename myVars{i}]);
    end
    for yr=1:701
            for j=5000:5500
                Nuptake(yr,j)=(154*agb(yr,j)/1000)*365; %Hill et al. 2019; reported in mg N m-2 d-1 but by dividing by 1000 and multiplying by 365, convert to g N m-2 yr-1
            end
                totalNuptake(yr,1)=sum(Nuptake(yr,:)); %sum across the marsh, for a total N uptake rate, reported in g N yr-1
                yearlyCaccumulation(yr,1)=(sum(bgb(yr,Marsh_edge(yr):5500))+sum(organic_dep_alloch(yr,Marsh_edge(yr):5500)))*0.4; % multiply by 0.4 to convert from OM to OC
                ElevAvg(yr)=mean(elev(yr,Marsh_edge(yr):realmarshend(yr)));
    end
    VASEAmatrix_control = zeros(701,5);
    VASEAmatrix_control(:,1)=msl;
    VASEAmatrix_control(:,2)=ElevAvg;
    VASEAmatrix_control(:,3)=resiliencymetrics(:,7);
    VASEAmatrix_control(:,4)=yearlyCaccumulation;
    VASEAmatrix_control(:,5)=totalNuptake;
    VASEAmatrix_control(:,6)=resiliencymetrics(:,11);
    VASEAmatrix_control(:,7)=totalmarshOMto5500;
    VASEAmatrix_control(:,8)=SLR;
    elev_control=elev;
    resiliencymetrics_control=resiliencymetrics;
    

    outputfilename=['CO' num2str(CO) '_RSLRA' num2str(RSLRA) '_Erosion' num2str(Erosion) '/']; 
    % %outputfilename=['CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Erosion' num2str(Erosion) '/'];
    if exist(outputfilename,'file') ~= 7
         mkdir(outputfilename) %If filename does not exist, create a folder to save the output variable to
    end
     VASEAmatrix=VASEAmatrix_control;
    % writematrix(VASEAmatrix_control,['CO' num2str(CO) '_RSLR' num2str(RSLR) '_Erosion' num2str(Erosion) '\VASEAmatrix_control.csv']);
    % writematrix(elev_control,['CO' num2str(CO) '_RSLR' num2str(RSLR) '_Erosion' num2str(Erosion) '\elev_control.csv']);
    % writematrix(VASEAmatrix_control);
     save([filename 'VASEAmatrix.mat'],'VASEAmatrix');
     save([outputfilename 'VASEAmatrix_control.mat'],'VASEAmatrix_control');
     save([outputfilename 'elev_control.mat'],'elev_control');
     save([outputfilename 'resiliencymetrics_control.mat'],'resiliencymetrics_control');
     clear("agb","bgb","ElevAvg","elev",'Marsh_edge',"msl","Nuptake","organic_dep_alloch","organic_dep_autoch","realmarshend","resiliencymetrics","totalmarshOMto5500","totalNuptake","VASEAmatrix_control","yearlyCaccumulation",'SLR',"VASEAmatrix")
%end
end
end
    end
end
% clear all

Nuptake=zeros(150,500);
totalNuptake=zeros(150,1);
yearlyCaccumulation=zeros(150,1);
ElevAvg=zeros(150,1);
for CO=10
%for runstartRSLR=3:4
    for RSLRA=[0.1592,0.1106,0.08023,0.03536,0.01881]
    for Erosion=[1.65e-09,2.82e-09,5.3e-09,7.8e-09,9.95e-09]
    for RSLR=3
            %for Erosion=2.984e-09 
    filename = ['Outputs_CO' num2str(CO) '_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR' num2str(RSLR) '_Sill1_SSCReduction20_silldelay551_TLP0_TLPthick0_TLPfreq0_TLPdelay551_Erosion' num2str(Erosion) '\'];
    %filename = ['Outputs_CO10_RSLRA0.221_spinupinitialRSLR1_runstartRSLR3_Sill0_SSCReduction0_TLP0_TLPthick0_TLPfreq0_Erosion2.988E-09\'];
    myVars={'bgb','agb','organic deposition','real marsh end','marsh edge','elev','meansealevel','resiliencymetrics','totalmarshOMto5500','SLR'}; 
    for i=1:10
        load([filename myVars{i}]);
    end
    for yr=1:701
            for j=5000:5500
                Nuptake(yr,j)=(154*agb(yr,j)/1000)*365; %Hill et al. 2019; reported in mg N m-2 d-1 but by dividing by 1000 and multiplying by 365, convert to g N m-2 yr-1
            end
                totalNuptake(yr,1)=sum(Nuptake(yr,:)); %sum across the marsh, for a total N uptake rate, reported in g N yr-1
                yearlyCaccumulation(yr,1)=(sum(bgb(yr,Marsh_edge(yr):5500))+sum(organic_dep_alloch(yr,Marsh_edge(yr):5500)))*0.4; % multiply by 0.4 to convert from OM to OC
                ElevAvg(yr)=mean(elev(yr,Marsh_edge(yr):realmarshend(yr)));
    end
    VASEAmatrix_sill = zeros(701,5);
    VASEAmatrix_sill(:,1)=msl;
    VASEAmatrix_sill(:,2)=ElevAvg;
    VASEAmatrix_sill(:,3)=resiliencymetrics(:,7);
    VASEAmatrix_sill(:,4)=yearlyCaccumulation;
    VASEAmatrix_sill(:,5)=totalNuptake;
    VASEAmatrix_sill(:,6)=resiliencymetrics(:,11);
    VASEAmatrix_sill(:,7)=totalmarshOMto5500;
    VASEAmatrix_sill(:,8)=SLR;
    elev_sill=elev;
    resiliencymetrics_sill=resiliencymetrics;

    %save(['CO' num2str(CO) '_RSLR' num2str(RSLR) '_Erosion' num2str(Erosion) '\VASEAmatrix_control.mat'],'VASEAmatrix_control');
    %save(['CO' num2str(CO) '_RSLR' num2str(RSLR) '_Erosion' num2str(Erosion) '\elev_control.mat'],'elev');
    %save(['CO' num2str(CO) '_RSLR' num2str(RSLR) '_Erosion' num2str(Erosion) '\resiliencymetrics_control.mat'],'resiliencymetrics');             
    outputfilename=['CO' num2str(CO) '_RSLRA' num2str(RSLRA) '_Erosion' num2str(Erosion) '/']; 
    %outputfilename=['CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Erosion' num2str(Erosion) '/'];
    if exist(outputfilename,'file') ~= 7
        mkdir(outputfilename) %If filename does not exist, create a folder to save the output variable to
    end
    VASEAmatrix=VASEAmatrix_sill;
%     writematrix(VASEAmatrix_sill,['CO' num2str(CO) '_RSLR' num2str(RSLR) '_Erosion' num2str(Erosion) '\VASEAmatrix_sill.csv']);
%     writematrix(elev_sill,['CO' num2str(CO) '_RSLR' num2str(RSLR) '_Erosion' num2str(Erosion) '\elev_sill.csv']);
    save([filename 'VASEAmatrix.mat'],'VASEAmatrix');
    save([outputfilename 'VASEAmatrix_sill.mat'],'VASEAmatrix_sill');
    save([outputfilename 'elev_sill.mat'],'elev_sill');
    save([outputfilename 'resiliencymetrics_sill.mat'],'resiliencymetrics_sill');
    clear("agb","bgb","ElevAvg","elev",'Marsh_edge',"msl","Nuptake","organic_dep_alloch","organic_dep_autoch","realmarshend",'SLR',"resiliencymetrics","totalmarshOMto5500","totalNuptake","VASEAmatrix_control","yearlyCaccumulation",'VASEAmatrix')
    end
    end
    end
    end
%clear all


Nuptake=zeros(150,500);
totalNuptake=zeros(150,1);
yearlyCaccumulation=zeros(150,1);
ElevAvg=zeros(150,1);
for CO=10
%for runstartRSLR=3:4
    for RSLRA=[0.1592,0.1106,0.08023,0.03536,0.01881]
        for Erosion=[1.65e-09,2.82e-09,5.3e-09,7.8e-09,9.95e-09]
%     for RSLR=3
%             for Erosion=2.984e-09 
    filename = ['Outputs_CO' num2str(CO) '_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR' num2str(RSLR) '_Sill0_SSCReduction0_silldelay551_TLP1_TLPthick0.05_TLPfreq15_TLPdelay551_Erosion' num2str(Erosion) '\'];
    %filename = ['Outputs_CO10_RSLRA0.221_spinupinitialRSLR1_runstartRSLR3_Sill0_SSCReduction0_TLP0_TLPthick0_TLPfreq0_Erosion2.988E-09\'];
    myVars={'bgb','agb','organic deposition','real marsh end','marsh edge','elev','meansealevel','resiliencymetrics','totalmarshOMto5500','SLR'}; 
    for i=1:10
        load([filename myVars{i}]);
    end
    for yr=1:701
            for j=5000:5500
                Nuptake(yr,j)=(154*agb(yr,j)/1000)*365; %Hill et al. 2019; reported in mg N m-2 d-1 but by dividing by 1000 and multiplying by 365, convert to g N m-2 yr-1
            end
                totalNuptake(yr,1)=sum(Nuptake(yr,:)); %sum across the marsh, for a total N uptake rate, reported in g N yr-1
                yearlyCaccumulation(yr,1)=(sum(bgb(yr,Marsh_edge(yr):5500))+sum(organic_dep_alloch(yr,Marsh_edge(yr):5500)))*0.4; % multiply by 0.4 to convert from OM to OC
                ElevAvg(yr)=mean(elev(yr,Marsh_edge(yr):realmarshend(yr)));
    end
    yearlyCaccumulation(551:15:701,1)=NaN;
    [yearlyCaccumulation,TF]=fillmissing(yearlyCaccumulation,'linear');
    VASEAmatrix_TLP = zeros(701,5);
    VASEAmatrix_TLP(:,1)=msl;
    VASEAmatrix_TLP(:,2)=ElevAvg;
    VASEAmatrix_TLP(:,3)=resiliencymetrics(:,7);
    VASEAmatrix_TLP(:,4)=yearlyCaccumulation;
    VASEAmatrix_TLP(:,5)=totalNuptake;
    VASEAmatrix_TLP(:,6)=resiliencymetrics(:,11);
    VASEAmatrix_TLP(:,7)=totalmarshOMto5500;
    VASEAmatrix_TLP(:,8)=SLR;
    elev_TLP=elev;
    resiliencymetrics_TLP=resiliencymetrics;
    

    %writematrix(VASEAmatrix_control,['CO' num2str(CO) '_RSLR' num2str(RSLR) '_Erosion' num2str(Erosion) '\VASEAmatrix_control.csv']);
    %writematrix(elev_control,['CO' num2str(CO) '_RSLR' num2str(RSLR) '_Erosion' num2str(Erosion) '\elev_control.csv']);
    %save(['CO' num2str(CO) '_RSLR' num2str(RSLR) '_Erosion' num2str(Erosion) '\VASEAmatrix_control.mat'],'VASEAmatrix_control');
    %save(['CO' num2str(CO) '_RSLR' num2str(RSLR) '_Erosion' num2str(Erosion) '\elev_control.mat'],'elev');
    %save(['CO' num2str(CO) '_RSLR' num2str(RSLR) '_Erosion' num2str(Erosion) '\resiliencymetrics_control.mat'],'resiliencymetrics');             
    %outputfilename=['CO10_RSLRA0.0707_spinupinitialRSLR1_runstartRSLR' num2str(runstartRSLR) '_Erosion1/']; 
    %outputfilename=['CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Erosion' num2str(Erosion) '/'];
    outputfilename=['CO' num2str(CO) '_RSLRA' num2str(RSLRA) '_Erosion' num2str(Erosion) '/']; 
    if exist(outputfilename,'file') ~= 7
        mkdir(outputfilename) %If filename does not exist, create a folder to save the output variable to
    end
%     writematrix(VASEAmatrix_TLP,['CO' num2str(CO) '_RSLR' num2str(RSLR) '_Erosion' num2str(Erosion) '\VASEAmatrix_TLP.csv']);
%     writematrix(elev_TLP,['CO' num2str(CO) '_RSLR' num2str(RSLR) '_Erosion' num2str(Erosion) '\elev_TLP.csv']);
    VASEAmatrix=VASEAmatrix_TLP;
    save([filename 'VASEAmatrix.mat'],'VASEAmatrix');
    save([outputfilename 'VASEAmatrix_TLP.mat'],'VASEAmatrix_TLP');
    save([outputfilename 'elev_TLP.mat'],'elev_TLP');
    save([outputfilename 'resiliencymetrics_TLP.mat'],'resiliencymetrics_TLP');
    clear("agb","bgb","ElevAvg","elev",'Marsh_edge',"msl","Nuptake","organic_dep_alloch","organic_dep_autoch","realmarshend",'SLR',"resiliencymetrics","totalmarshOMto5500","totalNuptake","VASEAmatrix_control","yearlyCaccumulation",'VASEAmatrix')
    end
    end
end
%clear all


Nuptake=zeros(150,500);
totalNuptake=zeros(150,1);
yearlyCaccumulation=zeros(150,1);
ElevAvg=zeros(150,1);
for CO=10
%for runstartRSLR=3:4
    for RSLRA = [0.1592,0.1106,0.08023,0.03536,0.01881]
    for Erosion=[1.65e-09,2.82e-09,5.3e-09,7.8e-09,9.95e-09]
        %for Erosion=2.984e-09
    for RSLR=3
    filename = ['Outputs_CO' num2str(CO) '_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR' num2str(RSLR) '_Sill1_SSCReduction20_silldelay551_TLP1_TLPthick0.05_TLPfreq15_TLPdelay551_Erosion' num2str(Erosion) '\'];
    %filename = ['Outputs_CO10_RSLRA0.221_spinupinitialRSLR1_runstartRSLR3_Sill0_SSCReduction0_TLP0_TLPthick0_TLPfreq0_Erosion2.988E-09\'];    
    myVars={'bgb','agb','organic deposition','real marsh end','marsh edge','elev','meansealevel','resiliencymetrics','totalmarshOMto5500','SLR'}; 
    for i=1:10
        load([filename myVars{i}]);
    end
    for yr=1:701
            for j=5000:5500
                Nuptake(yr,j)=(154*agb(yr,j)/1000)*365; %Hill et al. 2019; reported in mg N m-2 d-1 but by dividing by 1000 and multiplying by 365, convert to g N m-2 yr-1
            end
                totalNuptake(yr,1)=sum(Nuptake(yr,:)); %sum across the marsh, for a total N uptake rate, reported in g N yr-1
                yearlyCaccumulation(yr,1)=(sum(bgb(yr,Marsh_edge(yr):5500))+sum(organic_dep_alloch(yr,Marsh_edge(yr):5500)))*0.4; % multiply by 0.4 to convert from OM to OC
                ElevAvg(yr)=mean(elev(yr,Marsh_edge(yr):realmarshend(yr)));
    end
    yearlyCaccumulation(551:15:701,1)=NaN;
    [yearlyCaccumulation,TF]=fillmissing(yearlyCaccumulation,'linear');
    VASEAmatrix_TLPandsill = zeros(701,5);
    VASEAmatrix_TLPandsill(:,1)=msl;
    VASEAmatrix_TLPandsill(:,2)=ElevAvg;
    VASEAmatrix_TLPandsill(:,3)=resiliencymetrics(:,7);
    VASEAmatrix_TLPandsill(:,4)=yearlyCaccumulation;
    VASEAmatrix_TLPandsill(:,5)=totalNuptake;
    VASEAmatrix_TLPandsill(:,6)=resiliencymetrics(:,11);
    VASEAmatrix_TLPandsill(:,7)=totalmarshOMto5500;
    VASEAmatrix_TLPandsill(:,8)=SLR;
    elev_TLPandsill=elev;
    resiliencymetrics_TLPandsill=resiliencymetrics;
    
    outputfilename=['CO' num2str(CO) '_RSLRA' num2str(RSLRA) '_Erosion' num2str(Erosion) '/']; 
    %writematrix(VASEAmatrix_control,['CO' num2str(CO) '_RSLR' num2str(RSLR) '_Erosion' num2str(Erosion) '\VASEAmatrix_control.csv']);
    %writematrix(elev_control,['CO' num2str(CO) '_RSLR' num2str(RSLR) '_Erosion' num2str(Erosion) '\elev_control.csv']);
    %save(['CO' num2str(CO) '_RSLR' num2str(RSLR) '_Erosion' num2str(Erosion) '\VASEAmatrix_control.mat'],'VASEAmatrix_control');
    %save(['CO' num2str(CO) '_RSLR' num2str(RSLR) '_Erosion' num2str(Erosion) '\elev_control.mat'],'elev');
    %save(['CO' num2str(CO) '_RSLR' num2str(RSLR) '_Erosion' num2str(Erosion) '\resiliencymetrics_control.mat'],'resiliencymetrics');             
    %outputfilename=['CO10_RSLRA0.0707_spinupinitialRSLR1_runstartRSLR' num2str(runstartRSLR) '_Erosion1/']; 
    %outputfilename=['CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR3_Erosion' num2str(Erosion) '/'];
    if exist(outputfilename,'file') ~= 7
        mkdir(outputfilename) %If filename does not exist, create a folder to save the output variable to
    end
%     writematrix(VASEAmatrix_TLPandsill,['CO' num2str(CO) '_RSLR' num2str(RSLR) '_Erosion' num2str(Erosion) '\VASEAmatrix_TLPandsill.csv']);
%     writematrix(elev_TLPandsill,['CO' num2str(CO) '_RSLR' num2str(RSLR) '_Erosion' num2str(Erosion) '\elev_TLPandsill.csv']);
    VASEAmatrix=VASEAmatrix_TLPandsill;
    save([filename 'VASEAmatrix.mat'],'VASEAmatrix');
    save([outputfilename 'VASEAmatrix_TLPandsill.mat'],'VASEAmatrix_TLPandsill');
    save([outputfilename 'elev_TLPandsill.mat'],'elev_TLPandsill');
    save([outputfilename 'resiliencymetrics_TLPandsill.mat'],'resiliencymetrics_TLPandsill');
    clear("agb","bgb","ElevAvg","elev",'Marsh_edge',"msl","Nuptake","organic_dep_alloch","organic_dep_autoch","realmarshend",'SLR',"resiliencymetrics","totalmarshOMto5500","totalNuptake","VASEAmatrix_control","yearlyCaccumulation",'VASEAmatrix')
    end
    end
    end
end
%clear all
%%
save('CO10_RSLR4_Erosion4_TLP_CAR.mat',"CO10_RSLR4_Erosion4_TLP_CAR")
            % figure
            % plot(1:100,CO10_RSLR4_Erosion4_TLP_CAR(1,1:100))
            % xlabel('Model Run Year')
            % ylabel('C accumulation rate (g C yr^-^1)')


figure
plot(0:100,CO10_RSLR4_Erosion4_control)

figure
plot(1:100,yearlyCaccumulation(1,1:100))

plot(0:100,totalNuptake(551:701));
            save('totalNuptake.mat',"totalNuptake") %MBB do i need a save here? 

CO10_RSLR4_Erosion4_control_CAR=yearlyCaccumulation;
save('CO10_RSLR4_Erosion4_control_CAR.mat',"CO10_RSLR4_Erosion4_control_CAR")
%% new code for comparing sill delays
yearlyCaccumulation=zeros(100,1);
ElevAvg=zeros(100,1);
VASEAmatrix_sill=zeros(701,50);
for runstartRSLR=3
    for RSLRA=0.0707
    for delay=551:10:651
    filename = ['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR' num2str(runstartRSLR) '_Sill1_SSCReduction20_silldelay' num2str(delay) '_TLP0_TLPthick0_TLPfreq0_TLPdelay551_Erosion2.988E-09\'];
    myVars={'bgb','agb','organic deposition','real marsh end','marsh edge','elev','meansealevel','resiliencymetrics','totalmarshOMto5500','SLR'}; 
    for i=1:10
        load([filename myVars{i}]);
    end
    for yr=1:701
       yearlyCaccumulation(yr,1)=sum(bgb(yr,Marsh_edge(yr):5500))+sum(organic_dep_alloch(yr,Marsh_edge(yr):5500))*0.4; % multiply by 0.4 to convert from OM to OC
       ElevAvg(yr)=mean(elev(yr,Marsh_edge(yr):realmarshend(yr)));
    end
    VASEAmatrix(:,delay-550)=msl;
    VASEAmatrix(:,delay-549)=ElevAvg;
    VASEAmatrix(:,delay-548)=resiliencymetrics(:,7);
    VASEAmatrix(:,delay-547)=yearlyCaccumulation;
    VASEAmatrix(:,delay-546)=resiliencymetrics(:,11);
    VASEAmatrix(:,delay-545)=totalmarshOMto5500;
    VASEAmatrix(:,delay-544)=SLR;
    resiliencymetrics_sill((1000*(delay-551)/10)+1:(1000*(delay-551)/10)+701,:)=resiliencymetrics;
    elev_sill((1000*(delay-551)/10)+1:(1000*(delay-551)/10)+701,:)=elev;
    outputfilename=['CO10_RSLRA' num2str(RSLRA) '_rampedstartRSLR' num2str(runstartRSLR) '_Erosion1_delaysill/'];
    if exist(outputfilename,'file') ~= 7
        mkdir(outputfilename) %If filename does not exist, create a folder to save the output variable to
    end
    save([outputfilename 'VASEAmatrix_silldelay.mat'],'VASEAmatrix');
    save([outputfilename 'elev_silldelay.mat'],'elev_sill');
    save([outputfilename 'resiliencymetrics_silldelay.mat'],'resiliencymetrics_sill');
    %clear("agb","bgb","ElevAvg","elev",'Marsh_edge',"msl","organic_dep_alloch","organic_dep_autoch","realmarshend","resiliencymetrics","totalmarshOMto5500","VASEAmatrix","yearlyCaccumulation",'SLR')
    %end
    end
    end
end
%% new code for comparing TLP delays
yearlyCaccumulation=zeros(100,1);
ElevAvg=zeros(100,1);
VASEAmatrix=zeros(701,50);
for runstartRSLR=3
    for RSLRA=0.0707
    for delay=551:10:651
    filename = ['Outputs_CO10_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR' num2str(runstartRSLR) '_Sill0_SSCReduction0_silldelay551_TLP1_TLPthick0.05_TLPfreq15_TLPdelay' num2str(delay) '_Erosion2.988E-09\'];
    myVars={'bgb','agb','organic deposition','real marsh end','marsh edge','elev','meansealevel','resiliencymetrics','totalmarshOMto5500','SLR'}; 
    for i=1:10
        load([filename myVars{i}]);
    end
    for yr=1:701
       yearlyCaccumulation(yr,1)=sum(bgb(yr,Marsh_edge(yr):5500))+sum(organic_dep_alloch(yr,Marsh_edge(yr):5500))*0.4; % multiply by 0.4 to convert from OM to OC
       ElevAvg(yr)=mean(elev(yr,Marsh_edge(yr):realmarshend(yr)));
    end
    VASEAmatrix(:,delay-550)=msl;
    VASEAmatrix(:,delay-549)=ElevAvg;
    VASEAmatrix(:,delay-548)=resiliencymetrics(:,7);
    VASEAmatrix(:,delay-547)=yearlyCaccumulation;
    VASEAmatrix(:,delay-546)=resiliencymetrics(:,11);
    VASEAmatrix(:,delay-545)=totalmarshOMto5500;
    VASEAmatrix(:,delay-544)=SLR;
    resiliencymetrics_TLP((1000*(delay-551)/10)+1:(1000*(delay-551)/10)+701,:)=resiliencymetrics;
    elev_TLP((1000*(delay-551)/10)+1:(1000*(delay-551)/10)+701,:)=elev;
    outputfilename=['CO10_RSLRA' num2str(RSLRA) '_rampedstartRSLR' num2str(runstartRSLR) '_Erosion1_delayTLP/'];
    if exist(outputfilename,'file') ~= 7
        mkdir(outputfilename) %If filename does not exist, create a folder to save the output variable to
    end
    save([outputfilename 'VASEAmatrix_TLPdelay.mat'],'VASEAmatrix');
    save([outputfilename 'elev_TLPdelay.mat'],'elev_TLP');
    save([outputfilename 'resiliencymetrics_TLPdelay.mat'],'resiliencymetrics_TLP');
    %clear("agb","bgb","ElevAvg","elev",'Marsh_edge',"msl","organic_dep_alloch","organic_dep_autoch","realmarshend","resiliencymetrics","totalmarshOMto5500","VASEAmatrix","yearlyCaccumulation",'SLR')
    %end
    end
    end
end