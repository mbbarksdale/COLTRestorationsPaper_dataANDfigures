%%streamline code for saving for VASEA

function VASEAcodes(Restoration_strategy,CO,RSLR,RSLRA,Ke,SSCReduction,TLPthickness,TLPfrequency)

if nargin == 6
    TLPthickness=0;
    TLPfrequency=0;
else if nargin == 5
    SSCReduction=0;
    TLPthickness=0;
    TLPfrequency=0;    
end
end

if strcmp(Restoration_strategy, 'control') %aka no-restoration
    breakwater=0;
    TLP=0;
elseif strcmp(Restoration_strategy, 'breakwater')
    breakwater=1;
    TLP=0;
elseif strcmp(Restoration_strategy, 'TLP')
    breakwater=0;
    TLP=1; 
elseif strcmp(Restoration_strategy, 'breakwaterandTLP')
    breakwater=1;
    TLP=1;
else
    error('Error. Must correctly specify restoration strategy') %error message
end

Nuptake=zeros(150,500);
totalNuptake=zeros(150,1);
yearlyCaccumulation=zeros(150,1);
ElevAvg=zeros(150,1);
%for RSLRA=[0.1592,0.1106,0.08023,0.03536,0.01881]
%for Ke=[1.65e-09,2.82e-09,5.3e-09,7.8e-09,9.95e-09]

for Ke=Ke
    for RSLRA=RSLRA
        filename = ['Outputs_CO' num2str(CO) '_RSLRA' num2str(RSLRA) '_spinuprampedRSLR1_runstartRSLR' num2str(RSLR) '_Breakwater' num2str(breakwater) '_SSCReduction' num2str(SSCReduction) '_breakwaterdelay551_TLP' num2str(TLP) '_TLPthick' num2str(TLPthickness) '_TLPfreq' num2str(TLPfrequency) '_TLPdelay551_Erosion' num2str(Ke) '\'];
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
        if strcmp(Restoration_strategy, 'TLP') || strcmp(Restoration_strategy, 'breakwaterandTLP')
            yearlyCaccumulation(551:15:701,1)=NaN; %make years
            [yearlyCaccumulation,TF]=fillmissing(yearlyCaccumulation,'linear');
        end
        VASEAmatrix = zeros(701,5);
        VASEAmatrix(:,1)=msl;
        VASEAmatrix(:,2)=ElevAvg;
        VASEAmatrix(:,3)=resiliencymetrics(:,7);
        VASEAmatrix(:,4)=yearlyCaccumulation;
        VASEAmatrix(:,5)=totalNuptake;
        VASEAmatrix(:,6)=resiliencymetrics(:,11);
        VASEAmatrix(:,7)=totalmarshOMto5500;
        VASEAmatrix(:,8)=SLR;

        outputfilename=['CO' num2str(CO) '_RSLRA' num2str(RSLRA) '_Erosion' num2str(Ke) '/'];
        if exist(outputfilename,'file') ~= 7
            mkdir(outputfilename) %If filename does not exist, create a folder to save the output variable to
        end

        % commented code below for creating csv files of matrices
        % writematrix(VASEAmatrix,['CO' num2str(CO) '_RSLR' num2str(RSLR) '_Erosion' num2str(Erosion) '\VASEAmatrix_' Restoration_strategy '.csv']);
        % writematrix(elev,['CO' num2str(CO) '_RSLR' num2str(RSLR) '_Erosion' num2str(Erosion) '\elev_control' Restoration_strategy '.csv']);

        %save files
        save([filename 'VASEAmatrix.mat'],'VASEAmatrix'); %first save VASEAmatrix in the filename
        save([outputfilename 'VASEAmatrix_' Restoration_strategy '.mat'],'VASEAmatrix'); %then save all matrices in the outputfilename
        save([outputfilename 'elev_' Restoration_strategy '.mat'],'elev');
        save([outputfilename 'resiliencymetrics_' Restoration_strategy '.mat'],'resiliencymetrics');
        clear("agb","bgb","ElevAvg","elev",'Marsh_edge',"msl","Nuptake","organic_dep_alloch","organic_dep_autoch","realmarshend","resiliencymetrics","totalmarshOMto5500","totalNuptake","VASEAmatrix_control","yearlyCaccumulation",'SLR',"VASEAmatrix")
    end
end
end