%%% Project - Emotional Contagion - Pisa Pizzorusso 2024

%% This script is to generate the figures for the Manuscript
clear all
close all

%% Suppl. Figure - BarGraphs of all brain areas with statistics.
PathFile='...\Data'; %Set here the path where to find the raw data

%Define the type of data you want to plot
Type='cellden';

%Load the data
switch Type
    case 'cellden'
        Control_den=readtable(fullfile(PathFile,'noshock_density.csv'));
        Dem_den=readtable(fullfile(PathFile,'demonstrator_density.csv'));
        Obs_den=readtable(fullfile(PathFile,'observer_density.csv'));
        pVal=readtable(fullfile(PathFile,'pvalues_nonparamtest_density_198Areas_Fixed_20241021.csv'));

    case 'cellnum'
        Control_den=readtable(fullfile(PathFile,'noshock_n_cells.csv'));
        Dem_den=readtable(fullfile(PathFile,'demonstrator_n_cells.csv'));
        Obs_den=readtable(fullfile(PathFile,'observer_n_cells.csv'));
        pVal=readtable(fullfile(PathFile,'pvalues_nonparamtest_ncells_198Areas_Fixed_20241021.csv'));

    case 'energy'
        Control_den=readtable(fullfile(PathFile,'noshock_energy.csv'));
        Dem_den=readtable(fullfile(PathFile,'demonstrator_energy.csv'));
        Obs_den=readtable(fullfile(PathFile,'observer_energy.csv'));
        pVal=readtable(fullfile(PathFile,'pvalues_nonparamtest_energy_198Areas_Fixed_20241021.csv'));

end

%Remove the brain areas not included
iix=find(isnan(pVal.pval_noshock_vs_demonstrator)==1);

pVal(iix,:)=[];
Control_den(iix,:)=[];
Dem_den(iix,:)=[];
Obs_den(iix,:)=[];

%Now evaluate the median and 
MED=zeros(size(pVal,1),3);
Perc_l=zeros(size(pVal,1),3);
Perc_u=zeros(size(pVal,1),3);

for thisarea=1:size(pVal,1)
    MED(thisarea,3)=median(Control_den{thisarea,3:end});
    MED(thisarea,2)=median(Dem_den{thisarea,3:end});
    MED(thisarea,1)=median(Obs_den{thisarea,3:end});
 

    Perc_l(thisarea,3) = prctile(Control_den{thisarea,3:end},25);
    Perc_u(thisarea,3) = prctile(Control_den{thisarea,3:end},75);

    Perc_l(thisarea,2) = prctile(Dem_den{thisarea,3:end},25);
    Perc_u(thisarea,2) = prctile(Dem_den{thisarea,3:end},75);

    Perc_l(thisarea,1) = prctile(Obs_den{thisarea,3:end},25);
    Perc_u(thisarea,1) = prctile(Obs_den{thisarea,3:end},75);
end


%% Load Acronym
CSVTable=readtable('...Data\Resources\query.csv');

%"safe_name" is the correct variable to be used with ClearMap output to count cells
%"st_level" is the correct variable to select the areas you want to analyse
%"depth" is the correct variable to count cells
AllenTable=CSVTable(:,{'name','safe_name','id','acronym','depth','st_level'}); 

%%Now collect the acronym of areas
Areas=table;
for thisarea=1:size(pVal,1)
    ii=find(strcmp(AllenTable.safe_name,pVal.area{thisarea})==1);
    pVal.Acronym{thisarea}=AllenTable.acronym{ii};
    Areas(thisarea,:)=CSVTable(ii,:);
end


%% Plot data

%Define regions' order based on macro areas
ISO=1:1:23;
OLF=24:1:33;
HIP=34:1:42;
CSP=43:1:48;
STR=49:1:60;
PAL=61:1:68;
THL=69:1:104;
HYP=105:1:141;
MB=142:1:169;

%Plot ISO
C2=ISO;
figure ("Position",[100,100,1000,400])
hb=multibar_wError_prctl(MED(C2,:),Perc_l(C2,:),Perc_u(C2,:));
% set(gca, 'YScale', 'log')
set(gca,'XTick',1:1:length(C2),'xticklabel',pVal.Acronym(C2),'FontSize',10)
hb(1).FaceColor=[0.7031, 0.8281, 0.8711];
hb(2).FaceColor=[0.0625, 0.4375, 0.0625];
hb(3).FaceColor=[0.4961, 0.4961, 0.4961];
xtickangle(90)
ylabel('Cell Density','FontSize',10)
title('Isocortex')


%Plot OLF
C2=OLF;
figure ("Position",[100,100,1000,400])
hb=multibar_wError_prctl(MED(C2,:),Perc_l(C2,:),Perc_u(C2,:));
% set(gca, 'YScale', 'log')
set(gca,'XTick',1:1:length(C2),'xticklabel',pVal.Acronym(C2),'FontSize',10)
hb(1).FaceColor=[0.7031, 0.8281, 0.8711];
hb(2).FaceColor=[0.0625, 0.4375, 0.0625];
hb(3).FaceColor=[0.4961, 0.4961, 0.4961];
xtickangle(90)
ylabel('Cell Density','FontSize',10)
title('Olfactory Areas')

%Plot HIP
C2=HIP;
figure ("Position",[100,100,1000,400])
hb=multibar_wError_prctl(MED(C2,:),Perc_l(C2,:),Perc_u(C2,:));
% set(gca, 'YScale', 'log')
set(gca,'XTick',1:1:length(C2),'xticklabel',pVal.Acronym(C2),'FontSize',10)
hb(1).FaceColor=[0.7031, 0.8281, 0.8711];
hb(2).FaceColor=[0.0625, 0.4375, 0.0625];
hb(3).FaceColor=[0.4961, 0.4961, 0.4961];
xtickangle(90)
ylabel('Cell Density','FontSize',10)
title('Hippocampal Formation')


%Plot CSP
C2=CSP;
figure ("Position",[100,100,1000,400])
hb=multibar_wError_prctl(MED(C2,:),Perc_l(C2,:),Perc_u(C2,:));
% set(gca, 'YScale', 'log')
set(gca,'XTick',1:1:length(C2),'xticklabel',pVal.Acronym(C2),'FontSize',10)
hb(1).FaceColor=[0.7031, 0.8281, 0.8711];
hb(2).FaceColor=[0.0625, 0.4375, 0.0625];
hb(3).FaceColor=[0.4961, 0.4961, 0.4961];
xtickangle(90)
ylabel('Cell Density','FontSize',10)
title('Cortical Subplate')

%Plot STR
C2=STR;
figure ("Position",[100,100,1000,400])
hb=multibar_wError_prctl(MED(C2,:),Perc_l(C2,:),Perc_u(C2,:));
% set(gca, 'YScale', 'log')
set(gca,'XTick',1:1:length(C2),'xticklabel',pVal.Acronym(C2),'FontSize',10)
hb(1).FaceColor=[0.7031, 0.8281, 0.8711];
hb(2).FaceColor=[0.0625, 0.4375, 0.0625];
hb(3).FaceColor=[0.4961, 0.4961, 0.4961];
xtickangle(90)
ylabel('Cell Density','FontSize',10)
title('Striatum')

%Plot PAL
C2=PAL;
figure ("Position",[100,100,1000,400])
hb=multibar_wError_prctl(MED(C2,:),Perc_l(C2,:),Perc_u(C2,:));
% set(gca, 'YScale', 'log')
set(gca,'XTick',1:1:length(C2),'xticklabel',pVal.Acronym(C2),'FontSize',10)
hb(1).FaceColor=[0.7031, 0.8281, 0.8711];
hb(2).FaceColor=[0.0625, 0.4375, 0.0625];
hb(3).FaceColor=[0.4961, 0.4961, 0.4961];
xtickangle(90)
ylabel('Cell Density','FontSize',10)
title('Pallidum')

%Plot THL
C2=THL;
figure ("Position",[100,100,1000,400])
hb=multibar_wError_prctl(MED(C2,:),Perc_l(C2,:),Perc_u(C2,:));
% set(gca, 'YScale', 'log')
set(gca,'XTick',1:1:length(C2),'xticklabel',pVal.Acronym(C2),'FontSize',10)
hb(1).FaceColor=[0.7031, 0.8281, 0.8711];
hb(2).FaceColor=[0.0625, 0.4375, 0.0625];
hb(3).FaceColor=[0.4961, 0.4961, 0.4961];
xtickangle(90)
ylabel('Cell Density','FontSize',10)
title('Thalamus')

%Plot HYP
C2=HYP;
figure ("Position",[100,100,1000,400])
hb=multibar_wError_prctl(MED(C2,:),Perc_l(C2,:),Perc_u(C2,:));
% set(gca, 'YScale', 'log')
set(gca,'XTick',1:1:length(C2),'xticklabel',pVal.Acronym(C2),'FontSize',10)
hb(1).FaceColor=[0.7031, 0.8281, 0.8711];
hb(2).FaceColor=[0.0625, 0.4375, 0.0625];
hb(3).FaceColor=[0.4961, 0.4961, 0.4961];
xtickangle(90)
ylabel('Cell Density','FontSize',10)
title('Hypothalamus')


%Plot MB
C2=MB;
figure ("Position",[100,100,1000,400])
hb=multibar_wError_prctl(MED(C2,:),Perc_l(C2,:),Perc_u(C2,:));
% set(gca, 'YScale', 'log')
set(gca,'XTick',1:1:length(C2),'xticklabel',pVal.Acronym(C2),'FontSize',10)
hb(1).FaceColor=[0.7031, 0.8281, 0.8711];
hb(2).FaceColor=[0.0625, 0.4375, 0.0625];
hb(3).FaceColor=[0.4961, 0.4961, 0.4961];
xtickangle(90)
ylabel('Cell Density','FontSize',10)
title('Midbrain')
%% 
% 
% C2=find(pVal.pval_noshock_vs_demonstrator<=0.05|pVal.pval_noshock_vs_observer<=0.05|pVal.pval_demonstrator_vs_observer<=0.05);
% figure
% hb=multibar_wError_prctl(MED(C2,:),Perc_l(C2,:),Perc_u(C2,:));
% % set(gca, 'YScale', 'log')
% set(gca,'XTick',1:1:length(C2),'xticklabel',pVal.Acronym(C2),'FontSize',10)
% hb(1).FaceColor=[0.7031, 0.8281, 0.8711];
% hb(2).FaceColor=[0.0625, 0.4375, 0.0625];
% hb(3).FaceColor=[0.4961, 0.4961, 0.4961];
% xtickangle(90)
% ylabel('# of Cells','FontSize',10)
% title('Modulated by both Unfamiliar&Familiar')


%% Analyze data for Corr Coeff and Plot CC distribution

%Load the correlation Matrices

ttx=readtable(fullfile(PathFile,'corr_matrix_noshock_density.csv'));
Control_cc=table2array(ttx(:,2:end));
ttx=readtable(fullfile(PathFile,'corr_matrix_demonstrator_density.csv'));
Dem_cc=table2array(ttx(:,2:end));
ttx=readtable(fullfile(PathFile,'corr_matrix_observer_density.csv'));
Obs_cc=table2array(ttx(:,2:end));


% Flatten upper triangles of the correlation matrices
upper1 = Control_cc(triu(true(size(Control_cc)), 1));
upper2 = Dem_cc(triu(true(size(Dem_cc)), 1));
upper3 = Obs_cc(triu(true(size(Obs_cc)), 1));

Median_cc(1,1)=median(upper1);
Median_cc(2,1)=median(upper2);
Median_cc(3,1)=median(upper3);

% Combine the data and create a group label
data = [upper1; upper2; upper3];
group = [repmat(1, length(upper1), 1); repmat(2, length(upper2), 1); repmat(3, length(upper3), 1)];

% Perform the Kruskal-Wallis test
[p,h,sts] = kruskalwallis(data, group);
disp(p)

c=multcompare(sts);

% Perform K-S tests for pairwise comparisons
[h1, p1] = kstest2(upper1, upper2);
[h2, p2] = kstest2(upper1, upper3);
[h3, p3] = kstest2(upper2, upper3);

