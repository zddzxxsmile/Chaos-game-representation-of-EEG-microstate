% clear;
% clc;
% 
% % start EEGLAB to load all dependent paths
% eeglab
% %% set the path to the directory with the EEG files
% % change this path to the folder where the EEG files are saved
% EEGdir = 'F:\4447\';
% % retrieve a list of all EEG Files in EEGdir
% EEGFiles = dir([EEGdir '*.set']);
% %% 3.3 Data selection and aggregation
% %% 3.3.1 Loading datasets in EEGLAB
% seed=0;
% 
% load('exlabels4447');
% control=find(exlabels4447==0);
% psy=find(exlabels4447==1);
% 
% num=length(control);   
% mi=30;  
% for i=1:length(control)
% EEG = pop_loadset('filename',EEGFiles(control(i)).name,'filepath',EEGdir); 
% EEG=eeg_checkset(EEG);
% EEG = pop_select( EEG, 'channel',{'Fp1','Fp2','AF7','AF3','AF4','F7','F5','F3','F1','Fz','F2','F4','F6','F8','FT7','FC5','FC1','FC2','FC6','FT8','T7','C5','C3','C1','Cz','C2','C4','C6','T8','TP7','CP3','CP1','CP2','CP4','TP8','P7','P5','P3','Pz','P4','P6','P8','PO7','PO3','PO4','PO8','O1','Oz','O2'});
% EEG = eeg_checkset( EEG );
% EEG = pop_eegfiltnew(EEG, 'hicutoff',45);
% EEG = eeg_checkset( EEG );
% EEG = pop_reref( EEG, []); 
% EEG.event=[];
% EEG = pop_resample( EEG, 100);
% EEG = eeg_checkset( EEG );
% [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
% eeglab redraw % updates EEGLAB datasets
% end
% 
% save(strcat('pre332_control_145_100workspace_mi',num2str(mi),'_rng',num2str(seed)),'-v7.3');
% 
% clear;clc;
% for kk=1:100
% eeglab;
% load('pre332_control_145_100workspace_mi30_rng0.mat');
% load('M4_tem.mat');
% seed=kk;
% %%
% rng(seed);
% %% 3.3.2 Select data for microstate analysis
% [EEG, ALLEEG] = pop_micro_selectdata( EEG, ALLEEG, 'datatype', 'spontaneous',...
% 'avgref', 0, ...
% 'normalise', 1, ...
% 'MinPeakDist', 10, ...
% 'Npeaks', 1000, ...
% 'GFPthresh', 0, ...
% 'dataset_idx', 1:num);   
% 
% % store data in a new EEG structure
% [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
% eeglab redraw % updates EEGLAB datasets
% 
% %% 3.4 Microstate segmentation
% [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, num,'retrieve',num+1,'study',0);
% eeglab redraw
% 
% rng(seed);
% EEG = pop_micro_segment( EEG, 'algorithm', 'modkmeans', ...
% 'sorting', 'Global explained variance', ...
% 'Nmicrostates', 4, ...
% 'verbose', 1, ...
% 'normalise', 1, ...
% 'Nrepetitions', 100, ...
% 'max_iterations', 1000, ...
% 'threshold', 1e-06, ...
% 'fitmeas', 'CV',...
% 'optimised',1);
% [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
% 
% 
% %% 3.5 Review and select microstate segmentation
% %% 3.5.1 Plot microstate prototype topographies
% %figure;MicroPlotTopo( EEG, 'plot_range', [] );
% saveas(gcf, strcat('control4_rng',num2str(seed),'.jpg'));
% EEG = pop_saveset( EEG, 'filename',strcat('control_4_145_100workspace0_min',num2str(mi),'_rng',num2str(seed), '.set'), 'filepath','TEM');   
% end
% 
% %%
% clear;
% clc;
% eeglab
% EEGdir = 'F:\4447\TEM\';
% EEGFiles = dir([EEGdir 'control_4_145_100workspace0_min*.set']);
% load('M4_tem.mat');
% for i=1:length(EEGFiles)
%    
%    
%     EEG = pop_loadset('filename',EEGFiles(i).name,'filepath',EEGdir); 
%    
%    COR=abs(corr(M4_tem, EEG.microstate.prototypes));
%         
%    a=find(COR(3,:)== max(COR(3,:)));
%    temp(:,1)=EEG.microstate.prototypes(:,a);
%    a=find(COR(4,:)== max(COR(4,:)));
%     temp(:,2)=EEG.microstate.prototypes(:,a);
%      a=find(COR(2,:)== max(COR(2,:)));
%     temp(:,3)=EEG.microstate.prototypes(:,a);
%  a=find(COR(1,:)== max(COR(1,:)));
%     temp(:,4)=EEG.microstate.prototypes(:,a);
%     EEG.microstate.prototypes=temp;
%     EEG.microstate.Res.A_all{1,1}=temp;
%    EEG = pop_saveset( EEG, 'filename',strcat('sortABCD_',EEGFiles(i).name), 'filepath','TEM');
%     disp(i);
% 
% end
% 


%%
clear;clc;
for kk=1:100
eeglab;

load('pre332_145_100workspace_mi30_rng0.mat');
load('M4_tem.mat');
eeglab redraw;
seed=kk;
rng(seed);
EEG = pop_loadset('filename',strcat('sortABCD_control_4_145_100workspace0_min',num2str(mi),'_rng',num2str(seed), '.set'),'filepath','TEM');  
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
eeglab redraw;
EEG = pop_micro_selectNmicro( EEG,'Nmicro' ,4);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

%% Import microstate prototypes from other dataset to the datasets that should be back-fitted
for i = 1:length(EEGFiles)    
fprintf('Importing prototypes and backfitting for dataset %i\n',i)
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'retrieve',i,'study',0);
EEG = pop_micro_import_proto( EEG, ALLEEG, num+1);  
%% 3.6 Back-fit microstates on EEG
EEG = pop_micro_fit( EEG, 'polarity', 0 );
%% 3.7 Temporally smooth microstates labels
EEG = pop_micro_smooth( EEG, 'label_type', 'backfit', ...
'smooth_type', 'reject segments', ...
'minTime', mi, ...
'polarity', 0 );
%% 3.9 Calculate microstate statistics
EEG = pop_micro_stats( EEG, 'label_type', 'backfit', ...
'polarity', 0 );

GEV(i)=EEG.microstate.stats.GEVtotal;
Duration(i,:)=EEG.microstate.stats.Duration;
Coverage(i,:)=EEG.microstate.stats.Coverage;
Occurence(i,:)=EEG.microstate.stats.Occurence;
eee=EEG.microstate.stats.TP';
TP(i,:)=eee(:);
MspatCorr(i,:)=EEG.microstate.stats.MspatCorr;
MSClass=EEG.microstate.fit.labels;
MSC = reshape(MSClass, 1,  numel(MSClass));
ChangeIndex = find([0 diff(MSC)]);
DNA    = [MSC(1),MSC(ChangeIndex)];
DNA(find(DNA==0))=[];
ChangeIndex2 = find([0 diff(DNA)]);
DNA2    = [DNA(1),DNA(ChangeIndex2)];
SUB(i).DNA=DNA2;
SUB(i).name=ALLEEG(i).filename;
SUB(i).MS=MSC;
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
end

feature1 = [exlabels4447,GEV',Duration,Coverage, Occurence, TP,MspatCorr];  
writematrix(feature1, strcat('sortABCD_control_4_145_100workspace0_min',num2str(mi),'_rng',num2str(seed), '_feature_classic.txt'));

%%  microstate sequence for each subject
load('exlabels4447');
for j = 1:length(SUB)

bb=SUB(j).MS;

COR=abs(corr(M4_tem, EEG.microstate.prototypes));

for w=1:4
if find(COR(:,w)==max(COR(:,w)))==1
    bb(bb==w)='D';
elseif find(COR(:,w)==max(COR(:,w)))==2
    bb(bb==w)='C';
elseif find(COR(:,w)==max(COR(:,w)))==3
    bb(bb==w)='A';
elseif find(COR(:,w)==max(COR(:,w)))==4
    bb(bb==w)='B';
end
end

bb=char(bb);

    SUB(j).sequence=bb;

SUB(j).labels=exlabels4447(j);

end
 ee=rmfield(SUB, 'DNA');
 ee=rmfield(ee, 'MS');
 ee=rmfield(ee, 'name');
SUB4 = struct2table(ee);
writetable(SUB4, strcat('control_100RNA4447_145_4M_s',num2str(mi),'_rng',num2str(seed),'.txt'));
clear;

end

