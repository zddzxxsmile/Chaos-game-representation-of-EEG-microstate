
% clear;clc
% maindir = 'F:\ds003947-download\subjects';
% subdir =  dir( maindir );   
% for i = 1 : length( subdir )
%     if( isequal( subdir( i ).name, '.' ) ||  isequal( subdir( i ).name, '..' ) || ~subdir( i ).isdir )   
%         continue;
%     end
%     subdirpath = dir( [maindir,filesep, subdir( i ).name, filesep, 'eeg',filesep,'*'] );
%     copyfile(strcat(maindir,filesep, subdir( i ).name, filesep, 'eeg',filesep,'*'), maindir)
% end

clear;clc;
%%  
group1_dir = 'K:\ds003947-download\subjects';     
group1_files = dir([group1_dir, filesep, '*.vhdr']);
group1_loc = dir([group1_dir, filesep, '*.tsv']); 

tic;
for i=1:length(group1_files)
    subj_fn = group1_files(i).name;
    EEG = pop_loadbv(group1_dir, subj_fn);
    EEG = eeg_checkset( EEG );
    ss=round(EEG.srate);
    EEG = pop_editset(EEG, 'srate', [ss], 'run', []);
    EEG = pop_resample( EEG, 1000); 
    EEG = eeg_checkset( EEG );
    EEG=pop_chanedit(EEG, 'load',{strcat(group1_dir,filesep, group1_loc(i).name),'filetype','autodetect'},'lookup','D:\\eeglab2022.0\\plugins\\dipfit\\standard_BEM\\elec\\standard_1005.elc');
    EEG = eeg_checkset( EEG );
    EEG = pop_select( EEG,'nochannel',{'Iz' 'Misc' 'ECG' 'M2' 'FT9' 'FT10' 'T9' 'T10' 'TP9' 'TP10'});
    EEG = eeg_checkset( EEG );
    EEG = pop_saveset( EEG, 'filename',strcat(group1_files(i).name(1:end-4), '.set'), 'filepath',strcat(group1_dir, filesep, '_resam_remch'));  
end
toc;
   
%%
tic;
for i=1:length(group1_files)
    subj_fn = group1_files(i).name;
    EEG = pop_loadset('filename',strcat(subj_fn(1:end-4), '.set'), 'filepath', strcat(group1_dir, filesep, '_resam_remch'));
     EEG = pop_eegfiltnew(EEG, 'locutoff',0.1,'hicutoff',80);
     EEG = eeg_checkset( EEG );
     EEG = pop_eegfiltnew(EEG, 'locutoff',58,'hicutoff',62,'revfilt',1);
     EEG = eeg_checkset( EEG );
    EEG = pop_saveset( EEG, 'filename',strcat(group1_files(i).name(1:end-4), '.set'), 'filepath',strcat(group1_dir, filesep, '_filt_notch'));
end
toc;

%% EEG data was inspected manually to detect bad channels and segments
%% EEG data was saved to the '_preica' file
clc;clear all;close all
group1_dir = 'K:\ds003947-download\subjects';     
group1_dir1 = 'K:\ds003947-download\subjects\_preica';     
group1_files = dir([group1_dir1, filesep, '*.set']);  

for i=1:length(group1_files)
    subj_fn = group1_files(i).name;
    EEG = pop_loadset('filename',strcat(subj_fn(1:end-4), '.set'), 'filepath', strcat(group1_dir, filesep, '_preica'));
    EEG = pop_select( EEG, 'nochannel',{'VEOG'});
     EEG = pop_eegfiltnew(EEG, 'locutoff',1);
    EEG = pop_runica(EEG, 'icatype', 'runica', 'extended',1,'interrupt','on');
    EEG = eeg_checkset( EEG );
    EEG = pop_iclabel(EEG, 'default');
    EEG = eeg_checkset( EEG );
    EEG = pop_saveset( EEG, 'filename',strcat(group1_files(i).name(1:end-5), '.set'), 'filepath',strcat(group1_dir, filesep, '_ica'));    
end

%%  
clear; clc;
group1_dir = 'K:\ds003947-download\subjects';    
group1_dir2 = 'K:\ds003947-download\subjects\_ica'; 
group1_files = dir([group1_dir2, filesep, '*.set']);  
tic;
for i=1:length(group1_files)
    subj_fn = group1_files(i).name;
    EEG = pop_loadset('filename',strcat(subj_fn(1:end-4), '.set'), 'filepath', group1_dir2);
     EEG = pop_icflag(EEG, [0 0.2;0.7 1;0.7 1;0.7 1;0.7 1;0.7 1;0.8 1]);  % brain, muscle, eye, heart, line noise, channel noise, other
     EEG = eeg_checkset( EEG );
     EEG = pop_subcomp( EEG, [find(EEG.reject.gcompreject==1)], 0);
     EEG = eeg_checkset( EEG );
    EEG = pop_saveset( EEG, 'filename',strcat(group1_files(i).name(1:end-4), '.set'), 'filepath',strcat(group1_dir, filesep, '_rm_ica_label')); 
    clear artcomps info
end
toc;

%% The preprocessed data was uploaded.