
clear;clc;
b={};
SEQ={};


load('exlabels4447.mat')
con=find(exlabels4447==0);
psy=find(exlabels4447==1);

tic;
for i = 1:100
filename='psy_100RNA4447_145_4M_s30_rng';
a=readtable(strcat('.\M4_100_psy\',filename,num2str(i),'.txt'));
temp1=a(psy,:);
filename='control_100RNA4447_145_4M_s30_rng';
a=readtable(strcat('.\M4_100_control\',filename,num2str(i),'.txt'));
temp2=a(con,:);
c=[temp1;temp2];
% exporting sequence for computing FCGR using R
% writetable(c, strcat('M4_RNA_rng',num2str(i),'.txt'));
b{i}=table2struct(c);
for k=1:size(b{1},1)
SEQ{i,k} = b{1,i}(k).sequence;
end
end
toc;

i=size(SEQ,2);
k=size(SEQ,1);
rel_vec2=cell(k,i);
rel_vec=cell(k,i);
rel_o=cell(k,i);
rel_P1=cell(k,i);
rel_P12=cell(k,i);

tic;
Fs=100;
for i=1:size(SEQ,2)  % the i-th subject

for k=1:size(SEQ,1)  % the k-th randomization

real_seq=SEQ{k,i};
len=length(SEQ{k,i});
e=len;
[rel_vec{k,i},rel_o{k,i},rel_vec2{k,i}] = cgrDft(real_seq); % compute power spectrum
f=Fs*(0:round(e/2))/e;

%% Power spectrum for the time series Z
psd=rel_vec{k,i};
P2=abs(psd/e);
P1=P2(1:round(e/2+1)); 
P1(1:end-1)=2*P1(1:end-1);
rel_P1{k,i}=P1;

%% Power spectrum for the time series D
psd=rel_vec2{k,i};
P2=abs(psd/e);
P1=P2(1:round(e/2+1)); 
P1(1:end-1)=2*P1(1:end-1);
rel_P12{k,i}=P1;

%% Power spectrum indexes for the time series Z
AM(k,i)=mean(rel_P1{k,i});  
CF(k,i)=sum(f'.*rel_P1{k,i}')/sum(rel_P1{k,i}'); 
RMSF(k,i)=sqrt(sum((f.*f)'.*rel_P1{k,i}'/sum(rel_P1{k,i}'))); 
VF=sum(((f-CF(k,i)).^2)'.*rel_P1{k,i}')/sum(rel_P1{k,i}');
RVF(k,i)=sqrt(VF);

%% Power spectrum indexes for the time series D
AM2(k,i)=mean(rel_P12{k,i});  
CF2(k,i)=sum(f'.*rel_P12{k,i}')/sum(rel_P12{k,i}'); 
RMSF2(k,i)=sqrt(sum((f.*f)'.*rel_P12{k,i}'/sum(rel_P12{k,i}'))); 
VF2=sum(((f-CF2(k,i)).^2)'.*rel_P12{k,i}')/sum(rel_P12{k,i}'); 
RVF2(k,i)=sqrt(VF2);

%% Time domain characteristics for the time series D
rel_me(k,i)=mean(rel_o{k,i});
rel_med(k,i)=median(rel_o{k,i});
rel_st(k,i)=std(rel_o{k,i});
rel_rm(k,i)=rms(rel_o{k,i});   
disp(k);
end
disp(i);
end
toc;

control=82:142;
psy=1:81;

CF_mean=mean(CF);
RMSF_mean=mean(RMSF);
RVF_mean=mean(RVF);
AM_mean=mean(AM);
[~,p_AM,ci,stats]=ttest2(AM_mean(control),AM_mean(psy));
[~,p_CF,ci,stats]=ttest2(CF_mean(control),CF_mean(psy));
[~,p_RMSF,ci,stats]=ttest2(RMSF_mean(control),RMSF_mean(psy));
[~,p_RVF,ci,stats]=ttest2(RVF_mean(control),RVF_mean(psy));

CF2_mean=mean(CF2);
RMSF2_mean=mean(RMSF2);
RVF2_mean=mean(RVF2);
AM2_mean=mean(AM2);
[~,p_AM2,ci,stats]=ttest2(AM2_mean(control),AM2_mean(psy));
[~,p_CF2,ci,stats]=ttest2(CF2_mean(control),CF2_mean(psy));
[~,p_RMSF2,ci,stats]=ttest2(RMSF2_mean(control),RMSF2_mean(psy));
[~,p_RVF2,ci,stats]=ttest2(RVF2_mean(control),RVF2_mean(psy));

rel_med_mean=mean(rel_med);
rel_me_mean=mean(rel_me);
rel_rm_mean=mean(rel_rm);
rel_st_mean=mean(rel_st);
[~,p_med,ci,stats]=ttest2(rel_med_mean(control),rel_med_mean(psy));
[~,p_me,ci,stats]=ttest2(rel_me_mean(control),rel_me_mean(psy));
[~,p_st,ci,stats]=ttest2(rel_st_mean(control),rel_st_mean(psy));
[~,p_rm,ci,stats]=ttest2(rel_rm_mean(control),rel_rm_mean(psy));


a1=[AM_mean',CF_mean',RMSF_mean',RVF_mean'];
a12=[AM2_mean',CF2_mean',RMSF2_mean',RVF2_mean'];
ao=[rel_med_mean',rel_me_mean',rel_rm_mean',rel_st_mean'];

%% 
sss=1000; 

clear rel_P1_new rel_P12_new rel_o_new;
for i = 1:142
    for k=1:100  
    rel_P1_new(k,i,:)=Nspace(rel_P1{k,i},sss);
    rel_P12_new(k,i,:)=Nspace(rel_P12{k,i},sss);
   rel_o_new(k,i,:)=rel_o{k,i}(1:1000);  
    end
    

    disp(i);
end
     rel_P1_mean=squeeze(mean(rel_P1_new,1));
     rel_P12_mean=squeeze(mean(rel_P12_new,1));
     rel_o_mean=squeeze(mean(rel_o_new,1));      
    
rel_P1_mean_control=mean(rel_P1_mean(control,:),1);
rel_P1_mean_psy=mean(rel_P1_mean(psy,:),1);
rel_P12_mean_control=mean(rel_P12_mean(control,:),1);
rel_P12_mean_psy=mean(rel_P12_mean(psy,:),1);
rel_o_mean_control=mean(rel_o_mean(control,:),1);
rel_o_mean_psy=mean(rel_o_mean(psy,:),1);

save('results_sep_cgr_features','a1','a12','ao','rel_P1_mean','rel_P12_mean','rel_o_mean','rel_o_mean_control','rel_o_mean_psy','-v7.3');   