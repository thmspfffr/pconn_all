%% COMPUTE WHOLE BRAIN DFA IN SOURCE SPACE
% pconn_src_dfa

clear all

% --------------------------------------------------------
% VERSION 1
% --------------------------------------------------------
v         = 1;
v_rawdata = 6;
fsample   = 400;
SUBJLIST  = [3 4 5 6 7 8 9 10 11 12 13 15 16 17 19 20 21 22 23 24];
foi       = [2 4; 4 8; 8 12; 12 36];
i_fit     = [1 100];
i_calc    = [0.5 150];
gridsize  = 'cortex';
filt      = 'eloreta';
% --------------------------------------------------------
% VERSION 2
% --------------------------------------------------------
% v         = 2;
% v_rawdata = 6;
% fsample   = 400;
% SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24];
% foi       = [2 4; 4 8; 8 12; 12 36];
% i_fit     = [1 100];
% i_calc    = [0.5 150];
% gridsize  = 'cortex';
% filt      = 'lcmv';
% --------------------------------------------------------
% VERSION 3
% --------------------------------------------------------
% v         = 3;
% v_rawdata = 6;
% fsample   = 400;
% SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24];
% foi       = [2 4; 4 8; 8 12; 12 36];
% i_fit     = [1 100];
% i_calc    = [0.5 150];
% gridsize  = 'coarse';
% filt      = 'eloreta';
% --------------------------------------------------------
% VERSION 4
% --------------------------------------------------------
% v         = 4;
% v_rawdata = 6;
% fsample   = 400;
% SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24];
% foi_tmp   = [3:70];
% i_fit     = [3 100];
% i_calc    = [2 150];
% gridsize  = 'xcoarse';
% filt      = 'eloreta';
% --------------------------------------------------------

if v==4
  for i = 1 : length(foi_tmp)
    foi(i,:) = [foi_tmp(i)-1 foi_tmp(i)+1];
  end
end

restoredefaultpath

addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/
addpath /home/tpfeffer/Documents/MATLAB/fieldtrip-20130925/
addpath /home/tpfeffer/Documents/MATLAB/toolboxes/NBT-ReleaseNBTRC4a/

ft_defaults

indir   = '/home/tpfeffer/pconn/proc/src/';
outdir   = '/home/tpfeffer/pconn/proc/dfa/';
% mridir = '/home/gnolte/neuconn/mri_data/';
plotdir = '/home/tpfeffer/pconn/proc/plots/';
freq = 1;

run ~/Documents/MATLAB/toolboxes/NBT-NBTv0.5.3-alpha/installNBT.m
%%

for ifoi = 1 : length(foi)
  for m = 1 : 3
    for isubj = 17
      
      if ~exist(sprintf([outdir 'pconn_src_dfa_s%d_m%d_f%d_v%d_processing.txt'],isubj,m,ifoi,v))
        system(['touch ' outdir sprintf('pconn_src_dfa_s%d_m%d_f%d_v%d_processing.txt',isubj,m,ifoi,v)]);
      else
        continue
      end
      %
      disp(sprintf('Processing s%d m%d f%d ...', isubj,m,ifoi))
      
      for iblock = 1 : 2
        
        disp(sprintf('Loading MEG data ...'));
        
        if isubj > 3 && isubj ~= 17
          load(sprintf('/home/tpfeffer/pconn/proc/preproc/pconn_postproc_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,v_rawdata));
        else
          load(sprintf('/home/tpfeffer/pconn/proc/preproc/pconn_postproc_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,3));
        end
        if freq == 1
          clear data_hi
          [mydata,epleng] = megdata2mydata(data_low);
          clear data_low
        elseif freq == 2
          clear data_low
          [mydata,epleng] = megdata2mydata(data_hi);
          clear data_hi
        else
          error('Missing information on frequency!')
        end
        
        % ------------------------------------------
        % COMPUTE DFA IN SOURCE SPACE
        % ------------------------------------------
        
        clear L grid
        % load cross spectrum
        load(['~/pconn/proc/src/' sprintf('pconn_sens_cs_s%d_m%d_b%d_f%d_v%d.mat',isubj,m,iblock,freq,1)]);
        
        if strcmp(filt,'lcmv')
          if strcmp(gridsize,'xcoarse')
            load([indir sprintf('pconn_sa_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,2)]);
            [~,A1] = mkfilt_lcmv(sa.L_xcoarse,nanmean(cs(:,:,foi(ifoi,1):foi(ifoi,2)),3));
          elseif strcmp(gridsize,'coarse')
            load([indir sprintf('pconn_sa_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,1)]);
            [~,A1] = mkfilt_lcmv(sa.L_coarse,nanmean(cs(:,:,foi(ifoi,1):foi(ifoi,2)),3));
        	elseif strcmp(gridsize,'cortex')
            load([indir sprintf('pconn_sa_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,3)]);
            [~,A1] = mkfilt_lcmv(sa.L_coarse,nanmean(cs(:,:,foi(ifoi,1):foi(ifoi,2)),3));
           	A1 = getdipdir(nanmean(cs(:,:,foi(ifoi,1):foi(ifoi,2)),3),A1);           
          end

        elseif strcmp(filt,'eloreta')
          if strcmp(gridsize,'xcoarse')
            load([indir sprintf('pconn_sa_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,2)]);
           	A1 = mkfilt_eloreta_v2(sa.L_xcoarse);
            A1 = getdipdir(nanmean(cs(:,:,foi(ifoi,1):foi(ifoi,2)),3),A1);
          elseif strcmp(gridsize,'coarse')
            load([indir sprintf('pconn_sa_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,1)]);
            A1 = mkfilt_eloreta_v2(sa.L_coarse);
          	A1 = getdipdir(nanmean(cs(:,:,foi(ifoi,1):foi(ifoi,2)),3),A1);
          elseif strcmp(gridsize,'cortex')
            load([indir sprintf('pconn_sa_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,3)]);
            % watch out! L_coarse is correct even for cortex grid
           	A1 = mkfilt_eloreta_v2(sa.L_coarse);
          	A1 = getdipdir(nanmean(cs(:,:,foi(ifoi,1):foi(ifoi,2)),3),A1);
          end
        end

        siginfo = nbt_Info;
        siginfo.converted_sample_frequency = fsample;
        
        % compute bp-filtered signal
        ampenv = nbt_filter_fir(mydata,foi(ifoi,1),foi(ifoi,2),siginfo.converted_sample_frequency,2/foi(ifoi,1));
        
        clear mydata data_low data_hi
        % project bp-filtered signal into source space
        ampenv=ampenv*A1;

        for ivox = 1 : size(ampenv,2)
          disp(ivox)
          tmp=abs(hilbert(ampenv(:,ivox)));
          ampenv1(:,ivox) = resample(double(tmp),1,20);
          clear tmp
        end
        
        siginfo = nbt_Info;
        siginfo.converted_sample_frequency = fsample/20;
        
        save(sprintf([outdir 'pconn_src_ampenv_s%d_b%d_m%d_f%d_v%d.mat'],isubj,iblock,m,ifoi,v),'ampenv1','-v7.3');
        
        % compute DFA
        tmp  = nbt_doDFA(ampenv1, siginfo, i_fit,i_calc,0,0,0,[]);
        
        par.dfa(:,iblock)  = tmp.MarkerValues;
        par.var(:,iblock)  = nanvar(ampenv1);
        par.cvar(:,iblock) = nanstd(ampenv1)./nanmean(ampenv1);
             
        
      	clear ampenv ampenv1

        
      end
      
      save(sprintf([outdir 'pconn_src_dfa_s%d_m%d_f%d_v%d.mat'],isubj,m,ifoi,v),'par','-v7.3');
      clear expo tmp
      
    end
  end
end

error('STOP')

%% PLOT DFA
v = 1;

SUBJLIST    = [4 5 6 7 8 9 10 11 12 13 15 16 17 19 20 21 22 23 24];

clear exp_all exp_ord exp_all_f ord

ord   = pconn_randomization;

for ifoi = 1 : 4
  
  for m = 1 : 3
    cnt = 0;
    
    for isubj = SUBJLIST
      
      im = find(ord(isubj,:)==m);
          
      load(sprintf([outdir 'pconn_src_dfa_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,v));
      
      exp_all(:,m,isubj,ifoi) = nanmean(expo,2);
      
    end
  end
end

exp_all = exp_all(:,:,SUBJLIST,:);

ifoi = 3;
exp_all_f = squeeze(exp_all(:,:,:,ifoi));

col = [1 0.3 0; .5 .5 1; 0 1 1];


m1 = squeeze(nanmean(exp_all_f));
s  = squeeze(std(m1,[],2)./sqrt(size(m1,2)));
m2 = nanmean(m1,2);
u  = m2 + s;
l  = m2 - s;

h = figure;
set(h,'color','white')

% ERROR BARS
for i = 1 :3
  plot([i],m2(i),'o','color',col(i,:),'MarkerSize',10,'MarkerFaceColor',col(i,:)); hold on
  line([i-0.1 i+0.1],[u(i) u(i)],'color',col(i,:),'LineWidth',2);
  line([i-0.1 i+0.1],[l(i) l(i)],'color',col(i,:),'LineWidth',2)
  line([i i],[u(i) l(i)],'color',col(i,:),'LineWidth',4)
end

title(sprintf('DFA exponent for freq-range: %d - %d Hz',foi(ifoi,1),foi(ifoi,2)))
xlabel('1 = placebo; 2 = atomoxetine; 3 = donepezil');
ylabel('DFA exponent');

box off
xlim([0 4]); ylim([0.95*min(m2) 1.05*max(m2)])
%     xlim([0 4]); ylim([.60 .75])

[~,p] = ttest(m1(1,:),m1(2,:),'dim',2);
disp(sprintf('Placebo vs Atomoxetine: p = %.3f',p))
[~,p] = ttest(m1(1,:),m1(3,:),'dim',2);
disp(sprintf('Placebo vs Donepezil: p = %.3f',p))
[~,p] = ttest(m1(2,:),m1(3,:),'dim',2);
disp(sprintf('Atomoxetine vs Donepezil: p = %.3f',p))

saveas(gcf,sprintf([plotdir 'pconn_dfa_f%d_v%d.fig'],ifoi,v),'fig')

%%
clear exp_ord

ifoi = 3;
icond = 2;

load sa_meg_template;
grid = sa_meg_template.grid_coarse;

mri  = sa_meg_template.mri;
vc = sa_meg_template.vc;
%  GENERATE SOURCE PLOT

% figure('units','normalized','position',[0 0 .9 .9])

para                  = [];
% para.mydotmarkersize  = 40;
para.orientation      = 'coronal';
if icond == 1
  PARA.colorlimits      = [min(min(nanmean(exp_all_f(:,1,:),3))) max(nanmean(exp_all_f(:,1,:),3))];
end
% para.colorlimits = PARA.colorlimits;
para.colormaps        = {'jet'};

h = figure; hold on
set(h,'color','k');

% PLACEBO - ATOMOXETINE
para.colorlimits = [0 1];
showmri_transp_v3(mri,para,[grid a]);
showmri_transp_v3(mri,para,[grid squeeze(nanmean(exp_all_f(:,1,:),3))-squeeze(nanmean(exp_all_f(:,2,:),3))]);

% PLACEBO - DONEPEZIL
% showmri_transp_v3(mri,para,[grid squeeze(nanmean(exp_ord(:,1,:),3))-squeeze(nanmean(exp_ord(:,3,:),3))]);

% PLACEBO
% showmri_transp_v3(mri,para,[grid squeeze(nanmean(exp_ord(:,icond,:),3))]);
% ALL CONDITIONS
para.colorlimits = [0.5 0.7];
showmri_transp_v3(mri,para,[grid squeeze(nanmean(nanmean(exp_all_f,3),2))]);


set(gcf,'units','normalized','position',[0.1 0.1 0.8 0.8])
