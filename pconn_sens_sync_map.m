%% MAP KURAMOTO
% pconn_cnt_sens_sync_map

% Method by R Hardstone, discussed during meeting on 4th of March
% Takes band-pass filtered amplitude and time-resolved DFA as input
% and estimates E/I balance from that.

% --------------------------------------------------------
% VERSION 1
% --------------------------------------------------------
% v           = 1;
% v_rawdata   = 1;
% is_src      = 0;
% fsample     = 400;
% SUBJLIST    = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24];
% foi         = [2 4; 4 8; 8 12; 12 36];
% filt_ord    = 2;
% i_fit       = [1 100];
% i_calc      = [0.5 150];
% dfa_overlap = 0.5;
% --------------------------------------------------------
% VERSION 2 (with DFA)
% --------------------------------------------------------
v           = 2;
v_rawdata   = 1;
is_src      = 0;
fsample     = 400;
SUBJLIST    = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24];
foi         = [2 4; 4 8; 8 12; 12 36];
filt_ord    = 2;
i_fit       = [1 100];
i_calc      = [0.5 150];
dfa_overlap = 0.5;
% --------------------------------------------------------

addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/
addpath /home/tpfeffer/Documents/MATLAB/fieldtrip-20130925/

indir   = '/home/tpfeffer/pconn_cnt/proc/src/';
outdir   = '/home/tpfeffer/pconn_cnt/proc/dfa/';
% mridir = '/home/gnolte/neuconn/mri_data/';
plotdir = '/home/tpfeffer/pconn_cnt/proc/plots/';

run ~/Documents/MATLAB/toolboxes/NBT-NBTv0.5.3-alpha/installNBT.m
%%
for m = 1 : 3
  for isubj = SUBJLIST
    for ifoi = 1 : length(foi)
      
      if ~exist(sprintf([outdir 'pconn_cnt_sens_sync_map_s%d_m%d_f%d_v%d_processing.txt'],isubj,m,ifoi,v))
        system(['touch ' outdir sprintf('pconn_cnt_sens_sync_map_s%d_m%d_f%d_v%d_processing.txt',isubj,m,ifoi,v)]);
      else
        continue
      end
%       
      if isubj == 3
        load ~/pconn/matlab/pconn_sensorlabels.mat
      end
      
      for iblock = 1 : 2
        
        disp(sprintf('Processing MEG s%dm%df%d%b...',isubj,m,ifoi,iblock));
        
        load(sprintf('/home/tpfeffer/pconn_cnt/proc/preproc/pconn_cnt_postproc_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,v_rawdata));
                
        if isubj == 3
          [~,idx_lab]=intersect(data_low.label,lab);
          data_low.trial{1}=data_low.trial{1}(idx_lab,:);
          data_low.label = lab;
          save(['~/pconn/matlab/' sprintf('pconn_idxlab.mat')],'idx_lab','-v7.3');
        end
                
        clear data_hi
        
        [mydata,~] = megdata2mydata(data_low);
        
        siginfo = nbt_Info;
        siginfo.converted_sample_frequency = fsample;
        
        % compute bp-filtered signal
        
        tmp    = nbt_filter_fir(mydata,foi(ifoi,1),foi(ifoi,2),siginfo.converted_sample_frequency,filt_ord/foi(ifoi,1));
        
        pha    = angle(hilbert(tmp));
        amp    = abs(hilbert(tmp));
        
        idx = 1 :268;
        
        r     = abs(sum(exp(i*pha),2)/268);
        r_std = std(r);
           
        for ichan = 1 : 268
          
          ichan 
          
          c_r(ichan) = corr(amp(:,ichan),r);
%           r_chan(ichan)     = mean(abs(sum(exp(i*pha(:,idx(idx~=ichan))),2)/267))./mean(r);
%          	r_std_chan(ichan) = std(abs(sum(exp(i*pha(:,idx(idx~=ichan))),2)/267))./r_std;
% 
          
        end
%           ichan        
         
%         end
               
        save(sprintf([outdir 'pconn_cnt_sens_sync_map_s%d_b%d_m%d_f%d_v%d.mat'],isubj,iblock,m,ifoi,v),'c_r');
        
        clear r tmp pha dfa mydata
        
      end  
    end
  end
end


error('STOP')

%% COUNT REJECTED COMPONENTS
clear s x
v  = 3;

for isubj = SUBJLIST
  for m = 1 : 3
    im = find(ord(isubj,:)==m);
    for iblock = 1 : 2
      
        load(sprintf('~/pconn/proc/preproc/pconn_rejected_comps_s%d_m%d_b%d_f1_v%d.mat',isubj,im,iblock,v));
        
        s(iblock,m,isubj) = sum(rej_comp);
        
      	load(sprintf('~/pconn_cnt/proc/preproc/pconn_cnt_rejected_comps_s%d_m%d_b%d_f1_v%d.mat',isubj,im,iblock,1));

        x(iblock,m,isubj) = sum(rej_comp);
  
    end
  end
end

s= s(:,:,SUBJLIST);
s=squeeze(nanmean(s,1));
nanmean(s(:))


x= x(:,:,SUBJLIST);
x=squeeze(nanmean(x,1));
nanmean(x(:))

%% 
clear all_sync
ifoi = 2;

SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16  19 20 21 22 23 24];
outdir   = '/home/tpfeffer/pconn_cnt/proc/dfa/';
ord       = pconn_randomization;

for isubj = SUBJLIST
  for m = 1 : 3
    
    im = find(ord(isubj,:)==m);
    
    for iblock = 1 : 2
      
      load(sprintf([outdir 'pconn_cnt_sens_sync_map_s%d_b%d_m%d_f%d_v%d.mat'],isubj,iblock,im,ifoi,2));
      
      all_sync(:,iblock,m,isubj) = c_r; clear r_chan
      
    end
  end
end

all_sync = squeeze(nanmean(all_sync(:,:,:,SUBJLIST),2));

load /home/tpfeffer/pconn/proc/src/pconn_sa_s4_m1_b1_v1.mat

      
%%


d = nanmean(all_sync(:,2,:)-all_sync(:,1,:),3);
[t,p] = ttest(all_sync(:,2,:),all_sync(:,1,:),'dim',3);

pars = [];
pars.resolution = 300;
pars.linewidth = 9;
pars.scale = [min(d) max(d)];
pars.cbar = 0;
pars.markersize = 0;
d = d.*t;

showfield_colormap(d,sa.locs_2D,pars)