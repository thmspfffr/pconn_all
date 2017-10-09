%% CLUSTER BASED PERMUTATION TEST IN SOURCE SPACE
% all_src_clusterstat_power
% Performs cluster based permutation test on source space data


clear 

rng('shuffle','twister')

maxfoi = 5;
% v = 1: cortex, eloreta, v = 2: cortex, lcmv, v = 3: coarse, eloreta

% --------------------------------------------------------
% VERSION 1
% --------------------------------------------------------
v = 1;
SUBJLIST          = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
gridsize          = 'cortex';
para.minneigh     = 2;
para.method       = 'dependentT';
para.clusteralpha = 0.05;
para.alpha        = 0.025;
para.nperm        = 10000;
% --------------------------------------------------------

indir     = '/home/tpfeffer/pconn_cnt/proc/src/';
outdir    = '/home/tpfeffer/pconn_all/proc/';
plotdir   = '/home/tpfeffer/pconn_cnt/proc/plots/';

restoredefaultpath

addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/
addpath('/home/tpfeffer/Documents/MATLAB/fieldtrip-20160919/')
addpath ~/pconn/matlab/


ft_defaults


load sa_meg_template;

if strcmp(gridsize,'coarse')
  grid  = sa_meg_template.grid_coarse;
elseif strcmp(gridsize,'cortex')
  grid  = sa_meg_template.grid_cortex3000;
  g1 = sa_meg_template.grid_cortex3000;
  g2 = sa_meg_template.cortex10K.vc;
elseif strcmp(gridsize,'xcoarse')
  grid  = sa_meg_template.grid_xcoarse;
end

mri   = sa_meg_template.mri;
vc    = sa_meg_template.vc;
dd = .1;

n =  get_neighbours(grid);
para.neigh        = n;


%% PLOT DFA

ord   = pconn_randomization;

for ifoi = 1: maxfoi
  for isubj = SUBJLIST
    fprintf('Loading data s%d f%d ...\n',isubj,ifoi);
    for m = 1 : 3
      
      im = find(ord(isubj,:)==m);
      
      load(sprintf(['~/pconn_cnt/proc/src/pconn_cnt_src_pow_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,v));
      
      pow_all_cnt(:,m,isubj,ifoi)  = nanmean(par.pow,2); clear par
      %
      load(sprintf(['~/pconn/proc/src/pconn_src_pow_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,v));
      %
      pow_all_res(:,m,isubj,ifoi)  = nanmean(par.pow,2); clear par
      
      
      
    end
  end
end


pow_all_res  = pow_all_res(:,:,SUBJLIST,:);
pow_all_cnt  = pow_all_cnt(:,:,SUBJLIST,:);

%%
for ifoi = 1 : maxfoi
      
  fprintf('Computing stats v%d for and freq%d  ...\n',v,ifoi);
  
  contrasts = [2 1; 3 1; 2 3];
    
  par_all_res = pow_all_res(:,:,:,ifoi);
  par_all_cnt = pow_all_cnt(:,:,:,ifoi);

  for icontr = 1 : 3
    
    if ~exist(sprintf('~/pconn_all/proc/all_src_clusterstat_power_rest_c%d_f%d_v%d_processing.txt',icontr,ifoi,v))
      system(['touch ' '~/pconn_all/proc/' sprintf('all_src_clusterstat_power_rest_c%d_f%d_v%d_processing.txt',icontr,ifoi,v)]);
    else
      continue
    end
    
    z(:,:,1)=par_all_res(:,contrasts(icontr,1),:);
    z(:,:,2)=par_all_res(:,contrasts(icontr,2),:);
    
    stats = tp_clusterperm(z,para);
    
    save(sprintf('~/pconn_all/proc/all_src_clusterstat_power_rst_c%d_f%d_v%d.mat',icontr,ifoi,v),'stats');
    
  end
  
  clear z
  
  for icontr = 1 : 3
    
    if ~exist(sprintf('~/pconn_all/proc/all_src_clusterstat_power_task_c%d_f%d_v%d_processing.txt',icontr,ifoi,v))
      system(['touch ' '~/pconn_all/proc/' sprintf('all_src_clusterstat_power_task_c%d_f%d_v%d_processing.txt',icontr,ifoi,v)]);
    else
      continue
    end
    
    z(:,:,1)=par_all_cnt(:,contrasts(icontr,1),:);
    z(:,:,2)=par_all_cnt(:,contrasts(icontr,2),:);
    
    stats = tp_clusterperm(z,para);
    
    save(sprintf('~/pconn_all/proc/all_src_clusterstat_power_tsk_c%d_f%d_v%d.mat',icontr,ifoi,v),'stats');
    
  end
  
  % REST VS TASK
  if ~exist(sprintf('~/pconn_all/proc/all_src_clusterstat_power_tsk-rst_c%d_f%d_v%d_processing.txt',icontr,ifoi,v))
    system(['touch ' '~/pconn_all/proc/' sprintf('all_src_clusterstat_power_tsk-rst_c%d_f%d_v%d_processing.txt',icontr,ifoi,v)]);
  else
    continue
  end
  
  z(:,:,1)=nanmean(par_all_cnt,2);
  z(:,:,2)=nanmean(par_all_res,2);
  
  stats = tp_clusterperm(z,para);

  save(sprintf('~/pconn_all/proc/all_src_clusterstat_power_tsk-rst_c%d_f%d_v%d.mat',icontr,ifoi,v),'stats');

  for icontr = 1 : 3
    
    if ~exist(sprintf('~/pconn_all/proc/all_src_clusterstat_power_doublecontrast_c%d_f%d_v%d_processing.txt',icontr,ifoi,v))
      system(['touch ' '~/pconn_all/proc/' sprintf('all_src_clusterstat_power_doublecontrast_c%d_f%d_v%d_processing.txt',icontr,ifoi,v)]);
    else
      continue
    end
    
    z(:,:,1)=par_all_cnt(:,contrasts(icontr,1),:)-par_all_cnt(:,contrasts(icontr,2),:);
    z(:,:,2)=par_all_res(:,contrasts(icontr,1),:)-par_all_res(:,contrasts(icontr,2),:);
    
    stats = tp_clusterperm(z,para);
    
    save(sprintf('~/pconn_all/proc/all_src_clusterstat_power_diffdiff_c%d_f%d_v%d.mat',icontr,ifoi,v),'stats');
    
  end
  

end


