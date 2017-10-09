%% CLUSTER BASED PERMUTATION TEST IN SOURCE SPACE
% all_src_clusterstat
% Performs cluster based permutation test on source space data


clear 
rng('shuffle','twister')

FOI = 2;
% v = 1: cortex, eloreta, v = 2: cortex, lcmv, v = 3: coarse, eloreta

for v = [2]
  
  
  if v == 2
    % --------------------------------------------------------
    % VERSION 1
    % --------------------------------------------------------
    v_out     = v;
    v_res     = v_out;
    v_cnt     = v_out;
    SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
    gridsize = 'cortex';
    para.method       = 'dependentT';
    para.minneigh     = 2;
    para.clusteralpha = 0.05;
    para.alpha        = 0.025;
    para.nperm        = 10000;
    % --------------------------------------------------------
  elseif v == 22
    % --------------------------------------------------------
    % VERSION 2
    % --------------------------------------------------------
%     v_out     = v;
%     v_res     = 2;
%     v_cnt     = 2;
%     v = 2;
%     SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
%     gridsize = 'cortex';
%     para.method       = 'dependentT';
%      para.minneigh     = 2;
%      para.clusteralpha = 0.05;
%      para.alpha        = 0.025;
%      para.nperm        = 10000;
    % --------------------------------------------------------
  end
  
  indir     = '/home/tpfeffer/pconn_cnt/proc/src/';
  outdir    = '/home/tpfeffer/pconn_all/proc/';
  plotdir   = '/home/tpfeffer/pconn_cnt/proc/plots/';
  
allstr       = {'dfa';'cvar';'pow'};
  
  for istr = [3]
    
%     if ~exist(sprintf([outdir 'all_src_clusterstat_str%d_v%d_processing.txt'],istr,v_out))
%       system(['touch ' outdir sprintf('all_src_clusterstat_str%d_v%d_processing.txt',istr,v_out)]);
%     else
%       continue
%     end
%     
    fprintf('Computing stats for v%d istr%d ...\n',v,istr);
    
    str = allstr{istr};
    
    restoredefaultpath
    
    addpath /home/gnolte/meg_toolbox/toolbox/
    addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
    addpath /home/gnolte/meg_toolbox/toolbox_nightly/
    addpath /home/gnolte/meg_toolbox/meg/
    addpath('/home/tpfeffer/Documents/MATLAB/fieldtrip-20160919/')
    addpath /home/tpfeffer/Documents/MATLAB/toolboxes/NBT-ReleaseNBTRC4a/
    addpath ~/pconn/matlab/
    
    
    ft_defaults
    
    run ~/Documents/MATLAB/toolboxes/NBT-NBTv0.5.3-alpha/installNBT.m
    
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
    
    for ifoi = 1: max(FOI)
      for isubj = SUBJLIST
        fprintf('Loading data s%d f%d ...\n',isubj,ifoi);
        for m = 1 : 3
          
          im = find(ord(isubj,:)==m);
          
          load(sprintf(['~/pconn_cnt/proc/dfa/pconn_cnt_src_dfa_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,v));
          
          dfa_all_cnt(:,m,isubj,ifoi)  = nanmean(par.dfa,2);
          cvar_all_cnt(:,m,isubj,ifoi) = nanmean(par.cvar,2); 
          %
          load(sprintf(['~/pconn/proc/dfa/pconn_src_dfa_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,v));
          %
          dfa_all_res(:,m,isubj,ifoi)  = nanmean(par.dfa,2);
          cvar_all_res(:,m,isubj,ifoi) = nanmean(par.cvar,2); 
          
          load(sprintf(['~/pconn_cnt/proc/src/pconn_cnt_src_pow_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,1));
          pow_all_cnt(:,m,isubj,ifoi)  = nanmean(par.pow,2); clear par

          load(sprintf(['~/pconn/proc/src/pconn_src_pow_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,1));
          pow_all_res(:,m,isubj,ifoi)  = nanmean(par.pow,2); clear par

          
        end
      end
    end
    
    dfa_all_res  = dfa_all_res(:,:,SUBJLIST,:);
    cvar_all_res = cvar_all_res(:,:,SUBJLIST,:);
    pow_all_res  = pow_all_res(:,:,SUBJLIST,:);

    dfa_all_cnt  = dfa_all_cnt(:,:,SUBJLIST,:);
    cvar_all_cnt = cvar_all_cnt(:,:,SUBJLIST,:);
    pow_all_cnt  = pow_all_cnt(:,:,SUBJLIST,:);
 
    %%
    for ifoi = FOI
      
      fprintf('Computing stats v%d for %s and freq%d  ...\n',v,str,ifoi);

      contrasts = [2 1; 3 1; 2 3];

      if strcmp(str,'dfa')
        par_all_res = dfa_all_res(:,:,:,ifoi);
        par_all_cnt = dfa_all_cnt(:,:,:,ifoi);
      elseif strcmp(str,'cvar')
        par_all_res = cvar_all_res(:,:,:,ifoi);
        par_all_cnt = cvar_all_cnt(:,:,:,ifoi);
      elseif strcmp(str,'pow')
        par_all_res = pow_all_res(:,:,:,ifoi);
        par_all_cnt = pow_all_cnt(:,:,:,ifoi);
      end
      
      % POOLED ACROSS BEHAVIORAL CONDITONS (REST/TASK(
      for icontr = 1 : 3
        
        z(:,:,1)= (par_all_res(:,contrasts(icontr,1),:)+par_all_cnt(:,contrasts(icontr,1),:))./2;
        z(:,:,2)= (par_all_res(:,contrasts(icontr,2),:)+par_all_cnt(:,contrasts(icontr,2),:))./2;
        
        stats = tp_clusterperm(z,para);

        save(sprintf('~/pconn_all/proc/all_src_clusterstat_pooled_%s_c%d_f%d_v%d.mat',str,icontr,ifoi,v_out),'stats');
        
      end
      
%       for icontr = 1 : 3
%         
%         z(:,:,1)=par_all_res(:,contrasts(icontr,1),:);
%         z(:,:,2)=par_all_res(:,contrasts(icontr,2),:);
%         
%         stats = tp_clusterperm(z,para);
% 
%         save(sprintf('~/pconn_all/proc/all_src_clusterstat_rst_%s_c%d_f%d_v%d.mat',str,icontr,ifoi,v_out),'stats');
%         
%       end
%       
%       clear z
%       
%       for icontr = 1 : 3
%         
%         z(:,:,1)=par_all_cnt(:,contrasts(icontr,1),:);
%         z(:,:,2)=par_all_cnt(:,contrasts(icontr,2),:);
%         
% %         stats = clust_perm(z,para);
%         stats = tp_clusterperm(z,para);
% 
%         save(sprintf('~/pconn_all/proc/all_src_clusterstat_tsk_%s_c%d_f%d_v%d.mat',str,icontr,ifoi,v_out),'stats');
%         
%       end
%       
%       REST VS TASK
%       
%       z(:,:,1)=nanmean(par_all_cnt,2);
%       z(:,:,2)=nanmean(par_all_res,2);
%       
% %       stats = clust_perm(z,para);
%       stats = tp_clusterperm(z,para);
% 
%       save(sprintf('~/pconn_all/proc/all_src_clusterstat_tsk-rst_%s_c%d_f%d_v%d.mat',str,icontr,ifoi,v_out),'stats');
      
    end
    
    
    %% COMPUTE DOUBLE PHARMA CONTRAST
    
    
%     for ifoi = FOI
%       
%       n =  get_neighbours(grid);
%       para.neigh        = n;
%       contrasts = [2 1; 3 1; 2 3];
%       viewdir   = [-.5 -.5 .5; .5 -.5 .5; .5 .5 .5; -.5 .5 .5];
%       
%       if strcmp(str,'dfa')
%         par_all_res = dfa_all_res(:,:,:,ifoi);
%         par_all_cnt = dfa_all_cnt(:,:,:,ifoi);
%       elseif strcmp(str,'pow')
%         par_all_res = pow_all_res(:,:,:,ifoi);
%         par_all_cnt = pow_all_cnt(:,:,:,ifoi);
%       elseif strcmp(str,'cvar')
%         par_all_res = cvar_all_res(:,:,:,ifoi);
%         par_all_cnt = cvar_all_cnt(:,:,:,ifoi);
%       end
%       
%       for icontr = 1 : 3
%         
%         z(:,:,1)=par_all_cnt(:,contrasts(icontr,1),:)-par_all_cnt(:,contrasts(icontr,2),:);
%         z(:,:,2)=par_all_res(:,contrasts(icontr,1),:)-par_all_res(:,contrasts(icontr,2),:);
%         
% %         stats = clust_perm(z,para);
%        	stats = tp_clusterperm(z,para);
% 
%         save(sprintf('~/pconn_all/proc/all_src_clusterstat_diffdiff_%s_c%d_f%d_v%d.mat',str,icontr,ifoi,v_out),'stats');
%         
%       end
      
    end
  end
% end

