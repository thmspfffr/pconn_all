%% COMPUTE DFA SENSOR LEVEL TOPO STATISTICS
% all_sens_clusterstat
% Implement statistical test as described in Nichols & Holmes, 2001

% (1) Single threshold permutation test
% (2) Cluster-based permutation test

% tpfeffer | thms.pfffr@gmail.com | 05-05-15

clear

for v = [2]
  
  if v == 1
    % --------------------------------------------------------
    % VERSION 1
    % --------------------------------------------------------
    v_cnt             = 1;
    v_res             = 1;
    v_out             = 1;
    SUBJLIST          = [4 5 6 7 8 9 10 11 12 13 15 16 17 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
    CFG.clusteralpha  = 0.025;
    CFG.alpha         = 0.025;
  elseif v == 2
    % --------------------------------------------------------
    % VERSION 2
    % --------------------------------------------------------
    v_cnt             = 2;
    v_res             = 2;
    v_out             = 2;
    SUBJLIST          = [4 5 6 7 8 9 10 11 12 13 15 16 17 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
    CFG.clusteralpha  = 0.025;
    CFG.alpha         = 0.025;
  elseif v == 3
    % --------------------------------------------------------
    % VERSION 3
    % --------------------------------------------------------
    v_cnt             = 3;
    v_res             = 3;
    v_out             = 3;
    SUBJLIST          = [4 5 6 7 8 9 10 11 12 13 15 16 17 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
    CFG.clusteralpha  = 0.025;
    CFG.alpha         = 0.025;
  elseif v == 5
    % --------------------------------------------------------
    % VERSION 5
    % --------------------------------------------------------
    v_cnt             = 5;
    v_res             = 5;
    v_out             = 5;
    SUBJLIST          = [4 5 6 7 8 9 10 11 12 13 15 16 17 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
    CFG.clusteralpha  = 0.025;
    CFG.alpha         = 0.025;
  elseif v == 6
    % --------------------------------------------------------
    % VERSION 6
    % --------------------------------------------------------
    v_cnt             = 6;
    v_res             = 6;
    v_out             = 6;
    SUBJLIST          = [4 5 6 7 8 9 10 11 12 13 15 16 17 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
    CFG.clusteralpha  = 0.025;
    CFG.alpha         = 0.025;
  elseif v == 7
    % --------------------------------------------------------
    % VERSION 7
    % --------------------------------------------------------
    v_cnt             = 7;
    v_res             = 7;
    v_out             = 7;
    SUBJLIST          = [4 5 6 7 8 9 10 11 12 13 15 16 17 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
    CFG.clusteralpha  = 0.025;
    CFG.alpha         = 0.025;
    % --------------------------------------------------------
  elseif v == 8
    % --------------------------------------------------------
    % VERSION 8
    % --------------------------------------------------------
    v_cnt             = 8;
    v_res             = 8;
    v_out             = 8;
    SUBJLIST          = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24];
    CFG.clusteralpha  = 0.025;
    CFG.alpha         = 0.025;
    % --------------------------------------------------------
  elseif v == 9
    % --------------------------------------------------------
    % VERSION 9
    % --------------------------------------------------------
    v_cnt             = 9;
    v_res             = 9;
    v_out             = 9;
    SUBJLIST          = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24];
    CFG.clusteralpha  = 0.025;
    CFG.alpha         = 0.025;
    % --------------------------------------------------------
  end
  
  % addpath ~/Documents/MATLAB/fieldtrip-2016/fieldtrip-20160919/
  addpath /home/tpfeffer/Documents/MATLAB/fieldtrip-20130925/
  addpath /home/gnolte/meg_toolbox/toolbox/
  addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
  addpath /home/gnolte/meg_toolbox/toolbox_nightly/
  addpath /home/gnolte/meg_toolbox/meg/
  addpath /home/tpfeffer/pconn/matlab
  addpath /home/tpfeffer/pcbi/
  
  allstr = {'dfa';'amp';'cvar';'var'};
  
  outdir   = '/home/tpfeffer/pconn_all/proc/';
  ft_defaults
  
  %% CLUSTER-BASED PERMUTATION TEST
  % FOR WITHIN SUBJECTS DESIGNS
  
  for istr = 1 : length(allstr)
    
    str = allstr{istr};
    
    if ~exist(sprintf([outdir 'all_sens_clusterstat_s%d_v%d_processing.txt'],istr,v_out))
      system(['touch ' outdir sprintf('all_sens_clusterstat_s%d_v%d_processing.txt',istr,v_out)]);
    else
      continue
    end
    
    ord	= pconn_randomization;
    load ~/pconn/proc/pconn_label274.mat
    
%     load /home/tpfeffer/pconn/proc/src/pconn_sa_s4_m1_b1_v1.mat
    
    for ifoi = 1 : 4
      for isubj = SUBJLIST
        for m = 1 : 3
          
          im = find(ord(isubj,:)==m);
          
          tmp_cnt = pcbi_cnt(isubj);
          
          cnt(isubj,m) = nanmean(tmp_cnt((im*2-1):im*2));
          
          load(sprintf(['~/pconn_cnt/proc/dfa/' 'pconn_cnt_sens_dfa_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,v_cnt));
          
          dfa_cnt_all(:,isubj,m,ifoi)  = nanmean(par.dfa,2);
          var_cnt_all(:,isubj,m,ifoi)  = nanmean(par.var,2);
          cvar_cnt_all(:,isubj,m,ifoi) = nanmean(par.cvar,2);
          amp_cnt_all(:,isubj,m,ifoi)  = nanmean(par.amp,2);
          %       pow_cnt_all(:,isubj,m,ifoi)  = nanmean(par.pow,2); clear par
          
          load(sprintf(['~/pconn/proc/dfa/' 'pconn_sens_dfa_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,v_res));
          
          dfa_res_all(:,isubj,m,ifoi)  = nanmean(par.dfa,2);
          var_res_all(:,isubj,m,ifoi)  = nanmean(par.var,2);
          cvar_res_all(:,isubj,m,ifoi) = nanmean(par.cvar,2);
          amp_res_all(:,isubj,m,ifoi)  = nanmean(par.amp,2);
          
        end
      end
    end
    
    % cnt = cnt(SUBJLIST,:);
    dfa_cnt_all  = double(dfa_cnt_all(:,SUBJLIST,:,:));
    var_cnt_all  = double(var_cnt_all(:,SUBJLIST,:,:));
    cvar_cnt_all = double(cvar_cnt_all(:,SUBJLIST,:,:));
    amp_cnt_all  = double(amp_cnt_all(:,SUBJLIST,:,:));
    
    dfa_res_all  = double(dfa_res_all(:,SUBJLIST,:,:));
    var_res_all  = double(var_res_all(:,SUBJLIST,:,:));
    cvar_res_all = double(cvar_res_all(:,SUBJLIST,:,:));
    amp_res_all  = double(amp_res_all(:,SUBJLIST,:,:));
    
    %% COMPARE PHARACOLOGICAL CONDITIONS
    
    if strcmp(str,'dfa')
      par_cnt = dfa_cnt_all;
      par_res = dfa_res_all;
    elseif strcmp(str,'var')
      par_cnt = var_cnt_all;
      par_res = var_res_all;
    elseif strcmp(str,'cvar')
      par_cnt = cvar_cnt_all;
      par_res = cvar_res_all;
    elseif strcmp(str,'amp')
      par_cnt = amp_cnt_all;
      par_res = amp_res_all;
    end
    %
    load '/home/tpfeffer/pconn/proc/preproc/pconn_postproc_s4_m1_b1_v6.mat'
    data_low.label = lab274;
    
    clear data_hi
    %
    % -------------------------------------------------------
    % GET RID OF MISSING SENSORS
    % -------------------------------------------------------
    cfg             = [];
    cfg.method      = 'template';
    cfg.layout      = 'CTF275';
    n       = ft_prepare_neighbours(cfg);
    
    contrast = [2 1; 3 1; 2 3];
    
    for ifoi = 1 : 4
      
      for icontr = 1 : 3
        
        dat     = [par_cnt(:,:,contrast(icontr,1),ifoi) par_cnt(:,:,contrast(icontr,2),ifoi)];
        
        data_low.dimord                  = 'subj_chan_freq_time';
        data_low.time                   = [1 2];
        data_low.freq                   = [1 2];
        data_low.powspctrm(:,:,1:2,1:2) = repmat(permute(dat,[2 1]),[1 1 2 2]);
        
        cfg                  = [];
        cfg.channel          = 'all';
        cfg.latency          = [1 1];
        cfg.frequency        = [1 1];
        cfg.method           = 'montecarlo';
        cfg.statistic        = 'depsamplesT';
        cfg.computeprob      = 'yes';
        cfg.correctm         = 'cluster';
        cfg.clusteralpha     = CFG.clusteralpha;
        cfg.clusterstatistic = 'maxsum';
        cfg.clustertail      = 0; %1 = right
        cfg.tail             = 0; %1 = right
        cfg.alpha            = CFG.alpha;
        cfg.minnbchan        = 2;
        cfg.numrandomization = 10000;
        cfg.avgovertime      = 'yes';
        cfg.avgoverfreq      = 'yes';
        
        %     specifies with which sensors other sensors can form clusters
        cfg_neighb.method           = 'template';
        cfg_neighb.template         = 'CTF275_neighb.mat';
        cfg_neighb.feedback         = 'no';
        cfg.neighbours = n;
        
        n_subj = length(SUBJLIST);
        design = zeros(2,2*n_subj);
        design(1,:) = repmat(1:n_subj,1,2);
        design(2,:) = mod(floor([0:(2*n_subj-1)]/(n_subj/1)),2)+1;
        
        cfg.design   = design;
        cfg.uvar     = 1;
        cfg.ivar     = 2;
        
        [stats] = ft_freqstatistics(cfg, data_low);
        
        save(sprintf('~/pconn_all/proc/all_sens_tsk_%s_clusterstat_contr%d_f%d_v%d.mat',str,icontr,ifoi,v_out),'stats','SUBJLIST');
        
      end
    end
    
    % COMPARE PHARMA DURING REST
    
    for ifoi = 1 : 4
      
      for icontr = 1 : 3
        
        dat     = [par_res(:,:,contrast(icontr,1),ifoi) par_res(:,:,contrast(icontr,2),ifoi)];
        
        data_low.dimord                  = 'subj_chan_freq_time';
        data_low.time                   = [1 2];
        data_low.freq                   = [1 2];
        data_low.powspctrm(:,:,1:2,1:2) = repmat(permute(dat,[2 1]),[1 1 2 2]);
        
        cfg                  = [];
        cfg.channel          = 'all';
        cfg.latency          = [1 1];
        cfg.frequency        = [1 1];
        cfg.method           = 'montecarlo';
        cfg.statistic        = 'depsamplesT';
        cfg.computeprob      = 'yes';
        cfg.correctm         = 'cluster';
        cfg.clusteralpha     = CFG.clusteralpha;
        cfg.clusterstatistic = 'maxsum';
        cfg.clustertail      = 0; %1 = right
        cfg.tail             = 0; %1 = right
        cfg.alpha            = CFG.alpha;
        cfg.minnbchan        = 2;
        cfg.numrandomization = 10000;
        cfg.avgovertime      = 'yes';
        cfg.avgoverfreq      = 'yes';
        
        %     specifies with which sensors other sensors can form clusters
        cfg_neighb.method           = 'template';
        cfg_neighb.template         = 'CTF275_neighb.mat';
        cfg_neighb.feedback         = 'no';
        cfg.neighbours = n;
        
        n_subj = length(SUBJLIST);
        design = zeros(2,2*n_subj);
        design(1,:) = repmat(1:n_subj,1,2);
        design(2,:) = mod(floor([0:(2*n_subj-1)]/(n_subj/1)),2)+1;
        
        cfg.design   = design;
        cfg.uvar     = 1;
        cfg.ivar     = 2;
        
        [stats] = ft_freqstatistics(cfg, data_low);
        
        save(sprintf('~/pconn_all/proc/all_sens_res_%s_clusterstat_contr%d_f%d_v%d.mat',str,icontr,ifoi,v_out),'stats','SUBJLIST');
        
      end
    end
    
    % COMPARE TASK VS REST
    
    for ifoi = 1 : 4
      
      dat     = [nanmean(par_cnt(:,:,:,ifoi),3) nanmean(par_res(:,:,:,ifoi),3)];
      
      data_low.dimord                  = 'subj_chan_freq_time';
      data_low.time                   = [1 2];
      data_low.freq                   = [1 2];
      data_low.powspctrm(:,:,1:2,1:2) = repmat(permute(dat,[2 1]),[1 1 2 2]);
      
      cfg                  = [];
      cfg.channel          = 'all';
      cfg.latency          = [1 1];
      cfg.frequency        = [1 1];
      cfg.method           = 'montecarlo';
      cfg.statistic        = 'depsamplesT';
      cfg.computeprob      = 'yes';
      cfg.correctm         = 'cluster';
      cfg.clusteralpha     = CFG.clusteralpha;
      cfg.clusterstatistic = 'maxsum';
      cfg.clustertail      = 0; %1 = right
      cfg.tail             = 0; %1 = right
      cfg.alpha            = CFG.alpha;
      cfg.minnbchan        = 2;
      cfg.numrandomization = 10000;
      cfg.avgovertime      = 'yes';
      cfg.avgoverfreq      = 'yes';
      
      %   specifies with which sensors other sensors can form clusters
      cfg_neighb.method           = 'template';
      cfg_neighb.template         = 'CTF275_neighb.mat';
      cfg_neighb.feedback         = 'no';
      cfg.neighbours = n;
      
      n_subj = length(SUBJLIST);
      design = zeros(2,2*n_subj);
      design(1,:) = repmat(1:n_subj,1,2);
      design(2,:) = mod(floor([0:(2*n_subj-1)]/(n_subj/1)),2)+1;
      
      cfg.design   = design;
      cfg.uvar     = 1;
      cfg.ivar     = 2;
      
      [stats] = ft_freqstatistics(cfg, data_low);
      
      save(sprintf('~/pconn_all/proc/all_sens_%s_clusterstat_cnt-rest_f%d_v%d.mat',str,ifoi,v_out),'stats','SUBJLIST');
      
    end
    
    %% COMPARE ATOMOX(TASK-REST) VS. PLACEBO(TASK-REST)
    
    contrasts = [2 1; 3 1; 2 3];
    
    for ifoi = 1 : 4
      
      for icontr = 1 : 3
        
        dat     = [par_cnt(:,:,contrasts(icontr,1),ifoi)-par_cnt(:,:,contrasts(icontr,2),ifoi) par_res(:,:,contrasts(icontr,1),ifoi)-par_res(:,:,contrasts(icontr,2),ifoi)];
        
        data_low.dimord                  = 'subj_chan_freq_time';
        data_low.time                   = [1 2];
        data_low.freq                   = [1 2];
        data_low.powspctrm(:,:,1:2,1:2) = repmat(permute(dat,[2 1]),[1 1 2 2]);
        
        cfg                  = [];
        cfg.channel          = 'all';
        cfg.latency          = [1 1];
        cfg.frequency        = [1 1];
        cfg.method           = 'montecarlo';
        cfg.statistic        = 'depsamplesT';
        cfg.computeprob      = 'yes';
        cfg.correctm         = 'cluster';
        cfg.clusteralpha     = CFG.clusteralpha;
        cfg.clusterstatistic = 'maxsum';
        cfg.clustertail      = 0; %1 = right
        cfg.tail             = 0; %1 = right
        cfg.alpha            = CFG.alpha;
        cfg.minnbchan        = 2;
        cfg.numrandomization = 10000;
        cfg.avgovertime      = 'yes';
        cfg.avgoverfreq      = 'yes';
        
        % specifies with which sensors other sensors can form clusters
        cfg_neighb.method           = 'template';
        cfg_neighb.template         = 'CTF275_neighb.mat';
        cfg_neighb.feedback         = 'no';
        cfg.neighbours = n;
        
        n_subj = length(SUBJLIST);
        design = zeros(2,2*n_subj);
        design(1,:) = repmat(1:n_subj,1,2);
        design(2,:) = mod(floor([0:(2*n_subj-1)]/(n_subj/1)),2)+1;
        
        cfg.design   = design;
        cfg.uvar     = 1;
        cfg.ivar     = 2;
        
        [stats] = ft_freqstatistics(cfg, data_low);
        
        save(sprintf('~/pconn_all/proc/all_sens_%s_clusterstat_diffdiff_c%d_f%d_v%d.mat',str,icontr,ifoi,v_out),'stats','SUBJLIST');
        
      end
    end
    
  end
end
error('!')


%% BEHAVIOR

%
% ifoi = 5;
% contr = 1;
%
% clear im tmp dfa dfa_cnt_all dfa_res_all
%
% for isubj = SUBJLIST
%   for m = 1 : 3
%
%     im = find(ord(isubj,:)==m);
%
%     tmp_cnt = pcbi_cnt(isubj);
%
%     cnt(isubj,m) = nanmean(tmp_cnt((im*2-1):im*2));
%
%     for iblock = 1 : 2
%
%       load(sprintf([outdir 'pconn_cnt_sens_dfa_s%d_b%d_m%d_f%d_v%d.mat'],isubj,iblock,im,ifoi,v));
%
%       tmp(:,iblock) = dfa.MarkerValues; clear dfa
%
%     end
%
%     dfa_cnt_all(:,isubj,im) = nanmean(tmp,2); clear tmp
%
%     load(sprintf(['~/pconn/proc/dfa/' 'pconn_sens_dfa_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,v));
%     dfa_res_all(:,isubj,im) = nanmean(par.dfa,2);
%
%   end
% end
%
% load(sprintf('~/pconn_cnt/proc/dfa/pconn_cnt_sens_dfa_stat_clusterstat_contr%d_f%d_v%d.mat',icontr,ifoi,v));
%
% cnt = cnt(SUBJLIST,:);
% dfa_cnt_all = dfa_cnt_all(:,SUBJLIST,:);
% dfa_res_all = dfa_res_all(:,SUBJLIST,:);

%%

%
% clc
%
% % WHOLE BRAIN CORRELATION
%
% a     = squeeze(nanmean(dfa_cnt_all,1));
% slp   = pconn_regress(a(:)',cnt(:));
% [r,p] = corrcoef(a(:),cnt(:));
%
% figure; set(gcf,'color','white'); hold on
% line([min(a(:)) max(a(:))],[slp(2)*min(a(:))+slp(1) slp(2)*max(a(:))+slp(1)],'linewidth',7.5)
% scatter(a(:),cnt(:),200,'facebolor','r','markeredgecolor','w');
% set(gca,'TickDir','out','linewidth',3,'ticklength',[0.02 0.025],'fontsize',24);
% title(sprintf('Whole-brain: r = %.2f | p = %.2f',r(1,2),p(1,2)))
% xlabel('Difference DFA'); ylabel('Difference switches');
%
% axis([min(a(:))-0.1*min(a(:)) max(a(:))+0.1*max(a(:)) min(cnt(:))-2*min(cnt(:)) max(cnt(:))+0.2*max(cnt(:))]);
%
% print(gcf,'-depsc2',sprintf('~/pconn_cnt/plots/pconn_cnt_sens_behav_wb_f%d_v%d.eps',ifoi,v))
%
% if sum(stats.mask)==0
%   error('No significant clusters found!');
% end
%
% % CHANGE OF SWITCHES AFTER ATOMOX
%
% d_dfa = squeeze(nanmean(dfa_cnt_all(stats.mask,:,2),1)'-nanmean(dfa_cnt_all(stats.mask,:,1),1)');
% d_cnt = cnt(:,2)-cnt(:,1);
% slp   = pconn_regress(d_dfa',d_cnt);
%
% figure; set(gcf,'color','white'); hold on
% line([min(d_dfa) max(d_dfa)],[slp(2)*min(d_dfa)+slp(1) slp(2)*max(d_dfa)+slp(1)],'linewidth',7.5)
% scatter(d_dfa,d_cnt,200,'facebolor','r','markeredgecolor','w')
% set(gca,'TickDir','out','linewidth',3,'ticklength',[0.02 0.025]);
% axis([-0.05 0.15 -35 55]);
%
% xlabel('Difference DFA');ylabel('Difference switchs');
%
% [r,p]=corrcoef(d_dfa,d_cnt);
% title(sprintf('Cluster: r = %.2f | p = %.2f',r(1,2),p(1,2)))
%
% print(gcf,'-depsc2',sprintf('~/pconn_cnt/plots/pconn_cnt_sens_behav_pharmachange_f%d_c%d_v%d.eps',ifoi,icontr,v))
%
%
%
%









