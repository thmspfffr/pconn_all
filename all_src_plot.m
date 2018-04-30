%% COMPUTE WHOLE BRAIN DFA IN SOURCE SPACE
% all_src_plot

clear

% v = 1: cortex, eloreta, v = 2: cortex, lcmv, v = 3: coarse, eloreta
% str1 = dfa, str2 = amp, str3 = cvar, str4 = var
% freq in v = 2: [1 = 2 -12; 2 = 8 -12; 3 = 12 - 46; 4 = 54 - 100; 5 = 2 - 8]

foi = [6];
STR = [3];

addpath ~/Documents/MATLAB/cbrewer/cbrewer/
cmap1  = cbrewer('seq', 'YlOrRd', 150,'pchip'); cmap1 = cmap1(50:1:end,:);
cmap2  = cbrewer('seq', 'YlGnBu', 150,'pchip'); cmap2 = cmap2(end:-1:50,:);
cmap  = [cmap2; ones(50,3); cmap1];

for v = [2]
  %   v= 2;
  % --------------------------------------------------------
  % VERSION 1
  % --------------------------------------------------------
  v         = v;
  v_stat    = v;
  v_pow     = 1;
  SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
  % --------------------------------------------------------
  addpath ~/Documents/MATLAB/cbrewer/cbrewer/
  
  allstr    = {'dfa';'cvar';'pow'};
  gridsize  = 'cortex';
  contrasts = [2 1; 3 1; 2 3];
  
  outdir    = '/home/tpfeffer/pconn_all/proc/';
  
  addpath /home/gnolte/meg_toolbox/toolbox/
  addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
  addpath /home/gnolte/meg_toolbox/toolbox_nightly/
  addpath /home/gnolte/meg_toolbox/meg/
  addpath ~/pconn/matlab/
  addpath ~/Documents/MATLAB/Colormaps/Colormaps' (5)'/Colormaps/
  
  for istr = STR
    for ifoi = foi
      
      %       if ~exist(sprintf([outdir 'all_src_plot_s%d_f%d_v%d_processing.txt'],istr,ifoi,v))
      %         system(['touch ' outdir sprintf('all_src_plot_s%d_f%d_v%d_processing.txt',istr,ifoi,v)]);
      %       else
      %         continue
      %       end
      
      fprintf('Processing v%d s%d ...\n',v,istr);
      
      str = allstr{istr};
      
      % READ TEMPLATE MRI STUFF
      % ----------------
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
      dd    = .75;
      g1    = sa_meg_template.grid_cortex3000;
      % ----------------
      
      %% READ IN DATA
      fprintf('Reading data ...\n');
      ord   = pconn_randomization;
      
      for isubj = SUBJLIST
        %         fprintf('Reading data s%d ...\n',isubj);
        for m = 1 : 3
          
          im = find(ord(isubj,:)==m);
          
          if ~strcmp(str,'pow')
            load(sprintf(['~/pconn_cnt/proc/dfa/pconn_cnt_src_dfa_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,v));
          end
          if strcmp(str,'dfa')
            par_all_cnt(:,m,isubj)  = nanmean(par.dfa,2);
          elseif strcmp(str,'var')
            par_all_cnt(:,m,isubj)  = nanmean(par.var,2);
          elseif strcmp(str,'cvar')
            par_all_cnt(:,m,isubj) = nanmean(par.cvar,2);
          elseif strcmp(str,'amp')
            par_all_cnt(:,m,isubj)  = nanmean(par.amp,2); clear par
          elseif strcmp(str,'pow')
            load(sprintf(['~/pconn_cnt/proc/src/pconn_cnt_src_pow_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,v_pow));
            par_all_cnt(:,m,isubj)  = nanmean(par.pow,2);
          end
          %
          if ~strcmp(str,'pow')
            load(sprintf(['~/pconn/proc/dfa/pconn_src_dfa_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,v));
          end
          if strcmp(str,'dfa')
            par_all_res(:,m,isubj)  = nanmean(par.dfa,2);
          elseif strcmp(str,'var')
            par_all_res(:,m,isubj)  = nanmean(par.var,2);
          elseif strcmp(str,'cvar')
            par_all_res(:,m,isubj)  = nanmean(par.cvar,2);
          elseif strcmp(str,'amp')
            par_all_res(:,m,isubj)  = nanmean(par.amp,2);
          elseif strcmp(str,'pow')
            load(sprintf(['~/pconn/proc/src/pconn_src_pow_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,v_pow));
            par_all_res(:,m,isubj)  = nanmean(par.pow,2);
          end
          
          clear par
          
        end
      end
      
      par_all_res  = par_all_res(:,:,SUBJLIST);
      par_all_cnt  = par_all_cnt(:,:,SUBJLIST);
      
      fprintf('Reading data ... Done!\n');
      
      %% PHARMA COMPARISON
      bayes = 0;
      
      for mask = 0:1
        
        fprintf('Pharma comparison mask%d ...\n',mask);
        
        clear stats
        
        for icontr = 1:3
          
          % -------------------------------------------------------
          %   TASK
          % -------------------------------------------------------
          
          d = squeeze(nanmean(par_all_cnt(:,contrasts(icontr,1),:),3))-squeeze(nanmean(par_all_cnt(:,contrasts(icontr,2),:),3));
          %
          [~,~,~,tmp] = ttest(par_all_cnt(:,contrasts(icontr,1),:),par_all_cnt(:,contrasts(icontr,2),:),'dim',3);
          d = tmp.tstat; clear tmp
          
          if bayes
            fprintf('Computing Bayes factor...\n')
            for i = 1 : length(d)
              d(i) = t2smpbf(d(i),size(par_all_cnt,3),0.707);
            end
          end
          
          %       d(abs(d) < 2.5])=eps;
          
          if ~strcmp(str,'pow')
            if mask && exist(sprintf('~/pconn_all/proc/all_src_clusterstat_tsk_%s_c%d_f%d_v%d.mat',str,icontr,ifoi,v_stat))
              load(sprintf('~/pconn_all/proc/all_src_clusterstat_tsk_%s_c%d_f%d_v%d.mat',str,icontr,ifoi,v_stat))
              if isempty(stats.stat_pos) && isempty(stats.stat_neg)
                stats.mask = logical(zeros(size(stats.stat,1),1));
              end
              stats.mask = logical(stats.mask);
              if any(stats.mask)
                d(~stats.mask)=eps;
              else
                d(1:end) = eps;
              end
            end
          else
            if mask && exist(sprintf('~/pconn_all/proc/all_src_clusterstat_power_tsk_c%d_f%d_v%d.mat',icontr,ifoi,v_pow))
              load(sprintf('~/pconn_all/proc/all_src_clusterstat_power_tsk_c%d_f%d_v%d.mat',icontr,ifoi,v_pow))
              if isempty(stats.stat_pos) && isempty(stats.stat_neg)
                stats.mask = logical(zeros(size(stats.stat,1),1));
              end
              stats.mask = logical(stats.mask);
              if any(stats.mask)
                d(~stats.mask)=eps;
              else
                d(1:end) = eps;
              end
            end
          end
          par_interp = spatfiltergauss(d,g1,dd,g2);
          
          para = [];
          if bayes
            para.colorlimits = [0 2];
          else
            para.colorlimits = [-3 3];
          end
          
          % PLOT RESULTS
          para.filename = sprintf('~/pconn_all/plots/all_src_tsk_%s_mask%d_c%d_f%d_v%d.png',str,mask,icontr,ifoi,v_stat);
          tp_showsource(par_interp,cmap,sa_meg_template,para);
                    
          clear para r d stats
          
          % -------------------------------------------------------
          % RESTING STATE
          % -------------------------------------------------------
          
          d = squeeze(nanmean(par_all_res(:,contrasts(icontr,1),:),3))-squeeze(nanmean(par_all_res(:,contrasts(icontr,2),:),3));
          
          [~,~,~,tmp] = ttest(par_all_res(:,contrasts(icontr,1),:),par_all_res(:,contrasts(icontr,2),:),'dim',3);
          d = tmp.tstat; clear tmp
          
          if bayes
            fprintf('Computing Bayes factor...\n')
            for i = 1 : length(d)
              d(i) = t2smpbf(d(i),size(par_all_cnt,3),0.707);
            end
          end
          
          if ~strcmp(str,'pow')
            if mask && exist(sprintf('~/pconn_all/proc/all_src_clusterstat_rst_%s_c%d_f%d_v%d.mat',str,icontr,ifoi,v_stat))
              load(sprintf('~/pconn_all/proc/all_src_clusterstat_rst_%s_c%d_f%d_v%d.mat',str,icontr,ifoi,v_stat))
              if isempty(stats.stat_pos) && isempty(stats.stat_neg)
                stats.mask = logical(zeros(size(stats.stat,1),1));
              end
              stats.mask = logical(stats.mask);
              if any(stats.mask)
                d(~stats.mask)=eps;
              else
                d(1:end) = eps;
              end
            end
          else
            if mask && exist(sprintf('~/pconn_all/proc/all_src_clusterstat_power_rst_c%d_f%d_v%d.mat',icontr,ifoi,v_pow))
              load(sprintf('~/pconn_all/proc/all_src_clusterstat_power_rst_c%d_f%d_v%d.mat',icontr,ifoi,v_pow))
              if isempty(stats.stat_pos) && isempty(stats.stat_neg)
                stats.mask = logical(zeros(size(stats.stat,1),1));
              end
              stats.mask = logical(stats.mask);
              if any(stats.mask)
                d(~stats.mask)=eps;
              else
                d(1:end) = eps;
              end
            end
          end
          
          par_interp = spatfiltergauss(d,g1,dd,g2);
          
          para = [] ;
          % WHAT's WRONG HERE?
          %           para.colorlimits = [-1.96 1.96];
          if bayes
            para.colorlimits = [0 2];
          else
            para.colorlimits = [-3 3];
          end
          % PLOT RESULTS
          para.filename = sprintf('~/pconn_all/plots/all_src_rst_%s_mask%d_c%d_f%d_v%d.jpg',str,mask,icontr,ifoi,v_stat);
          tp_showsource(par_interp,cmap,sa_meg_template,para);
          
          clear stats
        end
        
      end
    end
  end
end

error('!')

%% SCATTER PLOTS (PLACEBO VS. DRUG)
str = 'var';
ifoi = 1;
v= 2;

ord   = pconn_randomization;

for isubj = SUBJLIST
  fprintf('Reading data s%d ...\n',isubj);
  for m = 1 : 3
    
    im = find(ord(isubj,:)==m);
    
    load(sprintf(['~/pconn_cnt/proc/dfa/pconn_cnt_src_dfa_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,v));
    
    if strcmp(str,'dfa')
      par_all_cnt(:,m,isubj)  = nanmean(par.dfa,2);
    elseif strcmp(str,'var')
      par_all_cnt(:,m,isubj)  = nanmean(par.var,2);
    elseif strcmp(str,'cvar')
      par_all_cnt(:,m,isubj) = nanmean(par.cvar,2);
    elseif strcmp(str,'amp')
      par_all_cnt(:,m,isubj)  = nanmean(par.amp,2); clear par
    end
    %
    load(sprintf(['~/pconn/proc/dfa/pconn_src_dfa_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,v));
    
    if strcmp(str,'dfa')
      par_all_res(:,m,isubj)  = nanmean(par.dfa,2);
    elseif strcmp(str,'var')
      par_all_res(:,m,isubj)  = nanmean(par.var,2);
    elseif strcmp(str,'cvar')
      par_all_res(:,m,isubj)  = nanmean(par.cvar,2);
    elseif strcmp(str,'amp')
      par_all_res(:,m,isubj)  = nanmean(par.amp,2);
    end
    
    clear par
    
  end
end


par_all_res  = par_all_res(:,:,SUBJLIST);
par_all_cnt  = par_all_cnt(:,:,SUBJLIST);

error('!')

%%  SHOW BLOCKS SEPARATELY


%% READ IN DATA



