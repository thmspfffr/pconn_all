%% COMPUTE WHOLE BRAIN DFA IN SOURCE SPACE
% all_src_plot

clear all

% v = 1: cortex, eloreta, v = 2: cortex, lcmv, v = 3: coarse, eloreta
% str1 = dfa, str2 = amp, str3 = cvar, str4 = var

STR = [1:4];
v = [2]

% --------------------------------------------------------
% VERSION 1
% --------------------------------------------------------
v         = v;
v_stat    = v;
v_rawdata = 2;
SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32];
% --------------------------------------------------------

allstr    = {'dfa';'amp';'cvar';'var'};
gridsize  = 'cortex';
contrasts = [2 1; 3 1; 2 3];
viewdir   = [0.001 -0 1; 0 -.5 0; 1 0 0; -.5 0 0; -.5 0 0;1 0 0 ];

outdir    = '/home/tpfeffer/pconn_all/proc/';
plotdir   = '/home/tpfeffer/pconn_all/plots/';

addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/
addpath ~/pconn/matlab/
addpath ~/Documents/MATLAB/Colormaps/Colormaps' (5)'/Colormaps/

for istr = STR
  
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
  %%
  
  ifoi = 1;
  
  str = 'amp'
  
  for isubj = SUBJLIST
    %         fprintf('Reading data s%d ...\n',isubj);
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
  
  %%
  icontr = 2;
  ifoi = 1;
  v_stat = 2;
  
  load(sprintf('~/pconn_all/proc/all_src_clusterstat_rst_%s_c%d_f%d_v%d.mat',str,icontr,ifoi,v_stat))

  idx1 = squeeze(nanmean(par_all_res(logical(stats.mask),1,:),1))>squeeze(nanmean(par_all_res(logical(stats.mask),2,:)));
  idx2 = ~idx1;
  
  figure; set(gcf,'color','white'); hold on
  
  scatter(squeeze(nanmean(par_all_res(logical(stats.mask),2,idx1),1)),squeeze(nanmean(par_all_res(logical(stats.mask),1,idx1),1)),200,'facecolor','k','markeredgecolor','w')
%   lsline
  scatter(squeeze(nanmean(par_all_res(logical(stats.mask),2,idx2),1)),squeeze(nanmean(par_all_res(logical(stats.mask),1,idx2),1)),200,'facecolor','r','markeredgecolor','w')

  
  rng = max([max(squeeze(nanmean(par_all_res(logical(stats.mask),1,:),1))) max(squeeze(nanmean(par_all_res(logical(stats.mask),2,:),1)))]);
  line([0 rng+0.1*rng],[0 rng+0.1*rng],'linestyle','--','color','k','linewidth',2)

  axis([0 rng+0.1*rng 0  rng+0.1*rng])
  tp_editplots;
  
  xlabel(sprintf('%s(PBO)',str),'fontweight','bold','fontsize',12)
  ylabel(sprintf('%s(DPZ)',str),'fontweight','bold','fontsize',12)

  print(gcf,'-depsc2',sprintf('~/pconn_all/plots/all_%s_scatter_ROI_f%d_v%d.eps',str,ifoi,v))
  
  
  
end
