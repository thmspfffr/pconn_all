%% COMPUTE WHOLE BRAIN DFA IN SOURCE SPACE
% all_src_plot
% Plots source representations individually for each block
% Also computes statistics (cluster-based permutation) on each block

clear

% v = 1: cortex, eloreta, v = 2: cortex, lcmv, v = 3: coarse, eloreta
% str1 = dfa, str2 = amp, str3 = cvar, str4 = var
% freq in v = 2: [1 = 2 -12; 2 = 8 -12; 3 = 12 - 46; 4 = 54 - 100; 5 = 2 - 8]

foi = [2];
STR = [1];

addpath ~/Documents/MATLAB/cbrewer/cbrewer/
cmap1  = cbrewer('seq', 'YlOrRd', 150,'pchip'); cmap1 = cmap1(50:1:end,:);
cmap2  = cbrewer('seq', 'YlGnBu', 150,'pchip'); cmap2 = cmap2(end:-1:50,:);
cmap  = [cmap2; ones(50,3); cmap1];

v = 2

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
    
    %       if ~exist(sprintf([outdir 'all_src_plot_indivblocks_s%d_f%d_v%d_processing.txt'],istr,ifoi,v))
    %         system(['touch ' outdir sprintf('all_src_plot_s%d_f%d_v%d_processing.txt',istr,ifoi,v)]);
    %       else
    %         continue
    %       end
    
    fprintf('Processing v%d s%d ...\n',v,istr);
    
    %     str = allstr{istr};
    
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
    
    %     allstr    = {'dfa';'cvar';'pow'};
    gridsize  = 'cortex';
    contrasts = [2 1; 3 1; 2 3];
    
    str = 'dfa';
    ifoi = 2;
    
    fprintf('Reading data ...\n');
    ord   = pconn_randomization;
    
    for isubj = SUBJLIST
      %         fprintf('Reading data s%d ...\n',isubj);
      for m = 1 : 3
        
        im = find(ord(isubj,:)==m);
        
        %         if ~strcmp(str,'pow')
                  load(sprintf(['~/pconn_cnt/proc/dfa/pconn_cnt_src_dfa_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,v));
        %         end
        %         if strcmp(str,'dfa')
        par_all_cnt(:,m,isubj,:)  = par.dfa;
        %         elseif strcmp(str,'var')
        %           par_all_cnt(:,m,isubj,:)  = par.var;
        %         elseif strcmp(str,'cvar')
        %           par_all_cnt(:,m,isubj,:)  = par.cvar;
        %         elseif strcmp(str,'amp')
        %           par_all_cnt(:,m,isubj,:)  = par.amp; clear par
        %         elseif strcmp(str,'pow')
        %           load(sprintf(['~/pconn_cnt/proc/src/pconn_cnt_src_pow_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,v_pow));
        %           par_all_cnt(:,m,isubj,:)  = par.pow;
      
      %         %
      %         if ~strcmp(str,'pow')
                load(sprintf(['~/pconn/proc/dfa/pconn_src_dfa_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,v));
      %         end
      %         if strcmp(str,'dfa')
      par_all_res(:,m,isubj,:)  = par.dfa;
      %         elseif strcmp(str,'var')
      %           par_all_res(:,m,isubj,:)  = par.var;
      %         elseif strcmp(str,'cvar')
      %           par_all_res(:,m,isubj,:)  = par.cvar;
      %         elseif strcmp(str,'amp')
      %           par_all_res(:,m,isubj,:)  = par.amp;
      %         elseif strcmp(str,'pow')
      %           load(sprintf(['~/pconn/proc/src/pconn_src_pow_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,v_pow));
      %           par_all_res(:,m,isubj,:)  = par.pow;
      %         end
      
      clear par
      end
    end
  end
end


par_all_res  = par_all_res(:,:,SUBJLIST,:);
par_all_cnt  = par_all_cnt(:,:,SUBJLIST,:);

fprintf('Reading data ... Done!\n');

%% COMPUTE STATISTICS


gridsize          = 'cortex';
para.method       = 'dependentT';
para.minneigh     = 2;
para.clusteralpha = 0.05;
para.alpha        = 0.025;
para.nperm        = 10000;
n                 =  get_neighbours(grid);
para.neigh        = n;

for icontr = 1 : 3
  for iblock = 1 : 2
    
    z(:,:,1)= par_all_res(:,contrasts(icontr,1),:,iblock);
    z(:,:,2)= par_all_res(:,contrasts(icontr,2),:,iblock);
    
    stats = tp_clusterperm(z,para);
    
    save(sprintf('~/pconn_all/proc/all_src_clusterstat_res_%s_block%d_c%d_f%d_v%d.mat',str,iblock,icontr,ifoi,v_out),'stats');
    
    z(:,:,1)= par_all_cnt(:,contrasts(icontr,1),:,iblock);
    z(:,:,2)= par_all_cnt(:,contrasts(icontr,2),:,iblock);
    
    stats = tp_clusterperm(z,para);
    
    save(sprintf('~/pconn_all/proc/all_src_clusterstat_tsk_%s_block%d_c%d_f%d_v%d.mat',str,iblock,icontr,ifoi,v_out),'stats');
    
  end
end

error('!')
%%
bayes = 0;
for mask = 0:1
  for icontr = 1 : 3
    for iblock = 1 : 2

      [~,~,~,s]=ttest(par_all_res(:,contrasts(icontr,1),:,iblock),par_all_res(:,contrasts(icontr,2),:,iblock),'dim',3)
      d = s.tstat;
      
      if mask == 1
        d(abs(d)<1.96)=0;
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
      tp_showsource(par_interp,cmap,sa_meg_template,para);

      print(gcf,'-djpeg100',sprintf('~/pconn_all/plots/all_src_rst_%s_block%d_mask%d_c%d_f%d_v%d.jpg',str,iblock,mask,icontr,ifoi,v_stat))

      [~,~,~,s]=ttest(par_all_cnt(:,contrasts(icontr,1),:,iblock),par_all_cnt(:,contrasts(icontr,2),:,iblock),'dim',3)
      d = s.tstat;

      if mask == 1
        d(abs(d)<1.96)=0;
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
      tp_showsource(par_interp,cmap,sa_meg_template,para);

      print(gcf,'-djpeg100',sprintf('~/pconn_all/plots/all_src_tsk_%s_block%d_mask%d_c%d_f%d_v%d.jpg',str,iblock,mask,icontr,ifoi,v_stat))

    end
  end
end

%% COMPUTE CORRELATIONS ACROSS SUBJECTS
clear d

icontr = 1;

  for isubj = 1 : 28
    
    d1 = nanmean(a(:,contrasts(icontr,1),isubj,1)-a(:,contrasts(icontr,2),isubj,1),3);
    d2 = nanmean(a(:,contrasts(icontr,1),isubj,2)-a(:,contrasts(icontr,2),isubj,2),3);
    
    [r(isubj)]=corr(d1,d2);
    
    %     d(isubj,1) = mean(d1);
    %     d(isubj,2) = mean(d2);
    
    %     fprintf('Correlation c%d: r = %.2f | p = %.3f\n',icontr,r,p)
  %   figure;
  %   scatter(d(:,1),d(:,2),40,'markerfacecolor','k','markeredgecolor','w')
  %   axis square
  %   lsline
  %   axis([-0.1 0.1 -0.1 0.1]);
  %   xlabel('\Delta(DFA) (Block1)')
  %   ylabel('\Delta(DFA) (Block2)')
  %   title(sprintf('Contrast %d',icontr))
  end
  
  
  

  %% second one: averaged maps correlations
  
  par_all_cnt(:,2,27,1) = par_all_cnt(:,2,27,2);

  a = (par_all_res+par_all_cnt)./2;

  
  d1 = nanmean(par_all_res(:,contrasts(icontr,1),:,1),3)-nanmean(par_all_res(:,contrasts(icontr,2),:,1),3);
  d2 = nanmean(par_all_res(:,contrasts(icontr,1),:,2),3)-nanmean(par_all_res(:,contrasts(icontr,2),:,2),3);


