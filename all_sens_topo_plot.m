%% PLOT SENSOR-LEVEL DFA RESULTS
% pconn_cnt_sens_dfa_plot

clear

v_cnt         = 22;
v_res         = 22;

SUBJLIST    = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24];

addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/
addpath /home/tpfeffer/pconn/matlab
addpath ~/pcbi/
addpath ~/Documents/MATLAB/cbrewer/cbrewer/

outdir   = '/home/tpfeffer/pconn_cnt/proc/dfa/';
plotdir = '/home/tpfeffer/pconn_all/plots/';


load /home/tpfeffer/pconn/proc/src/pconn_sa_s4_m1_b1_v1.mat

%%

ord       = pconn_randomization;


for ifoi = 1:4
  for isubj = SUBJLIST
    for m = 1 : 3
      
      im = find(ord(isubj,:)==m);
      
      tmp_cnt = pcbi_cnt(isubj);
      
      cnt(isubj,m) = nanmean(tmp_cnt((im*2-1):im*2));
      
      load(sprintf([outdir 'pconn_cnt_sens_dfa_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,v_cnt));
      
      dfa_cnt_all(:,isubj,m,ifoi)  = nanmean(par.dfa,2);
      var_cnt_all(:,isubj,m,ifoi)  = nanmean(par.var,2);
      cvar_cnt_all(:,isubj,m,ifoi) = nanmean(par.cvar,2);
      pow_cnt_all(:,isubj,m,ifoi)  = nanmean(par.pow,2); 
      amp_cnt_all(:,isubj,m,ifoi)  = nanmean(par.amp,2); clear par

      
      load(sprintf(['~/pconn/proc/dfa/' 'pconn_sens_dfa_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,v_res));
      
      dfa_res_all(:,isubj,m,ifoi)  = nanmean(par.dfa,2);
      var_res_all(:,isubj,m,ifoi)  = nanmean(par.var,2);
      cvar_res_all(:,isubj,m,ifoi) = nanmean(par.cvar,2);
      pow_res_all(:,isubj,m,ifoi)  = nanmean(par.pow,2);
      amp_res_all(:,isubj,m,ifoi)  = nanmean(par.amp,2); clear par
     
    end
  end
end

% cnt = cnt(SUBJLIST,:);
dfa_cnt_all  = double(dfa_cnt_all(:,SUBJLIST,:,:));
var_cnt_all  = double(var_cnt_all(:,SUBJLIST,:,:));
cvar_cnt_all = double(cvar_cnt_all(:,SUBJLIST,:,:));
pow_cnt_all  = double(pow_cnt_all(:,SUBJLIST,:,:));
amp_cnt_all  = double(amp_cnt_all(:,SUBJLIST,:,:));

dfa_res_all  = double(dfa_res_all(:,SUBJLIST,:,:));
var_res_all  = double(var_res_all(:,SUBJLIST,:,:));
cvar_res_all = double(cvar_res_all(:,SUBJLIST,:,:));
pow_res_all  = double(pow_res_all(:,SUBJLIST,:,:));
amp_res_all  = double(amp_res_all(:,SUBJLIST,:,:));

%% INDIV SUBJECTS
% SUBJ = [1:5; 6:10; 11:15; 16:19];

str = 'cvar';

for ifoi = 1 : 4
  
  if strcmp(str,'dfa')
    all_par_res = dfa_res_all(:,:,:,ifoi);
    all_par_cnt = dfa_cnt_all(:,:,:,ifoi);
  elseif strcmp(str,'var')
    all_par_res = var_res_all(:,:,:,ifoi);
    all_par_cnt = var_cnt_all(:,:,:,ifoi);
  elseif strcmp(str,'cvar')
    all_par_res = cvar_res_all(:,:,:,ifoi);
    all_par_cnt = cvar_cnt_all(:,:,:,ifoi);
  elseif strcmp(str,'pow')
    all_par_res = log10(pow_res_all(:,:,:,ifoi));
    all_par_cnt = log10(pow_cnt_all(:,:,:,ifoi));
  end
  
  for isubj = 1 : length(SUBJLIST)
    h=figure; set(h,'color','white'); %set(h,'Papertype','a4','visible','off')
    
    subj = SUBJLIST(isubj);
    
    subplot(1,3,1)
    
    par = nanmean(squeeze(all_par_res(:,isubj,:)),3);
    
    pars = [];
    pars.scale=[min(par(:)) max(par(:))];
    pars.cbar = 0;
    pars.markersize = 0;
    pars.linewidth = 4;
    pars.resolution = 300;
    showfield_colormap(nanmean(par,2),sa.locs_2D,pars);
    
    colormap(hot)
    
    subplot(1,3,2)
    
    par = nanmean(squeeze(all_par_cnt(:,isubj,:)),3);

    pars.cbar = 0;
    pars.markersize = 0;
    pars.scale=[min(par(:)) max(par(:))];
    pars.linewidth = 4;
    pars.resolution = 300;
    showfield_colormap(nanmean(par,2),sa.locs_2D,pars);
    
    subplot(1,3,3)
    
    par = nanmean(squeeze(all_par_cnt(:,isubj,:)-all_par_res(:,isubj,:)),3);
    
    pars.cbar = 0;
    pars.scale=[min(par(:)) max(par(:))];
    pars.markersize = 0;
    pars.linewidth = 4;
    pars.resolution = 300;
    showfield_colormap(nanmean(par,2),sa.locs_2D,pars);
    colormap(parula)
    print(gcf,'-djpeg',sprintf('~/pconn_cnt/plots/pconn_topo_all_%s_s%d_f%d_v%d.jpg',str,subj,ifoi,v_cnt))%%
    
    close
  end
end

error('!!')




%%
% ------------------------------------------------
% AVERAGE ACROSS CONDITIONS
% ------------------------------------------------
  
ifoi      = 1;
str       = 'cvar';
contrasts = [2 1; 3 1; 2 3];

if strcmp(str,'dfa')
  all_par_res = dfa_res_all(:,:,:,ifoi);
  all_par_cnt = dfa_cnt_all(:,:,:,ifoi);
elseif strcmp(str,'var')
  all_par_res = var_res_all(:,:,:,ifoi);
  all_par_cnt = var_cnt_all(:,:,:,ifoi);
elseif strcmp(str,'cvar')
  all_par_res = cvar_res_all(:,:,:,ifoi);
  all_par_cnt = cvar_cnt_all(:,:,:,ifoi);
elseif strcmp(str,'pow')
  all_par_res = log10(pow_res_all(:,:,:,ifoi));
  all_par_cnt = log10(pow_cnt_all(:,:,:,ifoi));
elseif strcmp(str,'amp')
  all_par_res = log10(amp_res_all(:,:,:,ifoi));
  all_par_cnt = log10(amp_cnt_all(:,:,:,ifoi));
end


figure; set(gcf,'color','white');

dfa_d1 = nanmean(nanmean(all_par_cnt,3),2);

pars.scale = [min(dfa_d1) max(dfa_d1)];
pars.markersize = 0;
pars.cbar = 0;
pars.linewidth = 9;
pars.resolution = 300;

showfield_colormap(dfa_d1,sa.locs_2D,pars);
colormap(hot)

print(gcf,'-djpeg',sprintf('~/pconn_all/plots/pconn_sens_topo_tsk_%s_avg_f%d_v%d.jpg',str,ifoi,v_cnt))

dfa_d1 = nanmean(nanmean(all_par_res,3),2);

pars.scale = [min(dfa_d1) max(dfa_d1)];

showfield_colormap(dfa_d1,sa.locs_2D,pars);
colormap(hot)

print(gcf,'-djpeg',sprintf('~/pconn_all/plots/pconn_sens_topo_rst_%s_avg_f%d_v%d.jpg',str,ifoi,v_cnt))

%
%% PHARMA COMPARISON & TASK VS REST
% -------------------------------------------

tval = 1;
ifoi      = 4;
str       = 'pow';
contrasts = [2 1; 3 1; 2 3];

if strcmp(str,'dfa')
  all_par_res = dfa_res_all(:,:,:,ifoi);
  all_par_cnt = dfa_cnt_all(:,:,:,ifoi);
elseif strcmp(str,'var')
  all_par_res = var_res_all(:,:,:,ifoi);
  all_par_cnt = var_cnt_all(:,:,:,ifoi);
elseif strcmp(str,'cvar')
  all_par_res = cvar_res_all(:,:,:,ifoi);
  all_par_cnt = cvar_cnt_all(:,:,:,ifoi);
elseif strcmp(str,'pow')
  all_par_res = log10(pow_res_all(:,:,:,ifoi));
  all_par_cnt = log10(pow_cnt_all(:,:,:,ifoi));
end

for icontr = 1 : 3
  
  pars = [];
  
  % ------------------------------------------------
  % PHARMA CONTRAST DURING TASK
  % ------------------------------------------------
  
  load(sprintf('~/pconn_all/proc/all_sens_tsk_%s_clusterstat_contr%d_f%d_v%d.mat',str,icontr,ifoi,v_cnt));
  senssel = find(stats.mask);
  
  if ~tval
    dfa_d1 = nanmean(all_par_cnt(:,:,contrasts(icontr,1))-all_par_cnt(:,:,contrasts(icontr,2)),2);
    r = max([min(abs(dfa_d1)) max(abs(dfa_d1))]);
    pars.scale=[-r r];
  else
    [~,~,~,s]=ttest(all_par_cnt(:,:,contrasts(icontr,1)),all_par_cnt(:,:,contrasts(icontr,2)),'dim',2)
    dfa_d1 = s.tstat;
    pars.scale = [-2.326 2.326];
  end
  
  figure; set(gcf,'color','white');
  
  if any(stats.mask)
    pars.markersize = 20;
    pars.markersel  = senssel;
  else
    pars.markersize = 0;
  end
  
  pars.cbar       = 0;
  pars.linewidth  = 9;
  pars.resolution = 300;
  
  showfield_colormap(dfa_d1,sa.locs_2D,pars);
  
  % add p-value information and colorlimits
  % ------------------------------------------------
  if isfield(stats,'negclusters') && ~isempty(stats.negclusters)
    m = stats.negclusters(1).prob;
    text(-1.05,-0.76,sprintf('- smallest neg. clust.: p = %.3f',m)); clear m
  else
    text(-1.05,-0.76,sprintf('- no neg. clust. found.'));
  end
  
  pars.markersize = 0;
  
  if isfield(stats,'posclusters') && ~isempty(stats.posclusters)
    m = stats.posclusters(1).prob;
    text(-1.05,-0.7,sprintf('- smallest pos. clust.: p = %.3f',m)); clear m
    
  else
    text(-1.05,-0.7,sprintf('- no pos. clust. found.'));
  end
  
  text(-1.05,-0.64,sprintf('- clim = [%.3f +%.3f]',pars.scale(1),pars.scale(2)));
  % ------------------------------------------------
  
  print(gcf,'-djpeg',sprintf('~/pconn_all/plots/pconn_topo_tsk_%s_c%d_f%d_v%d.jpg',str,icontr,ifoi,v_cnt))
  
  clear stats
  
  % ------------------------------------------------
  % PHARMA CONTRAST DURING REST
  % ------------------------------------------------
  
  load(sprintf('~/pconn_all/proc/all_sens_res_%s_clusterstat_contr%d_f%d_v%d.mat',str,icontr,ifoi,v_res));
  senssel = find(stats.mask);
  
  if ~tval
    dfa_d1 = nanmean(all_par_res(:,:,contrasts(icontr,1))-all_par_res(:,:,contrasts(icontr,2)),2);     
    r = max([min(abs(dfa_d1)) max(abs(dfa_d1))]);
    pars.scale=[-r r];
  else
    [~,~,~,s]=ttest(all_par_res(:,:,contrasts(icontr,1)),all_par_res(:,:,contrasts(icontr,2)),'dim',2)
    dfa_d1 = s.tstat;
    pars.scale = [-2.326 2.326];
  end
  
  figure; set(gcf,'color','white');

  
  if any(stats.mask)
    pars.markersize = 20;
    pars.markersel  = senssel;
  else
    pars.markersize = 0;
  end
  
  showfield_colormap(dfa_d1,sa.locs_2D,pars);
  
  % add p-value information and colorlimits
  % ------------------------------------------------
  if isfield(stats,'negclusters') && ~isempty(stats.negclusters)
    m = stats.negclusters(1).prob;
    text(-1.05,-0.76,sprintf('- smallest neg. clust.: p = %.3f',m)); clear m
  else
    text(-1.05,-0.76,sprintf('- no neg. clust. found.'));
  end
  
  if isfield(stats,'posclusters') && ~isempty(stats.posclusters)
    m = stats.posclusters(1).prob;
    text(-1.05,-0.7,sprintf('- smallest pos. clust.: p = %.3f',m)); clear m
  else
    text(-1.05,-0.7,sprintf('- no pos. clust. found.'));
  end
  
  text(-1.05,-0.64,sprintf('- clim = [%.3f +%.3f]',pars.scale(1),pars.scale(2)));
  % ------------------------------------------------
  
  print(gcf,'-djpeg',sprintf('~/pconn_all/plots/pconn_topo_rst_%s_c%d_f%d_v%d.jpg',str,icontr,ifoi,v_cnt))
  
  clear stats senssel
  close all
  
end

% ------------------------------------------------
% PLOT TASK VS REST
% ------------------------------------------------

figure; set(gcf,'color','white');

if ~tval
    dfa_d2 = nanmean(all_par_cnt(:,:,contrasts(icontr,1))-all_par_res(:,:,contrasts(icontr,2)),2);     
    r = max([min(abs(dfa_d1)) max(abs(dfa_d1))]);
    pars.scale=[-r r];
  else
    [~,~,~,s]=ttest(nanmean(all_par_cnt,3),nanmean(all_par_res,3),'dim',2)
    dfa_d2 = s.tstat;
    pars.scale = [-2.326 2.326];
  end

load(sprintf('~/pconn_all/proc/all_sens_%s_clusterstat_cnt-rest_f%d_v%d.mat',str,ifoi,v_cnt));
senssel = find(stats.mask);

if any(stats.mask)
  pars.markersize = 20;
  pars.markersel  = senssel;
else
  pars.markersize = 0;
end

showfield_colormap(dfa_d2,sa.locs_2D,pars);

% add p-value information and colorlimits
% ------------------------------------------------
if isfield(stats,'negclusters') && ~isempty(stats.negclusters)
  m = stats.negclusters(1).prob;
  text(-1.05,-0.76,sprintf('- smallest neg. clust.: p = %.3f',m)); clear m
else
  text(-1.05,-0.76,sprintf('- no neg. clust. found.'));
end

if isfield(stats,'posclusters') && ~isempty(stats.posclusters)
  m = stats.posclusters(1).prob;
  text(-1.05,-0.7,sprintf('- smallest pos. clust.: p = %.3f',m)); clear m
else
  text(-1.05,-0.7,sprintf('- no pos. clust. found.'));
end

text(-1.05,-0.64,sprintf('- clim = [%.3f +%.3f]',pars.scale(1),pars.scale(2)));
% ------------------------------------------------

print(gcf,'-djpeg',sprintf('~/pconn_all/plots/pconn_all_topo_%s_tsk-rst_f%d_v%d.jpg',str,ifoi,v_cnt))

%% PHARMA DIFFERENCE BETWEEN CONDITIONS 

tval = 1;
ifoi      = 4;
str       = 'pow';

if strcmp(str,'dfa')
  all_par_res = dfa_res_all(:,:,:,ifoi);
  all_par_cnt = dfa_cnt_all(:,:,:,ifoi);
elseif strcmp(str,'var')
  all_par_res = var_res_all(:,:,:,ifoi);
  all_par_cnt = var_cnt_all(:,:,:,ifoi);
elseif strcmp(str,'cvar')
  all_par_res = cvar_res_all(:,:,:,ifoi);
  all_par_cnt = cvar_cnt_all(:,:,:,ifoi);
elseif strcmp(str,'pow')
  all_par_res = log10(pow_res_all(:,:,:,ifoi));
  all_par_cnt = log10(pow_cnt_all(:,:,:,ifoi));
end


contrasts = [2 1; 3 1; 2 3];

for icontr = 1 : 3

  load(sprintf('~/pconn_all/proc/all_sens_%s_clusterstat_diffdiff_c%d_f%d_v%d.mat',str,icontr,ifoi,v_res));
  senssel = find(stats.mask);

  figure; set(gcf,'color','white');

  pars = [];
  pars.resolution = 300;
  pars.linewidth  = 9;
  pars.cbar       = 0;

  dfa_d1    = all_par_cnt(:,:,contrasts(icontr,1))-all_par_cnt(:,:,contrasts(icontr,2));
  dfa_d2    = all_par_res(:,:,contrasts(icontr,1))-all_par_res(:,:,contrasts(icontr,2));
  [~,~,~,d] = ttest(dfa_d1,dfa_d2,'dim',2);
  d         = d.tstat;

  pars.scale=[-2.326 2.326];



  if any(stats.mask)
    pars.markersize = 20;
    pars.markersel  = senssel;
  else
    pars.markersize = 0;
  end

  showfield_colormap(d,sa.locs_2D,pars);
  
  % add p-value information and colorlimits
  % ------------------------------------------------
  if isfield(stats,'negclusters') && ~isempty(stats.negclusters)
    m = stats.negclusters(1).prob;
    text(-1.05,-0.76,sprintf('- smallest neg. clust.: p = %.3f',m)); clear m
  else
    text(-1.05,-0.76,sprintf('- no neg. clust. found.'));
  end
  
  if isfield(stats,'posclusters') && ~isempty(stats.posclusters)
    m = stats.posclusters(1).prob;
    text(-1.05,-0.7,sprintf('- smallest pos. clust.: p = %.3f',m)); clear m
  else
    text(-1.05,-0.7,sprintf('- no pos. clust. found.'));
  end
  
  text(-1.05,-0.64,sprintf('- clim = [%.3f +%.3f]',pars.scale(1),pars.scale(2)));
  % ------------------------------------------------

  print(gcf,'-djpeg100',sprintf('~/pconn_all/plots/pconn_sens_clusterstat_%s_diffdiff_c%d_f%d_v%d.jpg',str,icontr,ifoi,v_cnt))
  clear stats
  
end

%% PLOT BAR GRAPHS AVERAGED ACROSS CHANNELS

ifoi      = 3;
str       = 'dfa';
contrasts = [2 1; 3 1; 2 3];

if strcmp(str,'dfa')
  all_par_res = dfa_res_all(:,:,:,ifoi);
  all_par_cnt = dfa_cnt_all(:,:,:,ifoi);
elseif strcmp(str,'var')
  all_par_res = var_res_all(:,:,:,ifoi);
  all_par_cnt = var_cnt_all(:,:,:,ifoi);
elseif strcmp(str,'cvar')
  all_par_res = cvar_res_all(:,:,:,ifoi);
  all_par_cnt = cvar_cnt_all(:,:,:,ifoi);
elseif strcmp(str,'pow')
  all_par_res = log10(pow_res_all(:,:,:,ifoi));
  all_par_cnt = log10(pow_cnt_all(:,:,:,ifoi));
end

par = squeeze(nanmean(all_par_res));


figure; set(gcf,'color','white'); hold on

m = mean(par);
s = std(par)/sqrt(size(par,1));

[~,p]=ttest(par(:,2),par(:,1));

plot(ones(1,size(par,1)),par(:,1),'rx','markersize',10);
plot(2*ones(1,size(par,1)),par(:,2),'bo','markersize',10);
plot(3*ones(1,size(par,1)),par(:,3),'md','markersize',10);

line([1.2 1.4],[m(1) m(1)],'linewidth',2,'color','r')
line([2.2 2.4],[m(2) m(2)],'linewidth',2,'color','b')
line([3.2 3.4],[m(3) m(3)],'linewidth',2,'color','m')

line([1.3 1.3],[m(1)-s(1) m(1)+s(1)],'linewidth',3,'color','r')
line([2.3 2.3],[m(2)-s(2) m(2)+s(2)],'linewidth',3,'color','b')
line([3.3 3.3],[m(3)-s(3) m(3)+s(3)],'linewidth',3,'color','m')

axis([-1 5 0.5 0.90]);

text(3,mean(m),sprintf('t-test: p = %.3f',p));

print(gcf,'-depsc2',sprintf('~/pconn_all/plots/pconn_topo_avgchanbar_%s_f%d_v%d.eps',str,ifoi,v_cnt))



%% COMPARE POWER IN SIGNIFICANT CLUSTERS

ifoi = 3;
str = 'dfa';

contrasts = [2 1; 3 1; 2 3];

if strcmp(str,'dfa')
  all_par_res = dfa_res_all(:,:,:,ifoi);
  all_par_cnt = dfa_cnt_all(:,:,:,ifoi);
elseif strcmp(str,'var')
  all_par_res = var_res_all(:,:,:,ifoi);
  all_par_cnt = var_cnt_all(:,:,:,ifoi);
elseif strcmp(str,'cvar')
  all_par_res = cvar_res_all(:,:,:,ifoi);
  all_par_cnt = cvar_cnt_all(:,:,:,ifoi);
elseif strcmp(str,'pow')
  all_par_res = log10(pow_res_all(:,:,:,ifoi));
  all_par_cnt = log10(pow_cnt_all(:,:,:,ifoi));
end

for icontr = 1 : 3
  
  
  load(sprintf('~/pconn_all/proc/all_sens_tsk_%s_clusterstat_contr%d_f%d_v%d.mat',str,icontr,ifoi,v_cnt));
  
  if any(stats.mask)
    
    senssel = find(stats.mask);
    
    figure; set(gcf,'color','white'); hold on;
    
    par = squeeze(nanmean(pow_cnt_all(senssel,:,:,ifoi),1));
    
    s = std(par,[],1)/sqrt(size(par,1));
    m = mean(par,1);
    
    figure; set(gcf,'color','white'); hold on;
    
    bar([1.0],m(1),.1,'facecolor','r','edgecolor','none');
    bar([1.2],m(2),.1,'facecolor','b','edgecolor','none');
    
    line([1 1],[m(1)-s(1) m(1)+s(1)],'linewidth',5)
    line([1.2 1.2],[m(2)-s(2) m(2)+s(2)],'linewidth',5)
    
    axis([0.8 1.4 min(max(m-s))-2*min(max(m-s)) max(max(m+s))+2*min(max(m-s))]);
    
    set(gca,'TickDir','out','XTick',[1 1.2],'XTickLabel',['P';'A']);
    
    [~,p] = ttest(par(:,contrasts(icontr,1)),par(:,contrasts(icontr,2)));
    
    text(1.1,max(max(m+s))+1.8*min(max(m-s)),sprintf('t-test: p = %.3f',p)); clear p
    
  else
    figure; set(gcf,'color','white'); hold on;
  end
  
  print(gcf,'-depsc2',sprintf('~/pconn_all/plots/pconn_clustersens_pow_%s_rst_c%d_f%d_v%d.eps',str,icontr,ifoi,v_cnt))
  
  clear stats m s par
  
  
  load(sprintf('~/pconn_all/proc/all_sens_res_%s_clusterstat_contr%d_f%d_v%d.mat',str,icontr,ifoi,v_cnt));
  
  if any(stats.mask)
    
    senssel = find(stats.mask);
    
    figure; set(gcf,'color','white'); hold on;
    
    par = squeeze(nanmean(pow_res_all(senssel,:,:,ifoi),1));
    
    s = std(par,[],1)/sqrt(size(par,1));
    m = mean(par,1);
    
    figure; set(gcf,'color','white'); hold on;
    
    bar([1.0],m(1),.1,'facecolor','r','edgecolor','none');
    bar([1.2],m(2),.1,'facecolor','b','edgecolor','none');
    
    line([1 1],[m(1)-s(1) m(1)+s(1)],'linewidth',5)
    line([1.2 1.2],[m(2)-s(2) m(2)+s(2)],'linewidth',5)
    
    axis([0.8 1.4 min(max(m-s))-2*min(max(m-s)) max(max(m+s))+2*min(max(m-s))]);
    
    set(gca,'TickDir','out','XTick',[1 1.2],'XTickLabel',['P';'A']);
    
    [~,p] = ttest(par(:,contrasts(icontr,1)),par(:,contrasts(icontr,2)));
    
    text(1.1,max(max(m+s))+1.8*min(max(m-s)),sprintf('t-test: p = %.3f',p)); clear p
    
  else
    figure; set(gcf,'color','white'); hold on;
  end
  
  print(gcf,'-depsc2',sprintf('~/pconn_all/plots/pconn_clustersens_pow_%s_rst_c%d_f%d_v%d.eps',str,icontr,ifoi,v_cnt))
  
  clear stats m s par
  
end

close all





















%% COMPARE MEASURES
% -------------------------------------

ifoi = 1;

for isubj = 1 : 18
  
  p_dfa  = nanmean(dfa_res_all(:,isubj,:,ifoi),3);
  p_var  = nanmean(var_res_all(:,isubj,:,ifoi),3);
  p_cvar = nanmean(cvar_res_all(:,isubj,:,ifoi),3);
  
  r_res(1,isubj) = corr(p_dfa,p_var);
  r_res(2,isubj) = corr(p_dfa,p_cvar);
  r_res(3,isubj) = corr(p_cvar,p_var);
  
  p_dfa  = nanmean(dfa_cnt_all(:,isubj,:,ifoi),3);
  p_var  = nanmean(var_cnt_all(:,isubj,:,ifoi),3);
  p_cvar = nanmean(cvar_cnt_all(:,isubj,:,ifoi),3);
  
  r_cnt(1,isubj) = corr(p_dfa,p_var);
  r_cnt(2,isubj) = corr(p_dfa,p_cvar);
  r_cnt(3,isubj) = corr(p_cvar,p_var);
  
  
  
end

for i = 1 : 3
  
  [~,p_res(i)] = ttest(r_res(:,i));
  [~,p_cnt(i)] = ttest(r_cnt(:,i));
  
end

%%
ifoi = 1;

clear r_res r_cnt

for isens = 1 : 268
  
  p_dfa  = nanmean(dfa_res_all(isens,:,:,ifoi),3);
  p_var  = nanmean(var_res_all(isens,:,:,ifoi),3);
  p_cvar = nanmean(cvar_res_all(isens,:,:,ifoi),3);
  
  r_res(1,isens) = corr(p_dfa',p_var');
  r_res(2,isens) = corr(p_dfa',p_cvar');
  r_res(3,isens) = corr(p_cvar',p_var');
  
  p_dfa  = nanmean(dfa_cnt_all(isens,:,:,ifoi),3);
  p_var  = nanmean(var_cnt_all(isens,:,:,ifoi),3);
  p_cvar = nanmean(cvar_cnt_all(isens,:,:,ifoi),3);
  
  r_cnt(1,isens) = corr(p_dfa',p_var');
  r_cnt(2,isens) = corr(p_dfa',p_cvar');
  r_cnt(3,isens) = corr(p_cvar',p_var');
  
end

figure; set(gcf,'color','white');
for i = 1 : 3
  
  subplot(3,2,i*2-1)
  
  pars.scale=[-1 1];
  pars.cbar = 0;
  pars.markersize = 0;
  pars.linewidth = 9;
  pars.resolution = 300;
  
  showfield_colormap(r_res(i,:),sa.locs_2D,pars);
  colormap(parula)
  
  subplot(3,2,i*2)
  
  showfield_colormap(r_cnt(i,:),sa.locs_2D,pars);
  colormap(parula)
  
end

%%