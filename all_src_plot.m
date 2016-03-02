%% COMPUTE WHOLE BRAIN DFA IN SOURCE SPACE
% pconn_src_dfa

clear all

% v = 1: cortex, eloreta, v = 2: cortex, lcmv, v = 3: coarse, eloreta

% --------------------------------------------------------
% VERSION 1
% --------------------------------------------------------
v         = 1;
v_stat    = 1;
v_rawdata = 2;
SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24];
% --------------------------------------------------------

if     v == 1 
  gridsize = 'cortex';
elseif v == 2
  gridsize = 'cortex';
elseif v == 3
  gridsize = 'coarse';
end

restoredefaultpath

addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/
addpath /home/tpfeffer/Documents/MATLAB/fieldtrip-20130925/
addpath /home/tpfeffer/Documents/MATLAB/toolboxes/NBT-ReleaseNBTRC4a/
addpath ~/pconn/matlab/
addpath ~/Documents/MATLAB/cbrewer/cbrewer/

ft_defaults

indir     = '/home/tpfeffer/pconn_cnt/proc/src/';
outdir    = '/home/tpfeffer/pconn_cnt/proc/dfa/';
plotdir   = '/home/tpfeffer/pconn_cnt/proc/plots/';

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


% idx = find(grid(:,1)>mean([-7.1200 7.1268]));

% grid(idx,:)=[]; g1 = grid;

% idx = find(g2(:,1)>mean([-7.1200 7.1268]));
% g2(idx,:) = [];

%% PLOT DFA
   
str = 'dfa'; 
ord   = pconn_randomization;

for ifoi = 3 : 3
  for isubj = SUBJLIST
    for m = 1 : 3

      im = find(ord(isubj,:)==m);

%       load(sprintf([outdir 'pconn_cnt_src_dfa_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,v));
%      
%       dfa_all_cnt(:,m,isubj,ifoi)  = nanmean(par.dfa,2);
%       var_all_cnt(:,m,isubj,ifoi)  = nanmean(par.var,2);
%      	cvar_all_cnt(:,m,isubj,ifoi) = nanmean(par.cvar,2); clear par
      try
      load(sprintf(['~/pconn/proc/dfa/pconn_src_dfa_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,v));
     
      dfa_all_res(:,m,isubj,ifoi)  = nanmean(par.dfa,2);
      var_all_res(:,m,isubj,ifoi)  = nanmean(par.var,2);
      cvar_all_res(:,m,isubj,ifoi) = nanmean(par.cvar,2);
      amp_all_res(:,m,isubj,ifoi)  = nanmean(par.amp,2); clear par
      catch me
        disp(sprintf('err s%d f%d m%d',isubj,ifoi,im))
      end

    end
  end
end

dfa_all_res  = dfa_all_res(:,:,SUBJLIST,:);
var_all_res  = var_all_res(:,:,SUBJLIST,:);
cvar_all_res = cvar_all_res(:,:,SUBJLIST,:);
amp_all_res  = amp_all_res(:,:,SUBJLIST,:);
% dfa_all_cnt = dfa_all_cnt(:,:,SUBJLIST,:);
% var_all_cnt  = var_all_cnt(:,:,SUBJLIST,:);
% cvar_all_cnt = cvar_all_cnt(:,:,SUBJLIST,:);

%% PHARMA COMPARISON

str = 'dfa';
ifoi = 3;
mask = 1;

contrasts = [2 1; 3 1; 2 3];
viewdir   = [-.5 -.5 .5; .5 -.5 .5; .5 .5 .5; -.5 .5 .5];

if strcmp(str,'dfa')
  par_all_res = dfa_all_res(:,:,:,ifoi);
%   par_all_cnt = dfa_all_cnt(:,:,:,ifoi);
elseif strcmp(str,'var')
  par_all_res = var_all_res(:,:,:,ifoi);
%   par_all_cnt = var_all_cnt(:,:,:,ifoi);
elseif strcmp(str,'cvar')
  par_all_res = cvar_all_res(:,:,:,ifoi);
%   par_all_cnt = cvar_all_cnt(:,:,:,ifoi);
elseif strcmp(str,'amp')
  par_all_res = amp_all_res(:,:,:,ifoi);
%   par_all_cnt = cvar_all_cnt(:,:,:,ifoi);
end

cmap = cbrewer('div', 'RdBu', 100,'pchip');% colormap(autumn)
cmap = cmap(end:-1:1,:);

for icontr = 1 : 3
  
%   mask = 0;
  % -------------------------------------------------------
  % TASK
  % -------------------------------------------------------
  
	d = squeeze(nanmean(par_all_cnt(:,contrasts(icontr,1),:),3))-squeeze(nanmean(par_all_cnt(:,contrasts(icontr,2),:),3));
% 
  [~,~,~,tmp] = ttest(par_all_cnt(:,contrasts(icontr,1),:),par_all_cnt(:,contrasts(icontr,2),:),'dim',3);
  d = tmp.tstat; clear tmp
  
  if mask && exist(sprintf('~/pconn_all/proc/all_src_clusterstat_tsk_%s_c%d_f%d_v%d.mat',str,icontr,ifoi,v_stat))
    load(sprintf('~/pconn_all/proc/all_src_clusterstat_tsk_%s_c%d_f%d_v%d.mat',str,icontr,ifoi,v_stat))
    d(~stats.mask)=eps;
    mask=1;
  end

  par_interp = spatfiltergauss(d,g1,dd,g2);

  figure; set(gcf,'color','white'); hold on;

  para = [] ;
  r = max(abs(min(d)),abs(max(d)));
  para.colorlimits = [-2.326 2.326];

  for iplot = 1 : 4

    subplot(2,2,iplot)
    para.myviewdir = viewdir(iplot,:);
    pconn_showsurface(sa_meg_template.cortex10K,para,par_interp)
    colormap(cmap)
    camlight headlight
axis tight
if iplot == 1 
  rotate(gca,[0 0 1],-90)
end
  end
%   
  print(gcf,'-djpeg100',sprintf('~/pconn_all/plots/pconn_src_tsk_%s_mask%d_c%d_f%d_v%d.jpg',str,mask,icontr,ifoi,v_stat)) 
  
  clear para r d stats
%   mask = 0;
  
  % -------------------------------------------------------
  % RESTING STATE
  % -------------------------------------------------------
    
  d = squeeze(nanmean(par_all_res(:,contrasts(icontr,1),:),3))-squeeze(nanmean(par_all_res(:,contrasts(icontr,2),:),3));
% 
  [~,~,~,tmp] = ttest(par_all_res(:,contrasts(icontr,1),:),par_all_res(:,contrasts(icontr,2),:),'dim',3);
  d = tmp.tstat; clear tmp
  
  if mask && exist(sprintf('~/pconn_all/proc/all_src_clusterstat_rst_%s_c%d_f%d_v%d.mat',str,icontr,ifoi,v_stat))
    load(sprintf('~/pconn_all/proc/all_src_clusterstat_rst_%s_c%d_f%d_v%d.mat',str,icontr,ifoi,v_stat))
    d(~stats.mask)=eps;
    mask=1;
  end
  
  par_interp = spatfiltergauss(d,g1,dd,g2);
   
  figure; set(gcf,'color','white'); hold on;

  para = [] ;
  r = max(abs(min(d)),abs(max(d)));
  para.colorlimits = [-2.326 2.326];

  for iplot = 1 : 4

    subplot(2,2,iplot)
    para.myviewdir = viewdir(iplot,:);
    pconn_showsurface(sa_meg_template.cortex10K,para,par_interp)
    colormap(cmap)
    camlight headlight
    
  end
  
  print(gcf,'-djpeg100',sprintf('~/pconn_all/plots/pconn_src_rst_%s_mask%d_c%d_f%d_v%d.jpg',str,mask,icontr,ifoi,v_stat)) 
  clear stats
end

% MEAN VALUES 
cmap = cbrewer('seq', 'PuRd', 100,'pchip');% colormap(autumn)

d = squeeze(nanmean(nanmean(par_all_res,3),2));

par_interp = spatfiltergauss(d,g1,dd,g2);

figure; set(gcf,'color','white'); hold on;

para = [] ;
para.colorlimits = [min(d) max(d)];

for iplot = 1 : 4
  
  subplot(2,2,iplot)
  para.myviewdir = viewdir(iplot,:);
  pconn_showsurface(sa_meg_template.cortex10K,para,par_interp)
  colormap(hot)
  camlight headlight
  
end

print(gcf,'-djpeg100',sprintf('~/pconn_all/plots/pconn_src_rst_avg_%s_c%d_f%d_v%d.jpg',str,icontr,ifoi,v))

%% COMPARISON TASK VS REST

for ifoi = 1 : 4
  
  str = 'var';

  viewdir   = [-.5 -.5 .5; .5 -.5 .5; .5 .5 .5; -.5 .5 .5];

  if strcmp(str,'dfa')
    par_all_res = dfa_all_res(:,:,:,ifoi);
    par_all_cnt = dfa_all_cnt(:,:,:,ifoi);
  elseif strcmp(str,'var')
    par_all_res = var_all_res(:,:,:,ifoi);
    par_all_cnt = var_all_cnt(:,:,:,ifoi);
  elseif strcmp(str,'cvar')
    par_all_res = cvar_all_res(:,:,:,ifoi);
    par_all_cnt = cvar_all_cnt(:,:,:,ifoi);
  end

  cmap = cbrewer('div', 'RdBu', 100,'pchip');% colormap(autumn)
  cmap = cmap(end:-1:1,:);

% -------------------------------------------------------
% TASK VS REST
% -------------------------------------------------------
  mask = 1;

  d = squeeze(nanmean(nanmean(par_all_cnt,3),2))-squeeze(nanmean(nanmean(par_all_res,3),2));
  [~,~,~,tmp] = ttest(nanmean(par_all_cnt,2),nanmean(par_all_res,2),'dim',3);
  d = tmp.tstat; clear tmp
  
  if mask && exist(sprintf('~/pconn_all/proc/all_src_clusterstat_tsk-rst_%s_c%d_f%d_v%d.mat',str,3,ifoi,v_stat))
    load(sprintf('~/pconn_all/proc/all_src_clusterstat_tsk-rst_%s_c%d_f%d_v%d.mat',str,3,ifoi,v_stat))
    d(~stats.mask)=eps;
    mask=1;
  end

  par_interp = spatfiltergauss(d,g1,dd,g2);

  figure; set(gcf,'color','white'); hold on;

  para = [] ;
    para.colorlimits = [-2.326 2.326];

  for iplot = 1 : 4

    subplot(2,2,iplot)
    para.myviewdir = viewdir(iplot,:);
    pconn_showsurface(sa_meg_template.cortex10K,para,par_interp)
    colormap(cmap)
    camlight headlight

  end
  %
  print(gcf,'-djpeg100',sprintf('~/pconn_all/plots/pconn_src_%s_tsk-rst_mask%d_f%d_v%d.jpg',str,mask,ifoi,v_stat))
end
clear para r d stats

%% COMPARE PHARMA ACROSS TASK/REST
mask = 0;
str  = 'dfa';

for ifoi = 1 : 4
  
  contrasts = [2 1; 3 1; 2 3];
  viewdir   = [-.5 -.5 .5; .5 -.5 .5; .5 .5 .5; -.5 .5 .5];

  if strcmp(str,'dfa')
    par_all_res = dfa_all_res(:,:,:,ifoi);
    par_all_cnt = dfa_all_cnt(:,:,:,ifoi);
  elseif strcmp(str,'var')
    par_all_res = var_all_res(:,:,:,ifoi);
    par_all_cnt = var_all_cnt(:,:,:,ifoi);
  elseif strcmp(str,'cvar')
    par_all_res = cvar_all_res(:,:,:,ifoi);
    par_all_cnt = cvar_all_cnt(:,:,:,ifoi);
  end

  cmap = cbrewer('div', 'RdBu', 100,'pchip');% colormap(autumn)
  cmap = cmap(end:-1:1,:);

  for icontr = 1 : 3

    d1 = squeeze(par_all_res(:,contrasts(icontr,1),:)-par_all_res(:,contrasts(icontr,2),:));

    d2 = squeeze(par_all_cnt(:,contrasts(icontr,1),:)-par_all_cnt(:,contrasts(icontr,2),:));

    [~,~,~,tmp] = ttest(d2,d1,'dim',2);

    d = tmp.tstat;

    if mask && exist(sprintf('~/pconn_all/proc/all_src_clusterstat_diffdiff_%s_c%d_f%d_v%d.mat',str,icontr,ifoi,v_stat))
      load(sprintf('~/pconn_all/proc/all_src_clusterstat_diffdiff_%s_c%d_f%d_v%d',str,icontr,ifoi,v_stat))
      d(~stats.mask)=eps;
      mask=1;
    end

    par_interp = spatfiltergauss(d,g1,dd,g2);

    figure; set(gcf,'color','white'); hold on;

    para = [] ;
    para.colorlimits = [-2.326 2.326];

    for iplot = 1 : 4

      subplot(2,2,iplot)
      para.myviewdir = viewdir(iplot,:);
      pconn_showsurface(sa_meg_template.cortex10K,para,par_interp)
      colormap(cmap)
      camlight headlight

    end

    print(gcf,'-djpeg',sprintf('~/pconn_all/plots/pconn_src_%s_tskrstpharm_mask%d_c%d_f%d_v%d.jpg',str,mask,icontr,ifoi,v))
%
  end
end


%% CORRELATE ALL MEASURES 

ifoi = 4;
str = 'dfa';

viewdir   = [-.5 -.5 .5; .5 -.5 .5; .5 .5 .5; -.5 .5 .5];

for isubj = 1 : length(SUBJLIST)
  
%   p1        = nanmean(dfa_all_res(:,:,isubj,ifoi),2);
%   p2        = nanmean(var_all_res(:,:,isubj,ifoi),2);
%   p3        = nanmean(cvar_all_res(:,:,isubj,ifoi),2);
%   
%   r1(isubj) = corr(p1,p2);
%   r2(isubj) = corr(p1,p3);
%   r3(isubj) = corr(p2,p3);
%   
  p1        = nanmean(dfa_all_cnt(:,:,isubj,ifoi),2);
  p2        = nanmean(var_all_cnt(:,:,isubj,ifoi),2);
  p3        = nanmean(cvar_all_cnt(:,:,isubj,ifoi),2);
  
  r1(isubj) = corr(p1,p2);
  r2(isubj) = corr(p1,p3);
  r3(isubj) = corr(p2,p3);
  
  
end

zerodist = zeros(1,length(SUBJLIST));

[~,p1] = ttest(r1,zerodist);
[~,p2] = ttest(r2,zerodist);
[~,p3] = ttest(r3,zerodist);


%% COMPARE MEASURES

ifoi = 4;

p_dfa  = nanmean(nanmean(dfa_all(:,:,:,ifoi),3),2);
p_var  = nanmean(nanmean(var_all(:,:,:,ifoi),3),2);
p_cvar = nanmean(nanmean(cvar_all(:,:,:,ifoi),3),2);

r(1) = corr(p_dfa,p_var);
r(2) = corr(p_dfa,p_cvar);
r(3) = corr(p_cvar,p_var);

figure; set(gcf,'color','white');

subplot(3,1,1)

scatter(p_dfa,p_var,20,'facecolor','k','markeredgecolor','w')
xlabel('DFA'); ylabel('Variance')
axis square
title(sprintf('r = %.2f',r(1)))

subplot(3,1,2)

scatter(p_dfa,p_cvar,20,'facecolor','k','markeredgecolor','w')
xlabel('DFA'); ylabel('Coef. of var.')
axis square
title(sprintf('r = %.2f',r(2)))

subplot(3,1,3)

scatter(p_var,p_cvar,20,'facecolor','k','markeredgecolor','w')
xlabel('Variance'); ylabel('Coef. of var.')
axis square
title(sprintf('r = %.2f',r(3)))

set(gcf,'Position',[50 50 800 1200])

print(gcf,'-djpeg100',sprintf('~/pconn_cnt/plots/pconn_cnt_src_dfavarcvarcorr_f%d_v%d.jpg',ifoi,v))

close 

p_dfa  = nanmean(nanmean(dfa_all_res(:,:,:,ifoi),3),2);
p_var  = nanmean(nanmean(var_all_res(:,:,:,ifoi),3),2);
p_cvar = nanmean(nanmean(cvar_all_res(:,:,:,ifoi),3),2);

r(1) = corr(p_dfa,p_var);
r(2) = corr(p_dfa,p_cvar);
r(3) = corr(p_cvar,p_var);

figure; set(gcf,'color','white');

subplot(3,1,1)

scatter(p_dfa,p_var,20,'facecolor','k','markeredgecolor','w')
xlabel('DFA'); ylabel('Variance')
axis square
title(sprintf('r = %.2f',r(1)))

subplot(3,1,2)

scatter(p_dfa,p_cvar,20,'facecolor','k','markeredgecolor','w')
xlabel('DFA'); ylabel('Coef. of var.')
axis square
title(sprintf('r = %.2f',r(2)))

subplot(3,1,3)

scatter(p_var,p_cvar,20,'facecolor','k','markeredgecolor','w')
xlabel('Variance'); ylabel('Coef. of var.')
axis square
title(sprintf('r = %.2f',r(3)))

set(gcf,'Position',[50 50 800 1200])

print(gcf,'-djpeg100',sprintf('~/pconn/proc/plots/pconn_src_dfavarcvarcorr_f%d_v%d.jpg',ifoi,v))

close 



%%  CORRELATION OF MEASURES ACROSS CONDITIONS

ifoi = 1;

p1_dfa  = nanmean(nanmean(dfa_all_res(:,:,:,ifoi),3),2);
p1_var  = nanmean(nanmean(var_all_res(:,:,:,ifoi),3),2);
p1_cvar = nanmean(nanmean(cvar_all_res(:,:,:,ifoi),3),2);

p2_dfa  = nanmean(nanmean(dfa_all(:,:,:,ifoi),3),2);
p2_var  = nanmean(nanmean(var_all(:,:,:,ifoi),3),2);
p2_cvar = nanmean(nanmean(cvar_all(:,:,:,ifoi),3),2);

r(1) = corr(p1_dfa,p2_dfa);
r(2) = corr(p1_var,p2_var);
r(3) = corr(p1_cvar,p2_cvar);

figure; set(gcf,'color','white');

subplot(3,1,1)

scatter(p1_dfa,p2_dfa,20,'facecolor','k','markeredgecolor','w')
xlabel('DFA (rest)'); ylabel('DFA (task)')
axis square
title(sprintf('r = %.2f',r(1)))

subplot(3,1,2)

scatter(log10(p1_var),log10(p2_var),20,'facecolor','k','markeredgecolor','w')
xlabel('Variance (rest)'); ylabel('Variance (task)')
axis square
title(sprintf('r = %.2f',r(2)))

subplot(3,1,3)

scatter(p1_cvar,p2_cvar,20,'facecolor','k','markeredgecolor','w')
xlabel('C. Variance (rest)'); ylabel('C. Variance (task)')
axis square
title(sprintf('r = %.2f',r(3)))

set(gcf,'Position',[50 50 800 1200])

print(gcf,'-djpeg100',sprintf('~/pconn/proc/plots/pconn_src_resttaskcorr_f%d_v%d.jpg',ifoi,v))

close 

%% 

%% SHOW SURFACE

g1 = sa_meg_template.grid_cortex3000;
g2 = sa_meg_template.cortex10K.vc;
dd = .01;
m2 = spatfiltergauss(m,g1,dd,g2);

para.colorlimits = [-0.03 0.03];

figure;  pconn_showsurface(a,[],z2(idx))

