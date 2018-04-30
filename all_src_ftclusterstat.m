clear
%% LOAD DATA
SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/
addpath('/home/tpfeffer/Documents/MATLAB/fieldtrip-20160919/')

v= 2;
FOI = 2;
addpath ~/pconn/matlab
addpath /home/gnolte/meg_toolbox/toolbox/
ord   = pconn_randomization;

for ifoi = FOI
  for isubj = SUBJLIST
    fprintf('Loading data s%d f%d ...\n',isubj,ifoi);
    for m = 1 : 3
      
      im = find(ord(isubj,:)==m);
      
      load(sprintf(['~/pconn_cnt/proc/dfa/pconn_cnt_src_dfa_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,v));
      
      dfa_all_cnt(:,m,isubj,ifoi)  = nanmean(par.dfa,2);
      %
      load(sprintf(['~/pconn/proc/dfa/pconn_src_dfa_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,v));
      %
      dfa_all_res(:,m,isubj,ifoi)  = nanmean(par.dfa,2);
      
      load(sprintf(['~/pconn_cnt/proc/src/pconn_cnt_src_pow_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,1));
      pow_all_cnt(:,m,isubj,ifoi)  = nanmean(par.pow,2); clear par
      
      load(sprintf(['~/pconn/proc/src/pconn_src_pow_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,1));
      pow_all_res(:,m,isubj,ifoi)  = nanmean(par.pow,2); clear par
      
      
    end
  end
end

dfa_all_res  = squeeze(dfa_all_res(:,:,SUBJLIST,FOI));
pow_all_res  = squeeze(pow_all_res(:,:,SUBJLIST,FOI));

dfa_all_cnt  = squeeze(dfa_all_cnt(:,:,SUBJLIST,FOI));
pow_all_cnt  = squeeze(pow_all_cnt(:,:,SUBJLIST,FOI));

%%
load sa_meg_template.mat
grid = sa_meg_template.grid_coarse;

[xrange,yrange,zrange] = grid2cube(grid); % doesnt work with cortical grids, only regular
xgrid = length(xrange);
ygrid = length(yrange);
zgrid = length(zrange);
n_voxel_cube = xgrid*ygrid*zgrid;


mni_cube = zeros(n_voxel_cube,3);
counter = 0;
for ix = xgrid:-1:1
  for iy = ygrid:-1:1
    for iz = zgrid:-1:1
      counter = counter + 1;
      mni_cube(counter,:) = [xrange(ix) yrange(iy) zrange(iz)];
    end
  end
end

for i =1  : size(mni_cube,1)
  vox(i) = any(sum(mni_cube(i,:) == round(grid),2) == 3);
end

n_subjects = 28;
cond1.dimord = 'pos_subj_freq_time';
cond1.dim = [xgrid,ygrid,zgrid];
cond1.pos = mni_cube; % das sind die mni coordinates aller voxel im cube (also n_voxel_cube X 3)
cond1.inside = false(n_voxel_cube,1);
cond1.inside(vox) = true; % vox sind die indices der voxel, die wirklich im grid existieren. alle anderen voxel sind 0 und damit als outside definiert
cond1.time = [1];
cond1.freq = [1];
cond1.pow = nan(n_voxel_cube,n_subjects,1,1);
cond2 = cond1;

cfg                  = [];
cfg.dim              = [xgrid,ygrid,zgrid];
cfg.latency          = 'all';
cfg.frequency        = 'all';
cfg.parameter        = 'pow';
cfg.method           = 'montecarlo'; %analytic montecarlo
cfg.statistic        = 'depsamplesT';
cfg.computeprob      = 'yes';
cfg.correctm         = 'cluster'; %cluster fdr no
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.clustertail      = 0; %1 = right
cfg.tail             = 0; %1 = right
cfg.alpha            = 0.05;
cfg.numrandomization = 10000;

cfg.avgovertime      = 'yes';
cfg.avgoverfreq      = 'yes';

design = zeros(2,2*n_subjects);
design(1,:) = repmat(1:n_subjects,1,2);
design(2,:) = mod(floor([0:(2*n_subjects-1)]/(n_subjects/1)),2)+1;
cfg.design   = design;
cfg.uvar     = 1;
cfg.ivar     = 2;
%%
cond1.pow(vox,:) = squeeze(dfa_all_res(:,2,:));
cond2.pow(vox,:) = squeeze(dfa_all_res(:,1,:));
[stats1] = ft_sourcestatistics(cfg, cond1, cond2);

cond1.pow(vox,:) = squeeze(dfa_all_cnt(:,2,:));
cond2.pow(vox,:) = squeeze(dfa_all_cnt(:,1,:));
[stats2] = ft_sourcestatistics(cfg, cond1, cond2);

cond1.pow(vox,:) = squeeze(dfa_all_cnt(:,3,:));
cond2.pow(vox,:) = squeeze(dfa_all_cnt(:,1,:));
[stats3] = ft_sourcestatistics(cfg, cond1, cond2);

cond1.pow(vox,:) = (squeeze(dfa_all_res(:,2,:))+squeeze(dfa_all_cnt(:,2,:)))./2;
cond2.pow(vox,:) = (squeeze(dfa_all_res(:,1,:))+squeeze(dfa_all_cnt(:,1,:)))./2;
[stats_maineffect1] = ft_sourcestatistics(cfg, cond1, cond2);

cond1.pow(vox,:) = (squeeze(dfa_all_res(:,3,:))+squeeze(dfa_all_cnt(:,3,:)))./2;
cond2.pow(vox,:) = (squeeze(dfa_all_res(:,1,:))+squeeze(dfa_all_cnt(:,1,:)))./2;
[stats_maineffect2] = ft_sourcestatistics(cfg, cond1, cond2);

cond1.pow(vox,:) = squeeze(dfa_all_cnt(:,1,:));
cond2.pow(vox,:) = squeeze(dfa_all_res(:,1,:));
[stats_tvr1] = ft_sourcestatistics(cfg, cond1, cond2);

cond1.pow(vox,:) = squeeze(dfa_all_cnt(:,2,:));
cond2.pow(vox,:) = squeeze(dfa_all_res(:,2,:));
[stats_tvr2] = ft_sourcestatistics(cfg, cond1, cond2);

cond1.pow(vox,:) = squeeze(dfa_all_cnt(:,3,:));
cond2.pow(vox,:) = squeeze(dfa_all_res(:,3,:));
[stats_tvr3] = ft_sourcestatistics(cfg, cond1, cond2);

cond1.pow(vox,:) = squeeze(dfa_all_cnt(:,2,:))-squeeze(dfa_all_cnt(:,1,:));
cond2.pow(vox,:) = squeeze(dfa_all_res(:,2,:))-squeeze(dfa_all_res(:,1,:));
[stats_context] = ft_sourcestatistics(cfg, cond1, cond2);

save(sprintf(['~/neighb/' 'all_src_ftclusterstat.mat']),'stats1','stats2','stats3','stats_maineffect1','stats_maineffect2','stats_tvr1','stats_tvr2','stats_tvr3','stats_context','cfg');

error('!')

%%
load sa_meg_template;

grid  = sa_meg_template.grid_coarse;
mri   = sa_meg_template.mri;
vc    = sa_meg_template.vc;
dd    = .75;
g1    = grid;
g2 = sa_meg_template.cortex10K.vc;


d = stats_tvr1.stat(vox).*stats_tvr1.mask(vox);
par_interp = spatfiltergauss(d,g1,dd,g2);

para = [];
para.colorlimits = [-3 3];

cmap1  = cbrewer('seq', 'YlOrRd', 150,'pchip'); cmap1 = cmap1(50:1:end,:);
cmap2  = cbrewer('seq', 'YlGnBu', 150,'pchip'); cmap2 = cmap2(end:-1:50,:);
cmap  = [cmap2; ones(50,3); cmap1];

% PLOT RESULTS
para.filename = sprintf('~/neighb/pconn_src_dfa_tvr_icontr1.png');
tp_showsource(par_interp,cmap,sa_meg_template,para);

d = stats_tvr2.stat(vox).*stats_tvr2.mask(vox);
par_interp = spatfiltergauss(d,g1,dd,g2);
para.filename = sprintf('~/neighb/pconn_src_dfa_tvr_icontr2.png');
tp_showsource(par_interp,cmap,sa_meg_template,para);

d = stats_tvr3.stat(vox).*stats_tvr3.mask(vox);
par_interp = spatfiltergauss(d,g1,dd,g2);
para.filename = sprintf('~/neighb/pconn_src_dfa_tvr_icontr3.png');
tp_showsource(par_interp,cmap,sa_meg_template,para);

% % CONJUNCTION
% para = [];
% cmap  = cbrewer('seq', 'YlGn', 150,'pchip'); cmap2 = cmap2(end:-1:50,:);
% para.colorlimits = [0 10];
% d = stats1.mask(vox)&stats2.mask(vox).*10;
% par_interp = spatfiltergauss(d,g1,dd,g2);
% para.filename = sprintf('~/neighb/pconn_src_dfa_conj_icontr1.png');
% tp_showsource(par_interp,cmap,sa_meg_template,para);



%%
% para.neigh = n;
% gridsize = 'cortex';
% para.method       = 'dependentT';
% para.minneigh     = 2;
% para.clusteralpha = 0.05;
% para.alpha        = 0.025;
% para.nperm        = 10000;


% %%
%
for i = 11 : 20
  
  figure; hold on; set(gcf,'color','w');
  set(gcf,'position',[100 100 1400 800])
  orient(gcf,'landscape')
  
  subplot(2,3,1)
  plot3(grid(i,1),grid(i,2),grid(i,3),'k.','markersize',15); hold on
  plot3(grid(find(n(i,:)),1),grid(find(n(i,:)),2),grid(find(n(i,:)),3),'r.','markersize',15); hold on
  plot3(grid(:,1),grid(:,2),grid(:,3),'k.','markersize',3); hold on
  view([1 1 4]); axis off tight
  
  subplot(2,3,2)
  plot3(grid(i,1),grid(i,2),grid(i,3),'k.','markersize',15); hold on
  plot3(grid(find(n(i,:)),1),grid(find(n(i,:)),2),grid(find(n(i,:)),3),'r.','markersize',15); hold on
  plot3(grid(:,1),grid(:,2),grid(:,3),'k.','markersize',3); hold on
  view([1 1 -4]); axis off tight
  
  subplot(2,3,3)
  plot3(grid(i,1),grid(i,2),grid(i,3),'k.','markersize',15); hold on
  plot3(grid(find(n(i,:)),1),grid(find(n(i,:)),2),grid(find(n(i,:)),3),'r.','markersize',15); hold on
  plot3(grid(:,1),grid(:,2),grid(:,3),'k.','markersize',3); hold on
  view([-1 -1 4]); axis off tight
  
  subplot(2,3,4)
  plot3(grid(i,1),grid(i,2),grid(i,3),'k.','markersize',15); hold on
  plot3(grid(find(n(i,:)),1),grid(find(n(i,:)),2),grid(find(n(i,:)),3),'r.','markersize',15); hold on
  plot3(grid(:,1),grid(:,2),grid(:,3),'k.','markersize',3); hold on
  view([1 -1 -4]); axis off tight
  
  subplot(2,3,4)
  plot3(grid(i,1),grid(i,2),grid(i,3),'k.','markersize',15); hold on
  plot3(grid(find(n(i,:)),1),grid(find(n(i,:)),2),grid(find(n(i,:)),3),'r.','markersize',15); hold on
  plot3(grid(:,1),grid(:,2),grid(:,3),'k.','markersize',3); hold on
  view([0 -1 -4]); axis off tight
  
  subplot(2,3,5)
  plot3(grid(i,1),grid(i,2),grid(i,3),'k.','markersize',15); hold on
  plot3(grid(find(n(i,:)),1),grid(find(n(i,:)),2),grid(find(n(i,:)),3),'r.','markersize',15); hold on
  plot3(grid(:,1),grid(:,2),grid(:,3),'k.','markersize',3); hold on
  view([0 90 0]); axis off tight
  
  subplot(2,3,6)
  plot3(grid(i,1),grid(i,2),grid(i,3),'k.','markersize',15); hold on
  plot3(grid(find(n(i,:)),1),grid(find(n(i,:)),2),grid(find(n(i,:)),3),'r.','markersize',15); hold on
  plot3(grid(:,1),grid(:,2),grid(:,3),'k.','markersize',3); hold on
  view([1 0 0]); axis off tight
  
  print(gcf,'-dpdf',sprintf('~/neighb/newgrid_n%d.pdf',i))
  close
  
end
%
%

%% NEW GRIDS
load sa_meg_template.mat
grid = sa_meg_template.grid_cortex3000;
gridsize = 'cortex';
para.method       = 'dependentT';
para.minneigh     = 2;
para.clusteralpha = 0.05;
para.alpha        = 0.05;
para.nperm        = 10000;
para.tail = 0;
para.paried = 1;
para.grid = grid;
g1    = grid;
dd    = .75;
para.gridtype = 'regular'
g2 = sa_meg_template.cortex10K.vc;

para.dim = max(grid)-min(grid)+1;

n =  get_neighbours_new(grid,2);
para.neigh        = n;

z(:,:,1) = squeeze(dfa_all_res(:,2,:));
z(:,:,2) = squeeze(dfa_all_res(:,1,:));

stats11 = tp_clusterperm(z,para)

z(:,:,1) = squeeze(dfa_all_cnt(:,1,:));
z(:,:,2) = squeeze(dfa_all_res(:,1,:));

stats1 = tp_clusterperm(z,para)


%%
addpath ~/Documents/MATLAB/cbrewer/cbrewer/
para = [];
para.colorlimits = [-3 3];
grid = sa_meg_template.grid_coarse;

cmap1  = cbrewer('seq', 'YlOrRd', 150,'pchip'); cmap1 = cmap1(50:1:end,:);
cmap2  = cbrewer('seq', 'YlGnBu', 150,'pchip'); cmap2 = cmap2(end:-1:50,:);
cmap  = [cmap2; ones(50,3); cmap1];

d = stats11.stat.*stats11.mask;
par_interp = spatfiltergauss(d,g1,dd,g2);
para.filename = sprintf('~/neighb/pconn_src_dfa_conj_icontr1.png');
tp_showsource(par_interp,cmap,sa_meg_template,para);
%

%%
n_subjects = 28;
clear cond1 cond2
cond1.dimord = 'pos_subj_freq_time';
% cond1.dim = [3000 1];
cond1.pos = grid; % das sind die mni coordinates aller voxel im cube (also n_voxel_cube X 3)
% cond1.inside = false(n_voxel_cube,1);
% cond1.inside(vox) = true; % vox sind die indices der voxel, die wirklich im grid existieren. alle anderen voxel sind 0 und damit als outside definiert
cond1.time = [1];
cond1.freq = [1];

% cond1.pow = nan(n_voxel_cube,n_subjects,1,1);
cond2 = cond1;


cfg                  = [];
% cfg.dim              = [3000 1];
% cfg.
cfg.latency          = 'all';
cfg.frequency        = 'all';
cfg.parameter        = 'pow';
cfg.method           = 'montecarlo'; %analytic montecarlo
cfg.statistic        = 'depsamplesT';
cfg.computeprob      = 'yes';
cfg.correctm         = 'cluster'; %cluster fdr no
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.clustertail      = 0; %1 = right
cfg.tail             = 0; %1 = right
cfg.alpha            = 0.025;
cfg.numrandomization = 10000;
cfg.connectivity = n;
cfg.minnbchan = 0;

cfg.avgovertime      = 'yes';
cfg.avgoverfreq      = 'yes';

design = zeros(2,2*n_subjects);
design(1,:) = repmat(1:n_subjects,1,2);
design(2,:) = mod(floor([0:(2*n_subjects-1)]/(n_subjects/1)),2)+1;
cfg.design   = design;
cfg.uvar     = 1;
cfg.ivar     = 2;
%

cond1.pow = squeeze(dfa_all_cnt(:,1,:));
cond2.pow = squeeze(dfa_all_res(:,1,:));
[stats_tvr1] = ft_sourcestatistics(cfg, cond1, cond2);

cond1.pow = squeeze(dfa_all_cnt(:,2,:));
cond2.pow = squeeze(dfa_all_res(:,2,:));
[stats_tvr2] = ft_sourcestatistics(cfg, cond1, cond2);

cond1.pow = squeeze(dfa_all_cnt(:,3,:));
cond2.pow = squeeze(dfa_all_res(:,3,:));
[stats_tvr3] = ft_sourcestatistics(cfg, cond1, cond2);

cond1.pow = squeeze(dfa_all_res(:,2,:));
cond2.pow = squeeze(dfa_all_res(:,1,:));
[stats1] = ft_sourcestatistics(cfg, cond1, cond2);

cond1.pow = squeeze(dfa_all_cnt(:,2,:));
cond2.pow = squeeze(dfa_all_cnt(:,1,:));
[stats2] = ft_sourcestatistics(cfg, cond1, cond2);

cond1.pow = (squeeze(dfa_all_cnt(:,2,:))+squeeze(dfa_all_res(:,2,:)))./2;
cond1.pow = (squeeze(dfa_all_cnt(:,1,:))+squeeze(dfa_all_res(:,1,:)))./2;
[stats_maineffect] = ft_sourcestatistics(cfg, cond1, cond2);

cond1.pow = (squeeze(dfa_all_cnt(:,2,:))+squeeze(dfa_all_res(:,2,:)))./2;
cond1.pow = (squeeze(dfa_all_cnt(:,1,:))+squeeze(dfa_all_res(:,1,:)))./2;
[stats_maineffect] = ft_sourcestatistics(cfg, cond1, cond2);

save(sprintf(['~/neighb/' 'all_src_ftclusterstat.mat']),'stats1','stats2','stats_maineffect','stats_tvr1','stats_tvr2','stats_tvr3','cfg');
%%

addpath ~/Documents/MATLAB/cbrewer/cbrewer/
para = [];
para.colorlimits = [-3 3];
grid = sa_meg_template.grid_cortex3000;

cmap1  = cbrewer('seq', 'YlOrRd', 150,'pchip'); cmap1 = cmap1(50:1:end,:);
cmap2  = cbrewer('seq', 'YlGnBu', 150,'pchip'); cmap2 = cmap2(end:-1:50,:);
cmap  = [cmap2; ones(50,3); cmap1];

d = stats_maineffect.stat.*stats_maineffect.mask;
par_interp = spatfiltergauss(d,g1,dd,g2);
para.filename = sprintf('~/neighb/pconn_src_dfa_main_icontr1.png');
tp_showsource(par_interp,cmap,sa_meg_template,para); close

d = stats1.stat.*stats1.mask;
par_interp = spatfiltergauss(d,g1,dd,g2);
para.filename = sprintf('~/neighb/pconn_src_dfa_icontr1.png');
tp_showsource(par_interp,cmap,sa_meg_template,para); close

d = stats2.stat.*stats2.mask;
par_interp = spatfiltergauss(d,g1,dd,g2);
para.filename = sprintf('~/neighb/pconn_src_dfa_icontr2.png');
tp_showsource(par_interp,cmap,sa_meg_template,para); close

d = stats_tvr1.stat.*stats_tvr1.mask;
par_interp = spatfiltergauss(d,g1,dd,g2);
para.filename = sprintf('~/neighb/pconn_src_dfa_tvr_icontr1.png');
tp_showsource(par_interp,cmap,sa_meg_template,para); close

d = stats_tvr2.stat.*stats_tvr2.mask;
par_interp = spatfiltergauss(d,g1,dd,g2);
para.filename = sprintf('~/neighb/pconn_src_dfa_tvr_icontr2.png');
tp_showsource(par_interp,cmap,sa_meg_template,para); close

d = stats_tvr3.stat.*stats_tvr3.mask;
par_interp = spatfiltergauss(d,g1,dd,g2);
para.filename = sprintf('~/neighb/pconn_src_dfa_tvr_icontr3.png');
tp_showsource(par_interp,cmap,sa_meg_template,para); close
%

s =stats1.mask&stats2.mask; 
s = spatfiltergauss(s,g1,dd,g2);

cmap  = cbrewer('seq', 'Greens', 200,'pchip'); cmap(1:50,:)=ones(50,3);
para = [] ;
para.colorlimits = [0 0.3];

para.filename = sprintf('~/neighb/pconn_src_dfa_conj.png');
tp_showsource(s,cmap,sa_meg_template,para); close

[~,p]=ttest(squeeze(mean(dfa_all_cnt(stats1.mask&stats2.mask,1,:))),squeeze(mean(dfa_all_res(stats1.mask&stats2.mask,1,:))))
[~,p]=ttest(squeeze(mean(dfa_all_cnt(stats1.mask&stats2.mask,2,:))),squeeze(mean(dfa_all_res(stats1.mask&stats2.mask,2,:))))
[~,p]=ttest(squeeze(mean(dfa_all_cnt(stats1.mask&stats2.mask,3,:))),squeeze(mean(dfa_all_res(stats1.mask&stats2.mask,3,:))))

