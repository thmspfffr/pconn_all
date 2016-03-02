%% CLUSTER BASED PERMUTATION TEST IN SOURCE SPACE

% Performs cluster based permutation test on source space data


clear all

% v = 1: cortex, eloreta, v = 2: cortex, lcmv, v = 3: coarse, eloreta

% --------------------------------------------------------
% VERSION 1
% --------------------------------------------------------
v         = 1;
v_rawdata = 2;
SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24];
v_out = 1;
% --------------------------------------------------------
% VERSION 2
% --------------------------------------------------------
% v         = 1;
% v_rawdata = 2;
% SUBJLIST  = [3 4 5 6 7 8 9 10 11 12 13 15 16 17 19 20 21 22 23 24];
% v_out     = 2;
% --------------------------------------------------------

str = 'cvar';

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

%% PLOT DFA

ord   = pconn_randomization;

for ifoi = 1 : 4
  for isubj = SUBJLIST
    for m = 1 : 3
      
      im = find(ord(isubj,:)==m);
      
      load(sprintf([outdir 'pconn_cnt_src_dfa_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,v));
      
      dfa_all_cnt(:,m,isubj,ifoi)  = nanmean(par.dfa,2);
      var_all_cnt(:,m,isubj,ifoi)  = nanmean(par.var,2);
      cvar_all_cnt(:,m,isubj,ifoi) = nanmean(par.cvar,2); clear par
      %
      load(sprintf(['~/pconn/proc/dfa/pconn_src_dfa_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,v));
      %
      dfa_all_res(:,m,isubj,ifoi)  = nanmean(par.dfa,2);
      var_all_res(:,m,isubj,ifoi)  = nanmean(par.var,2);
      cvar_all_res(:,m,isubj,ifoi) = nanmean(par.cvar,2); clear par
      
      
    end
  end
end

dfa_all_res  = dfa_all_res(:,:,SUBJLIST,:);
var_all_res  = var_all_res(:,:,SUBJLIST,:);
cvar_all_res = cvar_all_res(:,:,SUBJLIST,:);

dfa_all_cnt  = dfa_all_cnt(:,:,SUBJLIST,:);
var_all_cnt  = var_all_cnt(:,:,SUBJLIST,:);
cvar_all_cnt = cvar_all_cnt(:,:,SUBJLIST,:);

%%
for ifoi = 1 : 4
  % ifoi = 3;
  
  n =  get_neighbours(grid);
  
  para.nperm        = 10000;
  para.method       = 'dependentT';
  para.neigh        = n;
  para.minneigh     = 3;
  para.clusteralpha = 0.05;
  para.alpha        = 0.025;
  
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
  
  
  for icontr = 1 : 3
    
    z(:,:,1)=par_all_res(:,contrasts(icontr,1),:);
    z(:,:,2)=par_all_res(:,contrasts(icontr,2),:);
    
    stats = clust_perm(z,para);
    
    save(sprintf('~/pconn_all/proc/all_src_clusterstat_rst_%s_c%d_f%d_v%d.mat',str,icontr,ifoi,v_out),'stats');
    
  end
  
  clear z
  
  for icontr = 1 : 3
    
    z(:,:,1)=par_all_cnt(:,contrasts(icontr,1),:);
    z(:,:,2)=par_all_cnt(:,contrasts(icontr,2),:);
    
    stats = clust_perm(z,para);
    
    save(sprintf('~/pconn_all/proc/all_src_clusterstat_tsk_%s_c%d_f%d_v%d.mat',str,icontr,ifoi,v_out),'stats');
    
  end
  
  % REST VS TASK
  
  z(:,:,1)=nanmean(par_all_cnt,2);
  z(:,:,2)=nanmean(par_all_res,2);
  
  stats = clust_perm(z,para);
  
  save(sprintf('~/pconn_all/proc/all_src_clusterstat_tsk-rst_%s_c%d_f%d_v%d.mat',str,icontr,ifoi,v_out),'stats');
  
end


%% COMPUTE DOUBLE PHARMA CONTRAST


for ifoi = 1 : 4
  
  n =  get_neighbours(grid);
  
  para.nperm        = 10000;
  para.method       = 'dependentT';
  para.neigh        = n;
  para.minneigh     = 3;
  para.clusteralpha = 0.05;
  para.alpha = 0.025;
  
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
  
  for icontr = 1 : 3
    
    z(:,:,1)=par_all_cnt(:,contrasts(icontr,1),:)-par_all_cnt(:,contrasts(icontr,2),:);
    z(:,:,2)=par_all_res(:,contrasts(icontr,1),:)-par_all_res(:,contrasts(icontr,2),:);
        
    stats = clust_perm(z,para);
    
    save(sprintf('~/pconn_all/proc/all_src_clusterstat_diffdiff_%s_c%d_f%d_v%d.mat',str,icontr,ifoi,v_out),'stats');
    
  end
  
end
    














% SOURCE

% load sa_meg_template;
% 
% f=s.stat_pos{1}.freq;
% 
% f = [35];
% a = squeeze(nanmean(nanmean(z(:,contrats(icontr,1),:,f),3),4));
% b = squeeze(nanmean(nanmean(z(:,contrats(icontr,2),:,f),3),4));
% d = a-b;
% 
% m = zeros(size(z,1),1);
% m(s.stat_pos{1}.chan) = 1;
% 
% if strcmp(gridsize,'coarse')
%   grid  = sa_meg_template.grid_coarse;
% elseif strcmp(gridsize,'cortex')
%   grid  = sa_meg_template.grid_cortex3000;
% elseif strcmp(gridsize,'xcoarse')
%   grid  = sa_meg_template.grid_xcoarse;
% end
% 
% mri   = sa_meg_template.mri;
% vc    = sa_meg_template.vc;
% %  GENERATE SOURCE PLOT
% 
% para                  = [];
% para.orientation      = 'axial';
% para.colormaps        = {'jet'};
% % para.colorlimits = [0.5 0.7];
% h = figure; hold on
% set(h,'color',' k');
% 
% % ------------------------
% % PLACEBO - ATOMOXETINE
% % ------------------------
% r = max(abs(min(d)),abs(max(d)));
% para.colorlimits = [-r r];
% showmri_transp_v3(mri,para,[grid d]);
% 
% % print(gcf,'-djpeg100',sprintf('~/pconn/proc/plots/pconn_src_%s_raw_c%d_f%d_v%d.jpg',str,icontr,ifoi,v))
% 
% % d = d.*m;
% para.colorlimits = [-r r];
% showmri_transp_v3(mri,para,[grid d]);
% 
% % print(gcf,'-djpeg100',sprintf('~/pconn/proc/plots/pconn_src_%s_mask_c%d_f%d_v%d.jpg',str,icontr,ifoi,v))
% 
% %   end
% end
%%
%
% g1 = sa_meg_template.grid_cortex3000;
% g2 = sa_meg_template.cortex10K.vc;
% dd = .01;
% m2 = spatfiltergauss(m,g1,dd,g2);
%
% z2 = spatfiltergauss(d,g1,dd,g2);
% a.tri = sa_meg_template.cortex10K.tri;
%
% idx = [];
% idx = [idx; find(sa_meg_template.cortex10K.tri(:,1)>size(a.vc,1))];
% idx = [idx; find(sa_meg_template.cortex10K.tri(:,2)>size(a.vc,1))];
% idx = [idx; find(sa_meg_template.cortex10K.tri(:,3)>size(a.vc,1))];
%
% a.tri(unique(idx),:)=[];
%
% clear idx;
%
% idx = find(a.vc(:,1)>0);
%
% a.vc = a.vc(idx,:);
% a.normals = a.normals(idx,:);
% a.tri = delaunay(a.vc);
% % a.tri = a.tri(:,2:4);
% % a.tri = get_neighbours(a.vc);
% % z2(z2<0.0001) =0;
% para.colorlimits = [-0.03 0.03];
% % a.vc =
% figure;  pconn_showsurface(a,[],z2(idx))
%
%
% %%
%
% load ~/pconn_bttn/proc/pconn_bttn_mediandur_v4.mat
%
% d=nanmean(z(find(m),:,2),1)-nanmean(z(find(m),:,1),1)
% dd= dur(2,:)'-dur(1,:)''
%
% scatter(nanmean(nanmean(z,3),1),nanmean(dur))
%



