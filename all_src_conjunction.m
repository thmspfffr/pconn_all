%% COMPUTE WHOLE BRAIN DFA IN SOURCE SPACE
% pconn_src_dfa

clear all

% v = 1: cortex, eloreta, v = 2: cortex, lcmv, v = 3: coarse, eloreta

% --------------------------------------------------------
% VERSION 1
% --------------------------------------------------------
v         = 1;
v_stat    = 3;
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

% ft_defaults

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

%% PHARMA COMPARISON

str = 'dfa';
ifoi = 3;
mask = 1;
v_pup = 16;
v_stat = 1;

% cmap = cbrewer('div', 'RdBu', 100,'pchip');% colormap(autumn)
cmap = [1 1 1; 1 0 0; 1 0 0; 1 0 0; 0 0 1; 0 0 1; 0 0 1; 1 1 0; 1 1 0; 1 1 0;];

for icontr = 1 : 1
  
  % -------------------------------------------------------
  % RESTING STATE PHARMA
  % -------------------------------------------------------
    
  load(sprintf('~/pconn_all/proc/all_src_clusterstat_rst_%s_c%d_f%d_v%d.mat',str,icontr,ifoi,v_stat))
  s1 = stats.mask; clear stats
  s1 = round(spatfiltergauss(s1,g1,dd,g2));

  % -------------------------------------------------------
  % RESTING STATE PUPIL
  % -------------------------------------------------------
%   s1 = round(spatfiltergauss(s1,g1,dd,g2));
  load(sprintf('~/pconn_all/proc/all_src_pup_clusterstat_%s_f%d_v%d.mat',str,ifoi,v_pup))
  s2 = stats.mask; clear stats
  s2 = round(spatfiltergauss(s2,g1,dd,g2));

  s3 = zeros(size(s1,1),1);
  s3 = s1&s2;

  
   
  figure; set(gcf,'color','white'); hold on;

  para = [] ;
  para.fixedcolors{1} = [1 0 0];
  para.fixedcolors{2} = [0 0 1];
  para.fixedcolors{3} = [1 1 0]; 
%   para.colorlimits = [0 3];

  for iplot = 1 : 6

    subplot(3,2,iplot)
    
     para.myviewdir = viewdir(iplot,:);
     a = sa_meg_template.cortex10K;
     
     if iplot == 5
      	a=cutsurface(sa_meg_template.cortex10K,[mean(sa_meg_template.cortex10K.vc(:,1)) 0 0],[1 0 0]);
     elseif iplot == 6
      	a=cutsurface(sa_meg_template.cortex10K,[mean(sa_meg_template.cortex10K.vc(:,1)) 0 0],[-1 0 0]);
     end
     
    showsurface(a,para,nan(size(s1)),s1,s2,double(s3))
 
%     colormap(cmap)
    material dull
%     camlight headlight

  end

  print(gcf,'-djpeg100',sprintf('~/pconn_all/plots/pconn_src_conj_pharmpup_%s_c%d_f%d_v%d.jpg',str,icontr,ifoi,v_stat))

end

%% COMPARISON TASK VS REST

    load(sprintf('~/pconn_all/proc/all_src_clusterstat_tsk-rst_%s_c%d_f%d_v%d.mat',str,3,ifoi,v_stat))


    %%
    
viewdir = [0.001 -0 1; 0 -.5 0; 1 0 0; -.5 0 0; -.5 0 0;1 0 0 ];

for ivox = 1 : 3000
  
  d = squeeze(par_all_res(ivox,2,:)-par_all_res(ivox,1,:));
  [r(ivox) p(ivox)] = corr(d,d_behav1);
  
end

s1 = p'<0.05;

% -------------------------------------------------------
% RESTING STATE PHARMA
% -------------------------------------------------------

load(sprintf('~/pconn_all/proc/all_src_clusterstat_rst_%s_c%d_f%d_v%d.mat',str,icontr,ifoi,v_stat))

s2 = stats.mask; clear stats

s = s1&s2; s=double(s);
s(s1~=s2)=eps;

par_interp = spatfiltergauss(s,g1,dd,g2);

figure; set(gcf,'color','white'); hold on;

para = [] ;
para.colorlimits = [-1 1];

for iplot = 1 : 6
  
  subplot(3,2,iplot)
  
  para.myviewdir = viewdir(iplot,:);
  a = sa_meg_template.cortex10K;
  
  if iplot == 5
    a=cutsurface(sa_meg_template.cortex10K,[mean(sa_meg_template.cortex10K.vc(:,1)) 0 0],[1 0 0]);
  elseif iplot == 6
    a=cutsurface(sa_meg_template.cortex10K,[mean(sa_meg_template.cortex10K.vc(:,1)) 0 0],[-1 0 0]);
  end
  
  pconn_showsurface(a,para,par_interp)
  
  colormap(cmap)
  
  camlight headlight
  
end

print(gcf,'-djpeg100',sprintf('~/pconn_all/plots/pconn_src_conj_pharmbehav_%s_c%d_f%d_v%d.jpg',str,icontr,ifoi,v_stat))

