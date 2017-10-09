
clear

% --------------------------------------------------------
% VERSION 1
% --------------------------------------------------------
v         = 2;
SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
% --------------------------------------------------------

outdir = '/home/tpfeffer/pconn/proc/dfa/';

contr = [2 1; 3 1];


set(0,'Defaultfigurevisible','on')
tp_addpaths

% LOAD DATA
allstr       = {'dfa';'pow'};
all_behav    = {'count';'numb_switches';'avg'};%'mean_dur';'var_dur'};
all_cond     = {'rst';'tsk'};

cmap1  = cbrewer('seq', 'YlOrRd', 150,'pchip'); cmap1 = cmap1(50:1:end,:);
cmap2  = cbrewer('seq', 'YlGnBu', 150,'pchip'); cmap2 = cmap2(end:-1:50,:);
cmap  = [cmap2; ones(50,3); cmap1];

load sa_meg_template;

grid  = sa_meg_template.grid_cortex3000;
g1    = sa_meg_template.grid_cortex3000;
g2    = sa_meg_template.cortex10K.vc;
mri   = sa_meg_template.mri;
vc    = sa_meg_template.vc;
dd    = .75;
g1    = sa_meg_template.grid_cortex3000;
rh   = find(sa_meg_template.grid_cortex3000(:,1)>0.072);
%%
for ifoi = 1:4
  for istr = 1:2
    for ibehav = 1 : 3
      for icontr = 1 : 2
        for icond = 1 : 2
          for imeth = 1 : 2
            fprintf('Computing f%d s%d b%d c%d cond%d ...\n',ifoi,istr,ibehav,icontr,icond)

            load(sprintf('~/pconn_all/proc/all_behav_src_stats_pharm_pooled_str%d_beh%d_c%d_cond%d_meth%d_f%d_v%d.mat',istr,ibehav,icontr,icond,imeth,ifoi,v));
            
            r     = zeros(3000,1);
            r(rh) = s_pharm.corr;
            
            para = [];
            para.colorlimits  = [-0.3 0.3];
            par_interp        = spatfiltergauss(r,g1,dd,g2);
            
            tp_showsource(par_interp,cmap,sa_meg_template,para)
            
            print(gcf,'-djpeg100',sprintf('~/pconn_all/plots/all_behav_src_stats_pharm_map_pooled_mask0_str%d_beh%d_c%d_cond%d_meth%d_f%d_v%d.mat.jpg',istr,ibehav,icontr,icond,imeth,ifoi,v))
            
            if isfield(s_pharm,'mask')
              if any(s_pharm.mask)
                
                r     = zeros(3000,1);
                r(rh) = s_pharm.corr.*s_pharm.mask;
                
                para = [];
                para.colorlimits  = [-0.3 0.3];
                par_interp        = spatfiltergauss(r,g1,dd,g2);
                
                tp_showsource(par_interp,cmap,sa_meg_template,para)
                print(gcf,'-djpeg100',sprintf('~/pconn_all/plots/all_behav_src_stats_pharm_map_pooled_mask1_str%d_beh%d_c%d_cond%d_meth%d_f%d_v%d.mat.jpg',istr,ibehav,icontr,icond,imeth,ifoi,v))
                
              end
            end
            
            clear r s_pharm
            
            load(sprintf(['~/pconn_all/proc/' 'all_behav_src_stats_pbo_pooled_str%d_beh%d_c%d_cond%d_meth%d_f%d_v%d.mat'],istr,ibehav,icontr,icond,imeth,ifoi,v));
            
            r     = zeros(3000,1);
            r(rh) = s_pbo.corr;
            
            para = [];
            para.colorlimits  = [-0.3 0.3];
            par_interp        = spatfiltergauss(r,g1,dd,g2);
            
            tp_showsource(par_interp,cmap,sa_meg_template,para)
            
            print(gcf,'-djpeg100',sprintf('~/pconn_all/plots/all_behav_src_stats_pbo_map_pooled_mask0_str%d_beh%d_c%d_cond%d_meth%d_f%d_v%d.mat.jpg',istr,ibehav,icontr,icond,imeth,ifoi,v))
            
            if isfield(s_pbo,'mask')
              if any(s_pbo.mask)
                
                r     = zeros(3000,1);
                r(rh) = s_pbo.corr.*s_pbo.mask;
                
                para = [];
                para.colorlimits  = [-0.3 0.3];
                par_interp        = spatfiltergauss(r,g1,dd,g2);
                
                tp_showsource(par_interp,cmap,sa_meg_template,para)
                saveas(gcf,sprintf('~/pconn_all/plots/all_behav_src_stats_pbo_map_mask1_str%d_beh%d_c%d_cond%d_meth%d_f%d_v%d.mat.png',istr,ibehav,icontr,icond,imeth,ifoi,v),'png')

                print(gcf,'-djpeg100',sprintf('~/pconn_all/plots/all_behav_src_stats_pbo_map_pooled_mask1_str%d_beh%d_c%d_cond%d_meth%d_f%d_v%d.mat.jpg',istr,ibehav,icontr,icond,imeth,ifoi,v))
                
              end
            end
            
            
            clear s_pbo r
          end
          
          close all
          
        end
      end
    end
  end
end
%% AVERAGE ACROSS BEHAVIOR AND BRAIN DATA (COUNT + PRESS)

cmap1  = cbrewer('seq', 'YlOrRd', 150,'pchip'); cmap1 = cmap1(50:1:end,:);
cmap2  = cbrewer('seq', 'YlGnBu', 150,'pchip'); cmap2 = cmap2(end:-1:50,:);
cmap  = [cmap2; ones(50,3); cmap1];

load /home/tpfeffer/pconn_all/proc/all_src_clusterstat_pooled_cvar_c1_f2_v2.mat

para = [];
para.colorlimits  = [-3 3];
par_interp        = spatfiltergauss(stats.stat.*stats.mask,g1,dd,g2);

tp_showsource(par_interp,cmap,sa_meg_template,para)

print(gcf,'-djpeg100',sprintf('~/pconn_all/plots/all_src_cvar_pooled_c1_f2_v2.mat.jpg'))

load /home/tpfeffer/pconn_all/proc/all_src_clusterstat_pooled_cvar_c2_f2_v2.mat

para = [];
para.colorlimits  = [-3 3];
par_interp        = spatfiltergauss(stats.stat.*stats.mask,g1,dd,g2);

tp_showsource(par_interp,cmap,sa_meg_template,para)

print(gcf,'-djpeg100',sprintf('~/pconn_all/plots/all_src_cvar_pooled_c2_f2_v2.mat.jpg'))


%%
load /home/tpfeffer/pconn_all/proc/all_src_clusterstat_pooled_cvar_c1_f2_v2.mat

para = [];
para.colorlimits  = [-3 3];
par_interp        = spatfiltergauss(stats.stat.*stats.mask,g1,dd,g2);

tp_showsource(par_interp,cmap,sa_meg_template,para)

print(gcf,'-djpeg100',sprintf('~/pconn_all/plots/all_src_cvar_pooled_c1_f2_v2.mat.jpg'))

load /home/tpfeffer/pconn_all/proc/all_src_clusterstat_pooled_cvar_c2_f2_v2.mat

para = [];
para.colorlimits  = [-3 3];
par_interp        = spatfiltergauss(stats.stat.*stats.mask,g1,dd,g2);

tp_showsource(par_interp,cmap,sa_meg_template,para)

print(gcf,'-djpeg100',sprintf('~/pconn_all/plots/all_src_cvar_pooled_c2_f2_v2.mat.jpg'))




