%% all_behav_src_stats.m
% all_behav_src_stats_pooledhemispheres

% POOL THE HEMISPHERES!

% last change: 24th of november, 2016

% computes correlation between behavior (switch rate, mean duration etc.)
% and various brain signals (DFA, power, variance). corrects for multiple
% comparisons using a cluster-based permutation approach

clear 
restoredefaultpath
rng('shuffle','twister')

% --------------------------------------------------------
% VERSION 1
% --------------------------------------------------------
v         = 2;
v_stat = 2;
SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];

para.minneigh     = 3;
para.clusteralpha = 0.025;
para.alpha        = 0.025;
para.nperm        = 10000;
% --------------------------------------------------------
% VERSION 22
% --------------------------------------------------------
% v_stat         = 22;
% v = 2;
% SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
% 
% para.minneigh     = 2;
% para.clusteralpha = 0.05;
% para.alpha        = 0.025;
% para.nperm        = 10000;
% --------------------------------------------------------


outdir = '/home/tpfeffer/pconn/proc/dfa/';
indir  = '/home/tpfeffer/pconn/proc/preproc/';

addpath /home/tpfeffer/pconn/matlab/


addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/
addpath ~/pconn/matlab/

load sa_meg_template;
ord   = pconn_randomization;

contr = [2 1; 3 1];

%% LOAD DATA
allstr       = {'dfa';'pow'};
all_behav    = {'count';'numb_switches';'avg'};
all_cond     = {'rst';'tsk'};
all_meth     = {'correlation_pearson';'correlation_spearman'};
g1    = sa_meg_template.grid_cortex3000;

for ifoi = 1 : 4
  for istr = 1 : 2
    str         = allstr{istr};
    for icontr = 1 : 2
      for icond = 1 : 2
        cond        = all_cond{icond};
        for imeth = 1 : 2
          para.method = all_meth{imeth};
          for ibehav = 1 : 3

            if ~exist(sprintf(['~/pconn_all/proc/' 'all_behav_src_stats_pooled_str%d_beh%d_c%d_cond%d_meth%d_f%d_v%d_processing.txt'],istr,ibehav,icontr,icond,imeth,ifoi,v_stat))
              system(['touch ' '~/pconn_all/proc/' sprintf('all_behav_src_stats_pooled_str%d_beh%d_c%d_cond%d_meth%d_f%d_v%d_processing.txt',istr,ibehav,icontr,icond,imeth,ifoi,v_stat)]);
            else
              continue
            end        
%             
            fprintf('Computing f%d s%d b%d c%d cond%d ...\n',ifoi,istr,ibehav,icontr,icond)
            
            str_behav   = all_behav{ibehav};

            % READ IN BEHAVIORAL DATA
            if ~strcmp(str_behav,'avg')
              para.str_behav = str_behav;
              par_behav = pconn_read_behavioral_data(SUBJLIST,para);
            else
              for ib = 1 : 2
              	para.str_behav   = all_behav{ib};
                par_behav_tmp{ib} = pconn_read_behavioral_data(SUBJLIST,para);
              end
              
              str_behav = 'avg';
              par_behav = cat(3,par_behav_tmp{1},par_behav_tmp{2});
                
            end
            
            tmp=nanmean(par_behav,3)';
            [nanidx,~]=find(isnan(tmp));
            %% READ NEURAL DATA
            clear all_par
            
            for isubj = SUBJLIST
              fprintf('Reading subj %d ...\n',isubj);
              for m = 1 : 3
                
                im = find(ord(isubj,:)==m);
                if strcmp(cond,'rst')
                  if ~strcmp(str,'pow')
                    load(sprintf('~/pconn/proc/dfa/pconn_src_dfa_s%d_m%d_f%d_v%d.mat',isubj,im,ifoi,v));
                    eval(sprintf('all_par(:,:,m,isubj)  = par.%s;',str))
                  else
                    load(sprintf('~/pconn/proc/src/pconn_src_pow_s%d_m%d_f%d_v%d.mat',isubj,im,ifoi,1));
                    eval(sprintf('all_par(:,:,m,isubj)  = par.%s;',str))
                  end
                elseif strcmp(cond,'tsk')
                  if ~strcmp(str,'pow')
                    load(sprintf('~/pconn_cnt/proc/dfa/pconn_cnt_src_dfa_s%d_m%d_f%d_v%d.mat',isubj,im,ifoi,v));
                    eval(sprintf('all_par(:,:,m,isubj)  = par.%s;',str))
                  else
                    load(sprintf('~/pconn_cnt/proc/src/pconn_cnt_src_pow_s%d_m%d_f%d_v%d.mat',isubj,im,ifoi,1));
                    eval(sprintf('all_par(:,:,m,isubj)  = par.%s;',str))
                  end
                end
              end
            end
            
            all_par   = all_par(:,:,:,SUBJLIST);
            
            %% COMPUTE STATSISTICS
            
            load sa_meg_template.mat
            
            rh      = find(sa_meg_template.grid_cortex3000(:,1)>0.072);
            lh      = find(sa_meg_template.grid_cortex3000(:,1)<0.072);
            all_par = (all_par(rh,:,:,:)+all_par(lh,:,:,:))./2;

            all_brain     = squeeze(nanmean(all_par(:,:,[contr(icontr,1) contr(icontr,2)],:),2));
            all_brain     = permute(all_brain,[1 3 2]);
            par_behav_all = squeeze(nanmean(par_behav(contr(icontr,1),:,:),3))-squeeze(nanmean(par_behav(contr(icontr,2),:,:),3));
            
            if ~isempty(nanidx)
              all_brain(:,nanidx,:) = [];
              par_behav_all(nanidx) = [];
            end
              
            n                 = get_neighbours(g1(rh,:));
            para.neigh        = n;
            para.extvar       = par_behav_all;
            
            s_pharm           = tp_clusterperm(all_brain,para);
            
            s_pharm.behavior  = str_behav;
            s_pharm.contrast  = contr(icontr,:);
            s_pharm.metric    = str;
            
            save(sprintf(['~/pconn_all/proc/' 'all_behav_src_stats_pharm_pooled_str%d_beh%d_c%d_cond%d_meth%d_f%d_v%d.mat'],istr,ibehav,icontr,icond,imeth,ifoi,v_stat),'s_pharm','-v7.3');
            
            plac_brain = squeeze(nanmean(all_par(:,:,1,:),2));
            plac_behav = squeeze(nanmean(par_behav(1,:,:),3));
            
            para.extvar       = plac_behav;
            s_pbo             = tp_clusterperm(plac_brain,para);
            
            s_pbo.behavior    = str_behav;
            s_pbo.metric      = str;
            
            save(sprintf(['~/pconn_all/proc/' 'all_behav_src_stats_pbo_pooled_str%d_beh%d_c%d_cond%d_meth%d_f%d_v%d.mat'],istr,ibehav,icontr,icond,imeth,ifoi,v_stat),'s_pbo','-v7.3');
            clear s_pbo s_pharm
            
            
            
          end    
        end
      end
    end
  end
end


         