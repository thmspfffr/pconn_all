%% all_behav_src_stats.m
% all_behav_src_stats

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
SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];

para.minneigh     = 3;
para.clusteralpha = 0.025;
para.alpha        = 0.025;
para.nperm        = 10000;
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

%% 
ord       = pconn_randomization;
para.cond = 'rst';
para.str  = 'pow';
para.ifoi = 1;
all_par   = pconn_all_read_neural_data(SUBJLIST,ord,para);
icontr    = 2;
ifoi      = 1;

str_behav = 'count';

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

%%
ifoi = 2;
icontr = 1;
all_par   = pconn_all_read_neural_data(SUBJLIST,ord,para);

load(sprintf('/home/tpfeffer/pconn_all/proc/all_src_clusterstat_power_rst_c%d_f%d_v%d.mat',icontr,ifoi,1))
all_par = squeeze(nanmean(nanmean(all_par(logical(stats.mask),:,:,:),2),1));

d_behav = nanmean(par_behav(2,:,:),3)-nanmean(par_behav(1,:,:),3);
d_brain = all_par(2,:)-all_par(1,:);

            
          

plot(d_brain,d_behav,'.')
axis([min(d_brain) -min(d_brain) -max(d_behav) max(d_behav)]); axis square
lsline
[r,p]=corrcoef(d_brain,d_behav)







