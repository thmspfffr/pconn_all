
%% COMPUTE EXCITATION/INHIBITION BALANCE ACROSS SENSOR SPACE
% pconn_sens_sync

% Method by R Hardstone, discussed during meeting on 4th of March
% Takes band-pass filtered amplitude and time-resolved DFA as input
% and estimates E/I balance from that.


% --------------------------------------------------------
% VERSION 1
% --------------------------------------------------------
v         = 1;
SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24];
% --------------------------------------------------------


addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/
addpath /home/tpfeffer/Documents/MATLAB/fieldtrip-20130925/

indir   = '/home/tpfeffer/pconn/proc/src/';
outdir   = '/home/tpfeffer/pconn/proc/dfa/';
plotdir = '~/pconn_all/plots/';
addpath /home/tpfeffer/pconn/matlab/

run ~/Documents/MATLAB/toolboxes/NBT-NBTv0.5.3-alpha/installNBT.m
%%

clear res_std res_mean res_dfa cnt_std cnt_mean cnt_dfa


ord       = pconn_randomization;

ifoi = 4;

for m = 1 : 3
  for isubj = SUBJLIST
    	
    im = find(ord(isubj,:)==m);

    for iblock = 1 : 2

      load(sprintf('/home/tpfeffer/pconn/proc/dfa/pconn_sens_sync_s%d_b%d_m%d_f%d_v%d.mat',isubj,iblock,im,ifoi,v))
      
      res_std(iblock,m,isubj)  = std(r); 
    	res_mean(iblock,m,isubj) = mean(r);
      res_dfa(iblock,m,isubj)  = dfa.MarkerValues; 
      res_med(iblock,m,isubj) = median(abs(r-median(r))); clear r dfa
      
      load(sprintf('~/pconn_cnt/proc/dfa/pconn_cnt_sens_sync_s%d_b%d_m%d_f%d_v%d.mat',isubj,iblock,im,ifoi,v))
      
      cnt_std(iblock,m,isubj)  = std(r); 
    	cnt_mean(iblock,m,isubj) = mean(r);
      cnt_dfa(iblock,m,isubj)  = dfa.MarkerValues; 
      cnt_med(iblock,m,isubj) = median(abs(r-median(r))); clear r dfa

    end
  end
end

res_std  = squeeze(nanmean(res_std(:,:,SUBJLIST),1))';
res_mean = squeeze(nanmean(res_mean(:,:,SUBJLIST),1))';
res_dfa  = squeeze(nanmean(res_dfa(:,:,SUBJLIST),1))';
res_med  = squeeze(nanmean(res_med(:,:,SUBJLIST),1))';

cnt_std  = squeeze(nanmean(cnt_std(:,:,SUBJLIST),1))';
cnt_mean = squeeze(nanmean(cnt_mean(:,:,SUBJLIST),1))';
cnt_dfa  = squeeze(nanmean(cnt_dfa(:,:,SUBJLIST),1))';
cnt_med  = squeeze(nanmean(cnt_med(:,:,SUBJLIST),1))';



%%

str = 'std';

switch str
  case 'std'
    par_res  = res_std;
    par_cnt  = cnt_std;
  case 'mean'
    par_res  = res_mean;
    par_cnt  = cnt_mean;
  case 'dfa'
    par_res  = res_dfa;
    par_cnt  = cnt_dfa;
  case 'med'
    
end
    

% STD OF R
[~,p_std(1)]=ttest(par_res(:,2),par_res(:,1));
[~,p_std(2)]=ttest(par_res(:,3),par_res(:,1));
[~,p_std(3)]=ttest(par_res(:,3),par_res(:,2));

% BAR PLOTS

m_std    = nanmean(par_res);
sem_std  = nanstd(par_res)/sqrt(length(par_res));

figure; set(gcf,'color','white'); hold on; title(sprintf('STD(SYNC) RST: p1 = %.3f, p2 = %.3f, p3 = %.3f',p_std(1),p_std(2),p_std(3)))

lim = [min(m_std-sem_std)-0.1*min(m_std-sem_std) max(m_std+sem_std)+0.1*max(m_std+sem_std)];

bar([1 1.5 2],m_std,0.4); axis([0.5 2.5 lim]);
line([1 1],[m_std(1)-sem_std(1) m_std(1)+sem_std(1)],'linewidth',5)
line([1.5 1.5],[m_std(2)-sem_std(2) m_std(2)+sem_std(2)],'linewidth',5)
line([2 2],[m_std(3)-sem_std(3) m_std(3)+sem_std(3)],'linewidth',5)
set(gca,'TickDir','out')

print(gcf,'-depsc2',sprintf('~/pconn/proc/plots/pconn_sens_sync_rst_%s_f%d_v%d.eps',str,ifoi,v))

% STD OF R
[~,p_std(1)]=ttest(par_cnt(:,2),par_cnt(:,1));
[~,p_std(2)]=ttest(par_cnt(:,3),par_cnt(:,1));
[~,p_std(3)]=ttest(par_cnt(:,3),par_cnt(:,2));

% BAR PLOTS

m_std    = nanmean(par_cnt);
sem_std  = nanstd(par_cnt)/sqrt(length(par_cnt));

figure; set(gcf,'color','white'); hold on; title(sprintf('STD(SYNC) TSK: p1 = %.3f, p2 = %.3f, p3 = %.3f',p_std(1),p_std(2),p_std(3)))

lim = [min(m_std-sem_std)-0.1*min(m_std-sem_std) max(m_std+sem_std)+0.1*max(m_std+sem_std)];

bar([1 1.5 2],m_std,0.4); axis([0.5 2.5 lim]);
line([1 1],[m_std(1)-sem_std(1) m_std(1)+sem_std(1)],'linewidth',5)
line([1.5 1.5],[m_std(2)-sem_std(2) m_std(2)+sem_std(2)],'linewidth',5)
line([2 2],[m_std(3)-sem_std(3) m_std(3)+sem_std(3)],'linewidth',5)
set(gca,'TickDir','out')

print(gcf,'-depsc2',sprintf('~/pconn/proc/plots/pconn_sens_sync_tsk_%s_f%d_v%d.eps',str,ifoi,v))

[~,p]=ttest(nanmean(par_cnt,2),nanmean(par_res,2))
%%