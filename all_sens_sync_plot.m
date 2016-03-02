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

ifoi = 3;

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

res_std  = res_std(:,:,SUBJLIST);
res_mean = res_mean(:,:,SUBJLIST);
res_dfa  = res_dfa(:,:,SUBJLIST);
res_med  = res_med(:,:,SUBJLIST);

cnt_std  = cnt_std(:,:,SUBJLIST);
cnt_mean = cnt_mean(:,:,SUBJLIST);
cnt_dfa  = cnt_dfa(:,:,SUBJLIST);
cnt_med  = cnt_med(:,:,SUBJLIST);



%%


par_std  = [squeeze(nanmean(all_sync_std(:,1,:),1)) squeeze(nanmean(all_sync_std(:,2,:),1)) squeeze(nanmean(all_sync_std(:,3,:),1))];
par_mean = [squeeze(nanmean(all_sync_mean(:,1,:),1)) squeeze(nanmean(all_sync_mean(:,2,:),1)) squeeze(nanmean(all_sync_mean(:,3,:),1))];
par_dfa  = [squeeze(nanmean(all_sync_dfa(:,1,:),1)) squeeze(nanmean(all_sync_dfa(:,2,:),1)) squeeze(nanmean(all_sync_dfa(:,3,:),1))];

% STD OF R
[~,p_std(1)]=ttest(par_std(:,2),par_std(:,1));
[~,p_std(2)]=ttest(par_std(:,3),par_std(:,1));
[~,p_std(3)]=ttest(par_std(:,3),par_std(:,2));

% MEAN OF R
[~,p_mean(1)]=ttest(par_mean(:,2),par_mean(:,1));
[~,p_mean(2)]=ttest(par_mean(:,3),par_mean(:,1));
[~,p_mean(3)]=ttest(par_mean(:,3),par_mean(:,2));

% DFA OF R
[~,p_dfa(1)]=ttest(par_dfa(:,2),par_dfa(:,1));
[~,p_dfa(2)]=ttest(par_dfa(:,3),par_dfa(:,1));
[~,p_dfa(3)]=ttest(par_dfa(:,3),par_dfa(:,2));


% BAR PLOTS

m_std  = nanmean(par_std);
m_mean = nanmean(par_mean);
m_dfa  = nanmean(par_dfa);

sem_std  = nanstd(par_std)/sqrt(length(par_std));
sem_mean = nanstd(par_mean)/sqrt(length(par_std));
sem_dfa  = nanstd(par_dfa)/sqrt(length(par_std));

%

figure; set(gcf,'color','white'); hold on; title(sprintf('STD(SYNC): p1 = %.3f, p2 = %.3f, p3 = %.3f',p_std(1),p_std(2),p_std(3)))

lim = [min(m_std-sem_std)-0.1*min(m_std-sem_std) max(m_std+sem_std)+0.1*max(m_std+sem_std)];

bar([1 1.5 2],m_std,0.4); axis([0.5 2.5 lim]);
line([1 1],[m_std(1)-sem_std(1) m_std(1)+sem_std(1)],'linewidth',5)
line([1.5 1.5],[m_std(2)-sem_std(2) m_std(2)+sem_std(2)],'linewidth',5)
line([2 2],[m_std(3)-sem_std(3) m_std(3)+sem_std(3)],'linewidth',5)
set(gca,'TickDir','out')

print(gcf,'-depsc2',sprintf('~/pconn/proc/plots/pconn_sens_sync_std_f%d_v%d.eps',ifoi,v))

figure; set(gcf,'color','white'); hold on; title(sprintf('MEAN(SYNC): p1 = %.3f, p2 = %.3f, p3 = %.3f',p_mean(1),p_mean(2),p_mean(3)))

lim = [min(m_mean-sem_mean)-0.1*min(m_mean-sem_mean) max(m_mean+sem_mean)+0.1*max(m_mean+sem_mean)];

bar([1 1.5 2],m_mean,0.4); axis([0.5 2.5 lim]);
line([1 1],[m_mean(1)-sem_mean(1) m_mean(1)+sem_mean(1)],'linewidth',5)
line([1.5 1.5],[m_mean(2)-sem_mean(2) m_mean(2)+sem_mean(2)],'linewidth',5)
line([2 2],[m_mean(3)-sem_mean(3) m_mean(3)+sem_mean(3)],'linewidth',5)
set(gca,'TickDir','out')

print(gcf,'-depsc2',sprintf('~/pconn/proc/plots/pconn_sens_sync_mean_f%d_v%d.eps',ifoi,v))

figure; set(gcf,'color','white'); hold on; title(sprintf('DFA(SYNC): p1 = %.3f, p2 = %.3f, p3 = %.3f',p_dfa(1),p_dfa(2),p_dfa(3)))

lim = [min(m_dfa-sem_dfa)-0.1*min(m_dfa-sem_dfa) max(m_dfa+sem_dfa)+0.1*max(m_dfa+sem_dfa)];

bar([1 1.5 2],m_dfa,0.4); axis([0.5 2.5 lim]);
line([1 1],[m_dfa(1)-sem_dfa(1) m_dfa(1)+sem_dfa(1)],'linewidth',5)
line([1.5 1.5],[m_dfa(2)-sem_dfa(2) m_dfa(2)+sem_dfa(2)],'linewidth',5)
line([2 2],[m_dfa(3)-sem_dfa(3) m_dfa(3)+sem_dfa(3)],'linewidth',5)
set(gca,'TickDir','out')

print(gcf,'-depsc2',sprintf('~/pconn/proc/plots/pconn_sens_sync_dfa_f%d_v%d.eps',ifoi,v))



