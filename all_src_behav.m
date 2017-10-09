%% COMPUTE WHOLE BRAIN DFA IN SOURCE SPACE
% pconn_src_dfa

clear all

% v = 1: cortex, eloreta, v = 2: cortex, lcmv, v = 3: coarse, eloreta

% --------------------------------------------------------
% VERSION 1
% --------------------------------------------------------
v         = 2;
v_stat    = 2;
v_rawdata = 2;
SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
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
dd = 0.75;


% idx = find(grid(:,1)>mean([-7.1200 7.1268]));

% grid(idx,:)=[]; g1 = grid;

% idx = find(g2(:,1)>mean([-7.1200 7.1268]));
% g2(idx,:) = [];

%% PLOT DFA
   
str = 'dfa'; 
ord   = pconn_randomization;

for ifoi = 1 : 4
  for isubj = SUBJLIST
    for m = 1 : 3

      im = find(ord(isubj,:)==m);
% try
      load(sprintf([outdir 'pconn_cnt_src_dfa_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,v));
     
      dfa_all_cnt(:,m,isubj,ifoi)  = nanmean(par.dfa,2);
      var_all_cnt(:,m,isubj,ifoi)  = nanmean(par.var,2);
     	cvar_all_cnt(:,m,isubj,ifoi) = nanmean(par.cvar,2); 
      amp_all_cnt(:,m,isubj,ifoi)  = nanmean(par.amp,2); clear par
%       
     
      
      load(sprintf(['~/pconn/proc/dfa/pconn_src_dfa_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,v));
     
      dfa_all_res(:,m,isubj,ifoi)  = nanmean(par.dfa,2);
      var_all_res(:,m,isubj,ifoi)  = nanmean(par.var,2);
      cvar_all_res(:,m,isubj,ifoi) = nanmean(par.cvar,2);
      amp_all_res(:,m,isubj,ifoi)  = nanmean(par.amp,2); clear par
      
%       catch me
%         disp(sprintf('err s%d f%d m%d',isubj,ifoi,im))
%       end

    end
  end
end

dfa_all_res  = dfa_all_res(:,:,SUBJLIST,:);
var_all_res  = var_all_res(:,:,SUBJLIST,:);
cvar_all_res = cvar_all_res(:,:,SUBJLIST,:);
amp_all_res  = amp_all_res(:,:,SUBJLIST,:);
% 
dfa_all_cnt  = dfa_all_cnt(:,:,SUBJLIST,:);
var_all_cnt  = var_all_cnt(:,:,SUBJLIST,:);
cvar_all_cnt = cvar_all_cnt(:,:,SUBJLIST,:);
amp_all_cnt  = amp_all_cnt(:,:,SUBJLIST,:);

%% COMPUTE CORRELATION BETWEEN DFA DURING REST AND BEHAVIOR COUNT
% -------------------------------------------------------------------------

clear cnt cnt_all d_behav

addpath ~/pcbi

str = 'dfa';
ifoi = 3;

if strcmp(str,'dfa')
  par_all_res = dfa_all_res(:,:,:,ifoi);
  par_all_cnt = dfa_all_cnt(:,:,:,ifoi);
elseif strcmp(str,'var')
  par_all_res = var_all_res(:,:,:,ifoi).*10^32;
  par_all_cnt = var_all_cnt(:,:,:,ifoi).*10^32;
elseif strcmp(str,'cvar')
  par_all_res = cvar_all_res(:,:,:,ifoi);
  par_all_cnt = cvar_all_cnt(:,:,:,ifoi);
elseif strcmp(str,'amp')
  par_all_res = amp_all_res(:,:,:,ifoi);
  par_all_cnt = amp_all_cnt(:,:,:,ifoi);
end

% cmap = cbrewer('div', 'RdBu', 100,'pchip');% colormap(autumn)
% cmap = cmap(end:-1:1,:);

cnt=pcbi_cnt(1:34);

for isubj = 1:34
  for m = 1 : 3
    
    im = find(ord(isubj,:)==m);
    
    cnt_all(isubj,m) = nanmean(cnt(isubj,im*2-1:im*2));
    
  end
end

cnt_all = cnt_all(SUBJLIST,:);
 
d_behav1 = cnt_all(:,2)-cnt_all(:,1);

% PLOT BEHAVIOR

% figure; set (gcf,'color','w'); hold on;
% 
% plot(1:3,mean(cnt_all),'o','markersize',20)
% 
% line([1 1],[mean(cnt_all(:,1))-std(cnt_all(:,1))/sqrt(length(cnt_all)) mean(cnt_all(:,1))+std(cnt_all(:,1))/sqrt(length(cnt_all))],'color','k','linewidth',2)
% line([2 2],[mean(cnt_all(:,2))-std(cnt_all(:,2))/sqrt(length(cnt_all)) mean(cnt_all(:,2))+std(cnt_all(:,2))/sqrt(length(cnt_all))],'color','k','linewidth',2)
% line([3 3],[mean(cnt_all(:,3))-std(cnt_all(:,3))/sqrt(length(cnt_all)) mean(cnt_all(:,3))+std(cnt_all(:,3))/sqrt(length(cnt_all))],'color','k','linewidth',2)
% 
% axis([0.5 3.5 20 120]);
% 
% set(gca,'Tickdir','out');
% 
% print(gcf,'-depsc2',sprintf('~/pconn_all/plots/pconn_behav_cnt.eps'))

for icontr = 1 : 1
 
    % TASK
% 
    load(sprintf('~/pconn_all/proc/all_src_clusterstat_tsk_%s_c%d_f%d_v%d.mat',str,icontr,ifoi,v_stat))
    
    if any(stats.mask)
    
      stats.mask = logical(stats.mask);

      d=squeeze(nanmean(par_all_cnt(stats.mask,2,:)-par_all_cnt(stats.mask,1,:),1));

      figure; set(gcf,'color','white'); hold on

      s=pconn_regress(d',d_behav1);

      line([-0.04 0.12],[s(2)*-0.04+s(1) s(2)*0.12+s(1)],'linewidth',4,'color',[0.7 0.7 0.7])

      scatter(d,d_behav1,250,'markerfacecolor','k','markeredgecolor','w','linewidth',2);

      xlabel('DFA difference'); ylabel('Change in switch rate'); title(sprintf('Task: f%d',ifoi));
      
      [r,p]=corr(d,d_behav1)
      
      print(gcf,'-depsc2',sprintf('~/pconn_all/plots/pconn_src_behav_tsk_%s_f%d_vstat%d_v%d.eps',str,ifoi,v_stat,v))

    
    end
%     
    % RESTING STATE
    
   	load(sprintf('~/pconn_all/proc/all_src_clusterstat_rst_%s_c%d_f%d_v%d.mat',str,icontr,ifoi,v_stat))

    if any(stats.mask)
    

      stats.mask = logical(stats.mask);

      d=squeeze(nanmean(par_all_res(stats.mask,2,:)-par_all_res(stats.mask,1,:),1));

      figure; set(gcf,'color','white'); hold on

      s=pconn_regress(d',d_behav1);

      line([-0.06 0.14],[s(2)*-0.06+s(1) s(2)*0.14+s(1)],'linewidth',4,'color',[0.7 0.7 0.7])

      scatter(d,d_behav1,250,'markerfacecolor','k','markeredgecolor','w','linewidth',0.25);
      

      xlabel('DFA difference'); ylabel('Change in switch rate'); title(sprintf('Resting state: f%d',ifoi));

      [r,p]=corr(d,d_behav1)
      
      axis([-0.14 0.22 -70 70])
      
      box on; 
      
      print(gcf,'-depsc2',sprintf('~/pconn_all/plots/pconn_src_behav_rst_%s_f%d_vstat%d_v%d.eps',str,ifoi,v_stat,v))

  
    end
end

%% COMPUTE CORRELATIONM BETWEEN DFA-REST AND BEHAVIOR DURING BUTTON
% -------------------------------------------------------------------------

clear cnt d_behav d
str = 'dfa';

v_bttn = 4;
% SUBJLIST =  [3 4 5 6 7 8 9 10 11 12 13 15 16 17 19 20 21 22 23 24];

load ~/pconn_bttn/proc/pconn_bttn_hist_v4.mat

for isubj =  SUBJLIST
  for m = 1 : 3
    for iblock = 1 : 2
      
      try
        load(['~/pconn_bttn/proc/' sprintf('pconn_bttn_dur_s%d_b%d_m%d_v%d.mat',isubj,iblock,m,v_bttn)]) 
      catch 
        cnt(iblock,m,isubj) = nan;
          med_dur(iblock,m,isubj) = nan; 
          continue
      end
        
        if ~isempty(par.dur)
          cnt(iblock,m,isubj) = length(par.dur); 
          med_dur(iblock,m,isubj) = mean(par.dur); 
        else
          cnt(iblock,m,isubj) = nan;
          med_dur(iblock,m,isubj) = nan; 
        end
        clear par
     
      
    end    
  end
end

dur=squeeze(nanmean(med_dur(:,:,SUBJLIST)));
cnt=squeeze(nanmean(cnt(:,:,SUBJLIST),1));

%
% figure; set (gcf,'color','w'); hold on;
% cnt_all = cnt';
% 
% plot(1:3,mean(cnt_all),'o','markersize',20)
% line([1 1],[mean(cnt_all(:,1))-std(cnt_all(:,1))/sqrt(length(cnt_all)) mean(cnt_all(:,1))+std(cnt_all(:,1))/sqrt(length(cnt_all))],'color','k','linewidth',2)
% line([2 2],[mean(cnt_all(:,2))-std(cnt_all(:,2))/sqrt(length(cnt_all)) mean(cnt_all(:,2))+std(cnt_all(:,2))/sqrt(length(cnt_all))],'color','k','linewidth',2)
% line([3 3],[mean(cnt_all(:,3))-std(cnt_all(:,3))/sqrt(length(cnt_all)) mean(cnt_all(:,3))+std(cnt_all(:,3))/sqrt(length(cnt_all))],'color','k','linewidth',2)
% 
% axis([0.5 3.5 20 120]);
% 
% set(gca,'Tickdir','out');
% 
% print(gcf,'-depsc2',sprintf('~/pconn_all/plots/pconn_behav_bttn.eps'))

addpath ~/pcbi

ifoi = 3;

if strcmp(str,'dfa')
  par_all_res = dfa_all_res(:,:,:,ifoi);
%   par_all_cnt = dfa_all_cnt(:,:,:,ifoi);
elseif strcmp(str,'var')
  par_all_res = var_all_res(:,:,:,ifoi).*10^32;
%   par_all_cnt = var_all_cnt(:,:,:,ifoi).*10^32;
elseif strcmp(str,'cvar')
  par_all_res = cvar_all_res(:,:,:,ifoi);
%   par_all_cnt = cvar_all_cnt(:,:,:,ifoi);
elseif strcmp(str,'amp')
  par_all_res = amp_all_res(:,:,:,ifoi);
%   par_all_cnt = amp_all_cnt(:,:,:,ifoi);
end

cmap = cbrewer('div', 'RdBu', 100,'pchip');% colormap(autumn)
cmap = cmap(end:-1:1,:);

d_behav = cnt(2,:)-cnt(1,:);
icontr = 1

load(sprintf('~/pconn_all/proc/all_src_clusterstat_rst_%s_c%d_f%d_v%d.mat',str,icontr,ifoi,v_stat))

if any(stats.mask)
  
  
  stats.mask = logical(stats.mask);
  
  d=squeeze(nanmean(par_all_res(stats.mask,2,:)-par_all_res(stats.mask,1,:),1));
  
  figure; set(gcf,'color','white'); hold on
  
  s=pconn_regress(d',d_behav');
  
  line([-0.06 0.14],[s(2)*-0.06+s(1) s(2)*0.14+s(1)],'linewidth',4,'color',[0.7 0.7 0.7])
  
  scatter(d,d_behav,250,'markerfacecolor','k','markeredgecolor','w','linewidth',2);
  
  xlabel('DFA difference'); ylabel('Change in switch rate'); title(sprintf('Resting state: f%d',ifoi));
  
  [r,p]=corr(d,d_behav')
  
   axis([-0.14 0.22 -70 70])
  
  print(gcf,'-depsc2',sprintf('~/pconn_all/plots/pconn_src_behav_bttn_rst_%s_f%d_vstat%d_v%d.eps',str,ifoi,v_stat,v))
  
  
end

%% MEAN BEHAVIOR

cnt_all = (cnt_all'+cnt)./2;

addpath ~/pcbi

cmap = cbrewer('div', 'RdBu', 100,'pchip');% colormap(autumn)
cmap = cmap(end:-1:1,:);

d_behav = cnt_all(2,:)-cnt_all(1,:);
icontr = 1
%%
load(sprintf('~/pconn_all/proc/all_src_clusterstat_rst_%s_c%d_f%d_v%d.mat',str,icontr,ifoi,v_stat))

if any(stats.mask)
  
  
  stats.mask = logical(stats.mask);
  
  d=squeeze(nanmean(par_all_res(stats.mask,2,:)-par_all_res(stats.mask,1,:),1));
  
  figure; set(gcf,'color','white'); hold on
  
  s=pconn_regress(d',d_behav');
  
  line([-0.06 0.14],[s(2)*-0.06+s(1) s(2)*0.14+s(1)],'linewidth',4,'color',[0.7 0.7 0.7])
  
  scatter(d,d_behav,250,'markerfacecolor','k','markeredgecolor','w','linewidth',2);
  
  xlabel('DFA difference'); ylabel('Change in switch rate'); title(sprintf('Resting state: f%d',ifoi));
  
  [r,p]=corr(d,d_behav')
  
   axis([-0.14 0.22 -70 70])
  
  print(gcf,'-depsc2',sprintf('~/pconn_all/plots/pconn_src_behav_bttn_rst_avg_%s_f%d_vstat%d_v%d.eps',str,ifoi,v_stat,v))
  
  
end





















%% COMPUTE CORRELATION BETWEEN ORDER PARAMETER AND BEHAVIOR
% -------------------------------------------------------------------------

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

d = res_std(:,2)-res_std(:,1);

%%

d_behav = d_behav1;

figure; set(gcf,'color','white'); hold on
  
s=pconn_regress(d',d_behav);

line([-0.06 0.14],[s(2)*-0.06+s(1) s(2)*0.14+s(1)],'linewidth',4,'color',[0.7 0.7 0.7])

scatter(d,d_behav,250,'markerfacecolor','k','markeredgecolor','w','linewidth',2);

xlabel('DFA difference'); ylabel('Change in switch rate'); title(sprintf('Resting state: f%d',ifoi));

[r,p]=corr(d,d_behav)

axis([-0.14 0.22 -70 70])

print(gcf,'-djpeg100',sprintf('~/pconn_all/plots/pconn_sync_behav_rst_%s_f%d_v%d.jpg',str,ifoi,v))

%% REVERSE CORRELATE
viewdir = [0.001 -0 1; 0 -.5 0; 1 0 0; -.5 0 0; -.5 0 0;1 0 0 ];
for ivox = 1 : 3000
  d = squeeze(par_all_res(ivox,2,:)-par_all_res(ivox,1,:));
  
  [r(ivox) ~] = corr(d,d_behav1);
  
end

% alp = fdr(p,0.9)
r(p>0.05)=0;

 par_interp = spatfiltergauss(r',g1,dd,g2);

  figure; set(gcf,'color','white'); hold on;

  para = [] ;
  para.colorlimits = [-0.25 0.25];

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
  
  print(gcf,'-djpeg100',sprintf('~/pconn_all/plots/pconn_src_corr_%spharm_%s_c%d_f%d_v%d.jpg',str,1,ifoi,v_stat))
%% REVERSE CORRELATE RAW VALUES (NOT DIFFERENCES
clear r
viewdir = [0.001 -0 1; 0 -.5 0; 1 0 0; -.5 0 0; -.5 0 0;1 0 0 ];
for ivox = 1 : 3000
  
  d = squeeze(mean(par_all_res(ivox,1:3,:),2));
  
  [r(ivox) p(ivox)] = corr(d,mean(cnt_all(:,1:3),2));
  
end

% alp = fdr(p,0.9)
% r(p>0.05)=0;

 par_interp = spatfiltergauss(r',g1,dd,g2);

  figure; set(gcf,'color','white'); hold on;

  para = [] ;
  para.colorlimits = [-0.25 0.25];

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
  
  print(gcf,'-djpeg100',sprintf('~/pconn_all/plots/pconn_src_corr_%s_f%d_v%d.jpg',str,ifoi,v_stat))

  %% PERFORM CLUSTER BASED PERMUTATIn TEST ON CORRELATION
  
d         = squeeze(par_all_res(ivox,2,:)-par_all_res(ivox,1,:));
d_behav1  = cnt_all(:,2)-cnt_all(:,1);



  
  
%% IMAGE CORREALTION OF COUNT WITH DFA

for ivox = 1 : 3000
  
  d = squeeze(nanmean(par_all_res(ivox,2,:),2)-nanmean(par_all_res(ivox,1,:),2));
  
  [r(ivox) p(ivox)] = corr(d,d_behav1');
  
end

 par_interp = spatfiltergauss(-log10(p'),g1,dd,g2);

  figure; set(gcf,'color','white'); hold on;

  para = [] ;
  para.colorlimits = [-1.5 1.5];

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

%% COMPARE BLOCKS WITH HIGH COUNTS VS BLOCKS WITH FEWER COUNTS


%% BEHAVIORAL DATA
clear grid

behav{1}  = cnt_all;
behav{2}  = cnt';

m_bttn  = mean(behav{1},1);
m_cnt   = mean(behav{2},1);

s_bttn  = std(behav{1},1)./sqrt(length(behav{1}));
s_cnt   = std(behav{2},1)./sqrt(length(behav{2}));

figure; set(gcf,'color','white'); hold on;

plot(m_cnt(1),m_bttn(1),'r.','markersize',40);
plot(m_cnt(2),m_bttn(2),'b.','markersize',40);
plot(m_cnt(3),m_bttn(3),'m.','markersize',40);

legend('P','A','D')

line([m_cnt(1)-s_cnt(1) m_cnt(1)+s_cnt(1)],[m_bttn(1) m_bttn(1)],'color','r','linewidth',2)
line([m_cnt(2)-s_cnt(2) m_cnt(2)+s_cnt(2)],[m_bttn(2) m_bttn(2)],'color','b','linewidth',2)
line([m_cnt(3)-s_cnt(3) m_cnt(3)+s_cnt(3)],[m_bttn(3) m_bttn(3)],'color','m','linewidth',2)

line([m_cnt(1) m_cnt(1)],[m_bttn(1)-s_bttn(1) m_bttn(1)+s_bttn(1)],'color','r','linewidth',2)
line([m_cnt(2) m_cnt(2)],[m_bttn(2)-s_bttn(2) m_bttn(2)+s_bttn(2)],'color','b','linewidth',2)
line([m_cnt(3) m_cnt(3)],[m_bttn(3)-s_bttn(3) m_bttn(3)+s_bttn(3)],'color','m','linewidth',2)

axis([50 90 50 90]); box on; title('Switch rates')

set(gca,'tickdir','out'); xlabel('# of switches (count)'); ylabel('# of switches (button)');
grid  on; axis equal tight

print(gcf,'-djpeg100',sprintf('~/pconn_all/plots/pconn_scatter_behav_%s_f%d_v%d.jpg',str,ifoi,v))

%% 
b(:,1) = nanmean([behav{1}(:,1) behav{1}(:,1)],2);
b(:,2) = nanmean([behav{1}(:,2) behav{1}(:,2)],2);
b(:,3) = nanmean([behav{1}(:,3) behav{1}(:,3)],2);


%% COMPUTE CORRELATION BETWEEN DFA DURING REST AND BEHAVIOR COUNT
% -------------------------------------------------------------------------

para.cond = 'bth';

clear cnt cnt_all d_behav d_behav1

addpath ~/pcbi

str = 'dfa';
ifoi = 3;

if strcmp(str,'dfa')
  par_all_res = dfa_all_res(:,:,:,ifoi);
  par_all_cnt = dfa_all_cnt(:,:,:,ifoi);
elseif strcmp(str,'var')
  par_all_res = var_all_res(:,:,:,ifoi).*10^32;
  par_all_cnt = var_all_cnt(:,:,:,ifoi).*10^32;
elseif strcmp(str,'cvar')
  par_all_res = cvar_all_res(:,:,:,ifoi);
  par_all_cnt = cvar_all_cnt(:,:,:,ifoi);
elseif strcmp(str,'amp')
  par_all_res = amp_all_res(:,:,:,ifoi);
  par_all_cnt = amp_all_cnt(:,:,:,ifoi);
end


para.subj = SUBJLIST;

d_behav1 = pconn_getbehavior(para);

for icontr = 1 : 1
     
    % RESTING STATE
    
   	load(sprintf('~/pconn_all/proc/all_src_clusterstat_rst_%s_c%d_f%d_v%d.mat',str,icontr,ifoi,v_stat))

    if ~any(stats.mask)
      stats.mask = ~stats.mask;
    end
      
    
      stats.mask = logical(stats.mask);

      d=squeeze(nanmean(nanmean(par_all_res(stats.mask,1,:),2),1));

      figure; set(gcf,'color','white'); hold on

      s=pconn_regress(d',d_behav1);

      line([min(d) max(d)],[s(2)*min(d)+s(1) s(2)*max(d)+s(1)],'linewidth',4,'color',[0.7 0.7 0.7])

      scatter(d,d_behav1,250,'markerfacecolor','k','markeredgecolor','w','linewidth',0.25);
      

      xlabel('DFA difference'); ylabel('Change in switch rate'); title(sprintf('Resting state: f%d',ifoi));

      [r,p]=corr(d,d_behav1)
      
      axis([0.45 .95 -20 150])
      
      box on; 
      
      print(gcf,'-depsc2',sprintf('~/pconn_all/plots/pconn_src_behav_rst_nodiff_%s_%s_f%d_vstat%d_v%d.eps',str,para.cond,ifoi,v_stat,v))

  
%     end
end

%% VISUAL CORTEX

% V1
m_r = sa_meg_template.grid_cortex3000-repmat([2.8 -9.6 -0.6],[length(sa_meg_template.grid_cortex3000) 1]);
m_l = sa_meg_template.grid_cortex3000-repmat([-2.8 -9.6 -0.6],[length(sa_meg_template.grid_cortex3000) 1]);

% V5
m_r = sa_meg_template.grid_cortex3000-repmat([4.4 -8.4 1.0],[length(sa_meg_template.grid_cortex3000) 1]);
m_l = sa_meg_template.grid_cortex3000-repmat([-4.4 -8.4 1.0],[length(sa_meg_template.grid_cortex3000) 1]);


[~,i(1)]=min(sum(abs(m_r),2));
[~,i(2)]=min(sum(abs(m_l),2));

[~,p]=ttest(squeeze(nanmean(par_all_res(i,2,:),1)),squeeze(nanmean(par_all_res(i,1,:),1)));

d = squeeze(nanmean(par_all_res(i,2,:),1))-squeeze(nanmean(par_all_res(i,1,:),1));

[r p]=corr(squeeze(mean(nanmean(par_all_res(i,1:3,:),1),2)),mean(cnt_all,2))

z = zeros(3000,1);
z(i) = 100;

par_interp = spatfiltergauss(z,g1,dd,g2);

figure; set(gcf,'color','white'); hold on;

para = [] ;
para.colorlimits = [0 100];

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



