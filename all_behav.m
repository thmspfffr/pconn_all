%% pconn_behav_dfa.m
% all_behav_dfa

% last update: 30-09-2015

clear 
restoredefaultpath

% --------------------------------------------------------
% VERSION 1
% --------------------------------------------------------
v         = 1;
SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
% --------------------------------------------------------

outdir = '/home/tpfeffer/pconn/proc/dfa/';
indir  = '/home/tpfeffer/pconn/proc/preproc/';

tp_addpaths
  
load sa_meg_template;
 
g1    = sa_meg_template.grid_cortex3000;
g2    = sa_meg_template.cortex10K.vc;
mri   = sa_meg_template.mri;
vc    = sa_meg_template.vc;
dd    = .75;


addpath('/home/tpfeffer/Documents/MATLAB/fieldtrip-20160919/')
ft_defaults

%%

figure; hold on
for isubj =SUBJLIST
  for m = 1 : 3
    for ibl = 1 : 2
            
      if ~exist(sprintf([outdir 'pconn_behav_dfa_s%d_m%d_b%d_v%d_processing.txt'],isubj,m,ibl,v))
        system(['touch ' outdir sprintf('pconn_behav_dfa_s%d_m%d_b%d_v%d_processing.txt',isubj,m,ibl,v)]);
      else
        continue
      end

      disp(sprintf('Processing s%d b%d m%d',isubj,ibl,m));
%      	try
        load(['~/pconn_cnt/proc/preproc/' sprintf('pconn_cnt_preproc_data_button_s%d_m%d_b%d_v%d.mat',isubj,m,ibl,1)])
        clear data
        
        cfg = [];
        cfg = cfgs;
        cfg.channel     = {'UPPT002'};
        cfg.continuous  = 'yes';
        data.Fs         = 1200;
        data            = ft_preprocessing(cfg);
  
%         if max(diff(data.trial{1}))<3
%           warning('whaaat??')
%         else
          dat = diff(find(diff(data.trial{1})>0.5))
%         end
%          plot(dat); drawnow
%         par.acf = autocorr(dat,length(dat)-1);
        
%          [a,b]=hist(dat,15);
%          par = gamfit(a);
%          par = par(1);

        save([outdir sprintf('pconn_behav_dfa_s%d_m%d_b%d_v%d.mat',isubj,m,ibl,v)],'dat','-v7.3');

      clear data
    end
  end
end

error('!')

%% LOAD DATA

SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];

ord = pconn_randomization;
clear k b m_behav l_behav n_behav h_behav aa a_corr vv v_behav
addpath ~/pcbi
behav = pcbi_cnt(1:34);
% aa=zeros(2,)

for isubj = SUBJLIST
  for m = 1 : 3
  	
    im = find(ord(isubj,:)==m);

    d = dir(['~/pconn/proc/dfa/' sprintf('pconn_behav_dfa_s%d_m%d_b*_v%d.mat',isubj,im,1)]);
    b(m,isubj) = nanmean(behav(isubj,im*2-1:im*2));

    for iblock = 1:length(d)

      load(['~/pconn/proc/dfa/' d(iblock).name])
      
      if (isubj == 30 && m == 2 && iblock == 2) || (isubj == 26 && m == 1 && iblock == 1) || (isubj == 19 && m == 1 && iblock == 1)
        mm(iblock)   = nan;
        hh(iblock,:) = nan;
        ll(iblock)   = nan;
        nn(iblock)   = nan;

        aa(iblock,:) = nan;
        vv(iblock)   = nan;
        cv(iblock)   = nan;
      else
        mm(iblock)   = mean(dat)./400;
        ll(iblock)   = length(dat);
        nn(iblock)   = median(dat./400);

        aa(iblock,:) = autocorr(dat',9);
        vv(iblock)   = var(dat'./400);
        cv(iblock)   = var(dat'./400)./mean(dat);
      end
      
      clear dat
      
    end
    
    beh_mean_prcpt_dur(isubj,m)     = nanmean(mm);
    beh_switches_bttn(isubj,m)      = nanmean(ll);
    beh_median_prcpt_dur(isubj,m)   = nanmean(nn);
    beh_acf(:,m,isubj)              = nanmean(aa);
    beh_var_prct_dur(isubj,m)       = nanmean(vv);
    beh_cvar_prct_dur(isubj,m)      = nanmean(cv);

  end
end
   
beh_mean_prcpt_dur        = beh_mean_prcpt_dur(SUBJLIST,:);
beh_switches_bttn         = beh_switches_bttn(SUBJLIST,:);
beh_median_prcpt_dur      = beh_median_prcpt_dur(SUBJLIST,:);
beh_acf                   = beh_acf(:,:,SUBJLIST);
beh_var_prct_dur          = beh_var_prct_dur(SUBJLIST,:);
beh_cvar_prct_dur         = beh_cvar_prct_dur(SUBJLIST,:);
beh_switches_cnt          = b(:,SUBJLIST)';


%% PLOT HOW CHANGES IN CNT RELATE TO CHANGES IN BTTN
contrast = [1 2; 2 3];
 
for icontr = 1 : 2

  [r,p]=corr(beh_switches_bttn(:,contrast(icontr,2))-beh_switches_bttn(:,contrast(icontr,1)),beh_switches_cnt(:,contrast(icontr,2))-beh_switches_cnt(:,contrast(icontr,1)))
  slp = pconn_regress([beh_switches_bttn(:,contrast(icontr,2))-beh_switches_bttn(:,contrast(icontr,1))]',beh_switches_cnt(:,contrast(icontr,2))-beh_switches_cnt(:,contrast(icontr,1)));

  figure; set(gcf,'color','white'); hold on
       set(gcf,'Position',[800 20 800 500])

  line([-32 52],[-32*slp(2)+slp(1) 52*slp(2)+slp(1)],'linewidth',5,'color',[0.8 0.8 0.8])
  scatter(beh_switches_bttn(:,contrast(icontr,2))-beh_switches_bttn(:,contrast(icontr,1)),beh_switches_cnt(:,contrast(icontr,2))-beh_switches_cnt(:,contrast(icontr,1)),100,'facecolor','k','markeredgecolor','w')
  axis([-40 60 -45 20]);
  set(gca,'linewidth',2,'ticklength',[0.03 0.03],'tickdir','out');
  set(gca,'fontsize',12,'fontweight','bold');
axis square
  xlabel('\Delta (Count)','fontsize',14,'fontweight','bold')
  ylabel('\Delta (Button)','fontsize',14,'fontweight','bold')

  title('Relation between ATX changes')

  print(gcf,'-depsc2',sprintf('~/pconn_all/plots/behav/all_behav_gamma_v%d.eps',v))
end


%% HOW DOES PHARMA RELATE TO BEHAVIOR???
allstr    = {'dfa';'amp';'cvar'};
allcond      = {'rst';'tsk'};
v = 2;

addpath ~/Documents/MATLAB/cbrewer/cbrewer/
cmap  = cbrewer('div', 'RdBu', 500,'pchip');

cmap2  = cbrewer('seq', 'OrRd', 500,'pchip');
contrast = [2 1; 3 1; 2 3];

for icond = 1 : 1
  cond = allcond{icond};
  for ifoi = 1 : 1
    for icontr = 1 : 1
      for istr = 1 : 1
        
        str       = allstr{istr};
        
        clear dfa_all_res db1 db2 d stats
        
        fprintf('Reading data ...\n');
        ord   = pconn_randomization;
        
        for isubj = SUBJLIST
          isubj
          for m = 1 : 3
            
            im = find(ord(isubj,:)==m);
            load(sprintf(['~/pconn/proc/dfa/pconn_src_dfa_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,v));
            
            eval(sprintf('dfa_all_res(:,m,isubj)  = nanmean(par.%s,2);',str))
            
          end
        end
      end
    end
  end
end

%%
db2=b(:,contrast(icontr,1))-b(:,contrast(icontr,2));
db1=l_behav(:,contrast(icontr,1))-l_behav(:,contrast(icontr,2));

bhv_pbo = b(:,contrast(icontr,1));
dfa_pbo = squeeze(dfa_all_res(:,contrast(icontr,1),SUBJLIST));

d=squeeze(dfa_all_res(:,contrast(icontr,1),SUBJLIST)-dfa_all_res(:,contrast(icontr,2),SUBJLIST));

for i = 1 : 3000
  fprintf('voxel %d ...\n',i)
  [r(i), p(i)]=corr(d(i,:)',db1);
  t(i) = r(i)/sqrt((1-r(i)^2)/(length(db1)-2));
  [r2(i), p2(i)]=corr(d(i,:)',db2);
  t2(i) = r2(i)/sqrt((1-r2(i)^2)/(length(db1)-2));
  
  [r3(i), p3(i)]=corr(dfa_pbo(i,:)',bhv_pbo);
  t3(i) = r3(i)/sqrt((1-r3(i)^2)/(length(db1)-2));
  
end

%       par_interp = spatfiltergauss(r',g1,dd,g2);

para = [] ;
para.colorlimits = [-0.20 0.20];

par_interp = spatfiltergauss(r',g1,dd,g2);

% PLOT RESULTS
tp_showsource(par_interp,cmap,sa_meg_template,para);

print(gcf,'-djpeg100',sprintf('~/pconn_all/plots/all_src_corr_%s_behav_bttn_%s_c%d_f%d_v%d.jpg',str,cond,icontr,ifoi,v))

par_interp = spatfiltergauss(r2',g1,dd,g2);

% PLOT RESULTS
tp_showsource(par_interp,cmap,sa_meg_template,para);
print(gcf,'-djpeg100',sprintf('~/pconn_all/plots/all_src_corr_%s_behav_cnt_%s_c%d_f%d_v%d.jpg',str,cond,icontr,ifoi,v))

par_interp = spatfiltergauss(r3',g1,dd,g2);

% PLOT RESULTS
tp_showsource(par_interp,cmap,sa_meg_template,para);
print(gcf,'-djpeg100',sprintf('~/pconn_all/plots/all_src_corr_%s_behav_pbo_cnt_%s_c%d_f%d_v%d.jpg',str,cond,icontr,ifoi,v))

para = [] ;
para.colorlimits = [-1.96 1.96];

par_interp = spatfiltergauss(t',g1,dd,g2);

% PLOT RESULTS
tp_showsource(par_interp,cmap,sa_meg_template,para);

print(gcf,'-djpeg100',sprintf('~/pconn_all/plots/all_src_corr_%s_behav_bttn_pval_%s_c%d_f%d_v%d.jpg',str,cond,icontr,ifoi,v))

par_interp = spatfiltergauss(t2',g1,dd,g2);

% PLOT RESULTS
tp_showsource(par_interp,cmap,sa_meg_template,para);
print(gcf,'-djpeg100',sprintf('~/pconn_all/plots/all_src_corr_%s_behav_cnt_pval_%s_c%d_f%d_v%d.jpg',str,cond,icontr,ifoi,v))

par_interp = spatfiltergauss(t3',g1,dd,g2);

% PLOT RESULTS
tp_showsource(par_interp,cmap,sa_meg_template,para);
print(gcf,'-djpeg100',sprintf('~/pconn_all/plots/all_src_corr_%s_behav_cnt_pbo_%s_c%d_f%d_v%d.jpg',str,cond,icontr,ifoi,v))

clear r r2 p p2 d db1 db2

   
%%

icontr = 2;
contrast = [2 1; 3 1; 2 3];
      db1=m_behav(:,contrast(icontr,1))-m_behav(:,contrast(icontr,2));
      db2=b(:,contrast(icontr,1))-b(:,contrast(icontr,2));

str    = 'amp';
cond   = 'rst';
ifoi   = 4;


fprintf('Reading data ...\n');
ord   = pconn_randomization;

     if exist(sprintf('~/pconn_all/proc/all_src_clusterstat_%s_%s_c%d_f%d_v%d.mat',cond,str,icontr,ifoi,2))
        load(sprintf('~/pconn_all/proc/all_src_clusterstat_%s_%s_c%d_f%d_v%d.mat',cond,str,icontr,ifoi,2))
     else
        error('No clusters found!')
%         continue
      end

if strcmp(cond,'rst')
  for isubj = SUBJLIST
    isubj
    for m = 1 : 3

      im = find(ord(isubj,:)==m);
      load(sprintf(['~/pconn/proc/dfa/pconn_src_dfa_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,2));

      eval(sprintf('dfa_all_res(:,m,isubj)  = nanmean(par.%s,2);',str))

    end
  end
elseif strcmp(cond,'tsk')
  for isubj = SUBJLIST
    isubj
    for m = 1 : 3

      im = find(ord(isubj,:)==m);
      load(sprintf(['~/pconn_cnt/proc/dfa/pconn_cnt_src_dfa_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,2));

      eval(sprintf('dfa_all_res(:,m,isubj)  = nanmean(par.%s,2);',str))

    end
  end
end


d = squeeze(mean(dfa_all_res(logical(stats.mask),2,SUBJLIST),1)-mean(dfa_all_res(logical(stats.mask),1,SUBJLIST),1));
figure; set(gcf,'color','white'); hold on

set(gcf,'paperpositionmode','manual')

set(gcf,'Position',[800 20 1200 500])
set(gcf,'Papersize',[12 18])

subplot(1,2,1); 
scatter(db1,d,80,'facecolor','k','markeredgecolor','w');

set(gca,'linewidth',2,'ticklength',[0.03 0.03],'tickdir','out');
set(gca,'fontsize',12,'fontweight','bold');
% set(gca,'xtick',[2 3 4 5 6 7 8 9],'xticklabel',[4 8 16 32 64 128 256 512])
  
xlabel('\Delta (Behavior)','fontsize',14,'fontweight','bold')
ylabel('\Delta (DFA)','fontsize',14,'fontweight','bold')

[r,p]=corr(d,db1);

axis square; title(sprintf('BUTTON: p = %.2f',p)); lsline

subplot(1,2,2)

scatter(db2,d,80,'facecolor','k','markeredgecolor','w');

[r,p]=corr(d,db2);

set(gca,'linewidth',2,'ticklength',[0.03 0.03],'tickdir','out');
set(gca,'fontsize',12,'fontweight','bold');
% set(gca,'xtick',[2 3 4 5 6 7 8 9],'xticklabel',[4 8 16 32 64 128 256 512])
  
xlabel('\Delta (Behavior)','fontsize',14,'fontweight','bold')
ylabel('\Delta (DFA)','fontsize',14,'fontweight','bold')
  
axis square; title(sprintf('COUNT: p = %.2f',p)); lsline
% 
% 
% subplot(1,3,3)
% 
% scatter(db4,d,80,'facecolor','k','markeredgecolor','w');
% 
% set(gca,'linewidth',2,'ticklength',[0.03 0.03],'tickdir','out');
% set(gca,'fontsize',12,'fontweight','bold');
% % set(gca,'xtick',[2 3 4 5 6 7 8 9],'xticklabel',[4 8 16 32 64 128 256 512])
%   
% xlabel('\Delta (Behavior)','fontsize',14,'fontweight','bold')
% ylabel('\Delta (DFA)','fontsize',14,'fontweight','bold')
%   
% axis square; title('COUNT'); lsline

print(gcf,'-depsc2',sprintf('~/pconn_all/plots/behav/all_behav_change_f%d_%s_c%d_v%d.eps',ifoi,str,1,v))




%% BEHAVIORAL PLOTS
% VERSION FROM 22-10-16 in PDF

addpath ~/Documents/MATLAB/Colormaps/'Colormaps (5)'/Colormaps/
contrasts = [2 1; 3 1; 2 3];
for icontr = 1 : 2
  
col = [linspace(0,0.5,length(SUBJLIST)); linspace(0,0.5,length(SUBJLIST));linspace(0,0.5,length(SUBJLIST))]';

figure; set(gcf,'color','white'); hold on

for ip = 1 : 2
  subplot(1,3,ip); hold on

  if ip == 1 
    par = l_behav(:,[contrasts(icontr,2) contrasts(icontr,1)]);
  else
    par = b(:,[contrasts(icontr,2) contrasts(icontr,1)]);
  end
% plot(1*ones(size(m_behav,1),1),m_behav(:,3),'.','color',col{3},'markersize',20)

for i = 1 : length(SUBJLIST)
  plot(2,log2(par(i,1)),'.','color',col(i,:),'markersize',17)
  plot(3,log2(par(i,2)),'.','color',col(i,:),'markersize',17)
  plot(2:3,log2(par(i,1:2)),'color',col(i,:),'linewidth',1.5)
end

plot(1.75,log2(mean(par(:,1))),'.','color',[0.8 0.8 0.8],'markersize',30)
plot(3.25,log2(mean(par(:,2))),'.','color',[0.8 0.8 0.8],'markersize',30)

s1 = std(log2(par(:,1)))./sqrt(size(par,1));
s2 = std(log2(par(:,2)))./sqrt(size(par,1));

line([1.75 1.75],[log2(mean(par(:,1)))-s1 log2(mean(par(:,1)))+s1],'color',[0.6 0.6 0.6],'linewidth',3)
line([3.25 3.25],[log2(mean(par(:,2)))-s2 log2(mean(par(:,2)))+s2],'color',[0.6 0.6 0.6],'linewidth',3)

set(gca,'ytick',[2 3 4 5 6 7],'yticklabel',[4 8 16 32 64 128])
% plot([2*ones(1,size(m_behav,1)); 1*ones(1,size(m_behav,1))],m_behav(:,[1 3])','k')
if icontr == 1
  set(gca,'xtick',[2 3],'xticklabel',['PBO';'ATX'])
else
  set(gca,'xtick',[2 3],'xticklabel',['PBO';'DPZ'])
end

ylabel('# switches','fontsize',14)

[t,p] = ttest(par(:,1),par(:,2));

if ip == 1
  title(sprintf('Button: p = %.4f',p))
else
  title(sprintf('Count: p = %.4f',p))
end

tp_editplots
axis([1.5 3.5 3 7.5]);
end

% --------------------------
% AVERAGE
% --------------------------

subplot(1,3,3); hold on

par = (b(:,[contrasts(icontr,2) contrasts(icontr,1)])+l_behav(:,[contrasts(icontr,2) contrasts(icontr,1)]))./2;

for i = 1 : length(SUBJLIST)
  plot(2,log2(par(i,1)),'.','color',col(i,:),'markersize',17)
  plot(3,log2(par(i,2)),'.','color',col(i,:),'markersize',17)
  plot(2:3,log2(par(i,1:2)),'color',col(i,:),'linewidth',1.5)
end

plot(1.75,log2(mean(par(:,1))),'.','color',[0.8 0.8 0.8],'markersize',30)
plot(3.25,log2(mean(par(:,2))),'.','color',[0.8 0.8 0.8],'markersize',30)

s1 = std(log2(par(:,1)))./sqrt(size(par,1));
s2 = std(log2(par(:,2)))./sqrt(size(par,1));

line([1.75 1.75],[log2(mean(par(:,1)))-s1 log2(mean(par(:,1)))+s1],'color',[0.6 0.6 0.6],'linewidth',3)
line([3.25 3.25],[log2(mean(par(:,2)))-s2 log2(mean(par(:,2)))+s2],'color',[0.6 0.6 0.6],'linewidth',3)

set(gca,'ytick',[2 3 4 5 6 7],'yticklabel',[4 8 16 32 64 128])
% plot([2*ones(1,size(m_behav,1)); 1*ones(1,size(m_behav,1))],m_behav(:,[1 3])','k')
if icontr == 1
  set(gca,'xtick',[2 3],'xticklabel',['PBO';'ATX'])
else
  set(gca,'xtick',[2 3],'xticklabel',['PBO';'DPZ'])
end

ylabel('# switches','fontsize',14)

[t,p] = ttest(par(:,1),par(:,2));
title(sprintf('CNT+BTTN: p = %.4f',p))


tp_editplots
axis([1.5 3.5 3 7.5]);


print(gcf,'-depsc2',sprintf('~/pconn_all/plots/all_behav_bttn_cnt_c%d_v%d.eps',icontr,v))

end


