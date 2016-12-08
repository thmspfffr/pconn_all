%% COMPUTE POWER SPECTRUM FOR ALL SENSORS

% all_sens_powspec
clear

% --------------------------------------------------------
v         = 1;
fsample   = 400;
SUBJLIST 	= [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
f         = 2:0.2:200;
T         = 0.2;
planar    = 0;
% --------------------------------------------------------

restoredefaultpath

tp_addpaths();

addpath /home/tpfeffer/pconn/matlab
outdir   = '/home/tpfeffer/pconn/proc/powspec/';

%%
ord       = pconn_randomization;


for m = 1:3

  for isubj = SUBJLIST
    
  disp(sprintf('Processing subject %d ...',isubj));
  
   im = find(ord(isubj,:)==m);
%   
%     if ~exist(sprintf(['~/pconn_all/proc/' 'all_sens_powspec_s%d_m%d_v%d_processing.txt'],isubj,m,v))
%       system(['touch ' '~/pconn_all/proc/' sprintf('all_sens_powspec_s%d_m%d_v%d_processing.txt',isubj,m,v)]);
%     else
%       continue
%     end
   
   d = dir([outdir sprintf('pconn_powspec_dat_s%d_b*_m%d_v%d.mat',isubj,im,v)]);
   
    for iblock = 1 : 2
      
      if isubj >= 32
        load(sprintf('~/pconn/proc/preproc/pconn_chanidx_s%d_m%d.mat',31,3))
      elseif isubj < 4
        load(sprintf('~/pconn/proc/preproc/pconn_chanidx_s%d_m%d.mat',4,1))
      elseif isubj == 17
        load(sprintf('~/pconn/proc/preproc/pconn_chanidx_s%d_m%d.mat',16,1))
      else
        load(sprintf('~/pconn/proc/preproc/pconn_chanidx_s%d_m%d.mat',isubj,im))
      end
          
      disp(sprintf('Loading MEG data ...'));
      
      try
      load([outdir sprintf('pconn_powspec_dat_s%d_b%d_m%d_v%d.mat',isubj,iblock,im,v)]);
      
      if size(dat.label,1)~=274
        for ichan = 1 : 991
          pow_res(:,ichan,iblock) = pconn_sens_interp274(idx,squeeze(nanmean(dat.powspctrm(:,:,ichan),1)));
        end
      else
        pow_res(:,:,iblock) = squeeze(nanmean(dat.powspctrm,1));
      end
       catch me
        warning(sprintf('Data set (rest) not found: s%d b%d m%d',isubj,iblock,m))
        pow_res(:,:,iblock) = nan(274,991);
      end
      clear dat
      
      try
        load(['~/pconn_cnt/proc/powspec/' sprintf('pconn_cnt_powspec_dat_s%d_b%d_m%d_v%d.mat',isubj,iblock,im,v)]);

        if size(dat.label,1)~=274
          for ichan = 1 : 991
            pow_cnt(:,ichan,iblock) = pconn_sens_interp274(idx,squeeze(nanmean(dat.powspctrm(:,:,ichan),1)));
          end
        else
          pow_cnt(:,:,iblock) = squeeze(nanmean(dat.powspctrm,1));
        end
      catch me
        warning(sprintf('Data set (cnt) not found: s%d b%d m%d',isubj,iblock,m))
        pow_cnt(:,:,iblock) = nan(274,991);
      end
      
        figure; set(gcf,'color','w'); 
       title(sprintf('SUBJECT %d',isubj))
       subplot(2,2,1)
       
       plot(log10(f),squeeze(nanmean(pow_res(:,:,iblock))))
       
       set(gca,'xtick',[0 1 2],'xticklabel',[1 10 100],'tickdir','out')
       xlabel('Frequency [Hz]'); ylabel('Power');
       
       subplot(2,2,2)
       
       plot(log10(f),log10(squeeze(nanmean(pow_res(:,:,iblock)))))
       
       set(gca,'xtick',[0 1 2],'xticklabel',[1 10 100],'tickdir','out')
       xlabel('Frequency [Hz]'); ylabel('Log-Power');
       
       subplot(2,2,3)
      
       plot(log10(f),squeeze(nanmean(pow_cnt(:,:,iblock))))
       
       set(gca,'xtick',[0 1 2],'xticklabel',[1 10 100],'tickdir','out')
       xlabel('Frequency [Hz]'); ylabel('Power');

       
       subplot(2,2,4)
       
       plot(log10(f),log10(squeeze(nanmean(pow_cnt(:,:,iblock)))))
       
       set(gca,'xtick',[0 1 2],'xticklabel',[1 10 100],'tickdir','out')
       xlabel('Frequency [Hz]'); ylabel('Log-Power');
       
       print(gcf,'-djpeg100',sprintf('/home/tpfeffer/pconn_all/plots/all_sens_powspec_meanchan_s%d_m%d_b%d_v%d.jpg',isubj,m,iblock,v))
       
       % PLOT WITH ALL SENSORS VISIBLE
       
       figure; set(gcf,'color','w'); 
       
       subplot(2,2,1)
       
       plot(log10(f),pow_res(:,:,iblock))
       
       set(gca,'xtick',[0 1 2],'xticklabel',[1 10 100],'tickdir','out')
       xlabel('Frequency [Hz]'); ylabel('Power');
       
       subplot(2,2,2)
       
       plot(log10(f),log10(pow_res(:,:,iblock)))
       
       set(gca,'xtick',[0 1 2],'xticklabel',[1 10 100],'tickdir','out')
       xlabel('Frequency [Hz]'); ylabel('Log-Power');
       
       subplot(2,2,3)
      
       plot(log10(f),pow_cnt(:,:,iblock))
       
       set(gca,'xtick',[0 1 2],'xticklabel',[1 10 100],'tickdir','out')
       xlabel('Frequency [Hz]'); ylabel('Power');

       
       subplot(2,2,4)
       
       plot(log10(f),log10(pow_cnt(:,:,iblock)))
       
       set(gca,'xtick',[0 1 2],'xticklabel',[1 10 100],'tickdir','out')
       xlabel('Frequency [Hz]'); ylabel('Log-Power');
       
       print(gcf,'-djpeg100',sprintf('/home/tpfeffer/pconn_all/plots/all_sens_powspec_allchan_s%d_m%d_b%d_v%d.jpg',isubj,m,iblock,v))

       close all
      
    end
    
    pow_res = nanmean(pow_res,3);
    pow_cnt = nanmean(pow_cnt,3);
    
    save(sprintf('~/pconn_all/proc/pconn_all_pow_s%d_m%d_v%d.mat',isubj,m,v),'pow_res','pow_cnt','-v7.3');
    
    clear pow_res pow_cnt
    
  end
end



error('Stop here!')

%% PLOT EVERYTHING

for isubj = SUBJLIST
  for m = 1 : 3
    
    load(sprintf('~/pconn_all/proc/pconn_all_pow_s%d_m%d_v%d.mat',isubj,m,v));
    
    r_pow(:,isubj,m) = squeeze(nanmean(pow_res,1));
    t_pow(:,isubj,m) = squeeze(nanmean(pow_cnt,1));


  end
end

r_pow = r_pow(:,SUBJLIST,:);
t_pow = t_pow(:,SUBJLIST,:);

load(['~/pconn_all/proc/' sprintf('all_sens_cvar.mat')]);

%% PLOT POWER
clear res cnt

str = 'pow';

cmap = cbrewer('qual', 'Set1', 10,'pchip')

if strcmp(str,'pow')
  res     = squeeze(nanmean(r_pow,2));
  cnt     = squeeze(nanmean(t_pow,2));
  res_err = std(r_pow,[],2)/sqrt(size(r_pow,2));
  cnt_err = std(t_pow,[],2)/sqrt(size(t_pow,2));
  par_res = r_pow;
  par_cnt = t_pow;
elseif strcmp(str,'cvar')
  res     = squeeze(nanmean(nanmean(cvar_res_all,1),3));
  cnt     = squeeze(nanmean(nanmean(cvar_cnt_all,1),3));
  res_err = std(squeeze(nanmean(cvar_res_all)),[],3)/sqrt(size(r_pow,3));
  cnt_err = std(squeeze(nanmean(cvar_cnt_all)),[],3)/sqrt(size(t_pow,3));
  par_res = squeeze(nanmean(cvar_res_all));
  par_cnt = squeeze(nanmean(cvar_cnt_all));
end

res(find(f>46.5,1,'first'):find(f<53.5,1,'last'),:)=nan;
res(find(f>96.5,1,'first'):find(f<103.5,1,'last'),:)=nan;
res(find(f>146.5,1,'first'):find(f<153.5,1,'last'),:)=nan;

cnt(find(f>46.5,1,'first'):find(f<53.5,1,'last'),:)=nan;
cnt(find(f>96.5,1,'first'):find(f<103.5,1,'last'),:)=nan;
cnt(find(f>146.5,1,'first'):find(f<153.5,1,'last'),:)=nan;

res  = fixgaps(res,'pchip');
cnt  = fixgaps(cnt,'pchip');

maxf = find(f==100,1,'first');

res     = res(1:maxf,:);
res_err = res_err(1:maxf,:);
cnt     = cnt(1:maxf,:);
cnt_err = cnt_err(1:maxf,:);

para = [];
para.paired = 1;
% para.clusterthresh 

stats_r1 = tp_clustperm_1D(squeeze(par_res(1:maxf,:,1)),squeeze(par_res(1:maxf,:,2)),para);
if isempty(find(stats_r1.h))
  r_m1 = [];
else
  r_m1 = stats_r1.clustloc{find(stats_r1.h)};
end

stats_r2 = tp_clustperm_1D(squeeze(par_res(1:maxf,:,1)),squeeze(par_res(1:maxf,:,3)),para);
if isempty(find(stats_r2.h))
  r_m2 = [];
else
  r_m2 = stats_r2.clustloc{find(stats_r2.h)};
end

stats_t1 = tp_clustperm_1D(squeeze(par_cnt(1:maxf,:,1)),squeeze(par_cnt(1:maxf,:,2)),para);
if isempty(find(stats_t1.h))
  t_m1 = [];
else
  t_m1 = stats_t1.clustloc{find(stats_t1.h)};
end

stats_t2 = tp_clustperm_1D(squeeze(par_cnt(1:maxf,:,1)),squeeze(par_cnt(1:maxf,:,3)),para);
if isempty(find(stats_t2.h))
  t_m2 = [];
else
  t_m2 = stats_t2.clustloc{find(stats_t2.h)};
end
%% %% PLOT POWER
% ----------------------------------
% PLOT HERE
% ----------------------------------
cmap = cbrewer('qual', 'Set1', 10,'pchip');
cmap = [0.7 0.7 0.7; cmap(1,:); cmap(2,:)]
% cmap = cmap([1 2 6],:)

figure; set(gcf,'color','w'); hold on

subplot(1,2,1); hold on

shadedErrorBar(log2(f(1:maxf)),res(:,1),res_err(:,1),{'color',cmap(1,:),'linewidth',2}); alpha(0.4)
shadedErrorBar(log2(f(1:maxf)),res(:,2),res_err(:,2),{'color',cmap(2,:),'linewidth',2}); alpha(0.4)
shadedErrorBar(log2(f(1:maxf)),res(:,3),res_err(:,3),{'color',cmap(3,:),'linewidth',2}); alpha(0.4)

plot(log2(f(1:maxf)),res(:,1),'linewidth',2,'color',cmap(1,:));
plot(log2(f(1:maxf)),res(:,2),'linewidth',2,'color',cmap(2,:));
plot(log2(f(1:maxf)),res(:,3),'linewidth',2,'color',cmap(3,:));

set(gca,'xtick',[1 3 5 7],'xticklabel',[2 8  32  128],'tickdir','out')
set(gca,'linewidth',2,'ticklength',[0.03 0.03]);
if strcmp(str,'pow')
axis([ 0.5 7 0 5e-28])
elseif strcmp(str,'cvar')
axis([ 0.5 7 0 3e-28])
end

% legend('R-PBO','R-ATX','R-DPZ'); 
axis square

if ~isempty(r_m1); line([f(r_m1(1)) f(r_m1(end))],[5e-28 5e-28],'linewidth',3,'color','k'); end
if ~isempty(r_m2); line([log2(f(r_m2(1))) log2(f(r_m2(end)))],[5e-28 5e-28],'linewidth',3,'color','k'); end

% plot(log2(f(r_m1)),repmat(-90,length(r_m1),1),'.r','linewidth',3)
% plot(log2(f(r_m2)),repmat(-90.5,length(r_m2),1),'.b','linewidth',3)

set(gca,'fontsize',12,'fontweight','bold');
xlabel('Frequency [Hz]','fontsize',13,'fontweight','bold'); 
if strcmp('str','pow')
  ylabel('Power','fontsize',13,'fontweight','bold');
elseif strcmp(str,'cvar')
  ylabel('Coef. of variation','fontsize',13,'fontweight','bold');
end
% axis([0.5 7 -98 -90])

subplot(1,2,2); hold on

shadedErrorBar(log2(f(1:maxf)),cnt(:,1),cnt_err(:,1),{'color',cmap(1,:),'linewidth',2}); alpha(0.4)
shadedErrorBar(log2(f(1:maxf)),cnt(:,2),cnt_err(:,2),{'color',cmap(2,:),'linewidth',2}); alpha(0.4)
shadedErrorBar(log2(f(1:maxf)),cnt(:,3),cnt_err(:,3),{'color',cmap(3,:),'linewidth',2}); alpha(0.4)

plot(log2(f(1:maxf)),cnt(:,1),'linewidth',2,'color',cmap(1,:));
plot(log2(f(1:maxf)),cnt(:,2),'linewidth',2,'color',cmap(2,:));
plot(log2(f(1:maxf)),cnt(:,3),'linewidth',2,'color',cmap(3,:));


if ~isempty(t_m1); line([f(t_m1(1)) f(t_m1(end))],[-90 -90],'linewidth',3,'color','k'); end
if ~isempty(t_m2); line([log2(f(t_m2(1))) log2(f(t_m2(end)))],[-90 -90],'linewidth',3,'color','k'); end

set(gca,'xtick',[1 3 5 7],'xticklabel',[2 8  32  128],'tickdir','out')
set(gca,'linewidth',2,'ticklength',[0.03 0.03]);
if strcmp(str,'pow')
axis([ 0.5 7 0 5e-28])
elseif strcmp(str,'cvar')
axis([ 0.5 7 0 3e-28])
end
set(gca,'fontsize',12,'fontweight','bold'); axis square
xlabel('Frequency [Hz]','fontsize',13,'fontweight','bold'); 

print(gcf,'-depsc2',sprintf('~/pconn_all/plots/all_sens_%sspec_v%d.pdf',str,1))

%% 

col = {[0.7 0.7 0.7];[1 0.4 0];[0 0.5 1]}

[t1,p1]=ttest(squeeze(r_pow(:,:,1)),squeeze(r_pow(:,:,2)),'dim',2);
[t2,p2]=ttest(squeeze(r_pow(:,:,1)),squeeze(r_pow(:,:,3)),'dim',2);

figure; set(gcf,'color','w'); hold on
subplot(1,2,1); hold on

plot(log2(f),-log10(p1),'color',col{2},'linewidth',3)
plot(log2(f),-log10(p2),'color',col{3},'linewidth',3)

line([1 7],[1.3 1.3],'color','k','linestyle','--')

axis([1 7 0 3]); axis square
xlabel('Frequency [Hz]','fontsize',13,'fontweight','bold'); 
set(gca,'xtick',[1 3 5 7],'xticklabel',[2 8  32  128],'tickdir','out')
set(gca,'ytick',[0 1 2 3],'yticklabel',[1 0.1 0.01 0.001],'tickdir','out')

set(gca,'fontsize',12,'fontweight','bold'); 
set(gca,'linewidth',2,'ticklength',[0.03 0.03],'tickdir','out');

subplot(1,2,2); hold on


[t1,p1]=ttest(squeeze(t_pow(:,:,1)),squeeze(t_pow(:,:,2)),'dim',2);
[t2,p2]=ttest(squeeze(t_pow(:,:,1)),squeeze(t_pow(:,:,3)),'dim',2);

plot(log2(f),-log10(p1),'color',col{2},'linewidth',3)
plot(log2(f),-log10(p2),'color',col{3},'linewidth',3)

line([1 7],[1.3 1.3],'color','k','linestyle','--')
set(gca,'ytick',[0 1 2 3],'yticklabel',[1 0.1 0.01 0.001],'tickdir','out')

axis([1 7 0 3]); axis square
xlabel('Frequency [Hz]','fontsize',13,'fontweight','bold'); 
set(gca,'xtick',[1 3 5 7],'xticklabel',[2 8  32  128],'tickdir','out')

set(gca,'fontsize',12,'fontweight','bold'); 
set(gca,'linewidth',2,'ticklength',[0.03 0.03],'tickdir','out');

print(gcf,'-depsc2',sprintf('~/pconn_all/plots/all_sens_powspec_pval_v%d.eps',v))

error('Done!')

%% PLOT DFA OVER FREQ

v = 24;
ord = pconn_randomization;

for isubj = SUBJLIST
  for m = 1 : 3
  	im = find(ord(isubj,:)==m);

    for ifoi = 1 : 74
      
      load(sprintf(['~/pconn_cnt/proc/dfa/' 'pconn_cnt_sens_dfa_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,v));
      
      cnt_alldfa(isubj,m,ifoi) = nanmean(nanmean(par.dfa,2)); clear par
      
      load(sprintf(['~/pconn/proc/dfa/' 'pconn_sens_dfa_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,v));
      
      res_alldfa(isubj,m,ifoi) = nanmean(nanmean(par.dfa,2)); clear par
   
    end
  end
end

res_alldfa = res_alldfa(SUBJLIST,:,:);
cnt_alldfa = cnt_alldfa(SUBJLIST,:,:);

%%
maxf = 74;

par_res     = permute(res_alldfa,[3 1 2]);
par_cnt     = permute(cnt_alldfa,[3 1 2]);

res_err = std(par_res,[],2)/sqrt(size(res,2));
cnt_err = std(par_cnt,[],2)/sqrt(size(cnt,2));
  
para = [];
para.paired = 1;
% para.clusterthresh 

stats_r1 = tp_clustperm_1D(squeeze(par_res(1:maxf,:,1)),squeeze(par_res(1:maxf,:,2)),para);
if isempty(find(stats_r1.h))
  r_m1 = [];
else
  r_m1 = stats_r1.clustloc{find(stats_r1.h)};
end

stats_r2 = tp_clustperm_1D(squeeze(par_res(1:maxf,:,1)),squeeze(par_res(1:maxf,:,3)),para);
if isempty(find(stats_r2.h))
  r_m2 = [];
else
  r_m2 = stats_r2.clustloc{find(stats_r2.h)};
end

stats_t1 = tp_clustperm_1D(squeeze(par_cnt(1:maxf,:,1)),squeeze(par_cnt(1:maxf,:,2)),para);
if isempty(find(stats_t1.h))
  t_m1 = [];
else
  t_m1 = stats_t1.clustloc{find(stats_t1.h)};
end

stats_t2 = tp_clustperm_1D(squeeze(par_cnt(1:maxf,:,1)),squeeze(par_cnt(1:maxf,:,3)),para);
if isempty(find(stats_t2.h))
  t_m2 = [];
else
  t_m2 = stats_t2.clustloc{find(stats_t2.h)};
end

%% %% PLOT DFA SPECTRM 
% ----------------------------------
% PLOT HERE
% ----------------------------------
f=mean(foi,2);
res = squeeze(nanmean(par_res,2));
cnt = squeeze(nanmean(par_cnt,2));

cmap = cbrewer('qual', 'Set1', 10,'pchip');
cmap = [0.7 0.7 0.7; cmap(1,:); cmap(2,:)]
% cmap = cmap([1 2 6],:)

figure; set(gcf,'color','w'); hold on

subplot(1,2,1); hold on

shadedErrorBar(log2(f(1:maxf)),res(:,1),res_err(:,1),{'color',cmap(1,:),'linewidth',2}); alpha(0.4)
shadedErrorBar(log2(f(1:maxf)),res(:,2),res_err(:,2),{'color',cmap(2,:),'linewidth',2}); alpha(0.4)
shadedErrorBar(log2(f(1:maxf)),res(:,3),res_err(:,3),{'color',cmap(3,:),'linewidth',2}); alpha(0.4)

plot(log2(f(1:maxf)),res(:,1),'linewidth',2,'color',cmap(1,:));
plot(log2(f(1:maxf)),res(:,2),'linewidth',2,'color',cmap(2,:));
plot(log2(f(1:maxf)),res(:,3),'linewidth',2,'color',cmap(3,:));

set(gca,'xtick',[1 3 5 7],'xticklabel',[2 8  32  128],'tickdir','out')
set(gca,'linewidth',2,'ticklength',[0.03 0.03]);
axis([ 0.5 7 0.52 0.78])

axis square

if ~isempty(r_m1); line([f(r_m1(1)) f(r_m1(end))],[.7 .7],'linewidth',3,'color','k'); end
if ~isempty(r_m2); line([log2(f(r_m2(1))) log2(f(r_m2(end)))],[.72 .72],'linewidth',3,'color','k'); end

set(gca,'fontsize',12,'fontweight','bold');
xlabel('Frequency [Hz]','fontsize',13,'fontweight','bold'); 
ylabel('\beta','fontsize',13,'fontweight','bold');

subplot(1,2,2); hold on

shadedErrorBar(log2(f(1:maxf)),cnt(:,1),cnt_err(:,1),{'color',cmap(1,:),'linewidth',2}); alpha(0.4)
shadedErrorBar(log2(f(1:maxf)),cnt(:,2),cnt_err(:,2),{'color',cmap(2,:),'linewidth',2}); alpha(0.4)
shadedErrorBar(log2(f(1:maxf)),cnt(:,3),cnt_err(:,3),{'color',cmap(3,:),'linewidth',2}); alpha(0.4)

plot(log2(f(1:maxf)),cnt(:,1),'linewidth',2,'color',cmap(1,:));
plot(log2(f(1:maxf)),cnt(:,2),'linewidth',2,'color',cmap(2,:));
plot(log2(f(1:maxf)),cnt(:,3),'linewidth',2,'color',cmap(3,:));


if ~isempty(t_m1); line([t_m1(1) t_m1(end)],[.73 .73],'linewidth',3,'color','k'); end
if ~isempty(t_m2); line([t_m2(1) t_m2(end)],[.72 .72],'linewidth',3,'color','k'); end

set(gca,'xtick',[1 3 5 7],'xticklabel',[2 8  32  128],'tickdir','out')
set(gca,'linewidth',2,'ticklength',[0.03 0.03]);
axis([ 0.5 7 0.52 0.78])

set(gca,'fontsize',12,'fontweight','bold'); axis square
xlabel('Frequency [Hz]','fontsize',13,'fontweight','bold'); 

print(gcf,'-depsc2',sprintf('~/pconn_all/plots/all_sens_%sspec_v%d.pdf','dfa',1))


%%

f = mean(foi,2); 
col = cbrewer('qual', 'Set1', 10,'pchip');
col = [0.7 0.7 0.7; col(1,:); col(2,:)]

figure;set(gcf,'color','w'); hold on; set(gca,'tickdir','out')
alltitle = {'DFA'}

for iall = 1 : 2
  if iall == 1
    par = alldfa;
  elseif iall == 2
%     par = allcvar;
  elseif iall == 3
%     par = allcvar;
  elseif iall == 4
%     par = allamp;
  end
  
  subplot(1,2,1); hold on; title(alltitle{iall})

  for i = 1 : 3
    plot(log2(f),squeeze(nanmean(par(:,i,1:end),1)),'color',col(i,:),'linewidth',3)
  end

%   t = logical(squeeze(ttest(par(:,2,:),par(:,1,:),'dim',1)));
%   plot(log2(f(t)),log2(ones(sum(t),1)*min(par(:))),'*','color',col{2})
%   t = logical(squeeze(ttest(par(:,3,:),par(:,1,:),'dim',1)));
%   plot(log2(f(t)),log2(ones(sum(t),1)*min(par(:))),'*','color',col{3})
  
  set(gca,'xtick',[1  3  5  7],'xticklabel',[2  8  32  128],'tickdir','out')
  set(gca,'linewidth',2,'ticklength',[0.03 0.03]);
end

  print(gcf,'-depsc2',sprintf('~/pconn_all/plots/all_para_across_freq_v%d.eps',v))

figure;set(gcf,'color','w'); hold on;
alltitle = {'DFA'}

for iall = 1 : 1
  if iall == 1
    par = alldfa;
  elseif iall == 2
%     par = allvar;
  elseif iall == 3
%     par = allcvar;
  elseif iall == 4
%     par = allamp;
  end
  
  subplot(1,2,2); hold on; title(alltitle{iall})

  [t1,p1] = ttest(par(:,2,:),par(:,1,:),'dim',1);
  [t2,p2] = ttest(par(:,3,:),par(:,1,:),'dim',1);

      line([1 7],[1.3 1.3],'color','k','linestyle','--')

    plot(log2(f),-log10(squeeze(p1)),'color',col(2,:),'linewidth',3)
   plot(log2(f),-log10(squeeze(p2)),'color',col(3,:),'linewidth',3)
    axis([0 8 -0.2 3.5])
set(gca,'xtick',[1 3  5  7],'xticklabel',[2 8  32  128],'tickdir','out')
     
     set(gca,'linewidth',2,'ticklength',[0.03 0.03]);

end

  print(gcf,'-depsc2',sprintf('~/pconn_all/plots/all_para_across_freq_pval_v%d.eps',v))





