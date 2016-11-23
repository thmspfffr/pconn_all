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

addpath('/home/tpfeffer/Documents/MATLAB/fieldtrip-20160919/')

addpath /home/tpfeffer/pconn/matlab
ft_defaults

indir   = '/home/tpfeffer/pconn/proc/src/';
outdir   = '/home/tpfeffer/pconn/proc/powspec/';
plotdir = '/home/tpfeffer/pconn/proc/plots/';
addpath ~/pconn/matlab/
%%
ord       = pconn_randomization;


for m = 1:3

  for isubj = SUBJLIST
    
  disp(sprintf('Processing subject %d ...',isubj));
  
   im = find(ord(isubj,:)==m);
  
    if ~exist(sprintf(['~/pconn_all/proc/' 'all_sens_powspec_s%d_m%d_v%d_processing.txt'],isubj,m,v))
      system(['touch ' '~/pconn_all/proc/' sprintf('all_sens_powspec_s%d_m%d_v%d_processing.txt',isubj,m,v)]);
    else
      continue
    end
   
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

%%
clear res cnt

res     = log2(squeeze(nanmean(r_pow,2)));
cnt     = log2(squeeze(nanmean(t_pow,2)));
res_err = std(log2(r_pow),[],2)/sqrt(size(r_pow,2));
cnt_err = std(log2(t_pow),[],2)/sqrt(size(t_pow,2));

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

stats_r1 = tp_clustperm_1D(squeeze(r_pow(1:maxf,:,1)),squeeze(r_pow(1:maxf,:,2)),para);
if isempty(find(stats_r1.h))
  r_m1 = [];
else
  r_m1 = stats_r1.clustloc{find(stats_r1.h)};
end

stats_r2 = tp_clustperm_1D(squeeze(r_pow(1:maxf,:,1)),squeeze(r_pow(1:maxf,:,3)),para);
if isempty(find(stats_r2.h))
  r_m2 = [];
else
  r_m2 = stats_r2.clustloc{find(stats_r2.h)};
end

stats_t1 = tp_clustperm_1D(squeeze(t_pow(1:maxf,:,1)),squeeze(t_pow(1:maxf,:,2)),para);
if isempty(find(stats_t1.h))
  t_m1 = [];
else
  t_m1 = stats_t1.clustloc{find(stats_t1.h)};
end

stats_t2 = tp_clustperm_1D(squeeze(t_pow(1:maxf,:,1)),squeeze(t_pow(1:maxf,:,3)),para);
if isempty(find(stats_t2.h))
  t_m2 = [];
else
  t_m2 = stats_t1.clustloc{find(stats_t2.h)};
end
% ----------------------------------
% PLOT HERE
% ----------------------------------

figure; set(gcf,'color','w'); hold on

subplot(2,2,1); hold on

shadedErrorBar(log2(f(1:maxf)),res(:,1),res_err(:,1),{'color',cmap(1,:)}); %alpha(0.4)
shadedErrorBar(log2(f(1:maxf)),res(:,2),res_err(:,2),{'color',cmap(2,:)}); %alpha(0.4)
shadedErrorBar(log2(f(1:maxf)),res(:,3),res_err(:,3),{'color',cmap(3,:)}); %alpha(0.4)

set(gca,'xtick',[1 3 5 7],'xticklabel',[2 8  32  128],'tickdir','out')
set(gca,'linewidth',2,'ticklength',[0.03 0.03]);

% legend('R-PBO','R-ATX','R-DPZ'); 
axis square

if ~isempty(r_m1); line([f(r_m1(1)) f(r_m1(end))],[-90 -90],'linewidth',3,'color','k'); end
if ~isempty(r_m2); line([log2(f(r_m2(1))) log2(f(r_m2(end)))],[-90 -90],'linewidth',3,'color','k'); end

% plot(log2(f(r_m1)),repmat(-90,length(r_m1),1),'.r','linewidth',3)
% plot(log2(f(r_m2)),repmat(-90.5,length(r_m2),1),'.b','linewidth',3)

set(gca,'fontsize',12,'fontweight','bold');
xlabel('Frequency [Hz]','fontsize',13,'fontweight','bold'); 
ylabel('Log_{2} (Power)','fontsize',13,'fontweight','bold');
axis([0.5 7 -98 -90])

subplot(2,2,2); hold on

shadedErrorBar(log2(f(1:maxf)),cnt(:,1),cnt_err(:,1),{'color',cmap(1,:)}); %alpha(0.4)
shadedErrorBar(log2(f(1:maxf)),cnt(:,2),cnt_err(:,2),{'color',cmap(2,:)}); %alpha(0.4)
shadedErrorBar(log2(f(1:maxf)),cnt(:,3),cnt_err(:,3),{'color',cmap(3,:)}); %alpha(0.4)

if ~isempty(t_m1); line([f(t_m1(1)) f(t_m1(end))],[-90 -90],'linewidth',3,'color','k'); end
if ~isempty(t_m2); line([log2(f(t_m2(1))) log2(f(t_m2(end)))],[-90 -90],'linewidth',3,'color','k'); end

set(gca,'xtick',[1 3 5 7],'xticklabel',[2 8  32  128],'tickdir','out')
set(gca,'linewidth',2,'ticklength',[0.03 0.03]);

set(gca,'fontsize',12,'fontweight','bold'); axis square
xlabel('Frequency [Hz]','fontsize',13,'fontweight','bold'); 
% ylabel('Log_{2} (Power)','fontsize',13,'fontweight','bold');
axis([0.5 7 -98 -90])

% -----------------------------------
% PLOT DETAILS
% -----------------------------------

subplot(2,2,3); hold on

shadedErrorBar(log2(f(r_m2(1):r_m2(end))),cnt(:,1),cnt_err(:,1),{'color',cmap(1,:)}); %alpha(0.4)
shadedErrorBar(log2(f(r_m2(1):r_m2(end))),cnt(:,2),cnt_err(:,2),{'color',cmap(2,:)}); %alpha(0.4)
shadedErrorBar(log2(f(r_m2(1):r_m2(end))),cnt(:,3),cnt_err(:,3),{'color',cmap(3,:)}); %alpha(0.4)

if ~isempty(t_m1); line([f(t_m1(1)) f(t_m1(end))],[-90 -90],'linewidth',3,'color','k'); end
if ~isempty(t_m2); line([log2(f(t_m2(1))) log2(f(t_m2(end)))],[-90 -90],'linewidth',3,'color','k'); end

set(gca,'xtick',[1 3 5 7],'xticklabel',[2 8  32  128],'tickdir','out')
set(gca,'linewidth',2,'ticklength',[0.03 0.03]);

set(gca,'fontsize',12,'fontweight','bold'); axis square
xlabel('Frequency [Hz]','fontsize',13,'fontweight','bold'); 
% ylabel('Log_{2} (Power)','fontsize',13,'fontweight','bold');
axis([0.5 7 -98 -90])



print(gcf,'-depsc2',sprintf('~/pconn_all/plots/all_sens_powspec_v%d.eps',v))

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
%%
idx(3) = find(f>16.5&f<16.7)
idx(4) = find(f>16.7&f<16.9);


figure; set(gcf,'color','white'); hold on

p_res = squeeze(nanmean(nanmean(pow_res(:,1:length(f),:,:),1),4));
p_cnt = squeeze(nanmean(nanmean(pow_cnt(:,1:length(f),:,:),1),4));

m_pow_res = nanmean(p_res,2);
m_pow_cnt = nanmean(p_cnt,2);

s_pow_res = nanstd(p_res,[],2)/sqrt(size(p_res,2));
s_pow_cnt = nanstd(p_cnt,[],2)/sqrt(size(p_cnt,2));

% CUT OUT BS FILTER SEGMENTS
% m_pow_res(idx(1):idx(2),:)=NaN;
% m_pow_cnt(idx(1):idx(2),:)=NaN;
m_pow_res(idx(3):idx(4),:)=NaN;
m_pow_cnt(idx(3):idx(4),:)=NaN;

m_pow_res=fixgaps(m_pow_res,'pchip');
m_pow_cnt=fixgaps(m_pow_cnt,'pchip');

save('~/powerstuff.mat','m_pow_res','m_pow_cnt','s_pow_res','s_pow_res','f','idx');

% [t,p]= ttest(squeeze(nanmean(nanmean(pow_res(:,:,:,:),4),1)),squeeze(nanmean(nanmean(pow_cnt(:,1:length(f),:,:),4),1)),'dim',2);
% idx = find(p<0.01);

% plot(log10(2:0.2:max(f)),nanmean(m_pow_cnt(1:length(f),:),2),'linewidth',4);

x = log10(2:0.2:max(f));
X = [x,fliplr(x)];
Y = [m_pow_res'-s_pow_res',fliplr(m_pow_res'+s_pow_res')];
fill(X,Y,'k','edgecolor','none');

Y = [m_pow_cnt'-s_pow_cnt',fliplr(m_pow_cnt'+s_pow_cnt')];
fill(X,Y,'b','edgecolor','none');

% alpha(0.25)


plot(log10(2:0.2:max(f)),m_pow_res,'k','linewidth',4);
plot(log10(2:0.2:max(f)),m_pow_cnt,'b','linewidth',4);


set(gca,'tickdir','out','xtick',log10([2 4 8 16 32]),'xticklabel',[2 4 8 16 32]);

saveas(gcf,sprintf('~/pconn_all/plots/pconn_powspec_v%d.fig',v));
% print(gcf,'-deps',sprintf('~/pconn_all/plots/pconn_powspec_v%d.eps',v))

%%
set(gca,'TickDir','out')

set(gca,'XTick',log10([2 4 8 16 32 64]),'XTickLabel',[2 4 8 16 32 64])
xlabel('Frequency [Hz]');
ylabel('Log-Power');

plot(log10(2:0.2:98),-27*t(1:481))

plot(2:0.2:200,-log10(p))

%%

plot(log10(2:0.2:100),log10(m_pow_res(1:491,1)))
plot(log10(2:0.2:100),log10(m_pow_res(1:491,2)))
plot(log10(2:0.2:100),log10(m_pow_res(1:491,3)))

plot(log10(2:0.2:100),log10(m_pow_cnt(1:491,1)))
plot(log10(2:0.2:100),log10(m_pow_cnt(1:491,2)))
plot(log10(2:0.2:100),log10(m_pow_cnt(1:491,3)))

%% QUALITY CHECK OF SUBJECTS
im = 1;


isubj = 21;

m = ord(isubj,im);

iblock = 1;

load(sprintf('/home/tpfeffer/pconn/proc/preproc/pconn_postpostproc_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,1));

% cfg = [];
%       cfg.method    = 'summary';
%       tmp_data      = ft_rejectvisual(cfg,data);

figure; set(gcf,'color','white');
subplot(1,2,1); hold on
plot(zscore(var(data.trial{1},[],2))); tp_editplots
%%
idx = find(zscore(var(data.trial{1},[],2))>5);

segleng = 20;

nseg = floor(size(data.trial{1},2)/segleng);

for i = 1 : nseg
  v(i)=var(nanmean(data.trial{1}(:,(i-1)*segleng+1:i*segleng),1));
end
      
     v(i+1)=var(nanmean(data.trial{1}(:,i*segleng:end),1));

subplot(1,2,2); hold on
plot(zscore(v));

% iii=find(zscore(v)>20);

% data.trial{1}(:,(iii-1)*segleng+1:iii*segleng)=[];






