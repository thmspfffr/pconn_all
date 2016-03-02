%% COMPUTE POWER SPECTRUM FOR ALL SENSORS
% fits a straight line to the powerspectrum in order
% to identify peaks in the spectrum for subsequent
% DFA analysis.

% pconn_powspec

% last update: 18-02-2015, tpfeffer

% to be implemented:
% ***nothing left to implement***

% --------------------------------------------------------
% VERSION 6 - NO PLANAR TRANSFORMATION
% --------------------------------------------------------
v         = 1;
d         = 'data_low.trial{1}+data_hi.trial{1};';
v_rawdata = 6;
fsample   = 400;
SUBJLIST 	= [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24];
FOI       = 2:0.2:200;
T         = 0.2;
planar    = 0;
% --------------------------------------------------------

restoredefaultpath

addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/
addpath /home/tpfeffer/fieldtrip-20150318/
addpath /home/tpfeffer/Documents/MATLAB/toolboxes/NBT-ReleaseNBTRC4a/
addpath /home/tpfeffer/pconn/matlab
ft_defaults

indir   = '/home/tpfeffer/pconn/proc/src/';
outdir   = '/home/tpfeffer/pconn/proc/powspec/';
plotdir = '/home/tpfeffer/pconn/proc/plots/';

%%
ord       = pconn_randomization;

for m = 1 : 3

  for isubj = SUBJLIST  
 
  disp(sprintf('Processing subject %d ...',isubj));
  
   im = find(ord(isubj,:)==m);

    for iblock = 1 : 2
      

      disp(sprintf('Loading MEG data ...'));
      
      load([outdir sprintf('pconn_powspec_dat_s%d_b%d_m%d_v%d.mat',isubj,iblock,im,v)]);
      
      pow_res(:,:,iblock,isubj,m) = nanmean(dat.powspctrm,1); clear dat

      load(['~/pconn_cnt/proc/powspec/' sprintf('pconn_cnt_powspec_dat_s%d_b%d_m%d_v%d.mat',isubj,iblock,im,v)]);

    	pow_cnt(:,:,iblock,isubj,m) = nanmean(dat.powspctrm,1);

    end
  end
end

error('Stop here!')

%% PLOT EVERYTHING

f = 2:0.2:200;

idx(1) = find(f==48);
idx(2) = find(f==52);

pow_res = squeeze(nanmean(pow_res(:,:,:,SUBJLIST,:),3));

pow_cnt = squeeze(nanmean(pow_cnt(:,:,:,SUBJLIST,:),3));

%%

figure; set(gcf,'color','white'); hold on

m_pow_res = squeeze(nanmean(nanmean(pow_res,3),1));
m_pow_cnt = squeeze(nanmean(nanmean(pow_cnt,3),1));

% CUT OUT BS FILTER SEGMENTS
m_pow_res(idx(1):idx(2),:)=NaN;
m_pow_cnt(idx(1):idx(2),:)=NaN;

[t,p]= ttest(squeeze(nanmean(nanmean(pow_res,4),1)),squeeze(nanmean(nanmean(pow_cnt,4),1)),'dim',2);

plot(log10(2:0.2:98),log10(nanmean(m_pow_res(1:481,:),2)),'linewidth',4);
plot(log10(2:0.2:98),log10(nanmean(m_pow_cnt(1:481,:),2)),'linewidth',4);

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







