%% COMPUTE POWER SPECTRUM FOR ALL SENSORS
% fits a straight line to the powerspectrum in order
% to identify peaks in the spectrum for subsequent
% DFA analysis.

% pconn_all_powspec_norm

clear

% last update: 18-02-2015, tpfeffer

% to be implemented:
% ***nothing left to implement***

% --------------------------------------------------------
% VERSION 6 - NO PLANAR TRANSFORMATION
% --------------------------------------------------------
v         = 1;
fsample   = 400;
SUBJLIST 	= [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
FOI       = 2:0.2:200;
T         = 0.2;
planar    = 0;
% --------------------------------------------------------

restoredefaultpath

addpath('/home/tpfeffer/Documents/MATLAB/fieldtrip-20160919')
addpath('/home/tpfeffer/pconn/matlab')

ft_defaults

indir   = '/home/tpfeffer/pconn/proc/src/';
outdir   = '/home/tpfeffer/pconn/proc/powspec/';
plotdir = '/home/tpfeffer/pconn/proc/plots/';

%%
% ord   = pconn_randomization;

for m = 3 : 3
  for isubj = 33
    
% %     im = find(ord(isubj,:)==m);
%     if ~exist(sprintf([outdir 'pconn_all_powspec_norm_s%d_m%d_v%d_processing.txt'],isubj,m,v))
%       system(['touch ' outdir sprintf('pconn_all_powspec_norm_s%d_m%d_v%d_processing.txt',isubj,m,v)]);
%     else
%       continue
%     end
    
    disp(sprintf('Processing s%d m%d ...', isubj,m))
        
    for iblock = 1 : 1
            
      try 
        load(sprintf('/home/tpfeffer/pconn/proc/preproc/pconn_postpostproc_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,1));
      catch me
        powspec(:,:,1) = nan(274,length(FOI));
        continue
      end
      cfg = [];
      cfg.length = 1/T;
      cfg.overlap = 0;
      
      data      = ft_redefinetrial(cfg,data); clear data_low

      cfg         = [];      
      cfg.trials  = 1:size(data.trial,2)-1;
      data        = ft_redefinetrial(cfg,data);

      cfg             = [];
      cfg.method      = 'mtmfft';
      cfg.output      = 'pow';
      cfg.taper       = 'dpss';
      cfg.channel     = {'MEG'};
      cfg.foi         = FOI;
      cfg.keeptrials  = 'yes';
      cfg.pad = 'nextpow2';
      cfg.tapsmofrq   = 2;
      dat             = ft_freqanalysis(cfg, data); % run freqanalysis here
      
      powspec(:,:,1) = pconn_sens_interp274_alena(isubj,squeeze(mean(dat.powspctrm(:,:,:)))) %./ repmat(sum(squeeze(mean(dat.powspctrm(:,:,:))),2),[1 991]),m);
      clear dat data
      
      try 
        load(sprintf('/home/tpfeffer/pconn_cnt/proc/preproc/pconn_cnt_postpostproc_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,1));
      catch me
        powspec(:,:,2) = nan(274,length(FOI));
        save([outdir sprintf('pconn_all_powspec_norm_s%d_b%d_m%d_v%d.mat',isubj,iblock,m,v)],'powspec','-v7.3');

        continue
      end
      cfg = [];
      cfg.length = 1/T;
      cfg.overlap = 0;
      
      data      = ft_redefinetrial(cfg,data); clear data_low

      cfg         = [];      
      cfg.trials  = 1:size(data.trial,2)-1;
      data        = ft_redefinetrial(cfg,data);

      cfg             = [];
      cfg.method      = 'mtmfft';
      cfg.output      = 'pow';
      cfg.taper       = 'dpss';
      cfg.channel     = {'MEG'};
      cfg.pad = 'nextpow2';
      cfg.foi         = FOI;
      cfg.keeptrials  = 'yes';
      cfg.tapsmofrq   = 2;
      dat             = ft_freqanalysis(cfg, data); % ru
      
      powspec(:,:,2) = pconn_sens_interp274_alena(isubj,squeeze(mean(dat.powspctrm(:,:,:)))) %./ repmat(sum(squeeze(mean(dat.powspctrm(:,:,:))),2),[1 991]),m);;  
      
      save([outdir sprintf('pconn_all_powspec_norm_s%d_b%d_m%d_v%d.mat',isubj,iblock,m,v)],'powspec','-v7.3');
      clear data dat powspec 
      
    end
  end
end

error('Stop here!')

%% PLOT
ord = pconn_randomization;

for isubj = SUBJLIST
  isubj
  for m = 1 : 3
    for iblock = 1 : 2
      
      im = find(ord(isubj,:)==m);
       
        load([outdir sprintf('pconn_all_powspec_norm_s%d_b%d_m%d_v%d.mat',isubj,iblock,im,v)]);

        pow(isubj,m,:,:,:,iblock) = powspec;
 
      
    end
  end
end

pow = nanmean(pow(SUBJLIST,:,:,:,:,:),6);
%%
figure; hold on

plot(log(FOI(1:800)),log10(squeeze(nanmean(nanmean(pow(:,1,:,1:800,1),3)))))
plot(log(FOI(1:800)),log10(squeeze(nanmean(nanmean(pow(:,1,:,1:800,2),3)))))

%%
chanidx = startsWith(data.label,'MRO');
chanidx = chanidx + startsWith(data.label,'MLO')

% figure; hold on
load ~/pconn/proc/src/pconn_sa_s26_m2_b1_v4.mat

% plot(log(FOI(1:800)),log10(squeeze(nanmean(nanmean(pow(:,1,find(chanidx),1:800,1),3)))))
% plot(log(FOI(1:800)),log10(squeeze(nanmean(nanmean(pow(:,1,find(chanidx),1:800,2),3)))))
im = 2

diff = 100*(squeeze(nanmean(nanmean(pow(:,im,:,251:700,2),4)))-squeeze(nanmean(nanmean(pow(:,im,:,251:700,1),4)))) ./ squeeze(nanmean(nanmean(pow(:,im,:,251:700,1),4)))
[tt,p,~,tmp] = ttest(squeeze(nanmean(pow(:,im,:,251:700,2),4)),squeeze(nanmean(pow(:,im,:,251:700,1),4)),'alpha', 0.05)


t = tmp.tstat;
t(t < 2.5) = 0
pars = []
% diff=diff.*tt'
pars.linewidth = 9
figure; set(gcf,'color','white');
% subplot(1,3,1)
pars.scale = [-20 20]
% pars.mask = p<fdr1(p',0.05,0)
pars.mask = p<0.005
pconn_showfield(diff,sa.locs_2D,pars)
colormap(hot)
% pconn_showfield(squeeze(nanmean(diff,1)),sa.locs_2D,pars)
% print(gcf,'-depsc',sprintf('~/dfa_final/plots/figure4c.eps'))

% ttest(squeeze(nanmean(pow(:,1,find(chanidx),1:800,1),3)),squeeze(nanmean(pow(:,1,find(chanidx),1:800,2),3)),'dim',1)