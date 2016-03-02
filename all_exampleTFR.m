%% COMPUTE EXCITATION/INHIBITION BALANCE ACROSS SENSOR SPACE
% pconn_cnt_sens_planar_dfa

% --------------------------------------------------------
% VERSION 19
% --------------------------------------------------------
v         = 1;
v_rawdata = 1;

addpath /home/tpfeffer/Documents/MATLAB/fieldtrip-20130925/
addpath ~/Documents/MATLAB/

outdir   = '/home/tpfeffer/pconn/proc/dfa/';

addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/
%%
m = 1;
isubj = 4;
%       
if isubj == 3
  load ~/pconn/matlab/pconn_sensorlabels.mat
end
      
iblock = 1;

load(sprintf('/home/tpfeffer/pconn_cnt/proc/preproc/pconn_cnt_postproc_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,v_rawdata));

cfg = [];
cfg.method      = 'template';
cfg.template    = 'CTF275_neighb.mat';
cfg.neighbours  = ft_prepare_neighbours(cfg);
cfg.method      = 'sincos';
data_pl         = ft_megplanar(cfg,data_low);



cfg = [];
cfg.method = 'mtmconvol';
cfg.output = 'pow';
cfg.foi = 2:0.08:30;
cfg.taper = 'hanning';
cfg.tapsmofrq = 0.1*cfg.foi;
cfg.toi = 20:0.008:40;
cfg.t_ftimwin = 5./cfg.foi;
dat = ft_freqanalysis(cfg,data_low);

cfg = [];
cfg.bpfreq = [8 12];
cfg.bpfilter = 'yes';

raw = ft_preprocessing(cfg,data_low);

cfg = [];
cfg.bpfreq = [8 12];
cfg.bpfilter = 'yes';
cfg.hilbert = 'complex';

hilb = ft_preprocessing(cfg,data_low);


save([outdir sprintf('pconn_exampleTFR.mat')],'dat','-v7.3');

error('!')
%%

load([outdir sprintf('pconn_exampleTFR.mat')]);


a  = squeeze(nanmean(dat.powspctrm));
sc = [min(log10(a(:))) max(log10(a(:)))];
%%

figure; set(gcf,'color','white');

subplot(2,1,1);

imagesc(log10(squeeze(nanmean(dat.powspctrm(:,1:351,1:1250)))))
% 
set(gca,'TickDir','out');
set(gca,'YDir','normal');
set(gca,'YTick',[1 26 51 76 101 126 151 176 201 226 251 276 301 326 351],'YTickLabel',cfg.foi([1 26 51 76 101 126 151 176 201 226 251 276 301 326 351]))
set(gca,'XTick',[ 250  500  750  1000 ],'XTickLabel',[2 4 6 8 10])

ylabel('Frequency [Hz]');
xlabel('Time [s]')
colormap(jet)

%

t = round(cfg.toi(1:1250).*400);
% tmp = resample(tmp,1,2);

subplot(2,1,2)



all_t = 1 : length(t);
A = nanmean(abs(hilb.trial{1}(:,t)),1);
f = 10;
y = A.*sin(2*pi*f*all_t/125);
plot(all_t,y)

hold on

plot(all_t,(nanmean(abs(hilb.trial{1}(:,t)))))
ylabel('Amplitude [fT/cm2]');
xlabel('Time [s]')
axis tight
set(gca,'TickDir','out');
set(gca,'XTick',[1 250  500  750  1000 1250],'XTickLabel',[0 2 4 6 8 10])
