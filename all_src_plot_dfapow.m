%% COMPUTE WHOLE BRAIN DFA IN SOURCE SPACE
% pconn_src_dfa

clear all

% v = 1: cortex, eloreta, v = 2: cortex, lcmv, v = 3: coarse, eloreta

% --------------------------------------------------------
% VERSION 1
% --------------------------------------------------------
v         = 1;
v_stat    = 1;
v_rawdata = 2;
SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24];
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


%% PLOT DFA
   
ord   = pconn_randomization;

for ifoi = 1 : 5
  for isubj = SUBJLIST
    for m = 1 : 3

      im = find(ord(isubj,:)==m);

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
      
    end
  end
end

dfa_all_res  = dfa_all_res(:,:,SUBJLIST,:);
var_all_res  = var_all_res(:,:,SUBJLIST,:);
cvar_all_res = cvar_all_res(:,:,SUBJLIST,:);
amp_all_res  = amp_all_res(:,:,SUBJLIST,:);

dfa_all_cnt  = dfa_all_cnt(:,:,SUBJLIST,:);
var_all_cnt  = var_all_cnt(:,:,SUBJLIST,:);
cvar_all_cnt = cvar_all_cnt(:,:,SUBJLIST,:);
amp_all_cnt  = amp_all_cnt(:,:,SUBJLIST,:);

%% PHARMA COMPARISON

ifoi = 3;

if exist(sprintf('~/pconn_all/proc/all_src_clusterstat_tsk_%s_c%d_f%d_v%d.mat','dfa',1,ifoi,v_stat))
  load(sprintf('~/pconn_all/proc/all_src_clusterstat_tsk_%s_c%d_f%d_v%d.mat','dfa',1,ifoi,v_stat))
end

stats.mask=logical(stats.mask);

dfa_res = nanmean(dfa_all_res(stats.mask,:,:,ifoi),1);
dfa_cnt = nanmean(dfa_all_cnt(stats.mask,:,:,ifoi),1);

amp_res = nanmean(amp_all_res(stats.mask,:,:,ifoi),1);
amp_cnt = nanmean(amp_all_cnt(stats.mask,:,:,ifoi),1);


a=linspace(0,1,18);

col = [shuffle(a)' shuffle(a)' shuffle(a)'];
%%
figure_white(15,10); hold on

% subplot(1,2,1); hold on
title('Rest','Fontsize',14); 

for isubj = 1 : 18
  
  plot(squeeze(log10(amp_res(1,1,isubj))),squeeze(dfa_res(1,1,isubj)),'.','color','k','markersize',40)

  line([squeeze(log10(amp_res(1,1,isubj))) squeeze(log10(amp_res(1,2,isubj)))],[squeeze(dfa_res(1,1,isubj)) squeeze(dfa_res(1,2,isubj))],'color','k','linewidth',2)
  
  plot(squeeze(log10(amp_res(1,2,isubj))),squeeze(dfa_res(1,2,isubj)),'o','color','k','markersize',10,'markerfacecolor','w')

end

[r1,p1]=corr(squeeze(dfa_res(1,1,:)),squeeze(log10(amp_res(1,1,:))));
[r2,p2]=corr(squeeze(dfa_res(1,2,:)),squeeze(log10(amp_res(1,2,:))));
[r3,p3]=corr(squeeze(dfa_res(1,2,:))-squeeze(dfa_res(1,1,:)),squeeze(log10(amp_res(1,2,:)))-squeeze(log10(amp_res(1,1,:))));

axis square; box on; axis([-16.6 -15.8 0.48 0.82 ])
ylabel('DFA exponent','Fontsize',14); xlabel('Log_{10}(Amplitude)','Fontsize',14)
set(gca,'tickdir','out','ticklength',2*get(gca,'ticklength'))

text(-16.17,0.55,sprintf('r_{pbo}=%.2f; p_{pbo}=%.2f',r1,p1))
text(-16.17,0.52,sprintf('r_{atx}=%.2f; p_{atx}=%.2f',r2,p2))
text(-16.17,0.49,sprintf('r_{diff}=%.2f; p_{diff}=%.2f',r3,p3))

plot(-15.75,0.8,'k.','markersize',40,'clip','off');
text(-15.72,0.8,sprintf('Placebo'))
plot(-15.75,0.78,'ko','markersize',10,'clip','off');
text(-15.72,0.78,sprintf('Atomoxetine'))
set(gca,'ytick',[.5 .6 .7 .8],'yticklabel',[.5 .6 .7 .8])
set(gca,'xtick',[-16.6:0.2:-15.8],'xticklabel',[-16.6:0.2:-15.8])

dfa = [mean(squeeze(dfa_res(1,1,:))) mean(squeeze(dfa_res(1,2,:)))];
amp = [mean(squeeze(log10(amp_res(1,1,:)))) mean(squeeze(log10(amp_res(1,2,:))))];

% line([amp(1) amp(1)],[dfa(1) dfa(2)],'color','y','linewidth',5)
% line([dfa(1) dfa(1)],[amp(1) amp(2)],'color','y','linewidth',5)

% subplot(1,2,2); hold on
% title('Task','Fontsize',14); 
% 
% 
% for isubj = 1 : 18
%   
%   plot(squeeze(log10(amp_cnt(1,1,isubj))),squeeze(dfa_cnt(1,1,isubj)),'.','color','k','markersize',40)
% 
%   line([squeeze(log10(amp_cnt(1,1,isubj))) squeeze(log10(amp_cnt(1,2,isubj)))],[squeeze(dfa_cnt(1,1,isubj)) squeeze(dfa_cnt(1,2,isubj))],'color','k','linewidth',2)
%   
%   plot(squeeze(log10(amp_cnt(1,2,isubj))),squeeze(dfa_cnt(1,2,isubj)),'o','color','k','markersize',10,'markerfacecolor','w')
% 
% end
% 
% [r1,p1]=corr(squeeze(dfa_cnt(1,1,:)),squeeze(log10(amp_cnt(1,1,:))));
% [r2,p2]=corr(squeeze(dfa_cnt(1,2,:)),squeeze(log10(amp_cnt(1,2,:))));
% [r3,p3]=corr(squeeze(dfa_cnt(1,2,:))-squeeze(dfa_cnt(1,1,:)),squeeze(log10(amp_cnt(1,2,:)))-squeeze(log10(amp_cnt(1,1,:))));
% 
% axis square; box on; axis([-16.6 -15.8 0.48 0.82 ])
% ylabel('DFA exponent','Fontsize',14); xlabel('Log10(Amplitude)','Fontsize',14)
% set(gca,'tickdir','out','ticklength',2*get(gca,'ticklength'))
% 
% text(-16.17,0.55,sprintf('r_{pbo}=%.2f; p_{pbo}=%.2f',r1,p1))
% text(-16.17,0.52,sprintf('r_{atx}=%.2f; p_{atx}=%.2f',r2,p2))
% text(-16.17,0.49,sprintf('r_{diff}=%.2f; p_{diff}=%.2f',r3,p3))
% 
% plot(-15.75,0.8,'k.','markersize',40,'clip','off');
% text(-15.72,0.8,sprintf('Placebo'))
% plot(-15.75,0.78,'ko','markersize',10,'clip','off');
% text(-15.72,0.78,sprintf('Atomoxetine'))
% set(gca,'ytick',[.5 .6 .7 .8],'yticklabel',[.5 .6 .7 .8])
% set(gca,'xtick',[-16.6:0.2:-15.8],'xticklabel',[-16.6:0.2:-15.8])

set(gcf,'paperorientation','landscape')
set(gcf,'PaperPosition', [-2 10 30 10])

print(gcf,'-dpdf',sprintf('~/pconn_all/plots/pconn_src_dfapow_f%d_v%d.pdf',ifoi,v))

