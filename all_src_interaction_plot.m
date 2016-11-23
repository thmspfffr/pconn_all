
%% COMPUTE LOCAL DFA, AMPLITUDE AND PUPIL DILATION
% pconn_sens_pup_dfa

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


%% PLOT DFA
   
str = 'dfa'; 
ord   = pconn_randomization;

for ifoi = 1 : 4
  for isubj = SUBJLIST
    for m = 1 : 3

      im = find(ord(isubj,:)==m);
try
      load(sprintf([outdir 'pconn_cnt_src_dfa_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,v));
     
      dfa_all_cnt(:,m,isubj,ifoi)  = nanmean(par.dfa,2);
      var_all_cnt(:,m,isubj,ifoi)  = nanmean(par.var,2);
     	cvar_all_cnt(:,m,isubj,ifoi) = nanmean(par.cvar,2); 
      amp_all_cnt(:,m,isubj,ifoi)  = nanmean(par.amp,2); clear par
      
     
      
      load(sprintf(['~/pconn/proc/dfa/pconn_src_dfa_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,v));
     
      dfa_all_res(:,m,isubj,ifoi)  = nanmean(par.dfa,2);
      var_all_res(:,m,isubj,ifoi)  = nanmean(par.var,2);
      cvar_all_res(:,m,isubj,ifoi) = nanmean(par.cvar,2);
      amp_all_res(:,m,isubj,ifoi)  = nanmean(par.amp,2); clear par
      catch me
        disp(sprintf('err s%d f%d m%d',isubj,ifoi,im))
      end

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

%%
ifoi = 3

clear d

addpath ~/Documents/MATLAB/cbrewer/cbrewer/

cmap = cbrewer('seq', 'OrRd', 100,'pchip');% 

for isubj = 1 :  size(dfa_all_res,3)
  
  d(:,1) = squeeze(dfa_all_res(:,2,isubj,ifoi)-dfa_all_res(:,1,isubj,ifoi));
  d(:,2) = squeeze(cvar_all_res(:,2,isubj,ifoi)-cvar_all_res(:,1,isubj,ifoi));
  d(:,3) = squeeze(var_all_res(:,2,isubj,ifoi)-var_all_res(:,1,isubj,ifoi));
  d(:,4) = squeeze(amp_all_res(:,2,isubj,ifoi)-amp_all_res(:,1,isubj,ifoi));
  
  for i = 1 : 4
    for j = 1 : 4
      c_map(i,j,isubj) = corr(d(:,i),d(:,j));
    end
  end
  
end

h = figure; set(gcf,'color','white'); hold on

imagesc(nanmean(triu(nanmean(c_map,3)),3)); colormap(cmap)

% COMPUTE P-VALUES

for i = 1 : 4
  for j = 1 : 4
    [t(i,j) p(i,j)] = ttest(squeeze(c_map(i,j,:)),[],'alpha',0.05/6);
    if ~isnan(t(i,j)) && t(i,j) && i~=j
      if p(i,j) <.0001
%         plot(i,j,'*','markersize',10)
        text(i,j,'***')
      elseif p(i,j) < .001
        text(i,j,'**')
      elseif p(i,j) < .01
        text(i,j,'*')
      end
       
    end
  end
end

title(sprintf('Interaction: f%d',ifoi));

set(gca,'TickDir','out','XTick',[1 2 3 4],'XTickLabel',{'DFA';'CVAR';'VAR';'AMP'});
set(gca,'TickDir','out','YTick',[1 2 3 4],'YTickLabel',{'DFA';'CVAR';'VAR';'AMP'});

print(gcf,'-djpeg100',sprintf('~/pconn_all/plots/pconn_src_interaction_rst_f%d_v%d.jpg',ifoi,v))
 
%% TASK 
clear c_map clear d

for isubj = 1 :  size(dfa_all_res,3)
  
  d(:,1) = squeeze(dfa_all_cnt(:,2,isubj,ifoi)-dfa_all_cnt(:,1,isubj,ifoi));
  d(:,2) = squeeze(cvar_all_cnt(:,2,isubj,ifoi)-cvar_all_cnt(:,1,isubj,ifoi));
  d(:,3) = squeeze(var_all_cnt(:,2,isubj,ifoi)-var_all_cnt(:,1,isubj,ifoi));
  d(:,4) = squeeze(amp_all_cnt(:,2,isubj,ifoi)-amp_all_cnt(:,1,isubj,ifoi));
  
  for i = 1 : 4
    for j = 1 : 4
      c_map(i,j,isubj) = corr(d(:,i),d(:,j));
    end
  end
  
end

h = figure; set(gcf,'color','white'); hold on

imagesc(nanmean(triu(nanmean(c_map,3)),3)); colormap(cmap)

% COMPUTE P-VALUES

for i = 1 : 4
  for j = 1 : 4
    [t(i,j) p(i,j)] = ttest(squeeze(c_map(i,j,:)),[],'alpha',0.05/6);
     if p(i,j) <.0001
%         plot(i,j,'*','markersize',10)
        text(i,j,'***')
      elseif p(i,j) < .001
        text(i,j,'**')
      elseif p(i,j) < .01
        text(i,j,'*')
      end
  end
end

title(sprintf('Interaction: f%d',ifoi));

set(gca,'TickDir','out','XTick',[1 2 3 4],'XTickLabel',{'DFA';'CVAR';'VAR';'AMP'});
set(gca,'TickDir','out','YTick',[1 2 3 4],'YTickLabel',{'DFA';'CVAR';'VAR';'AMP'});

print(gcf,'-djpeg100',sprintf('~/pconn_all/plots/pconn_src_interaction_tsk_f%d_v%d.jpg',ifoi,v))
 
%%
d_res = squeeze(dfa_all_res(:,2,:)-dfa_all_res(:,1,:));
d_pup = squeeze(par_hi-par_lo);

for isubj = 1 : size(par_hi,2)
  
  r(isubj) = corr(d_res(:,isubj),d_pup(:,isubj));
  
end

%% COMPARE MEASURES
  
