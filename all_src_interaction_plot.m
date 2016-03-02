
%% COMPUTE LOCAL DFA, AMPLITUDE AND PUPIL DILATION
% pconn_sens_pup_dfa

clear
% --------------------------------------------------------
% SETTINGS
% --------------------------------------------------------
v               = 16;
v_rawdata       = 6;
is_src          = 0;
fsample         = 400;
SUBJLIST        = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24];
bin             = [50 50];

print_baseline  = 0;
% --------------------------------------------------------

addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/
addpath /home/tpfeffer/Documents/MATLAB/fieldtrip-20130925/

indir    = '~/pconn_all/proc/';
outdir   = '~/pconn_all/proc/';
plotdir  = '~/pconn_all/plots';

run ~/Documents/MATLAB/toolboxes/NBT-NBTv0.5.3-alpha/installNBT.m
%%
clear pup dfa amp par
cnt = 0;

addpath /home/tpfeffer/pconn/matlab
ord       = pconn_randomization;

s=0
for ifoi = 3
  
  for isubj = SUBJLIST
    for m = 1 : 3
      im = find(ord(isubj,:)==m);
      
      %     isubj
      d = dir(sprintf('~/pconn_all/proc/pconn_src_pup_s%d*m%d*f%d*v%d.mat',isubj,im,ifoi,v));
      s = s+length(d);
      if length(d)==0
        warning('!')
        fprintf('s%db%dm%d\n',isubj,iblock,im)
        
        continue
      end
      for iblock = 1 : length(d)
        clear par
        try
          load([outdir d(iblock).name])
        catch me
          %           warning('!')
          fprintf('Couldn''t load s%d m%d ...\n',isubj,m)
          continue
        end
        %
        idx_hi = par.pup_amp>prctile(nanmean(par.pup_amp,3),bin(1));
        idx_me = par.pup_amp>prctile(nanmean(par.pup_amp,3),bin(2))&par.pup_amp<prctile(nanmean(par.pup_amp,3),bin(1));
        idx_lo = par.pup_amp<prctile(nanmean(par.pup_amp,3),bin(2));
  
        pup_hi   = nanmean(par.pup_amp(idx_hi));
        pup_lo   = nanmean(par.pup_amp(idx_lo));
        pup_me   = nanmean(par.pup_amp(idx_me));
        
        if (pup_lo+pup_hi)/2<1000
          warning('!!!')
          fprintf('s%db%dm%d\n',isubj,iblock,im)
        end
        
        dfa_hi(:,iblock) = nanmean(par.meg_dfa(idx_hi,:),1);
        dfa_lo(:,iblock) = nanmean(par.meg_dfa(idx_lo,:),1);
        dfa_me(:,iblock) = nanmean(par.meg_dfa(idx_me,:),1);
        
        amp_hi(:,iblock) = nanmean(par.meg_amp(idx_hi,:));
        amp_lo(:,iblock) = nanmean(par.meg_amp(idx_lo,:));
        amp_me(:,iblock) = nanmean(par.meg_amp(idx_me,:));
        %
        cvar_hi(:,iblock) = nanmean(par.meg_cvar(idx_hi,:),1);
        cvar_lo(:,iblock) = nanmean(par.meg_cvar(idx_lo,:),1);
        
        %
        clear par
        
      end

      dfa_hi_all(:,m,isubj) = nanmean(dfa_hi,2);
      dfa_lo_all(:,m,isubj) = nanmean(dfa_lo,2);
      dfa_me_all(:,m,isubj) = nanmean(dfa_me,2);
      
      cvar_hi_all(:,m,isubj) = nanmean(cvar_hi,2);
      cvar_lo_all(:,m,isubj) = nanmean(cvar_lo,2);
      
      amp_hi_all(:,m,isubj) = nanmean(amp_hi,2);
      amp_lo_all(:,m,isubj) = nanmean(amp_lo,2);
      amp_me_all(:,m,isubj) = nanmean(amp_me,2);
      
      clear dfa_hi dfa_me dfa_lo pow_hi_ft pow_lo_ft pow_me_ft amp_hi amp_lo amp_me
      
    end
  end
end

%%
% end
dfa_hi_all = squeeze(nanmean(dfa_hi_all(:,:,SUBJLIST),2));
dfa_lo_all = squeeze(nanmean(dfa_lo_all(:,:,SUBJLIST),2));
dfa_me_all = squeeze(nanmean(dfa_me_all(:,:,SUBJLIST),2));

amp_hi_all = nanmean(amp_hi_all(:,:,SUBJLIST),2);
amp_lo_all = nanmean(amp_lo_all(:,:,SUBJLIST),2);
amp_me_all = nanmean(amp_me_all(:,:,SUBJLIST),2);

pup_amp_hi_all(:,isubj) = nanmean(nanmean(pup_hi));
pup_amp_lo_all(:,isubj) = nanmean(nanmean(pup_lo));
pup_amp_me_all(:,isubj) = nanmean(nanmean(pup_me));

cvar_hi_all = double(squeeze(nanmean(cvar_hi_all(:,:,SUBJLIST),2)));
cvar_lo_all = double(squeeze(nanmean(cvar_lo_all(:,:,SUBJLIST),2)));

%% COMPARE RESULTS TO PHARMA  
str = 'dfa'; 
v_res = 1;
ord   = pconn_randomization;

  switch str
    case 'amp'
      par_hi = squeeze(amp_hi_all);
      par_lo = squeeze(amp_lo_all);
    case 'dfa'
      par_hi = squeeze(dfa_hi_all);
      par_lo = squeeze(dfa_lo_all);
    case 'cvar'
      par_hi = squeeze(cvar_hi_all);
      par_lo = squeeze(cvar_lo_all);
  end

for ifoi = 3 : 3
  for isubj = SUBJLIST
    for m = 1 : 3

      im = find(ord(isubj,:)==m);

      load(sprintf(['~/pconn/proc/dfa/pconn_src_dfa_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,v_res));
     
      dfa_all_res(:,m,isubj,ifoi)  = nanmean(par.dfa,2);
      var_all_res(:,m,isubj,ifoi)  = nanmean(par.var,2);
      cvar_all_res(:,m,isubj,ifoi) = nanmean(par.cvar,2); clear par


    end
  end
end

dfa_all_res  = dfa_all_res(:,:,SUBJLIST,ifoi);
var_all_res  = var_all_res(:,:,SUBJLIST,ifoi);
cvar_all_res = cvar_all_res(:,:,SUBJLIST,ifoi);
%%

for isubj = 1 :  size(dfa_all_res,3)
  
  d(:,1) = squeeze(dfa_all_res(:,2,isubj)-dfa_all_res(:,1,isubj));
  d(:,2) = squeeze(cvar_all_res(:,2,isubj)-cvar_all_res(:,1,isubj));
  d(:,3) = squeeze(var_all_res(:,2,isubj)-var_all_res(:,1,isubj));
  
  c_map(1,1,isubj) = corr(d(:,1),d(:,1));
  c_map(1,2,isubj) = corr(d(:,1),d(:,2));
  c_map(1,3,isubj) = corr(d(:,1),d(:,3));
 
end
  

d_res = squeeze(dfa_all_res(:,2,:)-dfa_all_res(:,1,:));
d_pup = squeeze(par_hi-par_lo);

for isubj = 1 : size(par_hi,2)
  
  r(isubj) = corr(d_res(:,isubj),d_pup(:,isubj));
  
end
  
