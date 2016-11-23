%% COMPUTE WHOLE BRAIN DFA IN SOURCE SPACE
% all_src_dfa

clear all

% --------------------------------------------------------
% VERSION 1
% --------------------------------------------------------
v           = 1;
fsample     = 400;
SUBJLIST    = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33];
foi         = [2 4; 4 8; 8 12; 12 36; 60 120];
fit_interv  = [2 50];
calc_interv = [1.5 70];
gridsize    = 'cortex';
filt        = 'eloreta';
overlap     = 0.3;
para.detrend = 1;
dfa_overlap = 0.5;
v_pup = 1;
% --------------------------------------------------------
% VERSION 7
% --------------------------------------------------------
% v           = 7;
% fsample     = 400;
% SUBJLIST    = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33];
% foi         = [2 4; 4 8; 8 12; 12 36];
% fit_interv  = [2 30];
% calc_interv = [1.5 50];
% gridsize    = 'cortex';
% filt        = 'eloreta';
% dfa_overlap = 0.5;
% filt_ord    = 2;
% v_pup       = 1;
% para.detrend = 1;
% --------------------------------------------------------
% VERSION 11
% --------------------------------------------------------
% v           = 11;
% fsample     = 400;
% SUBJLIST    = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33];
% foi         = [2 4; 4 8; 8 12; 12 36];
% fit_interv  = [2 50];
% calc_interv = [1.5 70];
% gridsize    = 'cortex';
% filt        = 'eloreta';
% overlap     = 0.3;
% para.detrend = 0;
% dfa_overlap = 0.5;
% v_pup = 1;
% --------------------------------------------------------

% para.detrend = 1
% SUBJECT 20, m == 2, iblock == 1: no first trigger, perhaps mismatch
% between time courses. Consider excluding 1st block.

% SUBJECT 10, m == 3, iblock == 1: strange trigger data, excluded for now

% SUBJECT 32, m ==3, iblock == 1: No pupil data


restoredefaultpath

addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/
addpath /home/tpfeffer/Documents/MATLAB/fieldtrip-20130925/
addpath /home/tpfeffer/Documents/MATLAB/toolboxes/NBT-ReleaseNBTRC4a/
addpath ~/pconn/matlab

ft_defaults

indir   = '/home/tpfeffer/pconn/proc/src/';
outdir   = '/home/tpfeffer/pconn_all/proc/';
freq = 1;
run ~/Documents/MATLAB/toolboxes/NBT-NBTv0.5.3-alpha/installNBT.m
%%

for ifoi = 1 : length(foi)
  for m = 1 : 3
    for isubj = SUBJLIST
      
      if ~exist(sprintf([outdir 'pconn_src_pup_s%d_m%d_f%d_v%d_processing.txt'],isubj,m,ifoi,v))
        system(['touch ' outdir sprintf('pconn_src_pup_s%d_m%d_f%d_v%d_processing.txt',isubj,m,ifoi,v)]);
      else
        continue
      end
      %
      % ---------------------------------------------------
      % READ IN PUPIL
      % ---------------------------------------------------
      
      d=dir(['/home/tpfeffer/pconn/proc/pup/' sprintf('pconn_postproc_pupdfa_s%d_b*_m%d_v%d.mat',isubj,m,v_pup)]);
      
      if length(d)<1
        num_blocks(isubj,m) = 0;
        continue
      elseif length(d)==1 && str2double(d(1).name(end-10)) == 1
        num_blocks(isubj,m) = 1;
        bl = 1;
        warning(sprintf('only one block which is b%d',bl));
      elseif  length(d)==1 && str2double(d(1).name(end-10)) == 2
        num_blocks(isubj,m) = 2;
        bl = 2;
      else
        num_blocks(isubj,m) = 2;
        bl = 1;
      end
      
      disp(sprintf('Processing s%d m%d f%d ...', isubj,m,ifoi))
      
      
      for iblock = bl : num_blocks(isubj,m)
        
        disp(sprintf('Loading MEG data ...'));

        load(sprintf('~/pconn/proc/pup/pconn_postproc_pupdfa_s%d_b%d_m%d_v%d.mat',isubj,iblock,m,v_pup))

        % -----------------------------------------------------------------
        % INTERPOLATE BLINKS
        % -----------------------------------------------------------------
        x = abs(diff(zscore(pup)));
        [~,idx]=findpeaks(double(x>0.20),'MinPeakDistance',200);
        
        for iidx = 1 : size(idx,1)
          if idx(iidx)-50>0 && idx(iidx)+600<length(pup)
            pup(idx(iidx)-50:idx(iidx)+600)=NaN;
          elseif idx(iidx)-50<0
            pup(2:idx(iidx)+600)=NaN;
          elseif idx(iidx)+600>length(pup)
            pup(idx(iidx)-50:end)=NaN;
          end
        end
        
        clear idx x
        
        pup = fixgaps(pup,'pchip');
        
        if isnan(pup(1))
          tt = find(~isnan(pup),1,'first');
          pup(1) = pup(tt);
        elseif isnan(pup(end))
          tt = find(~isnan(pup),1,'last');
          pup(end) = pup(tt);
        end
        pup = fixgaps(pup,'pchip');
        % -----------------------------------------------------------------
        clear tt
        figure; set (gcf,'color','white')
subplot(1,2,1)
        plot(pup); title(sprintf('S%dM%dB%d',isubj,m,iblock))
        tp_editplots
%         print(gcf,'-djpeg100',sprintf('~/pconn_all/plots/all_src_pupil_timecourse_procproc_s%d_m%d_b%d_v%d.jpg',isubj,m,iblock,v))
                 
         % ------------------------------------------
        % COMPUTE DFA IN SOURCE SPACE
        % ------------------------------------------
        
        clear L grid
        % load cross spectrum
        load(['~/pconn/proc/src/' sprintf('pconn_sens_cs_s%d_m%d_b%d_f1_v%d.mat',isubj,m,iblock,3)]);
        
        if para.detrend 
          flp     = 0.01;           % lowpass frequency of filter
          fhi     = 5;           % highpass
          delt  	= 1/400;            % sampling interval
          k      	= 2;                  % 2nd order butterworth filter
          fnq    	= 1/(2*delt);       % Nyquist frequency
          Wn    	= [flp/fnq fhi/fnq]; % butterworth bandpass non-dimensional frequency
          [b,a]   = butter(k,Wn);  
          pup = filtfilt(b,a,pup);
          subplot(1,2,2)
          plot(pup); title(sprintf('S%dM%dB%d',isubj,m,iblock));tp_editplots
   
        end
        
      	print(gcf,'-djpeg100',sprintf('~/pconn_all/plots/all_src_pupil_timecourse_procproc_filt_s%d_m%d_b%d_v%d.jpg',isubj,m,iblock,v))

        
        if strcmp(filt,'lcmv')
          if strcmp(gridsize,'xcoarse')
            load([indir sprintf('pconn_sa_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,2)]);
            [~,A1] = mkfilt_lcmv(sa.L_xcoarse,nanmean(cs(:,:,foi(ifoi,1):foi(ifoi,2)),3));
          elseif strcmp(gridsize,'coarse')
            load([indir sprintf('pconn_sa_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,1)]);
            [~,A1] = mkfilt_lcmv(sa.L_coarse,nanmean(cs(:,:,foi(ifoi,1):foi(ifoi,2)),3));
          elseif strcmp(gridsize,'cortex')
            load([indir sprintf('pconn_sa_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,3)]);
            [~,A1] = mkfilt_lcmv(sa.L_coarse,nanmean(cs(:,:,foi(ifoi,1):foi(ifoi,2)),3));
            A1 = getdipdir(nanmean(cs(:,:,foi(ifoi,1):foi(ifoi,2)),3),A1);
          end
          
        elseif strcmp(filt,'eloreta')
          if strcmp(gridsize,'xcoarse')
            load([indir sprintf('pconn_sa_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,2)]);
            A1 = mkfilt_eloreta_v2(sa.L_xcoarse);
            A1 = getdipdir(nanmean(cs(:,:,foi(ifoi,1):foi(ifoi,2)),3),A1);
          elseif strcmp(gridsize,'coarse')
            load([indir sprintf('pconn_sa_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,1)]);
            A1 = mkfilt_eloreta_v2(sa.L_coarse);
            A1 = getdipdir(nanmean(cs(:,:,foi(ifoi,1):foi(ifoi,2)),3),A1);
          elseif strcmp(gridsize,'cortex')
            load([indir sprintf('pconn_sa_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,3)]);
            % watch out! L_coarse is correct even for cortex grid
            A1 = mkfilt_eloreta_v2(sa.L_coarse);
            A1 = getdipdir(nanmean(cs(:,:,foi(ifoi,1):foi(ifoi,2)),3),A1);
          end
        end
        
        clear cs sa
        
        siginfo = nbt_Info;
        siginfo.converted_sample_frequency = fsample;
        dat = dat';
        
        % compute bp-filtered signal
        ampenv = nbt_filter_fir(dat,foi(ifoi,1),foi(ifoi,2),siginfo.converted_sample_frequency,2/foi(ifoi,1));
        
        clear data_low data_hi dat
        % project bp-filtered signal into source space
        ampenv=ampenv*A1;
        
        eplen    = size(ampenv,1);
        seglen   = floor(calc_interv(2)*fsample);
        segshift = floor(overlap*seglen);
        nseg     = floor((eplen-seglen)/segshift+1);
        
        for itime = 1 : nseg
          
          timerange        = [(itime-1)*segshift+1 (itime-1)*segshift+seglen];
          pup_dat          = pup(timerange(1):timerange(2));
          tmp_dat          = abs(hilbert(ampenv(timerange(1):timerange(2),:)));
          
          pha              = angle(hilbert(ampenv(timerange(1):timerange(2),:)));
          r                = abs(sum(exp(i*pha),2)/size(pha,2)); clear pha
          r_std            = std(r);
          r_mean           = mean(r);
          
          fprintf('Computing DFA for seg%d ...\n',itime);
          
          tmp              = nbt_doDFA(tmp_dat, siginfo, fit_interv,calc_interv,dfa_overlap,0,0,[]);
          r_dfa            = nbt_doDFA(r, siginfo, fit_interv,calc_interv,dfa_overlap,0,0,[]);
          r_dfa            = r_dfa.MarkerValues;
          
%           par.meg_pow(itime,:)  = nanmean(p,2);
          par.meg_dfa(itime,:)  = tmp.MarkerValues;
          par.meg_amp(itime,:)  = nanmean(tmp_dat,1);
          par.pup_amp(itime,:)  = nanmean(pup_dat);
          par.pup_diff(itime,:) = diff(pup_dat);
          par.pup_var(itime,:)  = var(pup_dat);
          par.meg_var(itime,:)  = nanvar(tmp_dat);
          par.meg_cvar(itime,:) = nanstd(tmp_dat)./nanmean(tmp_dat);
          par.r_mean(itime,:)   = r_mean; 
          par.r_std(itime,:)    = r_std;
          par.r_dfa(itime,:)    = r_dfa;
          
          clear tmp pup_dat pup_dfa tmp_dat r_dfa r_mean r_std r 
          
        end
        
        clear mydata dat
        
        save(sprintf([outdir 'pconn_src_pup_s%d_b%d_m%d_f%d_v%d.mat'],isubj,iblock,m,ifoi,v),'par','-v7.3');
        
        clear par
        
      end
      
    end
  end
end

error('STOP')

