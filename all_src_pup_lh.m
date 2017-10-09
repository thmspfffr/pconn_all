%% COMPUTE WHOLE BRAIN DFA IN SOURCE SPACE
% all_src_pup_lh

clear all

% --------------------------------------------------------
% VERSION 1
% --------------------------------------------------------
v             = 1;
v_rawdata     = 6;
fsample       = 400;
SUBJLIST      = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24];
gridsize      = 'cortex';
filt          = 'eloreta';
para.detrend  = 1;
v_pup         = 10;
calc_interv(2) = 4000;
overlap = .5;
% --------------------------------------------------------

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
for m = 1 : 3
  for isubj = SUBJLIST
    
    if ~exist(sprintf([outdir 'pconn_src_pup_lh_s%d_m%d_v%d_processing.txt'],isubj,m,v))
      system(['touch ' outdir sprintf('pconn_src_pup_lh_s%d_m%d_v%d_processing.txt',isubj,m,v)]);
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
    disp(sprintf('Processing s%d m%d ...', isubj,m))
    
    
    for iblock = bl : num_blocks(isubj,m)
      
      disp(sprintf('Loading MEG data ...'));
      
      load(sprintf('~/pconn/proc/pup/pconn_postproc_pupdfa_s%d_b%d_m%d_v%d.mat',isubj,iblock,m,v_pup))
      
      % -----------------------------------------------------------------
      % INTERPOLATE BLINKS
      % -----------------------------------------------------------------
      x = abs(diff(zscore(pup.dil)));
      [~,idx]=findpeaks(double(x>0.20),'MinPeakDistance',200);
      
      for iidx = 1 : size(idx,1)
        if idx(iidx)-50>0 && idx(iidx)+600<length(pup.dil)
          pup.dil(idx(iidx)-50:idx(iidx)+600)=NaN;
        elseif idx(iidx)-50<0
          pup.dil(2:idx(iidx)+600)=NaN;
        elseif idx(iidx)+600>length(pup.dil)
          pup.dil(idx(iidx)-50:end)=NaN;
        end
      end
      
      clear idx x
      
      pup.dil = fixgaps(pup.dil,'pchip');
      
      if isnan(pup.dil(1))
        tt = find(~isnan(pup.dil),1,'first');
        pup.dil(1) = pup.dil(tt);
      elseif isnan(pup.dil(end))
        tt = find(~isnan(pup.dil),1,'last');
        pup.dil(end) = pup.dil(tt);
      end
      pup.dil = fixgaps(pup.dil,'pchip');
      % -----------------------------------------------------------------
      clear tt
      
      % ------------------------------------------
      % COMPUTE DFA IN SOURCE SPACE
      % ------------------------------------------
      
      clear L grid
      % load cross spectrum
      try
        load(['~/pconn/proc/src/' sprintf('pconn_sens_cs_s%d_m%d_b%d_f%d_v%d.mat',isubj,m,iblock,freq,1)]);
      catch me
        segleng    = 400;
        segshift   = 200;
        epleng     = size(dat,2);
        maxfreqbin = 50;
        cs         = data2cs_event(dat',segleng,segshift,epleng,maxfreqbin);
      end
      
      if para.detrend
        flp     = 0.01;           % lowpass frequency of filter
        fhi     = 5;           % highpass
        delt  	= 1/400;            % sampling interval
        k      	= 2;                  % 2nd order butterworth filter
        fnq    	= 1/(2*delt);       % Nyquist frequency
        Wn    	= [flp/fnq fhi/fnq]; % butterworth bandpass non-dimensional frequency
        [b,a]   = butter(k,Wn);
        pup.dil = filtfilt(b,a,pup.dil);
      end
      
      if strcmp(filt,'lcmv')
        if strcmp(gridsize,'xcoarse')
          load([indir sprintf('pconn_sa_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,2)]);
          [~,A1] = mkfilt_lcmv(sa.L_xcoarse,nanmean(cs(:,:,2:50),3));
          A1 = getdipdir(nanmean(cs(:,:,2:50),3),A1);
        elseif strcmp(gridsize,'coarse')
          load([indir sprintf('pconn_sa_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,1)]);
          [~,A1] = mkfilt_lcmv(sa.L_coarse,nanmean(cs(:,:,2:50),3));
          A1 = getdipdir(nanmean(cs(:,:,2:50),3),A1);
        elseif strcmp(gridsize,'cortex')
          load([indir sprintf('pconn_sa_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,3)]);
          [~,A1] = mkfilt_lcmv(sa.L_coarse,nanmean(cs(:,:,2:50),3));
          A1 = getdipdir(nanmean(cs(:,:,2:50),3),A1);
        end
        
      elseif strcmp(filt,'eloreta')
        if strcmp(gridsize,'xcoarse')
          load([indir sprintf('pconn_sa_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,2)]);
          A1 = mkfilt_eloreta_v2(sa.L_xcoarse);
          A1 = getdipdir(nanmean(cs(:,:,2:50),3),A1);
        elseif strcmp(gridsize,'coarse')
          load([indir sprintf('pconn_sa_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,1)]);
          A1 = mkfilt_eloreta_v2(sa.L_coarse);
          A1 = getdipdir(nanmean(cs(:,:,2:50),3),A1);
        elseif strcmp(gridsize,'cortex')
          load([indir sprintf('pconn_sa_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,3)]);
          % watch out! L_coarse is correct even for cortex grid
          A1 = mkfilt_eloreta_v2(sa.L_coarse);
          A1 = getdipdir(nanmean(cs(:,:,2:50),3),A1);
        end
      end
      
      clear cs sa
      
      dat = dat';
      
      flp     = 2;           % lowpass frequency of filter
      fhi     = 5;           % highpass
      delt  	= 1/400;            % sampling interval
      k      	= 4;                  % 2nd order butterworth filter
      fnq    	= 1/(2*delt);       % Nyquist frequency
      Wn    	= [flp/fnq fhi/fnq]; % butterworth bandpass non-dimensional frequency
      [b,a]   = butter(k,Wn);
      
      ampenv_lo = filtfilt(b,a,dat);
      ampenv_lo = ampenv_lo*A1;
      
      flp     = 20;           % lowpass frequency of filter
      fhi     = 40;           % highpass
      Wn    	= [flp/fnq fhi/fnq]; % butterworth bandpass non-dimensional frequency
      [b,a]   = butter(k,Wn);
      
      ampenv_hi = filtfilt(b,a,dat);
      ampenv_hi = ampenv_hi*A1;
      
      clear dat
      % project bp-filtered signal into source space
      
      eplen    = size(ampenv_lo,1);
      seglen   = floor(calc_interv(2));
      segshift = floor(overlap*seglen);
      nseg     = floor((eplen-seglen)/segshift+1);
      par.lh = zeros(3000,nseg);
      for itime = 1 : nseg
        itime
        timerange        = [(itime-1)*segshift+1 (itime-1)*segshift+seglen];
        pup_dat          = pup.dil(timerange(1):timerange(2));
        tmp_lo          = abs(hilbert(ampenv_lo(timerange(1):timerange(2),:)));
        tmp_hi          = abs(hilbert(ampenv_hi(timerange(1):timerange(2),:)));
        
        par.lh(:,itime)          = mean(tmp_lo)./mean(tmp_hi);
        par.pup(itime)         = mean(pup_dat);
   
      end
      
      clear mydata dat
      
      save(sprintf([outdir 'pconn_src_pup_lh_s%d_b%d_m%d_v%d.mat'],isubj,iblock,m,v),'par','-v7.3');
      
      clear par
      
    end
    
  end
end

error('STOP')

%%
cnt = 0;
for m = 1 : 3
  for isubj = SUBJLIST
    
    d = dir(sprintf([outdir 'pconn_src_pup_lh_s%d_b*_m%d_v%d.mat'],isubj,m,v));
    
    for iblock = 1 : length(d)
      
      cnt = cnt + 1;
      load([outdir d(iblock).name])
      
      for i = 1 : 3000
      t=corrcoef(par.lh(i,:),par.pup);
      
      r(i,cnt) = t(1,2);
      end
    end
  end
end
%%
    load sa_meg_template;
    
    tp_plot_surface(nanmean(r,2),'cortex',sa_meg_template)
    %
