%% PLOT SENSOR-LEVEL DFA RESULTS
% all_sens_amppowfreq

clear

v = 1;
SUBJLIST    = [3 4 5 6 7 8 9 10 11 12 13 15 16 17 19 20 21 22 23 24 25 26 27 28 29 30 31 32];
foi         = [2:2:178; 4:2:180]';
filt_ord    = 2;
v_rawdata  =6 ;
addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/

outdir   = '/home/tpfeffer/pconn_all/proc/';

addpath /home/tpfeffer/pconn/matlab
addpath ~/pcbi/

run ~/Documents/MATLAB/toolboxes/NBT-NBTv0.5.3-alpha/installNBT.m

%%
for m = 1 : 3
  for isubj = SUBJLIST
    for ifoi = 1 : length(foi)
          
      % %
      if isubj == 3
        load ~/pconn/matlab/pconn_sensorlabels.mat
      end
      
      for iblock = 1 : 2
        
        if ~exist(sprintf([outdir 'all_sens_amppowfreq_s%d_m%d_b%d_f%d_v%d_processing.txt'],isubj,m,iblock,ifoi,v))
          system(['touch ' outdir sprintf('all_sens_amppowfreq_s%d_m%d_b%d_f%d_v%d_processing.txt',isubj,m,iblock,ifoi,v)]);
        else
          continue
        end
        
        disp(sprintf('Processing MEG s%dm%df%d%b...',isubj,m,ifoi,iblock));
        
        if (isubj~=3 && isubj~=2 &&  isubj~=17) && isubj <= 24
          load(sprintf('/home/tpfeffer/pconn/proc/preproc/pconn_postproc_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,v_rawdata));
        else
          load(sprintf('/home/tpfeffer/pconn/proc/preproc/pconn_postproc_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,3));
        end
        
        if isubj == 3 || isubj == 2
          [~,idx_lab]=intersect(data_low.label,lab);
          data_low.trial{1}=data_low.trial{1}(idx_lab,:);
          data_low.label = lab;
          save(['~/pconn/matlab/' sprintf('pconn_idxlab.mat')],'idx_lab','-v7.3');
        end
        
        clear data_hi
        
        [mydata,epleng] = megdata2mydata(data_low);
        
        mydata = mydata(1:end-1000,:);
        
        siginfo = nbt_Info;
        siginfo.converted_sample_frequency = 400;
        
        % compute bp-filtered signal
        tmp    = single(nbt_filter_fir(mydata,foi(ifoi,1),foi(ifoi,2),siginfo.converted_sample_frequency,filt_ord/foi(ifoi,1)));
        ampenv = abs(hilbert(tmp)); clear tmp mydata
        
        for ichan = 1 : 268
          fprintf('Computing power chan%d ...\n',ichan)
          [p(ichan,:),f]=pwelch(double(ampenv(:,ichan)),[],50,logspace(-3,1,50),400);
        end
        
        save([outdir sprintf('all_sens_amppowfreq_s%d_m%d_b%d_f%d_v%d.mat',isubj,m,iblock,ifoi,v)],'p','-v7.3')
        
      end
    end
  end
end

error('!');

%% PLOT

for m = 1 : 3
  for isubj = SUBJLIST
    for ifoi = 1 : length(foi)
      for iblock = 1 : 2
        
        load([outdir sprintf('all_sens_amppowfreq_s%d_m%d_b%d_f%d_v%d.mat',isubj,m,iblock,ifoi,v)])
        allfreq(:,:,m,ifoi,isubj,iblock) = p;
        
      end
    end
  end
end

allfreq = nanmean(allfreq(:,:,:,:,SUBJLIST,:),6);
