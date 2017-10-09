%% ALL_SENS_GIM
% Computes time-resolved GIM and saves mean and std across freqs


% --------------------------------------------------------
% VERSION 1
% --------------------------------------------------------
v         = 1;
v_rawdata = 6;
fsample   = 400;
SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24];
% --------------------------------------------------------


addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/
addpath /home/tpfeffer/Documents/MATLAB/fieldtrip-20130925/

outdir   = '/home/tpfeffer/pconn_all/proc/';


%%

for m = 1 : 3
  for isubj = SUBJLIST
      
      if ~exist(sprintf([outdir 'all_sens_gim_s%d_m%d_v%d_processing.txt'],isubj,m,v))
        system(['touch ' outdir sprintf('all_sens_gim_s%d_m%d_v%d_processing.txt',isubj,m,v)]);
      else
        continue
      end
%       
      if isubj == 3
        load ~/pconn/matlab/pconn_sensorlabels.mat
      end
      
      for iblock = 1 : 2
        
        disp(sprintf('Processing MEG s%dm%df%d%b...',isubj,m,iblock));
        
        load(sprintf('/home/tpfeffer/pconn/proc/preproc/pconn_postproc_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,v_rawdata));
                                
        clear data_hi
        
        [mydata,~] = megdata2mydata(data_low);

        fs=400;

        mydata = mydata';
        % # channels, # samples
        [nch,ns]=size(mydata);

        % window length:
        wlength = 5*fs;

        nwin = floor(ns/wlength);

        segleng       = fs;
        segshift      = segleng/2;

        epleng        = wlength;
        df            = fs/segleng;
        maxfreqbin    = 40;
        para.segave   = 1;
        para.subave   = 1;
        para.zeropad  = 0;
        xf            = [0:maxfreqbin-1]*df;
        p             = 10;
        cstr          = zeros(nch,nch,length(xf),nwin);
        g             = zeros(length(xf),nwin);

        for iwin = 1 : nwin

            fprintf('Processing window %d/%d ...\n',iwin,nwin)
            dr=mydata(:,(iwin-1)*wlength+1:iwin*wlength);
            [cstr(:,:,:,iwin),~,~]=data2cs_event(dr',segleng,segshift,epleng,maxfreqbin,para);
            g(:,iwin)=proc_gim(squeeze(cstr(:,:,:,iwin)),p);    

        end
        
        gim.mean = mean(g,2);
        gim.freq = xf;
        gim.std  = std(g,[],2);
        
        save(sprintf(['~/pconn_all/proc/all_sens_gim_s%d_b%d_m%d_v%d.mat'],isubj,iblock,m,v),'gim');
        
      end
  end
end
        
%%

ord = pconn_randomization;

for m = 1 : 3
  for isubj = SUBJLIST
      for iblock = 1 : 2
        
        im = find(ord(isubj,:)==m);
        
       	load(sprintf('~/pconn_all/proc/all_sens_gim_s%d_b%d_m%d_v%d.mat',isubj,iblock,im,v));
      
        me(:,isubj,iblock,m) = gim.mean;
        st(:,isubj,iblock,m) = gim.std;
        
      end
  end
end

me = squeeze(nanmean(me(:,SUBJLIST,:,:),3));
st = squeeze(nanmean(st(:,SUBJLIST,:,:),3));

d=[mean(me(:,:,2))-mean(me(:,:,1))]








% f=10;
% figure;
% plot(g(:,:));
% set(gca, 'XLabel', 1:length(xf));
% set(gca, 'XTickLabel', xf);
% set(gca, 'YTick', 1:length(xf));
% set(gca, 'YTickLabel', xf);
% xlabel('Frequency', oplabel{:});
% set(gca,opgca{:});



