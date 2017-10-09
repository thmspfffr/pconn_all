%% COMPUTE WHOLE BRAIN DFA IN SOURCE SPACE
% all_src_pow_ROI

% requires execution of pconn_src_compute_common_filter.m first.
% this is where the common dipole orientation is determined.

clear

v     = 1;
FOI   = [2 4; 4 8; 8 12; 12 36];
smoo  = [1; 2; 4; 12];

SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 17 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];

restoredefaultpath

addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/
addpath('/home/tpfeffer/Documents/MATLAB/fieldtrip-20160919/')

ft_defaults

indir   = '/home/tpfeffer/pconn/proc/src/';
outdir   = '/home/tpfeffer/pconn/proc/src/';

%%

  for m = 1 : 3
    for isubj = SUBJLIST
      
      if ~exist(sprintf(['~/pconn_all/proc/pow/' 'all_src_pow_ROI_s%d_m%d_v%d_processing.txt'],isubj,m,v))
        system(['touch ' '~/pconn_all/proc/pow/' sprintf('all_src_pow_ROI_s%d_m%d_v%d_processing.txt',isubj,m,v)]);
      else
        continue
      end
      
      clear d
      
      disp(sprintf('Processing s%d m%d ...', isubj,m))
      
      d=dir(sprintf('/home/tpfeffer/pconn/proc/preproc/pconn_postpostproc_s%d_m%d_b*_v%d.mat',isubj,m,1));
      
      if length(d)==1
        if v < 10
          blocks = str2num(d(1).name(end-7));
        else
          blocks = str2num(d(1).name(end-8));
        end
      elseif length(d) == 2
        blocks = [1 2];
      end
            
      for iblock = blocks
        
        if length(blocks) == 1
          if iblock == 1
            par.pow(:,2) = nan(3000,1); 
          else
            par.pow(:,1) = nan(3000,1); 
          end
        end
        
        disp(sprintf('Processing MEG s%dm%d%b...',isubj,m,iblock));
        
        if length(d) == 2
          load(['~/pconn/proc/preproc/' d(iblock).name]);
        else
          load(['~/pconn/proc/preproc/' d(1).name]);
        end
        
        cfg = [];
        cfg.length =  1/0.2;
        cfg.overlap = 0;
        
        data      = ft_redefinetrial(cfg,data);
        
        cfg               = [];
        cfg.method        = 'mtmfft';
        cfg.output        = 'powandcsd';
        cfg.taper         = 'dpss';
        cfg.pad           = 'nextpow2';
        cfg.foi           = 2.^(1:0.125:7);
        cfg.keeptrials    = 'no';
        cfg.tapsmofrq     = 2;
        
        [~,  csd]         = ft_freqanalysis(cfg, data); clear data
        
        % watch out! L_coarse is correct even for cortex grid
        load([indir sprintf('pconn_sa_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,3)]);
        A = mkfilt_eloreta_v2(sa.L_coarse);
        
%         load(sprintf([outdir 'pconn_src_common_filter_s%d_f%d_v%d.mat'],isubj,ifoi,1));
  
          roi(1,:,1) = [11 -78 9]./10; % visual1
          roi(2,:,1) = [-11 -81 7]./10;
          roi(1,:,2) = [29 -92 2]./10; %visual2
          roi(2,:,2) = [-19 -92 2]./10;
          roi(1,:,3) = [50 -21 7]./10; % auditory
          roi(2,:,3) = [-52 -19 7]./10;
          roi(1,:,4) = [38 -18 45]./10; % prim motor
          roi(2,:,4) = [-36 -19 48]./10;
          roi(1,:,5) = [15 -33 48]./10; % sens assoc.
          roi(2,:,5) = [-14 -33 48]./10;
          roi(1,:,6) = [43 38 12]./10; % DLPFC (BA46)
          roi(2,:,6) = [-46 38 8]./10;
          roi(1,:,7) = [23 55 7]./10; % OFC (BA10) 
          roi(2,:,7) = [-23 55 4]./10; 
          roi(1,:,8) = [23 -60 61]./10; % precuneus (BA7)
          roi(2,:,8) = [-18 -61 55]./10;
          

        for iroi = 1 : 8
         
          [y(iroi,1),i(iroi,1)]=min(abs(sum(sa.grid_cortex3000-repmat(roi(1,:,iroi),[3000 1]),2)))
          [y(iroi,2),i(iroi,2)]=min(abs(sum(sa.grid_cortex3000-repmat(roi(2,:,iroi),[3000 1]),2)))
       
          for iff = 1 : size(csd,3)
            for ivox = 1:2

              Aloc       = squeeze(A(:,i(iroi,ivox),:));
              cs_src = Aloc'*csd(:,:,iff)*Aloc;
              [u s vv]   = svd(real(cs_src));
              A1(:,ivox) = Aloc*u(:,1);
              tmppow(ivox) =  diag(A1(:,ivox)'*real(csd(:,:,iff))*A1(:,ivox));

            end
            par.pow(iff,iroi) = mean(tmppow);
          end
        end
        
        clear csd A A1 
              save(sprintf(['~/pconn_all/proc/pow/' 'all_src_pow_ROI_s%d_b%d_m%d_v%d.mat'],isubj,iblock,m,v),'par','-v7.3');

      
      end
      
      clear par
      
    end
  end

error('STOP')

%% CLEAN NON PROCESSED FILES
ord           = pconn_randomization;
SUBJLIST = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
v = 1;

for iroi = 1 : 8
         
          [y(iroi,1),i(iroi,1)]=min(abs(sum(sa_meg_template.grid_cortex3000-repmat(roi(1,:,iroi),[3000 1]),2)))
          [y(iroi,2),i(iroi,2)]=min(abs(sum(sa_meg_template.grid_cortex3000-repmat(roi(2,:,iroi),[3000 1]),2)))
end

for m = 1 : 3
  for isubj = SUBJLIST
    im = find(ord(isubj,:)==m);
%     im = m;
    for iblock = 1 : 2
%       try
        load(sprintf(['~/pconn_all/proc/pow/' 'all_src_pow_ROI_s%d_b%d_m%d_v%d.mat'],isubj,iblock,im,v));
        p(:,:,m,isubj,iblock) = par.pow; clear par
%       catch me
%         p(:,:,m,isubj,iblock) = nan(49,8);
%       end
    end

    

      
  end
end

 
   
 
p = nanmean(p(:,:,:,SUBJLIST,:),5);

i{1} = find(strcmp(aalgrid.labels,'Precuneus_L')|strcmp(aalgrid.labels,'Precuneus_R'));
i{2} = find(strcmp(aalgrid.labels,'Cuneus_L')|strcmp(aalgrid.labels,'Cuneus_R'));
i{3} = find(strcmp(aalgrid.labels,'Frontal_Sup_Medial_L')|strcmp(aalgrid.labels,'Frontal_Sup_Medial_R'));
i{4} = find(strcmp(aalgrid.labels,'Postcentral_L')|strcmp(aalgrid.labels,'Postcentral_R'));
i{5} = find(strcmp(aalgrid.labels,'Calcarine_L')|strcmp(aalgrid.labels,'Calcarine_R'));
i{6} = find(strcmp(aalgrid.labels,'Postcentral_L')|strcmp(aalgrid.labels,'Postcentral_R'));
i{7} = find(strcmp(aalgrid.labels,'Frontal_Sup_L')|strcmp(aalgrid.labels,'Frontal_Sup_R'));
i{8} = find(strcmp(aalgrid.labels,'Frontal_Sup_Orb_L')|strcmp(aalgrid.labels,'Frontal_Sup_Orb_R'));


for ifoi = 1 : 4
  for isubj = SUBJLIST
    for m = 1 : 3
            im = find(ord(isubj,:)==m);
      load(sprintf(['~/pconn/proc/src/' 'pconn_src_pow_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,v));
      
      pow_all(:,isubj,m,ifoi) = nanmean(par.pow,2); clear par
    end
  end
end

pow_all = pow_all(:,SUBJLIST,:,:);

for iroi = 1: 8
  pow_roi(:,iroi,:,:) = nanmean((pow_all(i{iroi},:,:,:)),1);
end



%%
k=figure; set (gcf,'color','white'); hold on
str = {'Prim. vis.';'Vis. ass.';'Audit.';'Prim. mot.';'Sens. ass.';'DLPFC';'OFC';'Prec.'}
for iroi = 1 : 8
  h=subplot(2,4,iroi); hold on; get(h,'pos'); if iroi>4;set(h,'pos',[ans(1)+0 ans(2) ans(3)+0 ans(4)+0.3 ]);end
  title(sprintf('%s',str{iroi}))
  plot(log2(2.^(1:0.125:7)),nanmean(p(:,iroi,1,:),4),'color',[0.7 0.7 0.7],'linewidth',3)
  plot(log2(2.^(1:0.125:7)),nanmean(p(:,iroi,2,:),4),'color',[1 0.5 0],'linewidth',3)
  plot(log2(2.^(1:0.125:7)),nanmean(p(:,iroi,3,:),4),'color',[0 0.5 1],'linewidth',3)
  set(gca,'xtick',[1 3 5 7],'xticklabel',[2 8 32 128],'tickdir','out');
  legend('Pbo.','Atx.','Dpz.')
  axis square; 
end
set(k,'Position',[480 400 900 620]);

%% STATS
[~,p1,~,s1]=ttest(p(:,:,2,:),p(:,:,1,:),'dim',4)
[~,p2,~,s2]=ttest(p(:,:,3,:),p(:,:,1,:),'dim',4)

k=figure; set (gcf,'color','white'); hold on
str = {'Prim. vis.';'Vis. ass.';'Audit.';'Prim. mot.';'Sens. ass.';'DLPFC';'OFC';'Prec.'}
for iroi = 1 : 8
  h=subplot(2,4,iroi); hold on; get(h,'pos'); if iroi>4;set(h,'pos',[ans(1)+0 ans(2) ans(3)+0 ans(4)+0.3 ]);end
  title(sprintf('%s',str{iroi}))

  line([0.5 7.5],[-1.96 -1.96],'color','k');
  line([0.5 7.5],[1.96 1.96],'color','k')

    plot(log2(2.^(1:0.125:7)),s1.tstat(:,iroi),'color',[0.7 0.7 0.7],'linewidth',3)
  plot(log2(2.^(1:0.125:7)),s2.tstat(:,iroi),'color',[1 0.5 0],'linewidth',3)
  set(gca,'xtick',[1 3 5 7],'xticklabel',[2 8 32 128],'tickdir','out');
  if iroi == 8 ; legend('A-P','D-P'); end
  axis square;  axis([0.5 7.5 -3 3])
end
set(k,'Position',[480 400 900 620]);

%%
k=figure; set (gcf,'color','white'); hold on
str = {'Prim. vis.';'Vis. ass.';'Audit.';'Prim. mot.';'Sens. ass.';'DLPFC';'OFC';'Prec.'}
for iroi = 1 : 8
  h=subplot(2,4,iroi); hold on; get(h,'pos'); if iroi>4;set(h,'pos',[ans(1)+0 ans(2) ans(3)+0 ans(4)+0.3 ]);end
  title(sprintf('%s',str{iroi}))
  plot(1:4,squeeze(nanmean(pow_roi(:,iroi,1,:),1)),'color',[0.7 0.7 0.7],'linewidth',3)
  plot(1:4,squeeze(nanmean(pow_roi(:,iroi,2,:),1)),'color',[1 0.5 0],'linewidth',3)
  plot(1:4,squeeze(nanmean(pow_roi(:,iroi,3,:),1)),'color',[0 0.5 1],'linewidth',3)
  set(gca,'xtick',[1 2 3 4],'xticklabel',{'3';'6';'10';'24'},'tickdir','out');
  legend('Pbo.','Atx.','Dpz.')
  axis square; 
end
set(k,'Position',[480 400 900 620]);

%%
[~,p1,~,s1]=ttest(pow_roi(:,:,2,:),pow_roi(:,:,1,:),'dim',1)
[~,p2,~,s2]=ttest(pow_roi(:,:,3,:),pow_roi(:,:,1,:),'dim',1)

k=figure; set (gcf,'color','white'); hold on
str = {'Prim. vis.';'Vis. ass.';'Audit.';'Prim. mot.';'Sens. ass.';'DLPFC';'OFC';'Prec.'}
for iroi = 1 : 8
  h=subplot(2,4,iroi); hold on; get(h,'pos'); if iroi>4;set(h,'pos',[ans(1)+0 ans(2) ans(3)+0 ans(4)+0.3 ]);end
  title(sprintf('%s',str{iroi}))

  line([0.5 7.5],[-1.96 -1.96],'color','k');
  line([0.5 7.5],[1.96 1.96],'color','k')

    plot(1:4,squeeze(s1.tstat(:,iroi,:,:)),'color',[0.7 0.7 0.7],'linewidth',3)
  plot(1:4,squeeze(s2.tstat(:,iroi,:,:)),'color',[1 0.5 0],'linewidth',3)
%   set(gca,'xtick',[1 3 5 7],'xticklabel',[2 8 32 128],'tickdir','out');
  if iroi == 8 ; legend('A-P','D-P'); end
  axis square;  axis([0 5 -3 3])
end
set(k,'Position',[480 400 900 620]);

