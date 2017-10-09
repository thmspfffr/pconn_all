%% PLOT SENSOR-LEVEL DFA RESULTS
% all_sens_dfa_plot
clear


outdir    = '/home/tpfeffer/pconn_all/proc/';
plotdir   = '/home/tpfeffer/pconn_all/plots/';

tp_addpaths

allstr       = {'dfa';'cvar';'pow'};

% for v = [4]
v= 2;
if v<8
  SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
else
  SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 ];
end
%%
for istr = 1 : 3
  
  %     if ~exist(sprintf([outdir 'all_sens_dfa_plot_s%d_v%d_processing.txt'],istr,v))
  %       system(['touch ' outdir sprintf('all_sens_dfa_plot_s%d_v%d_processing.txt',istr,v)]);
  %     else
  %       continue
  %     end
  
  str = allstr{istr};
  
  v_cnt         = v;
  v_res         = v;
  v_clust       = v;
  
  addpath /home/gnolte/meg_toolbox/toolbox/
  addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
  addpath /home/gnolte/meg_toolbox/toolbox_nightly/
  addpath /home/gnolte/meg_toolbox/meg/
  
  load /home/tpfeffer/pconn/proc/src/pconn_sa_s26_m2_b1_v4.mat
  addpath /home/tpfeffer/pconn/matlab
  addpath ~/pcbi/
  
  ord           = pconn_randomization;
  
  for isubj = SUBJLIST
    
    fprintf('Processing v%d s%d isubj%d ...\n',v,istr,isubj);
    
    for m = 1 : 3
      
      if isubj >= 32
        load(sprintf('~/pconn/proc/preproc/pconn_chanidx_s%d_m%d.mat',31,3))
      elseif isubj < 4
        load(sprintf('~/pconn/proc/preproc/pconn_chanidx_s%d_m%d.mat',4,1))
      elseif isubj == 17
        load(sprintf('~/pconn/proc/preproc/pconn_chanidx_s%d_m%d.mat',16,1))
      else
        load(sprintf('~/pconn/proc/preproc/pconn_chanidx_s%d_m%d.mat',isubj,m))
      end
      
      im = find(ord(isubj,:)==m);
      
      %       if strcmp(str,'pow')
      cnt = 0; tmp_pow = zeros(274,991); tmp_cvar = zeros(274,991);
      
      for iblock = 1 : 2
        
        try
          cnt = cnt + 1;
          load(['~/pconn_cnt/proc/powspec/' sprintf('pconn_cnt_powspec_dat_s%d_b%d_m%d_v%d.mat',isubj,iblock,m,1)]);
        catch me
          warning('It is happening again! It is happening again!')
        end
        
        tmp_pow   = tmp_pow+pconn_sens_interp274(idx,squeeze(mean(dat.powspctrm,1)));
        tmp_cvar  = tmp_cvar+pconn_sens_interp274(idx,squeeze(var(dat.powspctrm)./mean(dat.powspctrm)));
        
      end
      
      pow_cnt_all(:,:,isubj,m) = tmp_pow./cnt;
      cvar_cnt_all(:,:,isubj,m) = tmp_cvar./cnt;
      
      cnt = 0; tmp_pow = zeros(274,991); tmp_cvar = zeros(274,991);
      
      for iblock = 1 : 2
        
        try
          cnt = cnt + 1;
          load(['~/pconn/proc/powspec/' sprintf('pconn_powspec_dat_s%d_b%d_m%d_v%d.mat',isubj,iblock,m,1)]);
        catch me
          warning('It is happening again! It is happening again!')
        end
        
        tmp_pow   = tmp_pow+pconn_sens_interp274(idx,squeeze(mean(dat.powspctrm,1)));
        tmp_cvar  = tmp_cvar+pconn_sens_interp274(idx,squeeze(var(dat.powspctrm)./mean(dat.powspctrm)));
        
      end
      
      pow_res_all(:,:,isubj,m)  = tmp_pow./cnt;
      cvar_res_all(:,:,isubj,m) = tmp_cvar./cnt;
      %       end
      
      for ifoi = 1:4
        
        
        tmp_cnt = pcbi_cnt(isubj);
        
        cnt(isubj,m) = nanmean(tmp_cnt((im*2-1):im*2));
        
        load(sprintf(['/home/tpfeffer/pconn_cnt/proc/dfa/' 'pconn_cnt_sens_dfa_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,v_cnt));
        
        dfa_cnt_all(:,isubj,m,ifoi)  = nanmean(par.dfa,2);
        %         cvar_cnt_all(:,isubj,m,ifoi) = nanmean(par.cvar,2);
        
        load(sprintf(['~/pconn/proc/dfa/' 'pconn_sens_dfa_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,v_res));
        
        dfa_res_all(:,isubj,m,ifoi)  = nanmean(par.dfa,2);
        %         cvar_res_all(:,isubj,m,ifoi) = nanmean(par.cvar,2);
        
        
        
      end
    end
  end
  
  dfa_cnt_all  = double(dfa_cnt_all(:,SUBJLIST,:,:));
  
  dfa_res_all  = double(dfa_res_all(:,SUBJLIST,:,:));
  
  cvar_cnt_all = double(cvar_cnt_all(:,:,SUBJLIST,:));
  %     pow_cnt_all =
  cvar_res_all = double(cvar_res_all(:,:,SUBJLIST,:));
  
  
  fprintf('Data loaded v%d s%d ...\n',v,istr); %error('!')
  %% INDIV SUBJECTS
  
  FOI_cvar = [2 4; 4 8; 8 12; 12 36];
  %     str = 'dfa';
  if strcmp(str,'dfa')
    all_par_res = dfa_res_all;
    all_par_cnt = dfa_cnt_all;
  elseif strcmp(str,'var')
    all_par_res = var_res_all;
    all_par_cnt = var_cnt_all;
  elseif strcmp(str,'cvar')
    all_par_res = cvar_res_all;
    all_par_cnt = cvar_cnt_all;
  elseif strcmp(str,'amp')
    all_par_res = zscore(amp_res_all);
    all_par_cnt = zscore(amp_cnt_all);
  end
  
  for nsubj = SUBJLIST
    
    fprintf('Processing v%d s%d: indiv s%d ...\n',v,istr,nsubj);
    
    h=figure; set(h,'color','white'); %set(h,'Papertype','a4','visible','off')
    
    isubj = find(SUBJLIST==nsubj);
    
    for i = 1 : 4
      subplot(2,4,i)
      
      par = nanmean(nanmean(nanmean(cvar_res_all(:,FOI_cvar(i,1):FOI_cvar(i,2),:,:),2),3),4);
      
      pars.cbar = 0;
      pars.markersize = 0;
      pars.linewidth = 4;
      pars.resolution = 300;
      showfield_colormap(par,sa.locs_2D,pars);
      
      subplot(2,4,4+i)
      
      par = nanmean(nanmean(nanmean(cvar_res_all(:,FOI_cvar(i,1):FOI_cvar(i,2),:,:),2),3),4);
      showfield_colormap(par,sa.locs_2D,pars);
  
      
    end
    colormap(inferno)
%     cmap = cbrewer('div', 'RdBu', 400,'pchip'); cmap = cmap(1:200,:); cmap=cmap(end:-1:1,:);
%     colormap(cmap);
    
    print(gcf,'-djpeg100',sprintf('~/pconn_all/plots/all_sens_%s_topo_indiv_tskrst_s%d_v%d.jpeg',str,nsubj,v_cnt))%%
    
    h=figure; set(h,'color','white'); %set(h,'Papertype','a4','visible','off')
    
    for i = 1 : 4
      
      subplot(1,4,i)
      
      pars.scale = [-0.1 0.1];
      par = nanmean(all_par_cnt(:,isubj,:,i),3)-nanmean(all_par_res(:,isubj,:,i),3);
      
      pars = [];
      pars.cbar = 0;
      pars.markersize = 0;
      pars.linewidth = 3;
      pars.resolution = 300;
      showfield_colormap(par,sa.locs_2D,pars);
      
    end
    
    print(gcf,'-djpeg100',sprintf('~/pconn_all/plots/all_sens_%s_topo_indiv_task-rest_s%d_v%d.jpeg',str,nsubj,v_cnt))%%
    close
  end
  
  
end

error('!')

%% PHARMA COMPARISON & TASK VS REST
% -------------------------------------------

for istr = 1 : 4
  
  str = allstr{istr};
  
  for mask = 0 : 1
    
    for ifoi = 1 : 4
      
      fprintf('Processing v%d ifoi%d s%d: pharma comparison ...\n',v,ifoi,istr);
      
      
      if strcmp(str,'dfa')
        all_par_res = dfa_res_all(:,:,:,ifoi);
        all_par_cnt = dfa_cnt_all(:,:,:,ifoi);
      elseif strcmp(str,'var')
        all_par_res = log10(var_res_all(:,:,:,ifoi));
        all_par_cnt = log10(var_cnt_all(:,:,:,ifoi));
      elseif strcmp(str,'cvar')
        all_par_res = cvar_res_all(:,:,:,ifoi);
        all_par_cnt = cvar_cnt_all(:,:,:,ifoi);
      elseif strcmp(str,'amp')
        all_par_res = log10(amp_res_all(:,:,:,ifoi));
        all_par_cnt = log10(amp_cnt_all(:,:,:,ifoi));
      end
      
      contrasts = [2 1; 3 1; 2 3];
      % dfa_d2 = nanmean(nanmean(all_par_cnt,3),2)-nanmean(nanmean(all_par_res,3),2);
      
      for icontr = 1 : 3
        
        load(sprintf('~/pconn_all/proc/all_sens_tsk_%s_clusterstat_contr%d_f%d_v%d.mat',str,icontr,ifoi,v_clust));
        %   senssel = find(stats.mask);
        
        % PHARMA CONTRAST DURING CNT
        if mask
          dfa_d1 = stats.mask.*(nanmean(all_par_cnt(:,:,contrasts(icontr,1))-all_par_cnt(:,:,contrasts(icontr,2)),2));
        else
          dfa_d1 = nanmean(all_par_cnt(:,:,contrasts(icontr,1))-all_par_cnt(:,:,contrasts(icontr,2)),2);
        end
        figure; set(gcf,'color','white');
        
        r = max([min(abs(dfa_d1)) max(abs(dfa_d1))]);
        pars.scale=[-r r];
        pars.scale=[-0.03 .03];
        
        pars.cbar = 0;
        pars.markersize = 0;
        pars.linewidth = 9;
        pars.resolution = 300;
        
        showfield_colormap(dfa_d1,sa.locs_2D,pars);
        %   colormap(parula)
        print(gcf,'-djpeg100',sprintf('~/pconn_all/plots/all_sens_%s_topo_tsk_c%d_f%d_mask%d_v%d.jpg',str,icontr,ifoi,mask,v_clust))
        
      end
      
      for icontr = 1 : 3
        
        clear stats
        
        load(sprintf('~/pconn_all/proc/all_sens_res_%s_clusterstat_contr%d_f%d_v%d.mat',str,icontr,ifoi,v_clust));
        %   senssel = find(stats.mask);
        
        % PHARMA CONTRAST DURING REST
        if mask
          dfa_d1 = stats.mask.*(nanmean(all_par_res(:,:,contrasts(icontr,1))-all_par_res(:,:,contrasts(icontr,2)),2));
        else
          dfa_d1 = (nanmean(all_par_res(:,:,contrasts(icontr,1))-all_par_res(:,:,contrasts(icontr,2)),2));
        end
        
        figure; set(gcf,'color','white');
        
        r = max([min(abs(dfa_d1)) max(abs(dfa_d1))]);
        pars.scale=[-r r];
        if strcmp(str,'dfa')
          pars.scale=[-0.03 .03];
        end
        
        pars.cbar = 0;
        pars.markersize = 0;
        pars.linewidth = 9;
        pars.resolution = 300;
        
        showfield_colormap(dfa_d1,sa.locs_2D,pars);
        %   colormap(parula)
        print(gcf,'-djpeg100',sprintf('~/pconn_all/plots/all_sens_%s_topo_rst_c%d_f%d_mask%d_v%d.jpg',str,icontr,ifoi,mask,v_clust))
        
      end
      
      for icontr = 1 : 3
        
        clear stats
        
        load(sprintf('~/pconn_all/proc/all_sens_res_%s_clusterstat_contr%d_f%d_v%d.mat',str,icontr,ifoi,v_clust));
        %   senssel = find(stats.mask);
        
        % PHARMA CONTRAST DURING REST
        dfa_d1 = squeeze(nanmean(all_par_res(:,:,icontr),2));
        
        
        figure; set(gcf,'color','white');
        
        r = [min(dfa_d1) max(dfa_d1)];
        %         pars.scale=[ r];
        if strcmp(str,'dfa')
          pars.scale=[0.55 0.7];
        end
        
        pars.cbar = 0;
        pars.markersize = 0;
        pars.linewidth = 9;
        pars.resolution = 300;
        
        showfield_colormap(dfa_d1,sa.locs_2D,pars);
        colormap(inferno)
        print(gcf,'-djpeg100',sprintf('~/pconn_all/plots/all_sens_%s_topo_rst_avg_f%d_mask%d_v%d.jpg',str,ifoi,mask,v_clust))
        
        
      end
      
      
    end
  end
end
% end

%% PLOT CVAR OVER FREQ

load(['~/pconn_all/proc/' sprintf('all_sens_cvar.mat')]);

c1=squeeze(nanmean(cvar_res_all(:,:,:,1),1));
c2=squeeze(nanmean(cvar_res_all(:,:,:,2),1));
c3=squeeze(nanmean(cvar_res_all(:,:,:,3),1));

para = [];
para.paired = 1;
para.nperm = 20000;
para.clusteralpha = 0.05;

s=tp_permtest2D(c1,c2,para);

%%
[h1,p1]=ttest(c1,c2,'dim',2,'alpha',0.025)
[h2,p2]=ttest(c1,c3,'dim',2,'alpha',0.025)

figure; set(gcf,'color','white'); hold on

plot(log2(dat.freq),squeeze(nanmean(c(:,:,1),1)),'linewidth',3)
plot(log2(dat.freq),squeeze(nanmean(c(:,:,2),1)),'linewidth',3)
plot(log2(dat.freq),squeeze(nanmean(c(:,:,3),1)),'linewidth',3)

plot(log2(dat.freq(find(h1))),1.3e-28*ones(sum(h1),1),'*')
plot(log2(dat.freq(find(h2))),1.2e-28*ones(sum(h2),1),'*')


