%% PLOT SENSOR-LEVEL DFA RESULTS
% all_sens_dfa_filt_plot
clear

    
outdir    = '/home/tpfeffer/pconn_all/proc/';
plotdir   = '/home/tpfeffer/pconn_all/plots/';


allstr = {'dfa';'amp';'cvar';'var'};

for v = [1 2 9]
  
  if v<8
    SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33];
  else
    SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 ];
  end
  
  for istr = 1 : 4
    
    if ~exist(sprintf([outdir 'all_sens_dfa_filt_plot_s%d_v%d_processing.txt'],istr,v))
      system(['touch ' outdir sprintf('all_sens_dfa_filt_plot_s%d_v%d_processing.txt',istr,v)]);
    else
      continue
    end
    
    fprintf('Processing v%d s%d ...\n',v,istr);
    
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
    
    for ifoi = 1:4
      for isubj = SUBJLIST
        for m = 1 : 3
          
          im = find(ord(isubj,:)==m);
          
          tmp_cnt = pcbi_cnt(isubj);
          
          cnt(isubj,m) = nanmean(tmp_cnt((im*2-1):im*2));
          
          load(sprintf(['/home/tpfeffer/pconn_cnt/proc/dfa/' 'pconn_cnt_sens_dfa_filt_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,v_cnt));
          
          dfa_cnt_all(:,isubj,m,ifoi)  = nanmean(par.dfa,2);
          var_cnt_all(:,isubj,m,ifoi)  = nanmean(par.var,2);
          cvar_cnt_all(:,isubj,m,ifoi) = nanmean(par.cvar,2);
          amp_cnt_all(:,isubj,m,ifoi)  = nanmean(par.amp,2);
          
          
          load(sprintf(['~/pconn/proc/dfa/' 'pconn_sens_dfa_filt_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,v_res));
          
          dfa_res_all(:,isubj,m,ifoi)  = nanmean(par.dfa,2);
          var_res_all(:,isubj,m,ifoi)  = nanmean(par.var,2);
          cvar_res_all(:,isubj,m,ifoi) = nanmean(par.cvar,2);
          amp_res_all(:,isubj,m,ifoi)  = nanmean(par.amp,2);
          
          
          
        end
      end
    end
    
    % cnt = cnt(SUBJLIST,:);
    dfa_cnt_all  = double(dfa_cnt_all(:,SUBJLIST,:,:));
    var_cnt_all  = double(var_cnt_all(:,SUBJLIST,:,:));
    cvar_cnt_all = double(cvar_cnt_all(:,SUBJLIST,:,:));
    amp_cnt_all  = double(amp_cnt_all(:,SUBJLIST,:,:));
    
    dfa_res_all  = double(dfa_res_all(:,SUBJLIST,:,:));
    var_res_all  = double(var_res_all(:,SUBJLIST,:,:));
    cvar_res_all = double(cvar_res_all(:,SUBJLIST,:,:));
    amp_res_all  = double(amp_res_all(:,SUBJLIST,:,:));
    
    fprintf('Data loaded v%d s%d ...\n',v,istr);
    %% INDIV SUBJECTS
    % SUBJ = [1:5; 6:10; 11:15; 16:19];
          
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
        all_par_res = log10(amp_res_all);
        all_par_cnt = log10(amp_cnt_all);
      end
      
      for nsubj = SUBJLIST
        
        fprintf('Processing v%d s%d: indiv s%d ...\n',v,istr,nsubj);
        
        h=figure; set(h,'color','white'); %set(h,'Papertype','a4','visible','off')
        
        pars.scale=[0.5 0.8];
        
        isubj = find(SUBJLIST==nsubj);
        
        subplot(2,4,1)
        
        par = nanmean(all_par_res(:,isubj,:,1),3);
        
        pars.cbar = 0;
        pars.markersize = 0;
        pars.linewidth = 4;
        pars.resolution = 300;
        showfield_colormap(par,sa.locs_2D,pars);
        % showfield_colormap(a,sa.locs_2D,pars);
        
        colormap(hot)
        
        subplot(2,4,2)
        
        par = nanmean(all_par_res(:,isubj,:,2),3);
        % figure; set(gcf,'color','white');
        % pars.scale=[0.5 0.75];
        pars.cbar = 0;
        pars.markersize = 0;
        pars.linewidth = 4;
        pars.resolution = 300;
        showfield_colormap(par,sa.locs_2D,pars);
        colormap(hot)
        
        subplot(2,4,3)
        
        par = nanmean(all_par_res(:,isubj,:,3),3);
        
        % figure; set(gcf,'color','white');
        pars.cbar = 0;
        pars.markersize = 0;
        pars.linewidth = 4;
        pars.resolution = 300;
        showfield_colormap(par,sa.locs_2D,pars);
        colormap(hot)
        
        subplot(2,4,4)        
        
        par = nanmean(all_par_res(:,isubj,:,4),3);
        
        pars = [];
        pars.cbar = 0;
        pars.markersize = 0;
        pars.linewidth = 4;
        pars.resolution = 300;
        showfield_colormap(par,sa.locs_2D,pars);
        colormap(hot)
        
        subplot(2,4,5)
        
        par = nanmean(all_par_cnt(:,isubj,:,1),3);
        % figure; set(gcf,'color','white');
        % pars.scale=[0.5 0.75];
        pars.cbar = 0;
        pars.markersize = 0;
        pars.linewidth = 4;
        pars.resolution = 300;
        showfield_colormap(par,sa.locs_2D,pars);
        colormap(hot)
        
        subplot(2,4,6)
        
        par = nanmean(all_par_cnt(:,isubj,:,2),3);
        
        % figure; set(gcf,'color','white');
        pars.cbar = 0;
        pars.markersize = 0;
        pars.linewidth = 4;
        pars.resolution = 300;
        showfield_colormap(par,sa.locs_2D,pars);
        colormap(hot)
        
        subplot(2,4,7)
        
        par = nanmean(all_par_cnt(:,isubj,:,3),3);
        
        % figure; set(gcf,'color','white');
        pars.cbar = 0;
        pars.markersize = 0;
        pars.linewidth = 4;
        pars.resolution = 300;
        showfield_colormap(par,sa.locs_2D,pars);
        colormap(hot)
        
        subplot(2,4,8)
        
        par = nanmean(all_par_cnt(:,isubj,:,4),3);
        
        % figure; set(gcf,'color','white');
        pars.cbar = 0;
        pars.markersize = 0;
        pars.linewidth = 4;
        pars.resolution = 300;
        showfield_colormap(par,sa.locs_2D,pars);
        colormap(hot)
        
        print(gcf,'-djpeg100',sprintf('~/pconn_all/plots/all_sens_filt_%s_topo_indiv_tskrst_s%d_v%d.jpeg',str,nsubj,v_cnt))%%
        
        h=figure; set(h,'color','white'); %set(h,'Papertype','a4','visible','off')
        
        subplot(1,4,1)
        
        pars.scale = [-0.1 0.1];
%         djpeg100
        par = nanmean(all_par_cnt(:,isubj,:,1),3)-nanmean(all_par_res(:,isubj,:,1),3);
        
        pars = [];
        pars.cbar = 0;
        pars.markersize = 0;
        pars.linewidth = 4;
        pars.resolution = 300;
        showfield_colormap(par,sa.locs_2D,pars);
        % showfield_colormap(a,sa.locs_2D,pars);
        
        % colormap(hot)
        %
        subplot(1,4,2)
        
        par = nanmean(all_par_cnt(:,isubj,:,2),3)-nanmean(all_par_res(:,isubj,:,2),3);
        % figure; set(gcf,'color','white');
        % pars.scale=[0.5 0.75];
        pars.cbar = 0;
        pars.markersize = 0;
        pars.linewidth = 4;
        pars.resolution = 300;
        showfield_colormap(par,sa.locs_2D,pars);
        colormap(hot)
        
        subplot(1,4,3)
        
        par = nanmean(all_par_cnt(:,isubj,:,3),3)-nanmean(all_par_res(:,isubj,:,3),3);
        
        % figure; set(gcf,'color','white');
        pars.cbar = 0;
        pars.markersize = 0;
        pars.linewidth = 4;
        pars.resolution = 300;
        showfield_colormap(par,sa.locs_2D,pars);
        % colormap(hot)
       
        subplot(1,4,4)

        par = nanmean(all_par_cnt(:,isubj,:,4),3)-nanmean(all_par_res(:,isubj,:,4),3);
        
        % figure; set(gcf,'color','white');
        pars.cbar = 0;
        pars.markersize = 0;
        pars.linewidth = 4;
        pars.resolution = 300;
        showfield_colormap(par,sa.locs_2D,pars);
        % colormap(hot)
        
        print(gcf,'-djpeg100',sprintf('~/pconn_all/plots/all_sens_filt_%s_topo_indiv_task-rest_s%d_v%d.jpeg',str,nsubj,v_cnt))%%
        % close
      end

    %% PHARMA COMPARISON & TASK VS REST
    % -------------------------------------------
    fprintf('Processing v%d s%d: pharma comparison ...\n',v,istr);

    for mask = 0 : 1
      
      for ifoi = 1 : 4
        
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
          
          load(sprintf('~/pconn_all/proc/all_sens_filt_tsk_%s_clusterstat_contr%d_f%d_v%d.mat',str,icontr,ifoi,v_clust));
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
          print(gcf,'-djpeg100',sprintf('~/pconn_all/plots/all_sens_filt_%s_topo_tsk_c%d_f%d_mask%d_v%d.jpg',str,icontr,ifoi,mask,v_clust))
          
        end
        
        for icontr = 1 : 3
          
          clear stats
          
          load(sprintf('~/pconn_all/proc/all_sens_filt_res_%s_clusterstat_contr%d_f%d_v%d.mat',str,icontr,ifoi,v_clust));
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
          pars.scale=[-0.035 .035];
          
          pars.cbar = 0;
          pars.markersize = 0;
          pars.linewidth = 9;
          pars.resolution = 300;
          
          showfield_colormap(dfa_d1,sa.locs_2D,pars);
          %   colormap(parula)
          print(gcf,'-djpeg100',sprintf('~/pconn_all/plots/all_sens_filt_%s_topo_rst_c%d_f%d_mask%d_v%d.jpg',str,icontr,ifoi,mask,v_clust))
          
        end
      end
    end
  end
end
% %% COMPARE MEASURES
% % -------------------------------------
%
% ifoi = 1;
%
% for isubj = 1 : 24
%
%   p_dfa  = nanmean(dfa_res_all(:,isubj,:,ifoi),3);
%   p_var  = nanmean(var_res_all(:,isubj,:,ifoi),3);
%   p_cvar = nanmean(cvar_res_all(:,isubj,:,ifoi),3);
%
%   r_res(1,isubj) = corr(p_dfa,p_var);
%   r_res(2,isubj) = corr(p_dfa,p_cvar);
%   r_res(3,isubj) = corr(p_cvar,p_var);
%
%   p_dfa  = nanmean(dfa_cnt_all(:,isubj,:,ifoi),3);
%   p_var  = nanmean(var_cnt_all(:,isubj,:,ifoi),3);
%   p_cvar = nanmean(cvar_cnt_all(:,isubj,:,ifoi),3);
%
%   r_cnt(1,isubj) = corr(p_dfa,p_var);
%   r_cnt(2,isubj) = corr(p_dfa,p_cvar);
%   r_cnt(3,isubj) = corr(p_cvar,p_var);
%
%
%
% end
%
% for i = 1 : 3
%
%   [~,p_res(i)] = ttest(r_res(:,i));
%   [~,p_cnt(i)] = ttest(r_cnt(:,i));
%
% end
%
% %%
% ifoi = 1;
%
% clear r_res r_cnt
%
% for isens = 1 : 274
%
%   p_dfa  = nanmean(dfa_res_all(isens,:,:,ifoi),3);
%   p_var  = nanmean(var_res_all(isens,:,:,ifoi),3);
%   p_cvar = nanmean(cvar_res_all(isens,:,:,ifoi),3);
%
%   r_res(1,isens) = corr(p_dfa',p_var');
%   r_res(2,isens) = corr(p_dfa',p_cvar');
%   r_res(3,isens) = corr(p_cvar',p_var');
%
%   p_dfa  = nanmean(dfa_cnt_all(isens,:,:,ifoi),3);
%   p_var  = nanmean(var_cnt_all(isens,:,:,ifoi),3);
%   p_cvar = nanmean(cvar_cnt_all(isens,:,:,ifoi),3);
%
%   r_cnt(1,isens) = corr(p_dfa',p_var');
%   r_cnt(2,isens) = corr(p_dfa',p_cvar');
%   r_cnt(3,isens) = corr(p_cvar',p_var');
%
% end
%
% figure; set(gcf,'color','white');
% for i = 1 : 3
%
%   subplot(3,2,i*2-1)
%
%   pars.scale=[-1 1];
%   pars.cbar = 0;
%   pars.markersize = 0;
%   pars.linewidth = 9;
%   pars.resolution = 300;
%
%   showfield_colormap(r_res(i,:),sa.locs_2D,pars);
%   %   colormap(parula)
%
%   subplot(3,2,i*2)
%
%   showfield_colormap(r_cnt(i,:),sa.locs_2D,pars);
%   %   colormap(parula)
%
% end
%
% %%