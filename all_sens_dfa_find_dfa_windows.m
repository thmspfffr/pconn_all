%% COMPUTE EXCITATION/INHIBITION BALANCE ACROSS SENSOR SPACE
% pconn_sens_dfa_find_dfa_windows

% Method by R Hardstone, discussed during meeting on 4th of March
% Takes band-pass filtered amplitude and time-resolved DFA as input
% and estimates E/I balance from that.

clear all

% aaa
% --------------------------------------------------------
% VERSION 1
% --------------------------------------------------------
v         = 2;
v_rawdata = 1;
is_src    = 0;
fsample   = 400;
SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32];
foi       = [8 12];
i_fit     = [3*ones(length(10:5:200),1) (10:5:200)'];
dfa_overlap = 0.5;
filt_ord = 2;
% --------------------------------------------------------

addpath ~/pconn/matlab/

addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/
addpath /home/tpfeffer/Documents/MATLAB/fieldtrip-20130925/

indir   = '/home/tpfeffer/pconn/proc/src/';
outdir   = '/home/tpfeffer/pconn_all/proc/';
% mridir = '/home/gnolte/neuconn/mri_data/';
plotdir = '/home/tpfeffer/pconn/proc/plots/';

run ~/Documents/MATLAB/toolboxes/NBT-NBTv0.5.3-alpha/installNBT.m
%%
for m = 1 : 3
  for isubj = SUBJLIST
    for ifoi = 1 : 1
      for ifit = 1 : length(i_fit)
        
        if ~exist(sprintf([outdir 'pconn_sens_dfa_find_dfa_windows_s%d_m%d_f%d_fit%d_v%d_processing.txt'],isubj,m,ifoi,ifit,v))
          system(['touch ' outdir sprintf('pconn_sens_dfa_find_dfa_windows_s%d_m%d_f%d_fit%d_v%d_processing.txt',isubj,m,ifoi,ifit,v)]);
        else
          continue
        end
        % %
        if isubj == 3
          load ~/pconn/matlab/pconn_sensorlabels.mat
        end
        
        d = dir(sprintf('/home/tpfeffer/pconn/proc/preproc/pconn_postpostproc_s%d_m%d_b*_v%d.mat',isubj,m,v_rawdata));
        
        for iblock = 1 : length(d)
          
          disp(sprintf('Processing MEG s%dm%df%d%b...',isubj,m,ifit,iblock));
          
          load(['/home/tpfeffer/pconn/proc/preproc/' d(iblock).name])
          
          if isubj == 3 || isubj == 2
            [~,idx_lab]=intersect(data.label,lab);
            data.trial{1}=data.trial{1}(idx_lab,:);
            data.label = lab;
            save(['~/pconn/matlab/' sprintf('pconn_idxlab.mat')],'idx_lab','-v7.3');
          end
          
          clear data_hi
          
          [mydata,epleng] = megdata2mydata(data);
          
          mydata = mydata(1:end-1000,:);
          
          siginfo = nbt_Info;
          siginfo.converted_sample_frequency = fsample;
          
          % compute bp-filtered signal
          tmp    = single(nbt_filter_fir(mydata,foi(ifoi,1),foi(ifoi,2),siginfo.converted_sample_frequency,filt_ord/foi(ifoi,1)));
          ampenv = abs(hilbert(tmp)); clear tmp mydata     
          dfa    = tp_dfa(ampenv, i_fit(ifit,:),fsample,dfa_overlap);
%           d=nbt_doDFA(ampenv(:,1), siginfo, i_fit(ifit,:), i_calc(ifit,:), .5,0,0,[]);
          clear data
          
          if isubj >= 32
            load(sprintf('~/pconn/proc/preproc/pconn_chanidx_s%d_m%d.mat',31,3))
          elseif isubj < 4
            load(sprintf('~/pconn/proc/preproc/pconn_chanidx_s%d_m%d.mat',4,1))
          elseif isubj == 17
            load(sprintf('~/pconn/proc/preproc/pconn_chanidx_s%d_m%d.mat',16,1))
          else
            load(sprintf('~/pconn/proc/preproc/pconn_chanidx_s%d_m%d.mat',isubj,m))
          end
          
          % INTERPOLATES THE RESULTS
          par.dfa(iblock)   = nanmean(pconn_sens_interp274(idx,dfa.exp));
          par.var(iblock)   = nanvar(pconn_sens_interp274(idx,dfa.exp));
          
          clear dfa dat ampenv p
        end
        
        save(sprintf([outdir 'pconn_sens_dfa_find_dfa_windows_s%d_m%d_f%d_fit%d_v%d.mat'],isubj,m,ifoi,ifit,v),'par','-v7.3');
        clear par
        
      end
    end
  end
end


error('STOP');

%%

SUBJLIST = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32]

for ifit = 1 : length(i_fit(:,2))
  for isubj = SUBJLIST
  load(sprintf('/home/tpfeffer/pconn_all/proc/pconn_sens_dfa_find_dfa_windows_s%d_m1_f1_fit%d_v1.mat',isubj,ifit))
  p(ifit,isubj) = nanmean(par.dfa);

  end
end

p = p(:,SUBJLIST);

for i = 1 : length(i_fit(:,2))-1
  
  [~,pval(i)]=ttest(p(i,:),p(i+1,:));
  
end

figure; set(gcf,'color','w'); hold on
% plot(nanmean(p(:,:,1),2),'color',[0.8 0.8 0.8],'linewidth',5);
% plot(nanmean(p(:,:,2),2),'color',[0.5 0.5 0.5],'linewidth',5);
subplot(1,3,1)
plot(i_fit(:,2),nanmean(p,2),'color',[0 0 0],'linewidth',3,'linestyle','-');
set(gca,'linewidth',2,'ticklength',[0.03 0.03],'tickdir','out');
set(gca,'fontsize',12,'fontweight','bold');

axis square


subplot(1,3,2)
plot(1:length(i_fit)-1,diff(nanmean(p,2)),'color',[0 0 0],'linewidth',3,'linestyle','-');
set(gca,'linewidth',2,'ticklength',[0.03 0.03],'tickdir','out');
set(gca,'fontsize',12,'fontweight','bold');

axis square

subplot(1,3,3)
plot(1:length(i_fit)-1,-log10(pval),'color',[0.5 0.5 0.5],'linewidth',3,'linestyle','-');
ylabel('Log_{2} (Power)','fontsize',13,'fontweight','bold');
set(gca,'fontsize',12,'fontweight','bold','tickdir','out');

axis square
set(gca,'linewidth',2,'ticklength',[0.03 0.03]);

print(gcf,'-depsc2',sprintf('~/pconn_all/plots/all_dfa_windowlength_v%d.eps',v))
%%
for ifit = 1 : 96
  
  x = rand(100000,1);

  [dfa,~,DFA_y_all] = nbt_doDFA(x, siginfo, i_fit(ifit,:),i_calc(ifit,:),dfa_overlap,0,0,[]);
  w(ifit,:)=[min(dfa.DFA_x) max(dfa.DFA_x)];
end 
          


%%
clear d

for j = 1 : 50
allwin = [10:5:100];

while 1
  xx = rand+0.5;
  if xx>= 0.5 && xx<1
    break
  end
end

  x = ffGn(4000,xx);

  for i = 1 : length(allwin)

    dfa = tp_dfa_test(x', [3 allwin(i)],4,0);
    d(i) = dfa.exp;

  end
  
  dd(j,:) = d-xx;
  
end
% end
% 
figure; set(gcf,'color','w'); hold on

plot(log10(dfa.win./400),log10(dfa.y{1}),'.')
plot(log10(dfa.win./400),0.5.*log10(dfa.win./400)+1.5)


%%



