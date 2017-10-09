%% Timescale: Pupil dilation analysis
% Loads preprocessed .mat files and plots processed data
% pconn_pup_preproc_092016

% takes input from pconn_pup_asc2mat (checked 10-2016)
% reads in only resting state data (*rest* indicated in raw filename)

clear


% -------------------------------------------------------------------------
% VERSION 01 - ignore saccades, *interpolate* blinks
% -------------------------------------------------------------------------
v = 1;
v_pup = 2;
% -------------------------------------------------------------------------

tp_addpaths
% restoredefaultpath

sampledir_cnt   = '/home/tpfeffer/pconn_cnt/proc/';
eventdir_cnt   = '/home/tpfeffer/pconn_cnt/proc/';
outdir      = '/home/tpfeffer/pconn_cnt/proc/';

sampledir_res   = '/home/tpfeffer/pconn/proc/pup/';
eventdir_res    = '/home/tpfeffer/pconn/proc/pup/';

% addpath /home/tpfeffer/Documents/MATLAB/fieldtrip-20130925/

SUBJLIST    = [4 5 6 7 9 10 11 12 13 15 16 17 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
addpath ~/pconn/matlab
% ft_defaults

%%
icond = 1;
ord   = pconn_randomization;

for m = 1:3
    for isubj = SUBJLIST
        im = find(ord(isubj,:)==m);
        
        if icond == 1
            smpdir    = dir(sprintf([sampledir_res '*samples*_s%d_*_m%d_*v%d.mat'],isubj,im,v_pup));
            evtdir    = dir(sprintf([eventdir_res '*events*_s%d_*_m%d_*v%d.mat'],isubj,im,v_pup));
        elseif icond == 2
            smpdir    = dir(sprintf([sampledir_cnt '*samples*_s%d_*_m%d_*v%d.mat'],isubj,im,v_pup));
            evtdir    = dir(sprintf([eventdir_cnt '*events*_s%d_*_m%d_*v%d.mat'],isubj,im,v_pup));
        end
        
        for iblock = 1 : length(evtdir)
            
            
            
            fprintf('Processing s%d b%d m%d ... \n',isubj,iblock,m)
            if icond == 1
                load([sampledir_res smpdir(iblock).name])
                load([eventdir_res evtdir(iblock).name])
                
            elseif icond == 2
                load([sampledir_cnt smpdir(iblock).name])
                load([eventdir_cnt evtdir(iblock).name])
            end
            if isubj == 26 && im == 1
                samples(:,2:4) = dat(:,1:3);
            end
            
%             if size(double(blinks),1)<=5
%                 b_dfa(isubj,m,iblock)=nan;
%             else
%                 tmp=tp_dfa(diff(double(blinks(:,1))),[3 size(blinks,1)-5],1,0.5,30);
%                 b_dfa(isubj,m,iblock) = tmp.exp;
%             end
            if size(samples,1) > 10000
                x = abs(diff(zscore(samples(:,4))));
                [~,idx]=findpeaks(double(x>0.20),'MinPeakDistance',200);
                allblinks(isubj,m,iblock) = length(idx);
                b(isubj,m,iblock)=length(blinks);
                
            else
                warning('not sufficient data...')
                allblinks(isubj,m,iblock) = nan;
                b(isubj,m,iblock)=nan;
                
            end
            
        end
    end
end

allblinks = nanmean(allblinks(SUBJLIST,:,:),3);
b = nanmean(b(SUBJLIST,:,:),3);
% b_dfa = nanmean(b_dfa(SUBJLIST,:,:),3);
error('Done')

%% PLOT


figure; set(gcf,'color','white'); hold on;
title('BLINKS')

bar(1.5,nanmean(b(:,1)),0.3,'FaceColor','r','EdgeColor','r'); ylabel('DFA-Exponent'); % ylim([0.65 0.85]);
bar(2,nanmean(b(:,2)),0.3,'FaceColor','b','EdgeColor','b'); ylabel('DFA-Exponent');   %ylim([0.65 0.85])
bar(2.5,nanmean(b(:,3)),0.3,'FaceColor','m','EdgeColor','m'); ylabel('DFA-Exponent');   %ylim([0.65 0.85])

lm = nanmean(b)-[nanstd(b)/sqrt(size(b,1))]; 
um = nanmean(b)+[nanstd(b)/sqrt(size(b,1))]; 

line([1.5 1.5],[lm(1) um(1)],'color',[0.6 0 0],'LineWidth',5);
line([2 2],[lm(2) um(2)],'color',[0 0 0.6],'LineWidth',5);
line([2.5 2.5],[lm(3) um(3)],'color',[0.6 0.1 0.6],'LineWidth',5);

ylabel('Blinks'); xlabel('Condition');
axis([1.3 2.8 -20 200]);

[~,p1,~,s1] = ttest(b(:,1),b(:,2)); p1
[~,p2,~,s2] = ttest(b(:,2),b(:,3)); p2
[~,p3,~,s3] = ttest(b(:,1),b(:,3)); p3

b1=t2smpbf(s1.tstat,28,28,1)
b2=t2smpbf(s2.tstat,28,28,1)
b3=t2smpbf(s3.tstat,28,28,1)

set(gca,'tickdir','out'); tp_editplots
print(gcf,'-depsc2',sprintf('~/pconn/proc/plots/pconn_blinks_cond%d_hb_v%d.eps',icond,v))


%%

%% HOW SHIT RELATES TO NEURAL STUFF

% --------------------------------------------------------
% VERSION 1
% --------------------------------------------------------
v_out     = 2;
v_res     = v_out;
v_cnt     = v_out;
SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];

outdir    = '/home/tpfeffer/dfa/proc/';

tp_addpaths

ord   = pconn_randomization;

for iifoi = 1:5
  for isubj = SUBJLIST
    fprintf('Loading data s%d f%d ...\n',isubj,iifoi);
    for m = 1 : 3
      
      im = find(ord(isubj,:)==m);
      
      load(sprintf(['~/pconn_cnt/proc/dfa/pconn_cnt_src_dfa_s%d_m%d_f%d_v%d.mat'],isubj,im,iifoi,v_out));
      
      dfa_all_cnt(:,m,isubj,iifoi)  = nanmean(par.dfa,2);
      %
      load(sprintf(['~/pconn/proc/dfa/pconn_src_dfa_s%d_m%d_f%d_v%d.mat'],isubj,im,iifoi,v_out));
      %
      dfa_all_res(:,m,isubj,iifoi)  = nanmean(par.dfa,2);
      
    end
  end
end

dfa_all_res  = dfa_all_res(:,:,SUBJLIST,:);
dfa_all_cnt  = dfa_all_cnt(:,:,SUBJLIST,:);

load(sprintf('~/pconn_all/proc/all_src_clusterstat_rst_%s_c%d_f%d_v%d.mat','dfa',1,2,2));

d_dfa_res = squeeze(nanmean(dfa_all_res(find(stats.mask),2,:,2),1)-nanmean(dfa_all_res(find(stats.mask),1,:,2),1))
d_dfa_cnt = squeeze(nanmean(dfa_all_cnt(find(stats.mask),2,:,2),1)-nanmean(dfa_all_cnt(find(stats.mask),1,:,2),1))

m_dfa_res = squeeze(nanmean(nanmean(dfa_all_res(:,:,:,2),2),1));
m_dfa_cnt = squeeze(nanmean(nanmean(dfa_all_cnt(:,:,:,2),2),1));


d_hrv     = b(:,2)-b(:,1);
% d_hrv_cnt = dfa_cnt(:,2)-dfa_cnt(:,1);

% d_hb =  hb_all(:,2)-hb_all(:,1);
% d_hb_cnt = hb_all_cnt(:,2)-hb_all_cnt(:,1);

% 
% m_hrv = nanmean(nanmean(dfa,2),3);
% m_hb  = nanmean(nanmean(hb,2),3);


%%
figure; set(gcf,'color','white');


subplot(2,2,1)
scatter(d_hrv,d_dfa_res,100,'markeredgecolor','w','markerfacecolor','k');% lsline
[r,p]=corr(d_hrv,d_dfa_res);
b = corrbf(r,28)

axis square; axis([-450 250 -0.14 0.14])
% print(gcf,'-depsc2',sprintf('~/pconn/proc/plots/pconn_hrv_scatter_hrv_res_v%d.eps',v))
set(gca,'tickdir','out'); tp_editplots

subplot(2,2,2)

% % figure; set(gcf,'color','white');
% scatter(d_hrv_cnt,d_dfa_cnt,100,'markeredgecolor','w','markerfacecolor','k'); lsline
% [r(2),p(2)]=corr(d_hrv_cnt,d_dfa_cnt);
% axis square; axis([-0.55 0.5 -0.15 0.15])
% set(gca,'tickdir','out'); tp_editplots
% 
% axis square

% print(gcf,'-depsc2',sprintf('~/pconn/proc/plots/pconn_hrv_scatter_hrv_cnt_v%d.eps',v))
% subplot(2,2,3)

% figure; set(gcf,'color','white');
% scatter(d_hb,d_dfa_res,100,'markeredgecolor','w','markerfacecolor','k'); lsline
% [r(3),p(3)]=corr(d_hb,d_dfa_res);
% axis square; axis([-5 30 -0.15 0.15])
% set(gca,'tickdir','out'); tp_editplots
% 
% axis square
% 
% % print(gcf,'-depsc2',sprintf('~/pconn/proc/plots/pconn_hrv_scatter_hb_res_v%d.eps',v))
% subplot(2,2,4)
% 
% % figure; set(gcf,'color','white');
% scatter(d_hb_cnt,d_dfa_cnt,100,'markeredgecolor','w','markerfacecolor','k'); lsline
% [r(4),p(4)]=corr(d_hb_cnt,d_dfa_cnt);
% axis square; axis square; axis([-5 30 -0.15 0.15])
% 
% set(gca,'tickdir','out'); tp_editplots

print(gcf,'-depsc2',sprintf('~/pconn/proc/plots/pconn_blinks_cond%d_scatter_v%d.eps',icond,v))
%%
figure; set(gcf,'color','white');

scatter(m_hrv,m_dfa_res,50,'markeredgecolor','w','markerfacecolor','k'); lsline
[r(1),p(1)]=corr(m_hrv,m_dfa_res);

figure; set(gcf,'color','white');

scatter(m_hb,m_dfa_res,50,'markeredgecolor','w','markerfacecolor','k'); lsline
[r(1),p(1)]=corr(m_hb,m_dfa_res);













%% CORR WITH BEHAVIOR

addpath ~/pcbi/

behav = pcbi_cnt(1:34);
clear bh

for isubj = SUBJLIST
    for m = 1 : 3
        
        im = find(ord(isubj,:)==m);
        
        bh(m,isubj) = nanmean(behav(isubj,im*2-1:im*2));
        
    end
end

all_behav	= bh(:,SUBJLIST)';
blinks    = nanmean(blinks,3);

%%
clear s_eng s_el

for m = 1:3
    for isubj = SUBJLIST
        
        im = find(ord(isubj,:)==m); isubj
        
        smpdir    = dir(sprintf([sampledir '*samples*_s%d_*_m%d_*v%d.mat'],isubj,im,v_pup));
        smpdir=smpdir(cell2mat({smpdir.bytes})>10000);
        
        for iblock = 1 : length(smpdir)
            
            
            load([sampledir smpdir(iblock).name])
            
            if isubj == 26 && im == 1
                samples(:,2:3) = dat(:,1:2);
            end
            
            if size(samples,1)>100000
                
                x=[samples(:,2) samples(:,3)];
                y=tp_detect_microsaccades(x,1000);
                s_eng(isubj,m,iblock) = length(y);
                
            else
                s_eng(isubj,m,iblock) = nan;
                s_el(isubj,m,iblock)=nan;
            end
            clear samples saccs x y
        end
        
        clear smpdir evtdir
    end
end

s_eng(s_eng==0)=nan;
s_eng = s_eng(SUBJLIST,:,:);


%% IMPORTANT: USE SS_ENG as estimated using tp_detect_microsaccades!!!
% ------------------------------
ss_eng = nanmean(s_eng,3);

% set(0,'defaultfigurepapertype','a4')

figure; set(gcf,'color','white');

set(gcf,'PaperPositionMode','manual','position',[400 400 650 850])
mksize = 130;

% PLOT 1
% ------------
subplot(3,2,1);

scatter(mean(blinks,2),nanmean(all_behav,2),mksize,'facecolor','k','markeredgecolor','w');

[r,p]=corr(mean(blinks,2),mean(all_behav,2))

axis([-100 500 -50 150])
lsline
title(sprintf('Blinks: r = %.2f | p = %.2f',r,p));xlabel('#Blinks'); ylabel('#Switches');

tp_editplots

% PLOT 2
% ------------
subplot(3,2,2);

scatter(mean(ss_eng,2),mean(all_behav,2),mksize,'facecolor','k','markeredgecolor','w');

axis([-300 1900 -50 150])
lsline
[r,p]=corr(mean(ss_eng,2),mean(all_behav,2))
title(sprintf('Sacc: r = %.2f | p = %.3f',r,p));xlabel('#Saccades'); ylabel('#Switches');

tp_editplots

% PLOT 3
% ------------
subplot(3,2,3); hold on

scatter(blinks(:,2)-blinks(:,1),all_behav(:,2)-all_behav(:,1),mksize,'facecolor','k','markeredgecolor','w');
[r,p]=corr(blinks(:,2)-blinks(:,1),all_behav(:,2)-all_behav(:,1))

axis([-600 600 -75 75])
lsline
title(sprintf('Blinks (ATX): r = %.2f | p = %.3f',r,p));xlabel('\Delta (#Blinks)'); ylabel('\Delta (#Switches)');

tp_editplots

% PLOT 4
% ------------
subplot(3,2,4); hold on

scatter(ss_eng(:,2)-ss_eng(:,1),all_behav(:,2)-all_behav(:,1),mksize,'facecolor','k','markeredgecolor','w');

axis([-900 900 -75 75])
lsline
[r,p]=corr(ss_eng(:,2)-ss_eng(:,1),all_behav(:,2)-all_behav(:,1))
title(sprintf('Sacc (ATX): r = %.2f | p = %.3f',r,p));xlabel('\Delta (#Saccades)'); ylabel('\Delta (#Switches)');

tp_editplots

% PLOT 5
% ------------
subplot(3,2,5); hold on

scatter(blinks(:,3)-blinks(:,1),all_behav(:,3)-all_behav(:,1),mksize,'facecolor','k','markeredgecolor','w');

[r,p]=corr(blinks(:,3)-blinks(:,1),all_behav(:,3)-all_behav(:,1))

axis([-600 600 -75 75])
lsline
title(sprintf('Blinks (DPZ): r = %.2f | p = %.3f',r,p));xlabel('\Delta (#Blinks)'); ylabel('\Delta (#Switches)');

tp_editplots

% PLOT 6
% ------------
subplot(3,2,6); hold on

scatter(ss_eng(:,3)-ss_eng(:,1),all_behav(:,3)-all_behav(:,1),mksize,'facecolor','k','markeredgecolor','w');

axis([-900 900 -75 75])
lsline
[r,p]=corr(ss_eng(:,3)-ss_eng(:,1),all_behav(:,3)-all_behav(:,1))
title(sprintf('Sacc (DPZ): r = %.2f | p = %.3f',r,p));xlabel('\Delta (#Saccades)'); ylabel('\Delta (#Switches)');

tp_editplots

print(gcf,'-depsc2',sprintf('~/pconn_all/plots/all_blinks_saccs_behav.eps'))

%% BAR PLOTS
figure; set(gcf,'color','w'); hold on

subplot(1,2,1); hold on

for i = 1 : length(SUBJLIST)
    scatter(1,ss_eng(i,3),100,'facecolor','k','markeredgecolor','w')
    
    scatter(2,ss_eng(i,1),100,'facecolor','k','markeredgecolor','w')
    scatter(3,ss_eng(i,2),100,'facecolor','k','markeredgecolor','w')
    plot(2:3,ss_eng(i,1:2),'color','k','linewidth',2)
    plot(1:2,ss_eng(i,[3 1]),'color','k','linewidth',2)
    
end

tp_editplots; axis([0.5 3.5 -200 2000]);

set(gca,'xtick',[1 2 3],'xticklabel',['DPZ';'PBO';'ATX']);
title('Saccades'); ylabel('# Saccades')

subplot(1,2,2); hold on

for i = 1 : length(SUBJLIST)
    scatter(1,bbb(i,3),100,'facecolor','k','markeredgecolor','w')
    
    scatter(2,bbb(i,1),100,'facecolor','k','markeredgecolor','w')
    scatter(3,bbb(i,2),100,'facecolor','k','markeredgecolor','w')
    plot(2:3,bbb(i,1:2),'color','k','linewidth',2)
    plot(1:2,bbb(i,[3 1]),'color','k','linewidth',2)
    
end

tp_editplots; axis([0.5 3.5 -50 650]);

set(gca,'xtick',[1 2 3],'xticklabel',['DPZ';'PBO';'ATX']);
title('Blinks'); ylabel('# Blinks')

print(gcf,'-depsc2',sprintf('~/pconn_all/plots/all_blinks_saccs_behav1.eps'))


%% PLOT FIXATION MAP

map = zeros(5000,5000);
ord   = pconn_randomization;

for m = 1:3
    
    im = find(ord(isubj,:)==m);
    
    for isubj = SUBJLIST
        
        isubj
        smpdir    = dir(sprintf([sampledir '*samples*_s%d_*_m%d_*v%d.mat'],isubj,im,v_pup));
        evtdir    = dir(sprintf([eventdir '*events*_s%d_*_m%d_*v%d.mat'],isubj,im,v_pup));
        
        for iblock = 1 : length(evtdir)
            
            load([sampledir smpdir(iblock).name])
            
            x=[samples(:,2) samples(:,3)];
            
            x(isnan(x(:,1)),:)=[];
            
            for i = 1 : length(x)
                
                map(round(x(i,1)+2000),round(x(i,2)+2000)) = map(round(x(1,1)+2000),round(x(1,2)+2000)) +1;
                
            end
        end
    end
end


%%

figure; set(gcf,'color','w');
scatter(nanmean(nanmean(s_el,3),2),nanmean(nanmean(s_eng,3),2),100,'facecolor','k','markeredgecolor','w')
tp_editplots
lsline
[r,p]=corr(nanmean(nanmean(s_el,3),2),nanmean(nanmean(s_eng,3),2))
xlabel('Saccades (Eye link)'); ylabel('Saccades (Engbert et al.)'); title(sprintf('r = %.2f | p = %.3f',r,p));


