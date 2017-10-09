%% READ IN DATA AND SAVE AS MAT
% pconn_preproc_read_data

clear all
restoredefaultpath

% -------------------------------------------------------------------------
% VERSION 4 - TASK
% -------------------------------------------------------------------------
denoise   = 'no';
v         = 1;
pad       = 1200/2; % 500 ms
rest      = 0;
task      = 1;
% -------------------------------------------------------------------------

indir1  = '/home/tpfeffer/pconn/rawdata/meg/';
outdir  = '/home/tpfeffer/pconn_cnt/proc/preproc/';

addpath /home/tpfeffer/pconn_cnt/matlab/
addpath('/home/tpfeffer/Documents/MATLAB/fieldtrip-20160919/')

% addpath('/home/tpfeffer/Documents/MATLAB/fieldtrip-20130925/')

ft_defaults

%%

SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
for m = 1:3
  for isubj = SUBJLIST
  
    indir = sprintf([indir1 'p%d/s%d/'],isubj,m);

    
    cont_dir = dir(indir);
    
    for i = 1 : length(cont_dir)
      tmp(i) = ~isempty(regexp(cont_dir(i).name,'(.ds)','once'));
    end
    
    ind = find(tmp>0);
    
    ibl(1:3) = 0;
    
    % no button presses registered for subject 10, m = 3
    if isubj == 10 && m == 3
      ind = [ind(3) ind(6)];
      cond = 3;
    end
    
    cnt = 0;
    
    for idir = 3:8
      
      tmp_dataset      = [cont_dir(idir).name];
     

        cfgs.datafile    = [indir '/' tmp_dataset '/' cont_dir(idir).name(1:end-2) 'meg4'];
        cfgs.headerfile  = [indir '/' tmp_dataset '/' cont_dir(idir).name(1:end-2) 'res4'];
      
      cfg = [];
      cfg.dataset    = [indir cont_dir(idir).name];
      
      if isubj == 22 && m == 3
        a=ft_read_event(cfgs.headerfile);
      else
        a=ft_read_event(cfg.dataset);
      end
      
      if ~(isubj ==10 && m==3)
        if size(cell2mat({a.value}),2)<10
            warning('Rest!');
            cond = 1;
            continue
        elseif sum(cell2mat({a.value})>50)>7
            warning('Task - button press!');
            cond = 2;
            continue
        elseif size(cell2mat({a.value}),2) > 10 && sum(cell2mat({a.value})>50)<7
          cnt = cnt + 1; 
          warning(sprintf('Task - counting - Block %d!',cnt));
            cond = 3;
            
            dat = ft_read_data(cfg.dataset,'chanindx',3); dat = squeeze(dat(1,:));
            [~,loc]=findpeaks(dat,'Minpeakprominence',30,'minpeakdistance',100);
            ts(isubj,m,cnt) = mean(diff(loc)/1200);
  %           continue
        end
      end
    end
  end
end
