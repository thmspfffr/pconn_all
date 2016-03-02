% pconn_preproc_pupdfa

% --------------------------------
% VERSION 10
% --------------------------------
v_in = 10;
v_out = 10;
v_rawdata = 6;
SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24];
dfa_overlap = 0.5;
fit_interv = [1 100];
calc_interv = [0.5 150];
% --------------------------------

outdir      = '/home/tpfeffer/pconn_all/proc/';
meg_in      = '/home/tpfeffer/pconn/proc/preproc/';

run ~/Documents/MATLAB/toolboxes/NBT-NBTv0.5.3-alpha/installNBT.m
%%

for isubj = SUBJLIST
  for iblock = 1 : 2
    for m = 1 : 3
      
      if ~exist([outdir sprintf('all_pup_dfa_s%d_b%d_m%d_v%d_processing.txt',isubj,iblock,m,v_out)])
        system(['touch ' outdir sprintf('all_pup_dfa_s%d_b%d_m%d_v%d_processing.txt',isubj,iblock,m,v_out)]);
      else
        continue
      end

      fprintf('Computing s%db%dm%d\n',isubj,iblock,m)
      try
        load(['/home/tpfeffer/pconn/proc/pup/' sprintf('pconn_pup_preproc_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,v_in)])
      catch me
        continue
      end
%       
      if ~isfield(pup,'dil')
        pup.dil = pup;
      end

      load(['~/pconn/proc/preproc/' sprintf('pconn_trig_resampled_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,3)])
     
      % get start and end triggers
      
      trigs(1) = find(trig.trial{1}>10,1,'first');
      trig.trial{1}(trigs:trigs+10000)=0;
      if ~isempty(find(trig.trial{1}>10,1,'first'))        
        trigs(2) = find(trig.trial{1}>10,1,'first');
      else
        warning('Second trigger not found!')
        continue
      end
      
      siginfo = nbt_Info;
      siginfo.converted_sample_frequency =  400;
      
      if size(pup.dil,1)>100000
        
        pup_dfa = nbt_doDFA(pup.dil, siginfo, fit_interv,calc_interv,dfa_overlap,0,0,[]);

      end
%       trigs = trigs - 20*400;
%       idx = idx(trigs(1):trigs(2));

%       while size(pup.dil,1)<size(idx,2)
%         idx = idx(2:end-1);
%       end
        
      
      % compare the two vectors again?
      
      
      save(sprintf([outdir 'all_pup_dfa_s%d_b%d_m%d_v%d.mat'],isubj,iblock,m,v_out),'pup_dfa');
      
 
    end
  end
end

%%

ord = pconn_randomization;


for isubj = SUBJLIST

  for iblock = 1 : 2
    for m = 1 : 3
      im = find(ord(isubj,:)==m);

      d=dir([outdir sprintf('all_pup_dfa_s%d_b*_m%d_v%d.mat',isubj,im,v_out)]);

      if length(d)<1
        num_blocks(isubj,m) = 0;
        continue
      elseif length(d)==1 && str2double(d(1).name(end-11)) == 1
        num_blocks(isubj,m) = 1;
        bl = 1;
        warning(sprintf('only one block which is b%d',bl));
      elseif  length(d)==1 && str2double(d(1).name(end-11)) == 2
        num_blocks(isubj,m) = 2;
        bl = 2;
      else
        num_blocks(isubj,m) = 2;
        bl = 1;
      end
        cnt = 0;
      for iblock = bl : num_blocks(isubj,m)
        cnt = cnt + 1;
        load([outdir d(cnt).name]);
        
        tmp(cnt) = pup_dfa.MarkerValues;
        
      end
      dfa(isubj,m) = nanmean(tmp);
    end
  end
end
        
        

dfa = dfa(SUBJLIST,:);
[x,y]=find(dfa==0);
dfa(x,:)=[];

