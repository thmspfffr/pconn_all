
function all_par = pconn_all_read_neural_data(SUBJLIST,ord,para)


for isubj = SUBJLIST
  fprintf('Reading subj %d ...\n',isubj);
  for m = 1 : 3
    
    im = find(ord(isubj,:)==m);
    if strcmp(para.cond,'rst')
      if ~strcmp(para.str,'pow')
        load(sprintf('~/pconn/proc/dfa/pconn_src_dfa_s%d_m%d_f%d_v%d.mat',isubj,im,para.ifoi,v));
        eval(sprintf('all_par(:,:,m,isubj)  = par.%s;',para.str))
      else
        load(sprintf('~/pconn/proc/src/pconn_src_pow_s%d_m%d_f%d_v%d.mat',isubj,im,para.ifoi,1));
        eval(sprintf('all_par(:,:,m,isubj)  = par.%s;',para.str))
      end
    elseif strcmp(para.cond,'tsk')
      if ~strcmp(para.str,'pow')
        load(sprintf('~/pconn_cnt/proc/dfa/pconn_cnt_src_dfa_s%d_m%d_f%d_v%d.mat',isubj,im,para.ifoi,v));
        eval(sprintf('all_par(:,:,m,isubj)  = par.%s;',para.str))
      else
        load(sprintf('~/pconn_cnt/proc/src/pconn_cnt_src_pow_s%d_m%d_f%d_v%d.mat',isubj,im,para.ifoi,1));
        eval(sprintf('all_par(:,:,m,isubj)  = par.%s;',para.str))
      end
    end
  end
end

all_par   = all_par(:,:,:,SUBJLIST);
