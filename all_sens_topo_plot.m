%% PLOT SENSOR-LEVEL DFA RESULTS
% pconn_cnt_sens_dfa_plot

clear

v_cnt         = 22;
v_res         = 22;

SUBJLIST    = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24];

addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/
addpath /home/tpfeffer/pconn/matlab
addpath ~/pcbi/
addpath ~/Documents/MATLAB/cbrewer/cbrewer/

outdir   = '/home/tpfeffer/pconn_cnt/proc/dfa/';
plotdir = '/home/tpfeffer/pconn_all/plots/';


load /home/tpfeffer/pconn/proc/src/pconn_sa_s4_m1_b1_v1.mat

%%

ord       = pconn_randomization;


for ifoi = 1:4
  for isubj = SUBJLIST
    for m = 1 : 3
      
      im = find(ord(isubj,:)==m);
      
      tmp_cnt = pcbi_cnt(isubj);
      
      cnt(isubj,m) = nanmean(tmp_cnt((im*2-1):im*2));
      
      load(sprintf([outdir 'pconn_cnt_sens_dfa_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,v_cnt));
      
      dfa_cnt_all(:,isubj,m,ifoi)  = nanmean(par.dfa,2);
      var_cnt_all(:,isubj,m,ifoi)  = nanmean(par.var,2);
      cvar_cnt_all(:,isubj,m,ifoi) = nanmean(par.cvar,2);
      pow_cnt_all(:,isubj,m,ifoi)  = nanmean(par.pow,2); 
      amp_cnt_all(:,isubj,m,ifoi)  = nanmean(par.amp,2); clear par

      
      load(sprintf(['~/pconn/proc/dfa/' 'pconn_sens_dfa_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,v_res));
      
      dfa_res_all(:,isubj,m,ifoi)  = nanmean(par.dfa,2);
      var_res_all(:,isubj,m,ifoi)  = nanmean(par.var,2);
      cvar_res_all(:,isubj,m,ifoi) = nanmean(par.cvar,2);
      pow_res_all(:,isubj,m,ifoi)  = nanmean(par.pow,2);
      amp_res_all(:,isubj,m,ifoi)  = nanmean(par.amp,2); clear par
     
    end
  end
end

% cnt = cnt(SUBJLIST,:);
dfa_cnt_all  = double(dfa_cnt_all(:,SUBJLIST,:,:));
var_cnt_all  = double(var_cnt_all(:,SUBJLIST,:,:));
cvar_cnt_all = double(cvar_cnt_all(:,SUBJLIST,:,:));
pow_cnt_all  = double(pow_cnt_all(:,SUBJLIST,:,:));
amp_cnt_all  = double(amp_cnt_all(:,SUBJLIST,:,:));

dfa_res_all  = double(dfa_res_all(:,SUBJLIST,:,:));
var_res_all  = double(var_res_all(:,SUBJLIST,:,:));
cvar_res_all = double(cvar_res_all(:,SUBJLIST,:,:));
pow_res_all  = double(pow_res_all(:,SUBJLIST,:,:));
amp_res_all  = double(amp_res_all(:,SUBJLIST,:,:));

%% INDIV SUBJECTS
% SUBJ = [1:5; 6:10; 11:15; 16:19];

str = 'dfa';
addpath ~/Documents/MATLAB/cbrewer/cbrewer/
cmap = cbrewer('div', 'RdBu', 100,'pchip');% colormap(autumn)
cmap = cmap(end:-1:1,:);

for ifoi = 1 : 4
  
  if strcmp(str,'dfa')
    all_par_res = dfa_res_all(:,:,:,ifoi);
    all_par_cnt = dfa_cnt_all(:,:,:,ifoi);
  elseif strcmp(str,'var')
    all_par_res = var_res_all(:,:,:,ifoi);
    all_par_cnt = var_cnt_all(:,:,:,ifoi);
  elseif strcmp(str,'cvar')
    all_par_res = cvar_res_all(:,:,:,ifoi);
    all_par_cnt = cvar_cnt_all(:,:,:,ifoi);
  elseif strcmp(str,'pow')
    all_par_res = log10(pow_res_all(:,:,:,ifoi));
    all_par_cnt = log10(pow_cnt_all(:,:,:,ifoi));
  end
  
  for isubj = 1 : length(SUBJLIST)
    h=figure; set(h,'color','white'); %set(h,'Papertype','a4','visible','off')
    
    subj = SUBJLIST(isubj);
    
    subplot(1,3,1)
    
    par = nanmean(squeeze(all_par_res(:,isubj,:)),3);
    
    pars = [];
    pars.scale=[min(par(:)) max(par(:))];
    pars.cbar = 0;
    pars.markersize = 0;
    pars.linewidth = 4;
    pars.resolution = 300;
    showfield_colormap(nanmean(par,2),sa.locs_2D,pars);
    
    colormap(hot)
    
    subplot(1,3,2)
    
    par = nanmean(squeeze(all_par_cnt(:,isubj,:)),3);

    pars.cbar = 0;
    pars.markersize = 0;
    pars.scale=[min(par(:)) max(par(:))];
    pars.linewidth = 4;
    pars.resolution = 300;
    showfield_colormap(nanmean(par,2),sa.locs_2D,pars);
    
    subplot(1,3,3)
    
    par = nanmean(squeeze(all_par_cnt(:,isubj,:)-all_par_res(:,isubj,:)),3);
    
    pars.cbar = 0;
    pars.scale=[min(par(:)) max(par(:))];
    pars.markersize = 0;
    pars.linewidth = 4;
    pars.resolution = 300;
    showfield_colormap(nanmean(par,2),sa.locs_2D,pars);
    colormap(parula)
    print(gcf,'-djpeg',sprintf('~/pconn_cnt/plots/pconn_topo_all_%s_s%d_f%d_v%d.jpg',str,subj,ifoi,v_cnt))%%
    
    close
  end
end

error('!!')




%%
% ------------------------------------------------
% AVERAGE ACROSS CONDITIONS
% ------------------------------------------------
  
ifoi      = 1;
str       = 'cvar';
contrasts = [2 1; 3 1; 2 3];

if strcmp(str,'dfa')
  all_par_res = dfa_res_all(:,:,:,ifoi);
  all_par_cnt = dfa_cnt_all(:,:,:,ifoi);
elseif strcmp(str,'var')
  all_par_res = var_res_all(:,:,:,ifoi);
  all_par_cnt = var_cnt_all(:,:,:,ifoi);
elseif strcmp(str,'cvar')
  all_par_res = cvar_res_all(:,:,:,ifoi);
  all_par_cnt = cvar_cnt_all(:,:,:,ifoi);
elseif strcmp(str,'pow')
  all_par_res = log10(pow_res_all(:,:,:,ifoi));
  all_par_cnt = log10(pow_cnt_all(:,:,:,ifoi));
elseif strcmp(str,'amp')
  all_par_res = log10(amp_res_all(:,:,:,ifoi));
  all_par_cnt = log10(amp_cnt_all(:,:,:,ifoi));
end


figure; set(gcf,'color','white');

dfa_d1 = nanmean(nanmean(all_par_cnt,3),2);

pars.scale = [min(dfa_d1) max(dfa_d1)];
pars.markersize = 0;
pars.cbar = 0;
pars.linewidth = 9;
pars.resolution = 300;

showfield_colormap(dfa_d1,sa.locs_2D,pars);
colormap(hot)

print(gcf,'-djpeg',sprintf('~/pconn_all/plots/pconn_sens_topo_tsk_%s_avg_f%d_v%d.jpg',str,ifoi,v_cnt))

dfa_d1 = nanmean(nanmean(all_par_res,3),2);

pars.scale = [min(dfa_d1) max(dfa_d1)]

showfield_colormap(dfa_d1,sa.locs_2D,pars);
colormap(hot)

print(gcf,'-djpeg',sprintf('~/pconn_all/plots/pconn_sens_topo_rst_%s_avg_f%d_v%d.jpg',str,ifoi,v_cnt))

%
%% PHARMA COMPARISON & TASK VS REST
% -------------------------------------------

tval = 0;
ifoi      = 4;
str       = 'dfa';
contrasts = [2 1; 3 1; 2 3];

if strcmp(str,'dfa')
  all_par_res = dfa_res_all(:,:,:,ifoi);
  all_par_cnt = dfa_cnt_all(:,:,:,ifoi);
elseif strcmp(str,'var')
  all_par_res = var_res_all(:,:,:,ifoi);
  all_par_cnt = var_cnt_all(:,:,:,ifoi);
elseif strcmp(str,'cvar')
  all_par_res = cvar_res_all(:,:,:,ifoi);
  all_par_cnt = cvar_cnt_all(:,:,:,ifoi);
elseif strcmp(str,'pow')
  all_par_res = log10(pow_res_all(:,:,:,ifoi));
  all_par_cnt = log10(pow_cnt_all(:,:,:,ifoi));
end

% cmap = cbrewer('div', 'RdBu', 100,'pchip');% colormap(autumn)
% cmap = cmap(end:-1:1,:);

cmap =   [0.941176474094391 0.941176474094391 0.941176474094391;0.935294091701508 0.935294091701508 0.935294091701508;0.929411768913269 0.929411768913269 0.929411768913269;0.923529386520386 0.923529386520386 0.923529386520386;0.917647063732147 0.917647063732147 0.917647063732147;0.911764681339264 0.911764681339264 0.911764681339264;0.905882358551025 0.905882358551025 0.905882358551025;0.899999976158142 0.899999976158142 0.899999976158142;0.867568910121918 0.852902054786682 0.877736866474152;0.835137844085693 0.805804193019867 0.855473756790161;0.802706778049469 0.758706271648407 0.833210647106171;0.770275771617889 0.711608350276947 0.81094753742218;0.737844705581665 0.664510428905487 0.78868442773819;0.705413639545441 0.617412567138672 0.766421318054199;0.672982573509216 0.570314645767212 0.744158208370209;0.647758424282074 0.533682942390442 0.72684246301651;0.622534275054932 0.497051239013672 0.709526658058167;0.597310066223145 0.460419565439224 0.692210912704468;0.572085916996002 0.423787862062454 0.674895167350769;0.54432600736618 0.383473604917526 0.655838668346405;0.516566097736359 0.343159347772598 0.636782169342041;0.488806158304214 0.30284509062767 0.617725670337677;0.461046248674393 0.262530833482742 0.598669171333313;0.433286309242249 0.222216591238976 0.579612672328949;0.405526399612427 0.181902334094048 0.560556173324585;0.38470646739006 0.151666641235352 0.546263813972473;0.363886535167694 0.121430955827236 0.531971454620361;0.343066573143005 0.0911952629685402 0.51767909526825;0.322246640920639 0.0609595701098442 0.503386676311493;0.301426708698273 0.0307238809764385 0.489094316959381;0.280606776475906 0.000488189951283857 0.474801957607269;0.293446481227875 0.00417250022292137 0.461163550615311;0.306286215782166 0.00785681046545506 0.447525143623352;0.319125920534134 0.0115411207079887 0.433886736631393;0.331965625286102 0.0152254309505224 0.420248329639435;0.344805359840393 0.0189097411930561 0.406609892845154;0.357645064592361 0.0225940514355898 0.392971485853195;0.370484799146652 0.0262783616781235 0.379333078861237;0.383324503898621 0.0299626719206572 0.365694671869278;0.396164208650589 0.0336469821631908 0.352056264877319;0.40900394320488 0.0373312942683697 0.338417857885361;0.421843647956848 0.0410156026482582 0.324779450893402;0.434683352708817 0.044699914753437 0.311141043901443;0.447523087263107 0.0483842231333256 0.297502636909485;0.460362792015076 0.0520685352385044 0.283864229917526;0.473202496767044 0.0557528436183929 0.270225793123245;0.486042231321335 0.0594371557235718 0.256587386131287;0.498881936073303 0.0631214678287506 0.242948994040489;0.511721670627594 0.0668057799339294 0.229310572147369;0.52456134557724 0.0704900845885277 0.215672165155411;0.537401080131531 0.0741743966937065 0.202033758163452;0.550240814685822 0.0778587087988853 0.188395351171494;0.563080549240112 0.0815430209040642 0.174756944179535;0.575920224189758 0.0852273255586624 0.161118522286415;0.588759958744049 0.0889116376638412 0.147480115294456;0.60159969329834 0.0925959497690201 0.133841708302498;0.614439368247986 0.0962802618741989 0.120203301310539;0.627279102802277 0.0999645665287971 0.106564886868;0.640118837356567 0.103648878633976 0.0929264798760414;0.652958512306213 0.107333190739155 0.0792880654335022;0.665798246860504 0.111017502844334 0.0656496584415436;0.678637981414795 0.114701807498932 0.0520112477242947;0.691477656364441 0.118386119604111 0.0383728370070457;0.704317390918732 0.12207043170929 0.024734428152442;0.70988517999649 0.127976059913635 0;0.714490175247192 0.141817703843117 0;0.719095170497894 0.155659362673759 0;0.723700165748596 0.169501006603241 0;0.728305160999298 0.183342665433884 0;0.73291015625 0.197184309363365 0;0.737515151500702 0.211025953292847 0;0.742120146751404 0.224867612123489 0;0.746725142002106 0.238709256052971 0;0.751330137252808 0.252550899982452 0;0.75593513250351 0.266392558813095 0;0.760540127754211 0.280234217643738 0;0.765145123004913 0.294075846672058 0;0.769750118255615 0.307917505502701 0;0.774355113506317 0.321759164333344 0;0.778960108757019 0.335600793361664 0;0.783565163612366 0.349442452192307 0;0.788170158863068 0.363284111022949 0;0.79277515411377 0.377125769853592 0;0.797380149364471 0.390967398881912 0;0.801985144615173 0.404809057712555 0;0.806590139865875 0.418650716543198 0;0.811195135116577 0.432492345571518 0;0.815800130367279 0.446334004402161 0;0.820405125617981 0.460175663232803 0;0.825010120868683 0.474017292261124 0;0.829615116119385 0.487858951091766 0;0.834220111370087 0.501700580120087 0;0.838825106620789 0.515542268753052 0;0.84343010187149 0.529383897781372 0;0.848035097122192 0.543225526809692 0;0.852640092372894 0.557067215442657 0;0.857245087623596 0.570908844470978 0;0.861850082874298 0.584750533103943 0;0.866455078125 0.598592162132263 0;0.871060073375702 0.612433791160583 0;0.875665068626404 0.626275479793549 0;0.880270063877106 0.640117108821869 0;0.884875059127808 0.653958737850189 0;0.88948005437851 0.667800426483154 0;0.894085049629211 0.681642055511475 0;0.898690044879913 0.695483684539795 0;0.903295040130615 0.70932537317276 0;0.907900035381317 0.72316700220108 0;0.912505030632019 0.737008631229401 0;0.917110025882721 0.750850319862366 0;0.921715021133423 0.764691948890686 0;0.926320016384125 0.778533577919006 0;0.930925071239471 0.792375266551971 0;0.935530066490173 0.806216895580292 0;0.940135061740875 0.820058524608612 0;0.944740056991577 0.833900213241577 0;0.949345052242279 0.847741842269897 0;0.953950047492981 0.861583530902863 0;0.958555042743683 0.875425159931183 0;0.963160037994385 0.889266788959503 0;0.967765033245087 0.903108477592468 0;0.972370028495789 0.916950106620789 0;0.97697502374649 0.930791735649109 0;0.981580018997192 0.944633424282074 0;0.986185014247894 0.958475053310394 0;0.990790009498596 0.972316682338715 0;0.995395004749298 0.98615837097168 0;1 1 0];

for icontr = 1 : 2
  
  pars = [];
  
  % --w----------------------------------------------
  % PHARMA CONTRAST DURING TASK
  % ------------------------------------------------
  
  load(sprintf('~/pconn_all/proc/all_sens_tsk_%s_clusterstat_contr%d_f%d_v%d.mat',str,icontr,ifoi,v_cnt));
  senssel = find(stats.mask);
  
  if ~tval
    dfa_d1 = nanmean(all_par_cnt(:,:,contrasts(icontr,1))-all_par_cnt(:,:,contrasts(icontr,2)),2);
    r = max([min(abs(dfa_d1)) max(abs(dfa_d1))]);
    pars.scale=[-.01 .01];
  else
    [~,~,~,s]=ttest(all_par_cnt(:,:,contrasts(icontr,1)),all_par_cnt(:,:,contrasts(icontr,2)),'dim',2)
    dfa_d1 = s.tstat;
    pars.scale = [-2.326 2.326];
  end
  
  figure; set(gcf,'color','white');
  
  if any(stats.mask)
    pars.markersize = 20;
    pars.markersel  = senssel;
  else
    pars.markersize = 0;
  end
  
  pars.cbar       = 0;
  pars.linewidth  = 9;
  pars.resolution = 300;
  
  showfield_colormap(dfa_d1,sa.locs_2D,pars);
  
  % add p-value information and colorlimits
  % ------------------------------------------------
  if isfield(stats,'negclusters') && ~isempty(stats.negclusters)
    m = stats.negclusters(1).prob;
    text(-1.05,-0.76,sprintf('- smallest neg. clust.: p = %.3f',m)); clear m
  else
    text(-1.05,-0.76,sprintf('- no neg. clust. found.'));
  end
  
  pars.markersize = 0;
  
  if isfield(stats,'posclusters') && ~isempty(stats.posclusters)
    m = stats.posclusters(1).prob;
    text(-1.05,-0.7,sprintf('- smallest pos. clust.: p = %.3f',m)); clear m
    
  else
    text(-1.05,-0.7,sprintf('- no pos. clust. found.'));
  end
  
  text(-1.05,-0.64,sprintf('- clim = [%.3f +%.3f]',pars.scale(1),pars.scale(2)));
  % ------------------------------------------------
  
  print(gcf,'-djpeg',sprintf('~/pconn_all/plots/pconn_topo_tsk_%s_c%d_f%d_v%d.jpg',str,icontr,ifoi,v_cnt))
  
  clear stats
  
  % ------------------------------------------------
  % PHARMA CONTRAST DURING REST
  % ------------------------------------------------
  
  load(sprintf('~/pconn_all/proc/all_sens_res_%s_clusterstat_contr%d_f%d_v%d.mat',str,icontr,ifoi,v_res));
  senssel = find(stats.mask);
  
  if ~tval
    dfa_d1 = nanmean(all_par_res(:,:,contrasts(icontr,1))-all_par_res(:,:,contrasts(icontr,2)),2);     
    r = max([min(abs(dfa_d1)) max(abs(dfa_d1))])*0.8;
    pars.scale=[-.033 .033];
  else
    [~,~,~,s]=ttest(all_par_res(:,:,contrasts(icontr,1)),all_par_res(:,:,contrasts(icontr,2)),'dim',2)
    dfa_d1 = s.tstat;
    pars.scale = [-2.326 2.326];
  end
  
  figure; set(gcf,'color','white');

  if any(stats.mask)
    pars.markersize = 20;
    pars.markersel  = senssel;
  else
    pars.markersize = 0;
  end
  
  showfield_colormap(dfa_d1,sa.locs_2D,pars);
  
  % add p-value information and colorlimits
  % ------------------------------------------------
  if isfield(stats,'negclusters') && ~isempty(stats.negclusters)
    m = stats.negclusters(1).prob;
    text(-1.05,-0.76,sprintf('- smallest neg. clust.: p = %.3f',m)); clear m
  else
    text(-1.05,-0.76,sprintf('- no neg. clust. found.'));
  end
  
  if isfield(stats,'posclusters') && ~isempty(stats.posclusters)
    m = stats.posclusters(1).prob;
    text(-1.05,-0.7,sprintf('- smallest pos. clust.: p = %.3f',m)); clear m
  else
    text(-1.05,-0.7,sprintf('- no pos. clust. found.'));
  end
  
  text(-1.05,-0.64,sprintf('- clim = [%.3f +%.3f]',pars.scale(1),pars.scale(2)));
  % ------------------------------------------------
  
  print(gcf,'-djpeg100',sprintf('~/pconn_all/plots/pconn_topo_rst_%s_c%d_f%d_v%d.jpg',str,icontr,ifoi,v_cnt))
  
  clear stats senssel
  close all
  
end

% ------------------------------------------------
% PLOT TASK VS REST
% ------------------------------------------------

figure; set(gcf,'color','white');

if ~tval
    dfa_d2 = nanmean(all_par_cnt(:,:,contrasts(icontr,1))-all_par_res(:,:,contrasts(icontr,2)),2);     
    r = max([min(abs(dfa_d1)) max(abs(dfa_d1))]);
    pars.scale=[-r r];
  else
    [~,~,~,s]=ttest(nanmean(all_par_cnt,3),nanmean(all_par_res,3),'dim',2)
    dfa_d2 = s.tstat;
    pars.scale = [-2.326 2.326];
  end

load(sprintf('~/pconn_all/proc/all_sens_%s_clusterstat_cnt-rest_f%d_v%d.mat',str,ifoi,v_cnt));
senssel = find(stats.mask);

if any(stats.mask)
  pars.markersize = 20;
  pars.markersel  = senssel;
else
  pars.markersize = 0;
end

showfield_colormap(dfa_d2,sa.locs_2D,pars);

% add p-value information and colorlimits
% ------------------------------------------------
if isfield(stats,'negclusters') && ~isempty(stats.negclusters)
  m = stats.negclusters(1).prob;
  text(-1.05,-0.76,sprintf('- smallest neg. clust.: p = %.3f',m)); clear m
else
  text(-1.05,-0.76,sprintf('- no neg. clust. found.'));
end

if isfield(stats,'posclusters') && ~isempty(stats.posclusters)
  m = stats.posclusters(1).prob;
  text(-1.05,-0.7,sprintf('- smallest pos. clust.: p = %.3f',m)); clear m
else
  text(-1.05,-0.7,sprintf('- no pos. clust. found.'));
end

text(-1.05,-0.64,sprintf('- clim = [%.3f +%.3f]',pars.scale(1),pars.scale(2)));
% ------------------------------------------------

print(gcf,'-djpeg',sprintf('~/pconn_all/plots/pconn_all_topo_%s_tsk-rst_f%d_v%d.jpg',str,ifoi,v_cnt))

%% PHARMA DIFFERENCE BETWEEN CONDITIONS 

tval = 1;
ifoi      = 4;
str       = 'pow';

if strcmp(str,'dfa')
  all_par_res = dfa_res_all(:,:,:,ifoi);
  all_par_cnt = dfa_cnt_all(:,:,:,ifoi);
elseif strcmp(str,'var')
  all_par_res = var_res_all(:,:,:,ifoi);
  all_par_cnt = var_cnt_all(:,:,:,ifoi);
elseif strcmp(str,'cvar')
  all_par_res = cvar_res_all(:,:,:,ifoi);
  all_par_cnt = cvar_cnt_all(:,:,:,ifoi);
elseif strcmp(str,'pow')
  all_par_res = log10(pow_res_all(:,:,:,ifoi));
  all_par_cnt = log10(pow_cnt_all(:,:,:,ifoi));
end


contrasts = [2 1; 3 1; 2 3];

for icontr = 1 : 3

  load(sprintf('~/pconn_all/proc/all_sens_%s_clusterstat_diffdiff_c%d_f%d_v%d.mat',str,icontr,ifoi,v_res));
  senssel = find(stats.mask);

  figure; set(gcf,'color','white');

  pars = [];
  pars.resolution = 300;
  pars.linewidth  = 9;
  pars.cbar       = 0;

  dfa_d1    = all_par_cnt(:,:,contrasts(icontr,1))-all_par_cnt(:,:,contrasts(icontr,2));
  dfa_d2    = all_par_res(:,:,contrasts(icontr,1))-all_par_res(:,:,contrasts(icontr,2));
  [~,~,~,d] = ttest(dfa_d1,dfa_d2,'dim',2);
  d         = d.tstat;

  pars.scale=[-2.326 2.326];



  if any(stats.mask)
    pars.markersize = 20;
    pars.markersel  = senssel;
  else
    pars.markersize = 0;
  end

  showfield_colormap(d,sa.locs_2D,pars);
  
  % add p-value information and colorlimits
  % ------------------------------------------------
  if isfield(stats,'negclusters') && ~isempty(stats.negclusters)
    m = stats.negclusters(1).prob;
    text(-1.05,-0.76,sprintf('- smallest neg. clust.: p = %.3f',m)); clear m
  else
    text(-1.05,-0.76,sprintf('- no neg. clust. found.'));
  end
  
  if isfield(stats,'posclusters') && ~isempty(stats.posclusters)
    m = stats.posclusters(1).prob;
    text(-1.05,-0.7,sprintf('- smallest pos. clust.: p = %.3f',m)); clear m
  else
    text(-1.05,-0.7,sprintf('- no pos. clust. found.'));
  end
  
  text(-1.05,-0.64,sprintf('- clim = [%.3f +%.3f]',pars.scale(1),pars.scale(2)));
  % ------------------------------------------------

  print(gcf,'-djpeg100',sprintf('~/pconn_all/plots/pconn_sens_clusterstat_%s_diffdiff_c%d_f%d_v%d.jpg',str,icontr,ifoi,v_cnt))
  clear stats
  
end

%% PLOT BAR GRAPHS AVERAGED ACROSS CHANNELS

ifoi      = 4;
str       = 'dfa';
contrasts = [2 1; 3 1; 2 3];

if strcmp(str,'dfa')
  all_par_res = dfa_res_all(:,:,:,ifoi);
  all_par_cnt = dfa_cnt_all(:,:,:,ifoi);
elseif strcmp(str,'var')
  all_par_res = var_res_all(:,:,:,ifoi);
  all_par_cnt = var_cnt_all(:,:,:,ifoi);
elseif strcmp(str,'cvar')
  all_par_res = cvar_res_all(:,:,:,ifoi);
  all_par_cnt = cvar_cnt_all(:,:,:,ifoi);
elseif strcmp(str,'pow')
  all_par_res = log10(pow_res_all(:,:,:,ifoi));
  all_par_cnt = log10(pow_cnt_all(:,:,:,ifoi));
end

par = squeeze(nanmean(all_par_res));


figure; set(gcf,'color','white'); hold on

m = mean(par);
s = std(par)/sqrt(size(par,1));

[~,p]=ttest(par(:,2),par(:,1));

plot(ones(1,size(par,1)),par(:,1),'rx','markersize',10);
plot(2*ones(1,size(par,1)),par(:,2),'bo','markersize',10);
plot(3*ones(1,size(par,1)),par(:,3),'md','markersize',10);

line([1.2 1.4],[m(1) m(1)],'linewidth',2,'color','r')
line([2.2 2.4],[m(2) m(2)],'linewidth',2,'color','b')
line([3.2 3.4],[m(3) m(3)],'linewidth',2,'color','m')

line([1.3 1.3],[m(1)-s(1) m(1)+s(1)],'linewidth',3,'color','r')
line([2.3 2.3],[m(2)-s(2) m(2)+s(2)],'linewidth',3,'color','b')
line([3.3 3.3],[m(3)-s(3) m(3)+s(3)],'linewidth',3,'color','m')

axis([-1 5 0.5 0.90]);

text(3,mean(m),sprintf('t-test: p = %.3f',p));

print(gcf,'-depsc2',sprintf('~/pconn_all/plots/pconn_topo_avgchanbar_%s_f%d_v%d.eps',str,ifoi,v_cnt))



%% COMPARE POWER IN SIGNIFICANT CLUSTERS

ifoi = 3;
str = 'dfa';

contrasts = [2 1; 3 1; 2 3];

if strcmp(str,'dfa')
  all_par_res = dfa_res_all(:,:,:,ifoi);
  all_par_cnt = dfa_cnt_all(:,:,:,ifoi);
elseif strcmp(str,'var')
  all_par_res = var_res_all(:,:,:,ifoi);
  all_par_cnt = var_cnt_all(:,:,:,ifoi);
elseif strcmp(str,'cvar')
  all_par_res = cvar_res_all(:,:,:,ifoi);
  all_par_cnt = cvar_cnt_all(:,:,:,ifoi);
elseif strcmp(str,'pow')
  all_par_res = log10(pow_res_all(:,:,:,ifoi));
  all_par_cnt = log10(pow_cnt_all(:,:,:,ifoi));
end

for icontr = 1 : 3
  
  
  load(sprintf('~/pconn_all/proc/all_sens_tsk_%s_clusterstat_contr%d_f%d_v%d.mat',str,icontr,ifoi,v_cnt));
  
  if any(stats.mask)
    
    senssel = find(stats.mask);
    
    figure; set(gcf,'color','white'); hold on;
    
    par = squeeze(nanmean(pow_cnt_all(senssel,:,:,ifoi),1));
    
    s = std(par,[],1)/sqrt(size(par,1));
    m = mean(par,1);
    
    figure; set(gcf,'color','white'); hold on;
    
    bar([1.0],m(1),.1,'facecolor','r','edgecolor','none');
    bar([1.2],m(2),.1,'facecolor','b','edgecolor','none');
    
    line([1 1],[m(1)-s(1) m(1)+s(1)],'linewidth',5)
    line([1.2 1.2],[m(2)-s(2) m(2)+s(2)],'linewidth',5)
    
    axis([0.8 1.4 min(max(m-s))-2*min(max(m-s)) max(max(m+s))+2*min(max(m-s))]);
    
    set(gca,'TickDir','out','XTick',[1 1.2],'XTickLabel',['P';'A']);
    
    [~,p] = ttest(par(:,contrasts(icontr,1)),par(:,contrasts(icontr,2)));
    
    text(1.1,max(max(m+s))+1.8*min(max(m-s)),sprintf('t-test: p = %.3f',p)); clear p
    
  else
    figure; set(gcf,'color','white'); hold on;
  end
  
  print(gcf,'-depsc2',sprintf('~/pconn_all/plots/pconn_clustersens_pow_%s_rst_c%d_f%d_v%d.eps',str,icontr,ifoi,v_cnt))
  
  clear stats m s par
  
  
  load(sprintf('~/pconn_all/proc/all_sens_res_%s_clusterstat_contr%d_f%d_v%d.mat',str,icontr,ifoi,v_cnt));
  
  if any(stats.mask)
    
    senssel = find(stats.mask);
    
    figure; set(gcf,'color','white'); hold on;
    
    par = squeeze(nanmean(pow_res_all(senssel,:,:,ifoi),1));
    
    s = std(par,[],1)/sqrt(size(par,1));
    m = mean(par,1);
    
    figure; set(gcf,'color','white'); hold on;
    
    bar([1.0],m(1),.1,'facecolor','r','edgecolor','none');
    bar([1.2],m(2),.1,'facecolor','b','edgecolor','none');
    
    line([1 1],[m(1)-s(1) m(1)+s(1)],'linewidth',5)
    line([1.2 1.2],[m(2)-s(2) m(2)+s(2)],'linewidth',5)
    
    axis([0.8 1.4 min(max(m-s))-2*min(max(m-s)) max(max(m+s))+2*min(max(m-s))]);
    
    set(gca,'TickDir','out','XTick',[1 1.2],'XTickLabel',['P';'A']);
    
    [~,p] = ttest(par(:,contrasts(icontr,1)),par(:,contrasts(icontr,2)));
    
    text(1.1,max(max(m+s))+1.8*min(max(m-s)),sprintf('t-test: p = %.3f',p)); clear p
    
  else
    figure; set(gcf,'color','white'); hold on;
  end
  
  print(gcf,'-depsc2',sprintf('~/pconn_all/plots/pconn_clustersens_pow_%s_rst_c%d_f%d_v%d.eps',str,icontr,ifoi,v_cnt))
  
  clear stats m s par
  
end

close all





















%% COMPARE MEASURES
% -------------------------------------

ifoi = 1;

for isubj = 1 : 18
  
  p_dfa  = nanmean(dfa_res_all(:,isubj,:,ifoi),3);
  p_var  = nanmean(var_res_all(:,isubj,:,ifoi),3);
  p_cvar = nanmean(cvar_res_all(:,isubj,:,ifoi),3);
  
  r_res(1,isubj) = corr(p_dfa,p_var);
  r_res(2,isubj) = corr(p_dfa,p_cvar);
  r_res(3,isubj) = corr(p_cvar,p_var);
  
  p_dfa  = nanmean(dfa_cnt_all(:,isubj,:,ifoi),3);
  p_var  = nanmean(var_cnt_all(:,isubj,:,ifoi),3);
  p_cvar = nanmean(cvar_cnt_all(:,isubj,:,ifoi),3);
  
  r_cnt(1,isubj) = corr(p_dfa,p_var);
  r_cnt(2,isubj) = corr(p_dfa,p_cvar);
  r_cnt(3,isubj) = corr(p_cvar,p_var);
  
  
  
end

for i = 1 : 3
  
  [~,p_res(i)] = ttest(r_res(:,i));
  [~,p_cnt(i)] = ttest(r_cnt(:,i));
  
end

%%
ifoi = 1;

clear r_res r_cnt

for isens = 1 : 268
  
  p_dfa  = nanmean(dfa_res_all(isens,:,:,ifoi),3);
  p_var  = nanmean(var_res_all(isens,:,:,ifoi),3);
  p_cvar = nanmean(cvar_res_all(isens,:,:,ifoi),3);
  
  r_res(1,isens) = corr(p_dfa',p_var');
  r_res(2,isens) = corr(p_dfa',p_cvar');
  r_res(3,isens) = corr(p_cvar',p_var');
  
  p_dfa  = nanmean(dfa_cnt_all(isens,:,:,ifoi),3);
  p_var  = nanmean(var_cnt_all(isens,:,:,ifoi),3);
  p_cvar = nanmean(cvar_cnt_all(isens,:,:,ifoi),3);
  
  r_cnt(1,isens) = corr(p_dfa',p_var');
  r_cnt(2,isens) = corr(p_dfa',p_cvar');
  r_cnt(3,isens) = corr(p_cvar',p_var');
  
end

figure; set(gcf,'color','white');
for i = 1 : 3
  
  subplot(3,2,i*2-1)
  
  pars.scale=[-1 1];
  pars.cbar = 0;
  pars.markersize = 0;
  pars.linewidth = 9;
  pars.resolution = 300;
  
  showfield_colormap(r_res(i,:),sa.locs_2D,pars);
  colormap(parula)
  
  subplot(3,2,i*2)
  
  showfield_colormap(r_cnt(i,:),sa.locs_2D,pars);
  colormap(parula)
  
end

%%