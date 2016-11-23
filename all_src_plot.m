%% COMPUTE WHOLE BRAIN DFA IN SOURCE SPACE
% all_src_plot

clear all

% v = 1: cortex, eloreta, v = 2: cortex, lcmv, v = 3: coarse, eloreta
% str1 = dfa, str2 = amp, str3 = cvar, str4 = var

foi = [5];
STR = [1:4];

for v = [2]
  
  % --------------------------------------------------------
  % VERSION 1
  % --------------------------------------------------------
  v         = v;
  v_stat    = v;
  v_rawdata = 2;
  SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32];
  % --------------------------------------------------------
  
  allstr    = {'dfa';'amp';'cvar';'var'};
  gridsize  = 'cortex';
  contrasts = [2 1; 3 1; 2 3];
  viewdir   = [0.001 -0 1; 0 -.5 0; 1 0 0; -.5 0 0; -.5 0 0;1 0 0 ];

  outdir    = '/home/tpfeffer/pconn_all/proc/';
  plotdir   = '/home/tpfeffer/pconn_all/plots/';
  
  addpath /home/gnolte/meg_toolbox/toolbox/
  addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
  addpath /home/gnolte/meg_toolbox/toolbox_nightly/
  addpath /home/gnolte/meg_toolbox/meg/
  addpath ~/pconn/matlab/
  addpath ~/Documents/MATLAB/Colormaps/Colormaps' (5)'/Colormaps/
  
    for istr = STR
      for ifoi = foi

%       if ~exist(sprintf([outdir 'all_src_plot_s%d_f%d_v%d_processing.txt'],istr,ifoi,v))
%         system(['touch ' outdir sprintf('all_src_plot_s%d_f%d_v%d_processing.txt',istr,ifoi,v)]);
%       else
%         continue
%       end
      
      fprintf('Processing v%d s%d ...\n',v,istr);
      
      str = allstr{istr};
       
      % READ TEMPLATE MRI STUFF
      % ----------------
      load sa_meg_template;
      
      if strcmp(gridsize,'coarse')
        grid  = sa_meg_template.grid_coarse;
      elseif strcmp(gridsize,'cortex')
        grid  = sa_meg_template.grid_cortex3000;
        g1 = sa_meg_template.grid_cortex3000;
        g2 = sa_meg_template.cortex10K.vc;
      elseif strcmp(gridsize,'xcoarse')
        grid  = sa_meg_template.grid_xcoarse;
      end
      
      mri   = sa_meg_template.mri;
      vc    = sa_meg_template.vc;
      dd    = .75;
      g1    = sa_meg_template.grid_cortex3000;
      % ----------------
      
      %% READ IN DATA
      fprintf('Reading data ...\n');
      ord   = pconn_randomization;

      for isubj = SUBJLIST
%         fprintf('Reading data s%d ...\n',isubj);
        for m = 1 : 3

          im = find(ord(isubj,:)==m);

          load(sprintf(['~/pconn_cnt/proc/dfa/pconn_cnt_src_dfa_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,v));

          if strcmp(str,'dfa')
            par_all_cnt(:,m,isubj)  = nanmean(par.dfa,2);
          elseif strcmp(str,'var')
            par_all_cnt(:,m,isubj)  = nanmean(par.var,2);
          elseif strcmp(str,'cvar')
            par_all_cnt(:,m,isubj) = nanmean(par.cvar,2);
          elseif strcmp(str,'amp')
            par_all_cnt(:,m,isubj)  = nanmean(par.amp,2); clear par
          end
          %
          load(sprintf(['~/pconn/proc/dfa/pconn_src_dfa_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,v));

          if strcmp(str,'dfa')
            par_all_res(:,m,isubj)  = nanmean(par.dfa,2);
          elseif strcmp(str,'var')
            par_all_res(:,m,isubj)  = nanmean(par.var,2);
          elseif strcmp(str,'cvar')
            par_all_res(:,m,isubj)  = nanmean(par.cvar,2);
          elseif strcmp(str,'amp')
            par_all_res(:,m,isubj)  = nanmean(par.amp,2);
          end

          clear par

        end
      end


      par_all_res  = par_all_res(:,:,SUBJLIST);
      par_all_cnt  = par_all_cnt(:,:,SUBJLIST);

      fprintf('Reading data ... Done!\n');

      %% PHARMA COMPARISON
      for mask = 0 : 1
        
         fprintf('Pharma comparison mask%d ...\n',mask);

    
        pl  = inferno(128);
        vi = viridis(128);
        cmap  = [vi(end:-1:10,:); pl(10:end,:)];
        
        len =  size(cmap,1)/2-15: size(cmap,1)/2+15;
        cmap(len,:) = repmat(1,[length(len) 3]);
        
        % cmap =   [0.941176474094391 0.941176474094391 0.941176474094391;0.935294091701508 0.935294091701508 0.935294091701508;0.929411768913269 0.929411768913269 0.929411768913269;0.923529386520386 0.923529386520386 0.923529386520386;0.917647063732147 0.917647063732147 0.917647063732147;0.911764681339264 0.911764681339264 0.911764681339264;0.905882358551025 0.905882358551025 0.905882358551025;0.899999976158142 0.899999976158142 0.899999976158142;0.867568910121918 0.852902054786682 0.877736866474152;0.835137844085693 0.805804193019867 0.855473756790161;0.802706778049469 0.758706271648407 0.833210647106171;0.770275771617889 0.711608350276947 0.81094753742218;0.737844705581665 0.664510428905487 0.78868442773819;0.705413639545441 0.617412567138672 0.766421318054199;0.672982573509216 0.570314645767212 0.744158208370209;0.647758424282074 0.533682942390442 0.72684246301651;0.622534275054932 0.497051239013672 0.709526658058167;0.597310066223145 0.460419565439224 0.692210912704468;0.572085916996002 0.423787862062454 0.674895167350769;0.54432600736618 0.383473604917526 0.655838668346405;0.516566097736359 0.343159347772598 0.636782169342041;0.488806158304214 0.30284509062767 0.617725670337677;0.461046248674393 0.262530833482742 0.598669171333313;0.433286309242249 0.222216591238976 0.579612672328949;0.405526399612427 0.181902334094048 0.560556173324585;0.38470646739006 0.151666641235352 0.546263813972473;0.363886535167694 0.121430955827236 0.531971454620361;0.343066573143005 0.0911952629685402 0.51767909526825;0.322246640920639 0.0609595701098442 0.503386676311493;0.301426708698273 0.0307238809764385 0.489094316959381;0.280606776475906 0.000488189951283857 0.474801957607269;0.293446481227875 0.00417250022292137 0.461163550615311;0.306286215782166 0.00785681046545506 0.447525143623352;0.319125920534134 0.0115411207079887 0.433886736631393;0.331965625286102 0.0152254309505224 0.420248329639435;0.344805359840393 0.0189097411930561 0.406609892845154;0.357645064592361 0.0225940514355898 0.392971485853195;0.370484799146652 0.0262783616781235 0.379333078861237;0.383324503898621 0.0299626719206572 0.365694671869278;0.396164208650589 0.0336469821631908 0.352056264877319;0.40900394320488 0.0373312942683697 0.338417857885361;0.421843647956848 0.0410156026482582 0.324779450893402;0.434683352708817 0.044699914753437 0.311141043901443;0.447523087263107 0.0483842231333256 0.297502636909485;0.460362792015076 0.0520685352385044 0.283864229917526;0.473202496767044 0.0557528436183929 0.270225793123245;0.486042231321335 0.0594371557235718 0.256587386131287;0.498881936073303 0.0631214678287506 0.242948994040489;0.511721670627594 0.0668057799339294 0.229310572147369;0.52456134557724 0.0704900845885277 0.215672165155411;0.537401080131531 0.0741743966937065 0.202033758163452;0.550240814685822 0.0778587087988853 0.188395351171494;0.563080549240112 0.0815430209040642 0.174756944179535;0.575920224189758 0.0852273255586624 0.161118522286415;0.588759958744049 0.0889116376638412 0.147480115294456;0.60159969329834 0.0925959497690201 0.133841708302498;0.614439368247986 0.0962802618741989 0.120203301310539;0.627279102802277 0.0999645665287971 0.106564886868;0.640118837356567 0.103648878633976 0.0929264798760414;0.652958512306213 0.107333190739155 0.0792880654335022;0.665798246860504 0.111017502844334 0.0656496584415436;0.678637981414795 0.114701807498932 0.0520112477242947;0.691477656364441 0.118386119604111 0.0383728370070457;0.704317390918732 0.12207043170929 0.024734428152442;0.70988517999649 0.127976059913635 0;0.714490175247192 0.141817703843117 0;0.719095170497894 0.155659362673759 0;0.723700165748596 0.169501006603241 0;0.728305160999298 0.183342665433884 0;0.73291015625 0.197184309363365 0;0.737515151500702 0.211025953292847 0;0.742120146751404 0.224867612123489 0;0.746725142002106 0.238709256052971 0;0.751330137252808 0.252550899982452 0;0.75593513250351 0.266392558813095 0;0.760540127754211 0.280234217643738 0;0.765145123004913 0.294075846672058 0;0.769750118255615 0.307917505502701 0;0.774355113506317 0.321759164333344 0;0.778960108757019 0.335600793361664 0;0.783565163612366 0.349442452192307 0;0.788170158863068 0.363284111022949 0;0.79277515411377 0.377125769853592 0;0.797380149364471 0.390967398881912 0;0.801985144615173 0.404809057712555 0;0.806590139865875 0.418650716543198 0;0.811195135116577 0.432492345571518 0;0.815800130367279 0.446334004402161 0;0.820405125617981 0.460175663232803 0;0.825010120868683 0.474017292261124 0;0.829615116119385 0.487858951091766 0;0.834220111370087 0.501700580120087 0;0.838825106620789 0.515542268753052 0;0.84343010187149 0.529383897781372 0;0.848035097122192 0.543225526809692 0;0.852640092372894 0.557067215442657 0;0.857245087623596 0.570908844470978 0;0.861850082874298 0.584750533103943 0;0.866455078125 0.598592162132263 0;0.871060073375702 0.612433791160583 0;0.875665068626404 0.626275479793549 0;0.880270063877106 0.640117108821869 0;0.884875059127808 0.653958737850189 0;0.88948005437851 0.667800426483154 0;0.894085049629211 0.681642055511475 0;0.898690044879913 0.695483684539795 0;0.903295040130615 0.70932537317276 0;0.907900035381317 0.72316700220108 0;0.912505030632019 0.737008631229401 0;0.917110025882721 0.750850319862366 0;0.921715021133423 0.764691948890686 0;0.926320016384125 0.778533577919006 0;0.930925071239471 0.792375266551971 0;0.935530066490173 0.806216895580292 0;0.940135061740875 0.820058524608612 0;0.944740056991577 0.833900213241577 0;0.949345052242279 0.847741842269897 0;0.953950047492981 0.861583530902863 0;0.958555042743683 0.875425159931183 0;0.963160037994385 0.889266788959503 0;0.967765033245087 0.903108477592468 0;0.972370028495789 0.916950106620789 0;0.97697502374649 0.930791735649109 0;0.981580018997192 0.944633424282074 0;0.986185014247894 0.958475053310394 0;0.990790009498596 0.972316682338715 0;0.995395004749298 0.98615837097168 0;1 1 0];
        clear stats
        
        for icontr = 1:2
          
          % -------------------------------------------------------
          %   TASK
          % -------------------------------------------------------
          
          d = squeeze(nanmean(par_all_cnt(:,contrasts(icontr,1),:),3))-squeeze(nanmean(par_all_cnt(:,contrasts(icontr,2),:),3));
          %
          [~,~,~,tmp] = ttest(par_all_cnt(:,contrasts(icontr,1),:),par_all_cnt(:,contrasts(icontr,2),:),'dim',3);
          d = tmp.tstat; clear tmp
          %       d(abs(d) < 2.5])=eps;
          
          if mask && exist(sprintf('~/pconn_all/proc/all_src_clusterstat_tsk_%s_c%d_f%d_v%d.mat',str,icontr,ifoi,v_stat))
            load(sprintf('~/pconn_all/proc/all_src_clusterstat_tsk_%s_c%d_f%d_v%d.mat',str,icontr,ifoi,v_stat))
            if isfield(stats,'mask')
              d(~stats.mask)=eps;
            else
              d(1:end) = eps;
            end
          end
          
          par_interp = spatfiltergauss(d,g1,dd,g2);
          
          para = [] ;
          para.colorlimits = [-1.96 1.96];
          
          % PLOT RESULTS
          tp_showsource(par_interp,cmap,sa_meg_template,para);
          
          print(gcf,'-djpeg100',sprintf('~/pconn_all/plots/all_src_tsk_%s_mask%d_c%d_f%d_v%d.jpg',str,mask,icontr,ifoi,v_stat))
          
          clear para r d stats
          
          % -------------------------------------------------------
          % RESTING STATE
          % -------------------------------------------------------
          
          d = squeeze(nanmean(par_all_res(:,contrasts(icontr,1),:),3))-squeeze(nanmean(par_all_res(:,contrasts(icontr,2),:),3));
          
          [~,~,~,tmp] = ttest(par_all_res(:,contrasts(icontr,1),:),par_all_res(:,contrasts(icontr,2),:),'dim',3);
          d = tmp.tstat; clear tmp
          
          if mask && exist(sprintf('~/pconn_all/proc/all_src_clusterstat_rst_%s_c%d_f%d_v%d.mat',str,icontr,ifoi,v_stat))
            load(sprintf('~/pconn_all/proc/all_src_clusterstat_rst_%s_c%d_f%d_v%d.mat',str,icontr,ifoi,v_stat))
            if isfield(stats,'mask')
              d(~stats.mask)=eps;
            else
              d(1:end) = eps;
            end
          end
          
          par_interp = spatfiltergauss(d,g1,dd,g2);
          
          para = [] ;
          para.colorlimits = [-1.96 1.96];
          % PLOT RESULTS
          tp_showsource(par_interp,cmap,sa_meg_template,para);

          print(gcf,'-djpeg100',sprintf('~/pconn_all/plots/all_src_rst_%s_mask%d_c%d_f%d_v%d.jpg',str,mask,icontr,ifoi,v_stat))
          clear stats
        end
        
      end
      end
    end
end

error('!')
      
%% SCATTER PLOTS (PLACEBO VS. DRUG)
str = 'var';
ifoi = 1;
v= 2;

ord   = pconn_randomization;

for isubj = SUBJLIST
  fprintf('Reading data s%d ...\n',isubj);
  for m = 1 : 3
    
    im = find(ord(isubj,:)==m);
    
    load(sprintf(['~/pconn_cnt/proc/dfa/pconn_cnt_src_dfa_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,v));
    
    if strcmp(str,'dfa')
      par_all_cnt(:,m,isubj)  = nanmean(par.dfa,2);
    elseif strcmp(str,'var')
      par_all_cnt(:,m,isubj)  = nanmean(par.var,2);
    elseif strcmp(str,'cvar')
      par_all_cnt(:,m,isubj) = nanmean(par.cvar,2);
    elseif strcmp(str,'amp')
      par_all_cnt(:,m,isubj)  = nanmean(par.amp,2); clear par
    end
    %
    load(sprintf(['~/pconn/proc/dfa/pconn_src_dfa_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,v));
    
    if strcmp(str,'dfa')
      par_all_res(:,m,isubj)  = nanmean(par.dfa,2);
    elseif strcmp(str,'var')
      par_all_res(:,m,isubj)  = nanmean(par.var,2);
    elseif strcmp(str,'cvar')
      par_all_res(:,m,isubj)  = nanmean(par.cvar,2);
    elseif strcmp(str,'amp')
      par_all_res(:,m,isubj)  = nanmean(par.amp,2);
    end
    
    clear par
    
  end
end


par_all_res  = par_all_res(:,:,SUBJLIST);
par_all_cnt  = par_all_cnt(:,:,SUBJLIST);

%%

figure; set(gcf,'color','white');

scatter(log10(nanmean(par_all_res(:,1,:),3)),log10(nanmean(par_all_res(:,3,:),3)),50,'facecolor','k','markeredgecolor','w')
tp_editplots
line([-33.8 -32.6],[-33.8 -32.6],'color','k')
      
      
      
      
      
      
      
      
      %% DATA FOR ALL CONDITIONS/SUBJECTS
%       
%       para = [] ;
%       para.colorlimits = [0.6 0.85];
%       cnt = 0;
%       para.myviewdir = viewdir(1,:);
%       a = sa_meg_template.cortex10K;
%                   cmap =   [0.941176474094391 0.941176474094391 0.941176474094391;0.935294091701508 0.935294091701508 0.935294091701508;0.929411768913269 0.929411768913269 0.929411768913269;0.923529386520386 0.923529386520386 0.923529386520386;0.917647063732147 0.917647063732147 0.917647063732147;0.911764681339264 0.911764681339264 0.911764681339264;0.905882358551025 0.905882358551025 0.905882358551025;0.899999976158142 0.899999976158142 0.899999976158142;0.867568910121918 0.852902054786682 0.877736866474152;0.835137844085693 0.805804193019867 0.855473756790161;0.802706778049469 0.758706271648407 0.833210647106171;0.770275771617889 0.711608350276947 0.81094753742218;0.737844705581665 0.664510428905487 0.78868442773819;0.705413639545441 0.617412567138672 0.766421318054199;0.672982573509216 0.570314645767212 0.744158208370209;0.647758424282074 0.533682942390442 0.72684246301651;0.622534275054932 0.497051239013672 0.709526658058167;0.597310066223145 0.460419565439224 0.692210912704468;0.572085916996002 0.423787862062454 0.674895167350769;0.54432600736618 0.383473604917526 0.655838668346405;0.516566097736359 0.343159347772598 0.636782169342041;0.488806158304214 0.30284509062767 0.617725670337677;0.461046248674393 0.262530833482742 0.598669171333313;0.433286309242249 0.222216591238976 0.579612672328949;0.405526399612427 0.181902334094048 0.560556173324585;0.38470646739006 0.151666641235352 0.546263813972473;0.363886535167694 0.121430955827236 0.531971454620361;0.343066573143005 0.0911952629685402 0.51767909526825;0.322246640920639 0.0609595701098442 0.503386676311493;0.301426708698273 0.0307238809764385 0.489094316959381;0.280606776475906 0.000488189951283857 0.474801957607269;0.293446481227875 0.00417250022292137 0.461163550615311;0.306286215782166 0.00785681046545506 0.447525143623352;0.319125920534134 0.0115411207079887 0.433886736631393;0.331965625286102 0.0152254309505224 0.420248329639435;0.344805359840393 0.0189097411930561 0.406609892845154;0.357645064592361 0.0225940514355898 0.392971485853195;0.370484799146652 0.0262783616781235 0.379333078861237;0.383324503898621 0.0299626719206572 0.365694671869278;0.396164208650589 0.0336469821631908 0.352056264877319;0.40900394320488 0.0373312942683697 0.338417857885361;0.421843647956848 0.0410156026482582 0.324779450893402;0.434683352708817 0.044699914753437 0.311141043901443;0.447523087263107 0.0483842231333256 0.297502636909485;0.460362792015076 0.0520685352385044 0.283864229917526;0.473202496767044 0.0557528436183929 0.270225793123245;0.486042231321335 0.0594371557235718 0.256587386131287;0.498881936073303 0.0631214678287506 0.242948994040489;0.511721670627594 0.0668057799339294 0.229310572147369;0.52456134557724 0.0704900845885277 0.215672165155411;0.537401080131531 0.0741743966937065 0.202033758163452;0.550240814685822 0.0778587087988853 0.188395351171494;0.563080549240112 0.0815430209040642 0.174756944179535;0.575920224189758 0.0852273255586624 0.161118522286415;0.588759958744049 0.0889116376638412 0.147480115294456;0.60159969329834 0.0925959497690201 0.133841708302498;0.614439368247986 0.0962802618741989 0.120203301310539;0.627279102802277 0.0999645665287971 0.106564886868;0.640118837356567 0.103648878633976 0.0929264798760414;0.652958512306213 0.107333190739155 0.0792880654335022;0.665798246860504 0.111017502844334 0.0656496584415436;0.678637981414795 0.114701807498932 0.0520112477242947;0.691477656364441 0.118386119604111 0.0383728370070457;0.704317390918732 0.12207043170929 0.024734428152442;0.70988517999649 0.127976059913635 0;0.714490175247192 0.141817703843117 0;0.719095170497894 0.155659362673759 0;0.723700165748596 0.169501006603241 0;0.728305160999298 0.183342665433884 0;0.73291015625 0.197184309363365 0;0.737515151500702 0.211025953292847 0;0.742120146751404 0.224867612123489 0;0.746725142002106 0.238709256052971 0;0.751330137252808 0.252550899982452 0;0.75593513250351 0.266392558813095 0;0.760540127754211 0.280234217643738 0;0.765145123004913 0.294075846672058 0;0.769750118255615 0.307917505502701 0;0.774355113506317 0.321759164333344 0;0.778960108757019 0.335600793361664 0;0.783565163612366 0.349442452192307 0;0.788170158863068 0.363284111022949 0;0.79277515411377 0.377125769853592 0;0.797380149364471 0.390967398881912 0;0.801985144615173 0.404809057712555 0;0.806590139865875 0.418650716543198 0;0.811195135116577 0.432492345571518 0;0.815800130367279 0.446334004402161 0;0.820405125617981 0.460175663232803 0;0.825010120868683 0.474017292261124 0;0.829615116119385 0.487858951091766 0;0.834220111370087 0.501700580120087 0;0.838825106620789 0.515542268753052 0;0.84343010187149 0.529383897781372 0;0.848035097122192 0.543225526809692 0;0.852640092372894 0.557067215442657 0;0.857245087623596 0.570908844470978 0;0.861850082874298 0.584750533103943 0;0.866455078125 0.598592162132263 0;0.871060073375702 0.612433791160583 0;0.875665068626404 0.626275479793549 0;0.880270063877106 0.640117108821869 0;0.884875059127808 0.653958737850189 0;0.88948005437851 0.667800426483154 0;0.894085049629211 0.681642055511475 0;0.898690044879913 0.695483684539795 0;0.903295040130615 0.70932537317276 0;0.907900035381317 0.72316700220108 0;0.912505030632019 0.737008631229401 0;0.917110025882721 0.750850319862366 0;0.921715021133423 0.764691948890686 0;0.926320016384125 0.778533577919006 0;0.930925071239471 0.792375266551971 0;0.935530066490173 0.806216895580292 0;0.940135061740875 0.820058524608612 0;0.944740056991577 0.833900213241577 0;0.949345052242279 0.847741842269897 0;0.953950047492981 0.861583530902863 0;0.958555042743683 0.875425159931183 0;0.963160037994385 0.889266788959503 0;0.967765033245087 0.903108477592468 0;0.972370028495789 0.916950106620789 0;0.97697502374649 0.930791735649109 0;0.981580018997192 0.944633424282074 0;0.986185014247894 0.958475053310394 0;0.990790009498596 0.972316682338715 0;0.995395004749298 0.98615837097168 0;1 1 0];
% 
%       for isubj = 1:length(SUBJLIST)
%         figure; set(gcf,'color','w');
% 
%         for i = 1 : 3
%           
%           subplot(1,3,i);
%           
%           par_interp = spatfiltergauss(par_all_res(:,i,isubj),g1,dd,g2);
%           pconn_showsurface(a,para,par_interp)
% 
%           colormap(cmap)
% %           drawnow
%         end
%        	print(gcf,'-djpeg100',sprintf('~/pconn_all/plots/all_src_rst_indiv_s%d_%s_mask%d_f%d_v%d.jpg',SUBJLIST(isubj),str,mask,ifoi,v_stat))
% 
%       end
        
%         camlight headlight
        
      
      
      
      %% COMPARISON TASK VS REST
      
%       for mask = 0 : 1
%         
%         fprintf('Context comparison mask%d ...\n',mask);
%         
%         % -------------------------------------------------------
%         % TASK VS REST
%         % -------------------------------------------------------
%         
%         d = squeeze(nanmean(nanmean(par_all_cnt,3),2))-squeeze(nanmean(nanmean(par_all_res,3),2));
%         [~,~,~,tmp] = ttest(nanmean(par_all_cnt,2),nanmean(par_all_res,2),'dim',3);
%         d = tmp.tstat; clear tmp
%         
%         if mask && exist(sprintf('~/pconn_all/proc/all_src_clusterstat_tsk-rst_%s_c%d_f%d_v%d.mat',str,3,ifoi,v_stat))
%           load(sprintf('~/pconn_all/proc/all_src_clusterstat_tsk-rst_%s_c%d_f%d_v%d.mat',str,3,ifoi,v_stat))
%           if isfield(stats,'mask')
%               d(~stats.mask)=eps;
%             else
%               d(1:end) = eps;
%             end
%         end
%         
%         par_interp = spatfiltergauss(d,g1,dd,g2);
%         
%         figure; set(gcf,'color','white'); hold on;
%         
%         para = [] ;
%         para.colorlimits = [-1.96 1.96];
%         
%         for iplot = 1 : 6
%           
%           subplot(3,2,iplot)
%           para.myviewdir = viewdir(iplot,:);
%           
%           a = sa_meg_template.cortex10K;
%           
%           if iplot == 5
%             a=cutsurface(sa_meg_template.cortex10K,[mean(sa_meg_template.cortex10K.vc(:,1)) 0 0],[1 0 0]);
%           elseif iplot == 6
%             a=cutsurface(sa_meg_template.cortex10K,[mean(sa_meg_template.cortex10K.vc(:,1)) 0 0],[-1 0 0]);
%           end
%           
%           pconn_showsurface(a,para,par_interp)
%           colormap(cmap)
%           camlight headlight
%           
%         end
%         set(gcf,'paperpositionmode','auto')
%         print(gcf,'-djpeg100',sprintf('~/pconn_all/plots/pconn_src_%s_tsk-rst_mask%d_f%d_v%d.jpg',str,mask,ifoi,v_stat))
%         clear para r d stats
%       end
      
      
      
      %% COMPARE PHARMA ACROSS TASK/REST
      
%       for mask = 0 : 1
        
%         fprintf('Pharma comparison context mask%d ...\n',mask);
% 
%         for icontr = 1 : 2
%           
%           d1 = squeeze(par_all_res(:,contrasts(icontr,1),:)-par_all_res(:,contrasts(icontr,2),:));
%           d2 = squeeze(par_all_cnt(:,contrasts(icontr,1),:)-par_all_cnt(:,contrasts(icontr,2),:));
%           
%           [~,~,~,tmp] = ttest(d2,d1,'dim',2);
%           
%           d = tmp.tstat;
%           
%           if mask && exist(sprintf('~/pconn_all/proc/all_src_clusterstat_diffdiff_%s_c%d_f%d_v%d.mat',str,icontr,ifoi,v_stat))
%             load(sprintf('~/pconn_all/proc/all_src_clusterstat_diffdiff_%s_c%d_f%d_v%d',str,icontr,ifoi,v_stat))
%             if isfield(stats,'mask')
%               d(~stats.mask)=eps;
%             else
%               d(1:end) = eps;
%             end
%           end
%           
%           par_interp = spatfiltergauss(d,g1,dd,g2);
%           
%           figure; set(gcf,'color','white'); hold on;
%           
%           para = [] ;
%           para.colorlimits = [-1.96 1.96];
%           
%           for iplot = 1 : 6
%             
%             subplot(3,2,iplot)
%             para.myviewdir = viewdir(iplot,:);
%             
%             a = sa_meg_template.cortex10K;
%             
%             if iplot == 5
%               a=cutsurface(sa_meg_template.cortex10K,[mean(sa_meg_template.cortex10K.vc(:,1)) 0 0],[1 0 0]);
%             elseif iplot == 6
%               a=cutsurface(sa_meg_template.cortex10K,[mean(sa_meg_template.cortex10K.vc(:,1)) 0 0],[-1 0 0]);
%             end
%             
%             pconn_showsurface(a,para,par_interp)
%             colormap(cmap)
%             camlight headlight
%             
%           end
%           set(gcf,'paperpositionmode','auto')
%           print(gcf,'-djpeg',sprintf('~/pconn_all/plots/pconn_src_%s_tskrstpharm_mask%d_c%d_f%d_v%d.jpg',str,mask,icontr,ifoi,v))
%           %
%         end
        

% end
% %% CORRELATE ALL MEASURES
%
% ifoi = 4;
% str = 'dfa';
%
% viewdir   = [-.5 -.5 .5; .5 -.5 .5; .5 .5 .5; -.5 .5 .5];
%
% for isubj = 1 : length(SUBJLIST)
%
% %   p1        = nanmean(dfa_all_res(:,:,isubj,ifoi),2);
% %   p2        = nanmean(var_all_res(:,:,isubj,ifoi),2);
% %   p3        = nanmean(cvar_all_res(:,:,isubj,ifoi),2);
% %
% %   r1(isubj) = corr(p1,p2);
% %   r2(isubj) = corr(p1,p3);
% %   r3(isubj) = corr(p2,p3);
% %
%   p1        = nanmean(dfa_all_cnt(:,:,isubj,ifoi),2);
%   p2        = nanmean(var_all_cnt(:,:,isubj,ifoi),2);
%   p3        = nanmean(cvar_all_cnt(:,:,isubj,ifoi),2);
%
%   r1(isubj) = corr(p1,p2);
%   r2(isubj) = corr(p1,p3);
%   r3(isubj) = corr(p2,p3);
%
%
% end
%
% zerodist = zeros(1,length(SUBJLIST));
%
% [~,p1] = ttest(r1,zerodist);
% [~,p2] = ttest(r2,zerodist);
% [~,p3] = ttest(r3,zerodist);
%
%
% %% COMPARE MEASURES
%
% ifoi = 4;
%
% p_dfa  = nanmean(nanmean(dfa_all(:,:,:,ifoi),3),2);
% p_var  = nanmean(nanmean(var_all(:,:,:,ifoi),3),2);
% p_cvar = nanmean(nanmean(cvar_all(:,:,:,ifoi),3),2);
%
% r(1) = corr(p_dfa,p_var);
% r(2) = corr(p_dfa,p_cvar);
% r(3) = corr(p_cvar,p_var);
%
% figure; set(gcf,'color','white');
%
% subplot(3,1,1)
%
% scatter(p_dfa,p_var,20,'facecolor','k','markeredgecolor','w')
% xlabel('DFA'); ylabel('Variance')
% axis square
% title(sprintf('r = %.2f',r(1)))
%
% subplot(3,1,2)
%
% scatter(p_dfa,p_cvar,20,'facecolor','k','markeredgecolor','w')
% xlabel('DFA'); ylabel('Coef. of var.')
% axis square
% title(sprintf('r = %.2f',r(2)))
%
% subplot(3,1,3)
%
% scatter(p_var,p_cvar,20,'facecolor','k','markeredgecolor','w')
% xlabel('Variance'); ylabel('Coef. of var.')
% axis square
% title(sprintf('r = %.2f',r(3)))
%
% set(gcf,'Position',[50 50 800 1200])
%
% print(gcf,'-djpeg100',sprintf('~/pconn_cnt/plots/pconn_cnt_src_dfavarcvarcorr_f%d_v%d.jpg',ifoi,v))
%
% close
%
% p_dfa  = nanmean(nanmean(dfa_all_res(:,:,:,ifoi),3),2);
% p_var  = nanmean(nanmean(var_all_res(:,:,:,ifoi),3),2);
% p_cvar = nanmean(nanmean(cvar_all_res(:,:,:,ifoi),3),2);
%
% r(1) = corr(p_dfa,p_var);
% r(2) = corr(p_dfa,p_cvar);
% r(3) = corr(p_cvar,p_var);
%
% figure; set(gcf,'color','white');
%
% subplot(3,1,1)
%
% scatter(p_dfa,p_var,20,'facecolor','k','markeredgecolor','w')
% xlabel('DFA'); ylabel('Variance')
% axis square
% title(sprintf('r = %.2f',r(1)))
%
% subplot(3,1,2)
%
% scatter(p_dfa,p_cvar,20,'facecolor','k','markeredgecolor','w')
% xlabel('DFA'); ylabel('Coef. of var.')
% axis square
% title(sprintf('r = %.2f',r(2)))
%
% subplot(3,1,3)
%
% scatter(p_var,p_cvar,20,'facecolor','k','markeredgecolor','w')
% xlabel('Variance'); ylabel('Coef. of var.')
% axis square
% title(sprintf('r = %.2f',r(3)))
%
% set(gcf,'Position',[50 50 800 1200])
%
% print(gcf,'-djpeg100',sprintf('~/pconn/proc/plots/pconn_src_dfavarcvarcorr_f%d_v%d.jpg',ifoi,v))
%
% close
%
%
%
% %%  CORRELATION OF MEASURES ACROSS CONDITIONS
%
% ifoi = 1;
%
% p1_dfa  = nanmean(nanmean(dfa_all_res(:,:,:,ifoi),3),2);
% p1_var  = nanmean(nanmean(var_all_res(:,:,:,ifoi),3),2);
% p1_cvar = nanmean(nanmean(cvar_all_res(:,:,:,ifoi),3),2);
%
% p2_dfa  = nanmean(nanmean(dfa_all(:,:,:,ifoi),3),2);
% p2_var  = nanmean(nanmean(var_all(:,:,:,ifoi),3),2);
% p2_cvar = nanmean(nanmean(cvar_all(:,:,:,ifoi),3),2);
%
% r(1) = corr(p1_dfa,p2_dfa);
% r(2) = corr(p1_var,p2_var);
% r(3) = corr(p1_cvar,p2_cvar);
%
% figure; set(gcf,'color','white');
%
% subplot(3,1,1)
%
% scatter(p1_dfa,p2_dfa,20,'facecolor','k','markeredgecolor','w')
% xlabel('DFA (rest)'); ylabel('DFA (task)')
% axis square
% title(sprintf('r = %.2f',r(1)))
%
% subplot(3,1,2)
%
% scatter(log10(p1_var),log10(p2_var),20,'facecolor','k','markeredgecolor','w')
% xlabel('Variance (rest)'); ylabel('Variance (task)')
% axis square
% title(sprintf('r = %.2f',r(2)))
%
% subplot(3,1,3)
%
% scatter(p1_cvar,p2_cvar,20,'facecolor','k','markeredgecolor','w')
% xlabel('C. Variance (rest)'); ylabel('C. Variance (task)')
% axis square
% title(sprintf('r = %.2f',r(3)))
%
% set(gcf,'Position',[50 50 800 1200])
%
% print(gcf,'-djpeg100',sprintf('~/pconn/proc/plots/pconn_src_resttaskcorr_f%d_v%d.jpg',ifoi,v))
%
% close
%
% %%
%
% %% SHOW SURFACE
%
% g1 = sa_meg_template.grid_cortex3000;
% g2 = sa_meg_template.cortex10K.vc;
% dd = .75;
% m2 = spatfiltergauss(m,g1,dd,g2);
%
% para.colorlimits = [-0.03 0.03];
%
% figure;  pconn_showsurface(a,[],z2(idx))
%
% %% COMPARE DIFF TASK WITH DIFF PHARM
% mask = 0;
% str  = 'dfa';
% str_cond = 'rst';
% clear p r d
%
% for ifoi = 3 : 3
%
%   contrasts = [2 1];
%   viewdir = [0.001 -0 1; 0 -.5 0; 1 0 0; -.5 0 0; -.5 0 0;1 0 0 ];
%
%   if strcmp(str,'dfa')
%     par_all_res = dfa_all_res(:,:,:,ifoi);
%     par_all_cnt = dfa_all_cnt(:,:,:,ifoi);
%   elseif strcmp(str,'var')
%     par_all_res = var_all_res(:,:,:,ifoi);
%     par_all_cnt = var_all_cnt(:,:,:,ifoi);
%   elseif strcmp(str,'cvar')
%     par_all_res = cvar_all_res(:,:,:,ifoi);
%     par_all_cnt = cvar_all_cnt(:,:,:,ifoi);
%   end
%
%   cmap = cbrewer('div', 'RdBu', 100,'pchip');% colormap(autumn)
%   cmap = cmap(end:-1:1,:);
%
%   for icontr = 1 : 1
%
%     if strcmp(str_cond,'rst')
%       d1 = squeeze(par_all_res(:,contrasts(icontr,1),:)-par_all_res(:,contrasts(icontr,2),:));
%     else
%       d1 = squeeze(par_all_cnt(:,contrasts(icontr,1),:)-par_all_cnt(:,contrasts(icontr,2),:));
%     end
%
%     d2 = squeeze(nanmean(par_all_cnt(:,1,:)-par_all_res(:,1,:),2));
%
%     for isubj = 1 : size(par_all_res,3)
%
%       d(isubj,:) = [mean(d1(:,isubj)) mean(d2(:,isubj))];
%       r(isubj)   = corr(d1(:,isubj),d2(:,isubj));
%
%     end
%
%     [tr tp] = corr(d(:,1),d(:,2))
%
%     slp = pconn_regress(d(:,1)',d(:,2));
%
%     figure; set(gcf,'color','white'); hold on;
%
%     scatter(d(:,1),d(:,2),250,'markerfacecolor','k','markeredgecolor','w');
%
%     xlabel('Change DFA (ATX - PBO)'); ylabel('Change DFA (Tsk - Rst)');
%
%     line([-0.075 0.075],[slp(2)*(-0.075)+slp(1) slp(2)*(0.075)+slp(1)],'color','k','linewidth',3)
%
%     axis([-0.15 0.15 -0.25 0.15]); box on; set(gca,'TickDir','out')
%
%     print(gcf,'-deps',sprintf('~/pconn_all/plots/pconn_src_%s_%s_tskrstpharm_scatter_c%d_f%d_v%d.eps',str,str_cond,icontr,ifoi,v))
%
%     clear d r
%
%     for ivox = 1 : size(d1,1)
%
%       [r(ivox) p(ivox)] = corr(d1(ivox,:)',d2(ivox,:)');
%
%     end
%
% %     alp = fdr(p,0.05);
%
%     para = [] ;
%
%     para.colorlimits = [-.7 .7];
%
% %     r(p>alp) = eps;
%
% %
%     par_interp = spatfiltergauss(r',g1,dd,g2);
% %
%     figure; set(gcf,'color','white'); hold on;
%
% %     para.colorlimits = [min(r(r~=0)) max(r(r~=0))];
%
%     for iplot = 1 : 6
%
%       subplot(3,2,iplot)
%       para.myviewdir = viewdir(iplot,:);
%
%       a = sa_meg_template.cortex10K;
%
%       if iplot == 5
%         a=cutsurface(sa_meg_template.cortex10K,[mean(sa_meg_template.cortex10K.vc(:,1)) 0 0],[1 0 0]);
%       elseif iplot == 6
%         a=cutsurface(sa_meg_template.cortex10K,[mean(sa_meg_template.cortex10K.vc(:,1)) 0 0],[-1 0 0]);
%       end
%
%       pconn_showsurface(a,para,par_interp)
%       colormap(cmap)
%       camlight headlight
%
%     end
% %
%     print(gcf,'-djpeg100',sprintf('~/pconn_all/plots/pconn_src_%s_%s_tskrstpharm_mask%d_c%d_f%d_v%d.jpg',str,str_cond,mask,icontr,ifoi,v))
% % %
%   end
% end
%
% s = p<alp;
%
% %% PLOT AMPLIUDE
%
% str = 'amp';
% ifoi = 3;
% mask = 1;
%
% contrasts = [2 1; 3 1; 2 3];
% viewdir = [0.001 -0 1; 0 -.5 0; 1 0 0; -.5 0 0; -.5 0 0;1 0 0 ];
%
% if strcmp(str,'dfa')
%   par_all_res = dfa_all_res(:,:,:,ifoi);
%   par_all_cnt = dfa_all_cnt(:,:,:,ifoi);
% elseif strcmp(str,'var')
%   par_all_res = var_all_res(:,:,:,ifoi);
%   par_all_cnt = var_all_cnt(:,:,:,ifoi);
% elseif strcmp(str,'cvar')
%   par_all_res = cvar_all_res(:,:,:,ifoi);
%   par_all_cnt = cvar_all_cnt(:,:,:,ifoi);
% elseif strcmp(str,'amp')
%   par_all_res = amp_all_res(:,:,:,ifoi);
%   par_all_cnt = amp_all_cnt(:,:,:,ifoi);
% end
%
% load(sprintf('~/pconn_all/proc/all_src_clusterstat_rst_%s_c%d_f%d_v%d.mat','dfa',icontr,ifoi,v_stat))
%
% par = squeeze(nanmean(par_all_res(find(stats.mask),:,:),1));
%
% m1 = nanmean(par(1,:),2);
% m2 = nanmean(par(2,:),2);
%
% s1 = nanstd(par(1,:),[],2)/sqrt(size(par,2));
% s2 = nanstd(par(2,:),[],2)/sqrt(size(par,2));
%
% [~,p] = ttest(par(1,:),par(2,:));
%
% figure; set(gcf,'color','white'); hold on;
%
% bar(1, m1,.02,'facecolor','r','edgecolor','none');
% bar(1.05, m2,.02,'facecolor','b','edgecolor','none');
%
% line([1 1 ],[m1-s1 m1+s1]);
% line([1.05 1.05],[m2-s2 m2+s2]);
%
% axis([0.95 1.10 0 max([m1 m2])+max([m1 m2])*0.2])
%
% box on; set(gca,'Tickdir','out')
%
% % print(gcf,'-depsc2',sprintf('~/pconn_all/plots/pconn_src_dfaamp_f%d_v%d.eps,ifoi,v))
%
% %% PLOT AAL
% a=aalgrid.mask
% par_interp = spatfiltergauss(a,g1,dd,g2)
% para.colorlimits =[min(par_interp) max(par_interp)]
% for iplot = 1 : 6
%
%     subplot(3,2,iplot)
%
%      para.myviewdir = viewdir(iplot,:);
%      a = sa_meg_template.cortex10K;
%
%      if iplot == 5
%       	a=cutsurface(sa_meg_template.cortex10K,[mean(sa_meg_template.cortex10K.vc(:,1)) 0 0],[1 0 0]);
%      elseif iplot == 6
%       	a=cutsurface(sa_meg_template.cortex10K,[mean(sa_meg_template.cortex10K.vc(:,1)) 0 0],[-1 0 0]);
%      end
%
%     pconn_showsurface(a,para,par_interp); colormap(parula)
% end
%
%
% %% PLOT INDIVIDUAL SUBJECTS
%
% str = 'amp';
%
% for ifoi = 1 : 5
% % ifoi = 1;
%   mask = 0;
%
%   contrasts = [2 1; 3 1; 2 3];
%   viewdir = [0.001 -0 1; 0 -.5 0; 1 0 0; -.5 0 0; -.5 0 0;1 0 0 ];
%
%   try
%   if strcmp(str,'dfa')
%     par_all_res = dfa_all_res(:,:,:,ifoi);
%     par_all_cnt = dfa_all_cnt(:,:,:,ifoi);
%   elseif strcmp(str,'var')
%     par_all_res = var_all_res(:,:,:,ifoi);
%     par_all_cnt = var_all_cnt(:,:,:,ifoi);
%   elseif strcmp(str,'cvar')
%     par_all_res = cvar_all_res(:,:,:,ifoi);
%     par_all_cnt = cvar_all_cnt(:,:,:,ifoi);
%   elseif strcmp(str,'amp')
%     par_all_res = amp_all_res(:,:,:,ifoi);
%     par_all_cnt = amp_all_cnt(:,:,:,ifoi);
%   end
%   catch me
%     warning(me.message)
%   end
%
%   cmap = cbrewer('div', 'RdBu', 100,'pchip');% colormap(autumn)
%   cmap = cmap(end:-1:1,:);
%
%   for icontr = 1 : 2
%
%     figure; set(gcf,'color','white'); hold on;
%   %   -------------------------------------------------------
%   %   TASK
%     % -------------------------------------------------------
%   	for isubj = 1:length(SUBJLIST)
%
%       d = squeeze(nanmean(par_all_cnt(:,contrasts(icontr,1),isubj),3))-squeeze(nanmean(par_all_cnt(:,contrasts(icontr,2),isubj),3));
%       par_interp = spatfiltergauss(d,g1,dd,g2);
%
%
%       para = [] ;
%       r1 = min(min(nanmean(par_all_cnt(:,contrasts(icontr,1),:),3)-nanmean(par_all_cnt(:,contrasts(icontr,2),:),3)));
%       r2 = max(max(nanmean(par_all_cnt(:,contrasts(icontr,1),:),3)-nanmean(par_all_cnt(:,contrasts(icontr,2),:),3)));
%       para.colorlimits = [r1 r2];
%
%       subplot(4,5,isubj)
%       para.myviewdir = viewdir(1,:);
%       a = sa_meg_template.cortex10K;
%
%       pconn_showsurface(a,para,par_interp)
%       colormap(cmap)
%       camlight headlight
%
%     end
%
%     print(gcf,'-djpeg100',sprintf('~/pconn_all/plots/pconn_src_tsk_indivsubj_%s_mask%d_c%d_f%d_v%d.jpg',str,mask,icontr,ifoi,v_stat))
%   %
%     clear para r d stats
%
%     % -------------------------------------------------------
%     % RESTING STATE
%     % -------------------------------------------------------
%
%     figure; set(gcf,'color','white'); hold on;
%   %   -------------------------------------------------------
%   %   TASK
%     % -------------------------------------------------------
%     for isubj = 1:length(SUBJLIST)
%
%       d = squeeze(nanmean(par_all_res(:,contrasts(icontr,1),isubj),3))-squeeze(nanmean(par_all_res(:,contrasts(icontr,2),isubj),3));
%       par_interp = spatfiltergauss(d,g1,dd,g2);
%
%
%       para = [] ;
%       r1 = min(min(nanmean(par_all_res(:,contrasts(icontr,1),:),3)-nanmean(par_all_res(:,contrasts(icontr,2),:),3)));
%       r2 = max(max(nanmean(par_all_res(:,contrasts(icontr,1),:),3)-nanmean(par_all_res(:,contrasts(icontr,2),:),3)));
%       para.colorlimits = [r1 r2];
%
%       subplot(4,5,isubj)
%       para.myviewdir = viewdir(1,:);
%       a = sa_meg_template.cortex10K;
%
%       pconn_showsurface(a,para,par_interp)
%       colormap(cmap)
%       camlight headlight
%
%     end
%
%     print(gcf,'-djpeg100',sprintf('~/pconn_all/plots/pconn_src_res_indivsubj_%s_mask%d_c%d_f%d_v%d.jpg',str,mask,icontr,ifoi,v_stat))
%
%     clear para r d stats
%
%   end
%
% end