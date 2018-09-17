% Model Information
% 
% Calculate the abundance of Heat Producing Elements in crust from LITHO1.0
% Scott A. Wipperfurth (
% updated September 2018
% 
% see http://igppweb.ucsd.edu/~gabi/crust1.html for more info
% 
% The 8 crustal layers:
% ====================
% 1) water             
% 2) ice               
% 3) upper sediments   (VP, VS, Litho1.rho not defined in all cells) 
% 4) middle sediments  "
% 5) lower sediments   "
% 6) upper crystalline crusts
% 7) middle crystalline crust
% 8) lower crystalline crust
% 
% a ninth layer gives V_Pn, V_Sn and Litho1.rho below the Moho. The values
% are associated with LLNL model G3Cv3 on continents and a thermal
% model in the oceans.
% 
% The model is defined from 89.5 to -89.5 deg latitude and -179.5 to 179.5 deg
% longitude (64,800 total cells (= 360*180).
% 
% 
% --- REFERENCES CITED IN CODE: ---
% 
% Arevalo, R., McDonough, W.F., Stracke, A., Willbold, M., Ireland, T.J., 
%     Walker, R.J., 2013. Simplified mantle architecture and distribution of 
%     radiogenic power: Mantle Architecture and Radiogenic Power. Geochemistry, 
%     Geophysics, Geosystems 14, 2265–2285.
%     https://doi.org/10.1002/ggge.20152
% 
% Chambat, F., Ricard, Y., Valette, B., 2010. Flattening of the Earth: 
%     further from hydrostaticity than previously estimated: Hydrostatic flattening.
%      Geophysical Journal International 183, 727–732. 
%     https://doi.org/10.1111/j.1365-246X.2010.04771.x
% 
% Dziewonski, A.M., Anderson, D.L., 1981. Preliminary reference Earth model.
%     Physics of the Earth and Planetary Interiors 25, 297–356.
%     https://doi.org/10.1016/0031-9201(81)90046-7
% 
% Huang, Y., Chubakov, V., Mantovani, F., Rudnick, R.L., McDonough, W.F., 2013.
%     A reference Earth model for the heat-producing elements and associated 
%     geoneutrino flux. Geochem. Geophys. Geosystems 14, 2003–2029.
%     https://doi.org/10.1002/ggge.20129
% 
% Laske, G., Masters, G., Ma, Z., Pasyanos, M.E., 2013. Update on CRUST1.0
%     - A 1-degree Global Model of Earth’s Crust, in: Geophysical Research 
%     Abstracts, EGU2013-2658. Presented at the EGU.
% 
% Olugboji, T.M., Lekic, V., McDonough, W., 2017. A statistical assessment of
%      seismic models of the U.S. continental crust using Bayesian inversion of 
%     ambient noise surface wave dispersion data. Tectonics 2017TC004468. 
%     https://doi.org/10.1002/2017TC004468
% 
% Pasyanos, M.E., Masters, T.G., Laske, G., Ma, Z., 2014. LITHO1.0: An 
%     updated crust and lithospheric model of the Earth: LITHO1.0. 
%     Journal of Geophysical Research: Solid Earth 119, 2153–2173. 
%     https://doi.org/10.1002/2013JB010626
% 
% Plank, T., 2014. The Chemical Composition of Subducting Sediments, 
%     in: Treatise on Geochemistry. Elsevier, pp. 607–629.
% 
% Rudnick, R.L., Gao, S., 2014. Composition of the Continental Crust, 
%     in: Treatise on Geochemistry. Elsevier, pp. 1–51. 
% 
% Ruedas, T., 2017. Radioactive heat production of six geologically important
%     nuclides. Geochem. Geophys. Geosystems 18, 3530–3541. 
%     https://doi.org/10.1002/2017GC006997
% 
% White, W.M., Klein, E.M., 2014. Composition of the Oceanic Crust,
%     in: Treatise on Geochemistry. Elsevier, pp. 457–496. 
%     https://doi.org/10.1016/B978-0-08-095975-7.00315-6
% 
% Yoder, C.F., 1995. Astrometric and Geodetic Properties of Earth and the 
%     Solar System, in: Global Earth Physics. American Geophysical Union, pp.
%      1–31. https://doi.org/10.1029/RF001p0001
% 
% 
% OTHER: 
% NASA Earth Fact Sheet: https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html
% 
% LITHO1.0 Download Page: https://igppweb.ucsd.edu/~gabi/crust1.html




for dets = 2 %1:6
   for methods = 1 %1 = H13, 2 = Bivariate
       
    clearvars -except dets methods 

    
    
%% Set Cluster Information
%{

% Get a handle to the cluster
c = parcluster;

% Get the job ID
id = j.ID;


% Get state of job
state = j.state



% Submit a batch pool job using 8 workers for 16 simulations
j = c.batch(@parallel_example, 1, {}, ‘Pool’, 8);

ClusterInfo.setWallTime('50')
ClusterInfo.setMemUsage('4000mb')
ClusterInfo.UserDefinedOptions('--ntasks=19','--share','--account=schmerr-prj-hi','--mail-type=ALL',...
    '--job-name=geonu_test');

addpath('/lustre/swipp/code/Functions') %add path to function location

numCores = 7; 

j = c.batch

%}

%myCluster=parcluster('local'); myCluster.NumWorkers=numCores; parpool(myCluster,numCores)
    
    maxNumCompThreads(8);
    
%% 1) ---- Define Model Space ----
    % Set maximum number of cores to use
    
fprintf('    Starting time = %s \n',datestr(now,'mmmm dd, yyyy HH:MM AM')) %print time of loop start

tic; %clearvars -except testing iterations
%cd(fileparts(matlab.desktop.editor.getActiveFilename));% moves to current folder

MASTER.StartTime = datestr(now,'mmmm dd, yyyy HH:MM AM');
addpath('/lustre/swipp/code/Functions') %add path to function location


   % gcp %print cluster information
% -- Define Possible Detectors --
    % Includes Name, longitude, latitude, depth (m), detector efficiency,
    % protons, nf_botleft_lonlat
    detectors = {137.31, 36.43, 1000, 0.7, 5.98*10^31,[134,34];... %KamLAND - Araki etal. (2005) 
                 13.57, 42.45, 1400, 0.842, 9.76E+30,[10,40];...   %Borexino - Wikipedia LNGS
                 -81.201, 46.475, 2092, 0.8, 5.76E+31,[-84,44];...  %SNO+ - Chen (2006), Andringa etal. (2016)
                 112.518, 22.118, 700, 0.8, 1.29E+33,[109,20];...   %JUNO - An etal. (2016) JUNO Yellowbook
                 101.71, 28.15, 2400, 0.8, 2.16E+32,[99,26];...    %Jinping - Beacom etal. (2016) letter of intent
                 -156.32, 19.72, 0, 0.8, 0,[-160,18];         %Hawaii
                 }; 
    detectors = cell2table(detectors);
    detectors.Properties.RowNames = {'KamLAND','Borexino','SNO','JUNO','Jinping','Hawaii'};
    detectors.Properties.VariableNames = {'Longitude','Latitude','Depth','Efficiency','Protons','nf_botleft_lonlat'};

    % - Choose Detector -
    det = detectors(dets,:); 

    
    iter = 10; 
    simple2.meth = methods; %1 = H13 method, 2 = bivariate 
    %det = detectors(all_det,:); 
    %det = 0;
    %abundance = 1; % 1 = don't calculate abundance again
    
    MASTER.iter = iter; 
    MASTER.detector = det; 
    if simple2.meth ==1
        MASTER.method = 'Huang et al. 2013'; %record method to variable with run information
    elseif simple2.meth ==2
        MASTER.method = 'Bivariate Analysis ';
    end

    
    
    % -- Calculate geoneutrino flux? --
    simple2.calcFlux = true; 
    
    % -- Calculate mantle signal? --
    calcMantle = true; 
    
    % -- Decide Mantle structure --
    if calcMantle == true
        % Two Options for mantle: Homogenous (no enriched or depleted
        % mantle) or layered (enriched mantle = 19% of mantle)
        calcMantleLayered = true; 
        MASTER.mantle.calc = 'true';
        
        if calcMantleLayered == true
            MASTER.mantle.Style = 'Layered';
        else
            MASTER.mantle.Style = 'Homogenous';
        end
    end
    simple2.calcMantleLayered = calcMantleLayered; 
   
     % -- Near Field Method --
     calcNearFieldTraditional = true;

    
%% -- Load Model Data --   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -- Load Bivariate Data (Vp vs SiO2 and SiO2 vs U,Th,K for granulite/amphibolite) --
    simple2.bivar = load('BivarData_04042018.mat'); 


% - Define energies we will calculate flux at -
    % this has to be here so we can pre-define variables
    E_int = 75; %(keV) Energy bin size
    simple2.energy(:,1) = 1806+E_int/2:E_int:3276-E_int/2; %(keV) 3272 = max but this is easier
    MASTER.energy.bin = simple2.energy; 
    MASTER.energy.binsize = '100 KeV';


% -- Load LITHO1.0 Data -- (Pasyanos et al. 2014)
    load('LITHO1_BaseData_OndrejReFormat.mat')
    MASTER.model = 'LITHO1.0';
    
    
% -- Load Preliminary Reference Earth Model -- (Dziewonski and Anderson, 1981)    
    load PREM.mat

% ----------  Logically index cc, oc, stable, and Archean crust  --------
    %{ 
    Assign a value of 1 if the cell has a crust type that we designate as cc,
    oc, stable, or Archean crust.  Designation as part of one of these 4 broad 
    categories is performed by us based on the names given by LITHO1.0
    (actually CRUST1.0). All other values are 0. 

    Oceanic Crust: 
    26 = A1 normal oceanic
    27 = A0 oceans 3 Myrs and younger
    28 = B- melt affected o.c. and oceanic plateaus
    31 = V1 Inactive Ridge, Alpha ridge
    36 = Y3 Caspian Sea Oceanic

    Stable Crust: (Archean + Proterozoic crust)
    3 = Archean (Antarctica) 
    4 = early Archean
    5 = late archean
    6 = early/mid proterozoic
    7 = early mid proterozoic (Antarctica, slow)
    8 = late proterozoic
    9 = slow late proterozoic

    Archean Crust: 
    3 = Archean (Antarctica) 
    4 = early Archean
    5 = late archean
    %}
    oc_types = [26,27,28,31,36];  %Crust1.0 data types
    oc = any(bsxfun(@eq,oc_types,Litho1.type(:,3)),2); clear oc_types
    cc = ~oc; 

    stable_types = [3,4,5,6,7,8,9];
    stable = any(bsxfun(@eq,stable_types,Litho1.type(:,3)),2); clear stable_types

    archean_types = [3,4,5];
    archean = any(bsxfun(@eq,archean_types,Litho1.type(:,3)),2); clear archean_types
   
    
   
 clear data_load x detectors selection p poolsize 

%% 2) ---- Pre-allocate Variables ----
% -- Extract layer info from original file & create new structures --
%{
  --Variables--
    latitude
    longitude
    elevation (m) (orignally in km)
    layer thickness (m) (orignally in km)
    density (kg/m3) (orignally in g/m3)
    Vp (km/s)
    Vs (km/s)
%}


% -- Allocate empty subfields to each layer (to be filled in later) --
    n = {'aU';'aTh';'aK';'hf'};
    x = zeros(length(Litho1.latlon),3); % 3 columns for mean, + uncertainty, - unc.
    for i = 1:length(n)
    ice(:).(n{i})   = deal(x);
    water(:).(n{i}) = deal(x);
    UC(:).(n{i})    = deal(x);
    MC(:).(n{i})    = deal(x);
    LC(:).(n{i})    = deal(x);
    LM(:).(n{i})    = deal(x);
    s1(:).(n{i})    = deal(x);
    s2(:).(n{i})    = deal(x);
    s3(:).(n{i})    = deal(x);
    end           
   
    % - Allocate empty subfields within the 'mass' subfield
    n = {'mass';'U';'Th';'K'};
    for i = 1:length(n)
    ice.mass(:).(n{i})   = deal(x);
    water.mass(:).(n{i}) = deal(x);
    UC.mass(:).(n{i})    = deal(x);
    MC.mass(:).(n{i})    = deal(x);
    LC.mass(:).(n{i})    = deal(x);
    LM.mass(:).(n{i})    = deal(x);
    s1.mass(:).(n{i})    = deal(x);
    s2.mass(:).(n{i})    = deal(x);
    s3.mass(:).(n{i})    = deal(x);
    end           
    
    
     % - Allocate empty subfields within the 'flux' subfield
    n = {'U238','Th232'};
    for i = 1:length(n)
    ice.flux(:).(n{i})   = deal(x);
    water.flux(:).(n{i}) = deal(x);
    UC.flux(:).(n{i})    = deal(x);
    MC.flux(:).(n{i})    = deal(x);
    LC.flux(:).(n{i})    = deal(x);
    LM.flux(:).(n{i})    = deal(x);
    s1.flux(:).(n{i})    = deal(x);
    s2.flux(:).(n{i})    = deal(x);
    s3.flux(:).(n{i})    = deal(x);
    end       

% -- Predefine variables used in "parfor" loop --
    % These variables have values added to them in every iteration and must
    % exist before the loop starts. See notes below for why length = 18 & 36
     
    %- sums variables
    [s1_cc_sums,s2_cc_sums,s3_cc_sums,UC_cc_sums,MC_cc_sums,LC_cc_sums,LM_cc_sums,...
        s1_nf_sums,s2_nf_sums,s3_nf_sums,UC_nf_sums,MC_nf_sums,LC_nf_sums,LM_nf_sums]...
        = deal(zeros(iter,9));   
    
    [s1_oc_sums,s2_oc_sums,s3_oc_sums,UC_oc_sums,MC_oc_sums,LC_oc_sums,LM_oc_sums]...
        = deal(zeros(iter,9));     
    
    [s1_cc_flux_sums,s2_cc_flux_sums,s3_cc_flux_sums,UC_cc_flux_sums,...
        MC_cc_flux_sums,LC_cc_flux_sums,LM_cc_flux_sums,s1_nf_flux_sums,...
        s2_nf_flux_sums,s3_nf_flux_sums,UC_nf_flux_sums,MC_nf_flux_sums,...
        LC_nf_flux_sums,LM_nf_flux_sums]...
        = deal(zeros(iter,length(simple2.energy(:,1))*2));   
    
    [s1_oc_flux_sums,s2_oc_flux_sums,s3_oc_flux_sums,UC_oc_flux_sums,...
        MC_oc_flux_sums,LC_oc_flux_sums,LM_oc_flux_sums]...
        = deal(zeros(iter,length(simple2.energy(:,1))*2));     
    
clear n x i

%% 3) ---- Define variable correlations ----
% ----------------- HUANG ET AL. 2013 METHOD -----------------------
% -- Define random uniform distribution for correlation --
    %{ 
    To provide correlation we will create a variable of random numbers which
    will define where on a distribution we select a value. E.g. both random
    numbers will be from the high-end of a distribution rather than
    random. THIS HAS TO BE OUT OF LOOP, otherwise correlation will not be
    perfect, which is unjustified.  We need different correlation for each
    layer, although within a layer the cells are correlated (i.e. use the 
    same correlation variable). The function "rand_n" is a local function
    defined at the end of the script.  This function outputs "0" if iter
    =1, which causes the code to use the central value for the
    distribution. 
    Need seperate correlation for:  
    1) thickness
    2) Vp and density (same correlation value)
    3) abundance


    %}
     
    n = {'s1';'s2';'s3';'UC';'MC';'LC';'LM'};
    x = rand_n(iter); %for thickness correlation, which is same for all layers
        
    for i = 1:length(n)    
    cor(:).(n{i}).vp = rand_n(iter); %"rand_n" is local function, see end of script
    cor(:).(n{i}).abund = rand_n(iter);
    cor(:).(n{i}).thick = x;
    end

        % Abundances endmember (for Huang et al. 2013 method)
    cor.MC.end.abund = rand_n(iter); % for felsic and mafic endmembers in MC (Huang et al. method)
    cor.LC.end.abund = rand_n(iter); % for felsic and mafic endmembers in LC (Huang et al. method)
    
    cor.MC.bivar.abund = rand_n(iter); % for abundance of U, Th, K for each SiO2 bin
    cor.LC.bivar.abund = rand_n(iter); 
    
        % Vp endmember (for Huang et al. 2013 method)
    cor.MC.end.vp = rand_n(iter); % for felsic and mafic endmembers in MC (Huang et al. method)
    cor.LC.end.vp = rand_n(iter); % for felsic and mafic endmembers in LC (Huang et al. method)
        
        % Abundance endmember (for BSE and mantle mass balance)
    cor.bse = rand_n(iter); 
    
% --------------------- BIVARIATE METHOD --------------------
% -- Calculate correlation for bi-variate analysis -- 
    %{ 
       The bivariate analysis (Vp-SiO2) does not have a set PDF so this needs to be
       done carefully to take into account the distribution (which isnt
       gaussian or log-normal). The following code populates a variable
       with the index to an SiO2 bin, where the frequency of occurance for
       each bin is equal to the probability of each bin from the sample
       set.
    %}

    % Granulite
    for i = 1:length(simple2.bivar.seis.gran.prob(1,:)) % loop through Vp bins
        P = simple2.bivar.seis.gran.prob(:,i); %probability of each SiO2 bin (row) for each Vp bin (column)

        x = cumsum([0 P(:).'/sum(P(:))]); %cumulative sum of # of samples for givin Vp bin
        x(end) = 1e3*eps + x(end); % sets end value slightly larger than 1 (why?)
        % The following code creates a matrix of length = iter of randomly
        % selected (based on probability from P) SiO2 bins. The syntax
        % "[~,~,X] = histcounts"returns the bin edges (really SiO2 bin) 
        % which we will use later for indexing SiO2 fit parameters. 
        [~,~,cor.LC.bivar.sio2(:,i)] = histcounts(rand(iter,1),x); % provides length = iter of bins to 
    end
    %cor.bivar.gran = sort(cor.bivar.gran); 
    
    % Amphibolite
    for i = 1:length(simple2.bivar.seis.amp.prob(1,:))
        P = simple2.bivar.seis.amp.prob(:,i); 

        x = cumsum([0 P(:).'/sum(P(:))]);
        x(end) = 1e3*eps + x(end);
        [~,~,cor.MC.bivar.sio2(:,i)] = histcounts(rand(iter,1),x);% sorted so there is correlation, ok since no other correlation is sorted
    end
    %cor.bivar.amp = sort(cor.bivar.amp); 
    clear P x    
    
   
    %load correlation_7-3-18.mat % THIS IS ONLY TO CMPARE GRIDDING EVERY
    %ITERATION AND GRIDDING ONLY ONCE. 
    
%% 4) ---- Define abundance and Vp information ----

% -- Define K ratio of K2O mass -- (it is first used in this section)
    K.r = (2*39.0983)/((2*39.0983)+(15.999)); %ratio of K/K20 mass = ~.83
    K.b = 0.00011959;  % mass fraction of K that is K40 (0.000117 is atomic frac)
    
    
% -- Assign HPE abundances to oc, UC, LM, and sediments --

    % -- Oceanic Crust abundances -- (White and Klein 2014, tbl 8; sigma = 30% (Huang etal. 2013 pg 2013)
    % We need to do this in two steps because assignment of multiple columns
    % to certain rows in a structure is a pain (i.e. 'oc' rows)...
    oc_aU = zeros(length(oc),3);  %need to pre-define so result is 64800 long (otherwise wont be)
    oc_aTh = zeros(length(oc),3); 
    oc_aK = zeros(length(oc),3); 
    x  = repmat([0.07 .021 .021]*10^-6,sum(oc),1); oc_aU(oc,:) = x; %create temp variable oc_a*
    x  = repmat([0.210 .063 .063]*10^-6,sum(oc),1); oc_aTh(oc,:) = x; 
    x  = repmat([0.0716 .0215 .0215]*10^-2,sum(oc),1); oc_aK(oc,:) = x;
    
    
    % - Assign OC stats to UC, MC, and LC -
    [UC.aU,  MC.aU,  LC.aU]  = deal(oc_aU);  %(kg/kg)
    [UC.aTh, MC.aTh, LC.aTh] = deal(oc_aTh); %(kg/kg)
    [UC.aK,  MC.aK,  LC.aK]  = deal(oc_aK);  %(kg/kg K (not K40))
    clear oc_aU oc_aTh oc_aK
    
    
    % -- UC abundances -- (Rudnick and Gao 2014, tbl 3)
    UC.aU(cc,:)  = repmat([2.7 0.6 0.6]*10^-6,sum(cc),1);    %(kg/kg)
    UC.aTh(cc,:) = repmat([10.5 1.0 1.0]*10^-6,sum(cc),1);   %(kg/kg)
    UC.aK(cc,:)  = repmat([2.32 0.19 0.19]*10^-2,sum(cc),1); %(kg/kg K, not K40)
    
    
    %--  LM (continent) peridotite abundance -- (Huang etal. 2013, tbl 5)
    LM.aU(cc,:)  = repmat([0.033 0.049 0.020]*10^-6,sum(cc),1); %(kg/kg)
    LM.aTh(cc,:) = repmat([0.150 0.277 0.097]*10^-6,sum(cc),1); %(kg/kg)
    LM.aK(cc,:)  = repmat([0.038 0.052 0.022]*10^-2,sum(cc),1); %(kg/kg K, not K40)
    
    LM.aU(oc,:)  = repmat([0.008 0 0]*10^-6,sum(oc),1); %(kg/kg)
    LM.aTh(oc,:) = repmat([0.022 0 0]*10^-6,sum(oc),1); %(kg/kg)
    LM.aK(oc,:)  = repmat([0.015 0 0]*10^-2,sum(oc),1); %(kg/kg K, not K40)
    
 
    % -- Sediment abundances -- (Plank 2014 GLOSS II model, tbl 2)
    sed.U  = repmat([1.73 0.09 0.09]*10^-6,length(UC.lat),1); % temp variable
    sed.Th = repmat([8.10 0.59 0.59]*10^-6,length(UC.lat),1); 
    sed.K  = repmat([2.21 0.14 0.14]*K.r*10^-2,length(UC.lat),1); 
    
    [s1.aU,  s2.aU,  s3.aU]    = deal(sed.U);  %(kg/kg)
    [s1.aTh, s2.aTh, s3.aTh]   = deal(sed.Th); %(kg/kg)
    [s1.aK,  s2.aK,  s3.aK]    = deal(sed.K);  %(kg/kg K, not K40)
    clear sed
    
    
    % -- Mass balance for U, Th, and K in mantle -- (from Arevalo et al. 2013)
    % Define elemental abundance
    bse.aU = randist(20.3,2.03,iter,cor.bse)*10^-9;  % kg/kg U
    bse.aTh = randist(79.5,7.95,iter,cor.bse)*10^-9; % kg/kg Th
    bse.aK40 = randist(240,48,iter,cor.bse)*K.b*10^-6; % kg/kg K40    
    
        
% -- Define Granulite and Amphibolite abundance endmembers -- 
% This can be outside loop since there is lateral correlation
    ppm = 10^-6; 
    wt = 10^-2 * 0.000112 * 0.83; %10^-2 to bring to kg/kg, .000112 to bring K to K40, and 0.83 for K20 to K 
      
    % - Middle crust (amphibolite) (log normal distribution)
    am.f.U  = logdist(1.37,1.03,0.59,iter,cor.MC.end.abund)*ppm;    am.m.U  = logdist(0.37,0.39,0.19,iter,cor.MC.end.abund)*ppm; %U
    am.f.Th = logdist(8.27,8.12,4.10,iter,cor.MC.end.abund)*ppm;    am.m.Th = logdist(0.58,0.57,0.29,iter,cor.MC.end.abund)*ppm; %Th
    am.f.K  = logdist(2.89,1.81,1.11,iter,cor.MC.end.abund)*wt;     am.m.K  = logdist(0.50,0.41,0.23,iter,cor.MC.end.abund)*wt; %K
    % - Lower crust (granulite)
    gr.f.U  = logdist(0.42,0.41,0.21,iter,cor.LC.end.abund)*ppm;    gr.m.U  = logdist(0.10,0.14,0.06,iter,cor.LC.end.abund)*ppm;
    gr.f.Th = logdist(3.87,7.35,2.54,iter,cor.LC.end.abund)*ppm;    gr.m.Th = logdist(0.30,0.46,0.18,iter,cor.LC.end.abund)*ppm;
    gr.f.K  = logdist(2.71,2.05,1.17,iter,cor.LC.end.abund)*wt;     gr.m.K  = logdist(0.39,0.31,0.17,iter,cor.LC.end.abund)*wt;    
    
    simple2.am = am; 
    simple2.gr = gr; 
    
% -- Define Granulite and Amphibolite Vp endmembers --  
% This can be outside loop since there is lateral correlation

    % - Middle crust (amphibolite)
    Vp.MC.f = randist(6.34,0.16,iter,cor.MC.end.vp); %(km/s) (normal dist.)
    Vp.MC.m = randist(6.98,0.20,iter,cor.MC.end.vp); %(km/s) (normal dist.)
    %- Lower crust (granulite)
    Vp.LC.f = randist(6.52,0.19,iter,cor.LC.end.vp);  %(km/s) (normal dist.)
    Vp.LC.m = randist(7.21,0.20,iter,cor.LC.end.vp);  %(km/s) (normal dist.) 
        
%% Define Isotope information (hp, mass, molar ratio)   
% -- Define heat Production from kg of isotope -- (Ruedas 2017)
    simple2.hp.U238 = 94.946*10^-6;  %(Watt/kg)
    simple2.hp.U235 = 568.402*10^-6; %(Watt/kg)
    simple2.hp.U    = 98.314*10^-6;  %(Watt/kg)
    simple2.hp.Th   = 26.368*10^-6;  %(Watt/kg)
    simple2.hp.K40  = 28.761*10^-6;  %(Watt/kg) K40, not K 
    
 
% -- Define isotope parameters -- 
    % Molar ratios from NNDC.BNL.gov Chart of Nuclides
    % Isotope masses from http://atom.kaeri.re.kr Table of Nuclides   
    
    % - Molar Ratio - 
    simple2.iso.molar.U238  = 0.992742; 
    simple2.iso.molar.U235  = 1 - simple2.iso.molar.U238;
    simple2.iso.molar.Th232 = 1.00; 
    simple2.iso.molar.K40   = 0.000117; % mass ratio = 0.000112

    % - Isotope Mass - (amu)
    simple2.iso.amu.U238  = 238.0507826; 
    simple2.iso.amu.U235  = 235.0439231;
    simple2.iso.amu.Th232 = 232.0380504; 
    simple2.iso.amu.K40   = 39.9639987; 
    simple2.iso.amu.amuKg = 1.660539040*10^-27;

    % - Avogadros Number - (atom/mole)
    simple2.iso.avgd = 6.022140857*10^23; % http://physics.nist.gov/cgi-bin/cuu/Value?na    

    % - Decay constant - %(/s)
    simple2.iso.dc.U238  = log(2)/(4.4683e9*3.154e7);  %(/s)
    simple2.iso.dc.U235  = log(2)/(0.70348e9*3.154e7); %(/s) 
    simple2.iso.dc.Th232 = log(2)/(14.00e9*3.154e7);   %(/s)
    simple2.iso.dc.K40   = log(2)/(1.253e9*3.154e7);   %(/s)

 
    
%% Define Antineutrino Parameters (oscillation, spectrum, cross-section)
% -- Define geoneutrino oscillation parameters   --
    % -- Values from Review of Particle Physics 2018, C. Patrignani et al. (Particle Data Group), Chin. Phys. C, 40, 100001 (2016) and 2017 update. 
    %    <http://pdg.lbl.gov/> (really from Capozzi et al. 2017)
    % No correlation among parameters (or at least it is not known) (source = Beda Roskovec)
    x = (asind(sqrt(0.0240)) - asind(sqrt(0.0215)))/3; % uncertainty (symetrical)
    the13               = randist(asind(sqrt(0.0215)),x,1);        %(degree) theta 13 center, +, -
    x = (asind(sqrt(0.354)) - asind(sqrt(0.297)))/3; %postive uncertainty (1 sigma)
    y = (asind(sqrt(0.297)) - asind(sqrt(0.250)))/3; % negative uncertainty
    the12               = logdist(asind(sqrt(0.297)),x,y,1);        %(degree) theta 12 
    simple2.pee.delm21  = logdist(7.37,0.59,0.44,1)*10^-5;  %(eV^2) delta-mass 21 
    simple2.pee.delm32n = logdist(2.56,0.13,0.11,1)*10^-3;  %(eV^2) delta-mass 32 normal hierachy 
    simple2.pee.delm32i = logdist(-2.54,0.12,0.12,1)*10^-3;  %(eV^2) delta-mass 32 inverted hierachy 
    weak                = asind(sqrt(0.23122));                     %(degree) weak mixing angle from Canas et al. 2016 Physics Letter B
    
    
% -- Calculate Oscillation probability constants -- (these don't involve
    % distance or energy so we calculate now)
    simple2.pee.p1 = cosd(the13(1)).^4 .* sind(2.*the12(1)).^2;
    simple2.pee.p2 = sind(2*the13(1)).^2 .* cosd(2.*the12(1)).^2;
    simple2.pee.p3 = sind(2*the13(1)).^2;        
        
        
% -- Load and bin geonuetrino spectrum -- (Enomoto et al. 2007)
%{
File structure: 
1 energy (keV) (1 keV bins)
2 energy (MeV)
3 energy (pJ)
4 238U frequency (sums to 6)
5 235U frequency (sums to 4)
6 232Th frequency (sums to 4)
7 40K frequency (sums to 0.893)
8 187Re 
9 176Lu
10 138La
11 115Inwww
12 113Cd
13 87Rb
%}
%[enomoto,~,~] = xlsread('Enomoto2007_AntineutrinoSpectrum.xlsx,2,'A2:M4387'); % 
load 'Enomoto2007_AntineutrinoSpectrum.mat'

    %- Re-bin enomoto data around energies (original bin size was 1 KeV)-
        for i = 1:length(simple2.energy)
            top = simple2.energy(i)+0.5*E_int;   % top indice
            bot = simple2.energy(i)-0.5*E_int+1; % bottom indice
        simple2.eno.U238(i,1)  = sum(enomoto(bot:top,4));  
        simple2.eno.U235(i,1)  = sum(enomoto(bot:top,5));
        simple2.eno.Th232(i,1) = sum(enomoto(bot:top,6));  
        simple2.eno.K40(i,1)   = sum(enomoto(bot:top,7));  
        end
    
    % Calc. number of neutrinos per decay chain (for Sramek flux equation)
        simple2.eno.num.U238  = sum(enomoto(:,4)); 
        simple2.eno.num.U235  = sum(enomoto(:,5)); 
        simple2.eno.num.Th232 = sum(enomoto(:,6)); 
        simple2.eno.num.K40   = sum(enomoto(:,7)); 

    % - Calculate antineutrino interaction cross-section - (see Dye 2012, eqn 15)
        prot = 0.511; %(MeV/c)
        E = simple2.energy/1000; % convert to (MeV)
        x = sind(weak)^2; % from wikipedia "Weak Mixing Angle"
        Tmax = E./(1+prot./(2*E));
        simple2.sigma = 0.43.*(x.^2.*Tmax+(x+1).^2.*E./(3.*(1-(1-Tmax./E).^3))... 
                - x.*(x+1).*prot.*Tmax.^2./(2.*E.^2)).*(10^-44); %(cm^2)
        simple2.sigma = simple2.sigma/1e4; %(m^2)
        

% -- Define Bin-edges for "Flux vs Distance" plot -- (m)
        simple2.centers = logspace(2,7.1,100);

        [s1_distCount_sums,s2_distCount_sums,s3_distCount_sums,UC_distCount_sums,...
            MC_distCount_sums,LC_distCount_sums,LM_distCount_sums,man_distCount_sums]...
            = deal(zeros(length(logspace(2,7,100)),length(simple2.energy(:,1))*2));       
    
    
clear prot Tmax x E enomoto top bot E_int

%% 5) ---- Perform non-loop calculations ----

% -- Find radius at detector location -- (m) (radius Earth at point 
if size(det,2)>1
    x = knnsearch([Litho1.latlon(:,1),Litho1.latlon(:,2)],[det{1,1},det{1,2}]);
    det{1,3} = Litho1.r(x) - det{1,3}; clear x 
    det.Properties.VariableNames = {'Longitude','Latitude','Radius','Efficiency','Proton','NearField'};
end
    
% ---- Find near and far-field tiles ----
    % near-field crust is the 4x6 degree tiles around detector (so 24 tiles)

if calcNearFieldTraditional == true
    
    % -- Near Field = 4x6 degrees (rectangle; traditional way)
        z = det{1,6}; 
        x = z(1)+0.5:1:z(1)+6-0.5; % longitude centers 
        y = z(2)+0.5:1:z(2)+4-0.5; % latitude centers
        [X,Y] = meshgrid(x,y);
        temp = horzcat(X(:),Y(:));
        

        nearField.idx = knnsearch(Litho1.latlon,temp); %indeces of nearField cells
        nearField.logic = ismember(Litho1.latlon,Litho1.latlon(nearField.idx,:),'rows'); %logical index
%{
        figure
        coasts(gca,'k'); axis tight equal; hold on;
        scatter(Litho1.latlon(nearField,1),Litho1.latlon(nearField,2),50,'r','x')
        scatter(det{1,1},det{1,2},'b');
        for i = 1:length(nearField)
        rectangle('Position',[Litho1.latlon(nearField(i),1)-0.5 Litho1.latlon(nearField(i),2)-0.5 1 1])
        end
        axis([z(1)-1 z(1)+7 z(2)-1 z(2)+5])
        pos(5)
%}
        
    else

    % -- Near Field = closest 24 cells (not rectangle)
        % - closest cell
        temp = Litho1.latlon(knnsearch([UC.lat UC.lon],[det{1,2} det{1,1}]),:); %lon lat

        % Find closest 24 cells
        distance = dis(Litho1.latlon(:,1),Litho1.latlon(:,2),det{1,1:3},Litho1.r,0); 
        x = sort(distance); 
        nearField.idx = knnsearch(distance,x(1:(4*6))); %indeces of nearField cells
        nearField.logic = ismember(Litho1.latlon,Litho1.latlon(nearField.idx,:),'rows'); %logical index

        %nearField = knnsearch(distance,distance(distance<500000)); 
%{
        figure
        coasts(gca,'k'); axis tight equal; hold on;
        scatter(Litho1.latlon(nearField,1),Litho1.latlon(nearField,2),50,'r','x')
        scatter(det{1,1},det{1,2},'b');
        for i = 1:length(nearField)
        rectangle('Position',[Litho1.latlon(nearField(i),1)-0.5 Litho1.latlon(nearField(i),2)-0.5 1 1])
        end
        %axis([temp(1)-6 temp(1)+6 temp(2)-4 temp(2)+4])
        pos(5)
        %}
end


% ---- Calculate Geo-flux parameters to finish calculation to get TNU ---- 

iso = simple2.iso; 
eno = simple2.eno;

TNU.U238 = 1/(iso.amu.U238*iso.amu.amuKg*iso.molar.U238) .*eno.U238' /sum(eno.U238)...
    *6 * iso.dc.U238 /10^4 /(7.67*10^4); % nu/cm2/s

TNU.Th232 = 1/(iso.amu.Th232*iso.amu.amuKg) .* eno.Th232' /sum(eno.Th232)...
    *4 * iso.dc.Th232 /10^4 /(2.48*10^5); % nu/cm2/s


clear iso eno temp distance z x y X Y 


%% 6) ---- Lithosphere Monte Carlo  ----
%{
The monte carlo is used to calculate the mass distribution and heat production
for each layer. For the middle and lower crust, the abundance of U, Th, and
K is calculated from comparison of seismic velocity and petrologic
constraints, then the heat production is calculated. 
%}
%disp('Starting Monte carlo')
rng shuffle %reseed random number generator based on time (DO NOT PUT INTO LOOP, time intensive and maybe logically not correct)

tic
fprintf('Lithosphere Monte Carlo...')
nf_temp = nearField.logic; % Needed otherwise "nearField" becomes broadcast variable

parfor n = 1:50; %length(Litho1.latlon) % n = a specific cell (out of 64,800)
%percent_done(n,length(Litho1.latlon),5)

temp_P = zeros(iter,1); %temporary pressure, reset every new cell (otherwise it will continue summing between cells
    if cc(n) == 1 % cc(n) = 1 is continental crust, 0 is oceanic crust
        %fprintf('cc %d \n',n)
    % -- Input data into "Huang13" function --
%{ 
    "Huang13" function outputs data in matrix form (w/ length = "iter") as 
    "layer"_output_sums and statistics form (median, + uncertainty, - uncertainty)
    as "layer"_output_stat. The output_sums matrix is used to add together
    those parameters every iteration of "n". In this way we can sum the mass
    of the crust or mass of each element in each layer. 
    
    "output_sums" 
    Col 
     1      mass (kg)
     2      mass of U (kg)
     3      mass of Th (kg) 
     4      mass of K40 (kg)
     5      heat flow (W/m^2)
     6      total heat production (W)
     7      heat production from U (W)
     8      heat production from Th (W)
     9      heat production from K (W)
    
    "ouput_stat"
    Col
    1-3     mass (kg; median, + uncertainty, - uncertainty)
    4-6     mass of U (kg)
    7-9     mass of Th (kg) 
    10-12   mass of K40 (kg)
    13-15   heat production (W)
    16-18   heat flow (W/m^2)
    19-21   abundance of U (kg/kg)
    22-24   abundance of Th (kg/kg)
    25-27   abundance of K40 (kg/kg)
    28-30   "f" fraction of felsic (only exists for MC and LC)
    31-33   temperature in center of layer (only exists for MC and LC)
    34-36   # of times repeated (only exists for MC and LC)
     

      - "Water" and "Ice" are commented out because they contain no U, Th or K --
    
         % -  Water  (water; layer 1)    
               [water_output_stat, water_output_sums] = Huang13(n,iter,'water',water,simple);
               water_sums = water_sums + water_output_sums;    
    
         % -  Ice  (ice; layer 2)    
               [ice_output_stat, ice_output_sums] = Huang13(n,iter,'ice',ice,simple);
               ice_sums = ice_sums + ice_output_sums;    
      %}
  
         % -  Upper Sediment  (s1; layer 3)    
               [s1_out(n), s1_output_sums,s1_geoResponse(n,:),s1_flux(n),s1_flux_sums,P,s1_distCount]...
                   = Huang13_cluster(n,iter,'s1',s1,simple2,cor.s1,Litho1.r(n),det,temp_P);
              s1_cc_sums = s1_cc_sums + s1_output_sums;
              s1_cc_flux_sums = s1_cc_flux_sums + s1_flux_sums;
              temp_P = temp_P + P; % add up pressure
              
            
         % -  Middle Sediment  (s2; layer 4)    
               [s2_out(n), s2_output_sums,s2_geoResponse(n,:),s2_flux(n),s2_flux_sums,P,s2_distCount]...
                   = Huang13_cluster(n,iter,'s2',s2,simple2,cor.s2,Litho1.r(n),det,temp_P);
              s2_cc_sums = s2_cc_sums + s2_output_sums;  
              s2_cc_flux_sums = s2_cc_flux_sums + s2_flux_sums;
              temp_P = temp_P + P; % add up pressure
      
         % -  Lower Sediment  (s3; layer 5)    
               [s3_out(n), s3_output_sums,s3_geoResponse(n,:),s3_flux(n),s3_flux_sums,P,s3_distCount]...
                   = Huang13_cluster(n,iter,'s3',s3,simple2,cor.s3,Litho1.r(n),det,temp_P);
              s3_cc_sums = s3_cc_sums + s3_output_sums;               
              s3_cc_flux_sums = s3_cc_flux_sums + s3_flux_sums;
              temp_P = temp_P + P; % add up pressure
           
         % - Upper Crust (UC; layer 6)      
               [UC_out(n), UC_output_sums,UC_geoResponse(n,:),UC_flux(n),UC_flux_sums,P,UC_distCount]...
                   = Huang13_cluster(n,iter,'UC',UC,simple2,cor.UC,Litho1.r(n),det,temp_P);
              UC_cc_sums = UC_cc_sums + UC_output_sums; 
              UC_cc_flux_sums = UC_cc_flux_sums + UC_flux_sums;
              temp_P = temp_P + P; % add up pressure
              
         % - Middle Crust (MC; layer 7)      
               [MC_out(n), MC_output_sums,MC_geoResponse(n,:),MC_flux(n),MC_flux_sums,P,MC_distCount]...
                   = Huang13_cluster(n,iter,'MC',MC,simple2,cor.MC,Litho1.r(n),det,temp_P);
              MC_cc_sums = MC_cc_sums + MC_output_sums; 
              MC_cc_flux_sums = MC_cc_flux_sums + MC_flux_sums;
              temp_P = temp_P + P; % add up pressure
              
         % - Lower Crust (LC; layer 8)      
               [LC_out(n), LC_output_sums,LC_geoResponse(n,:),LC_flux(n),LC_flux_sums,P,LC_distCount]...
                   = Huang13_cluster(n,iter,'LC',LC,simple2,cor.LC,Litho1.r(n),det,temp_P);
              LC_cc_sums = LC_cc_sums + LC_output_sums;   
              LC_cc_flux_sums = LC_cc_flux_sums + LC_flux_sums;
              temp_P = temp_P + P; %moho(n).pressure = stat(P); % Pressure at base of LC (i.e. moho)
  
         % - Lithospheric Mantle (LM; layer 9)      
               [LM_out(n), LM_output_sums,LM_geoResponse(n,:),LM_flux(n),LM_flux_sums,P,LM_distCount]...
                   = Huang13_cluster(n,iter,'LM',LM,simple2,cor.LM,Litho1.r(n),det,temp_P);
               LM_cc_sums = LM_cc_sums + LM_output_sums;  
               LM_cc_flux_sums = LM_cc_flux_sums + LM_flux_sums; 
               %LAB(n).pressure = stat(temp_P + P); % pressure at base of lithosphere (i.e. LAB)
     
	
    else  % ----  CALCULATE ATTRIBUTES FOR OCEANIC CRUST --------
           %fprintf('   oc %d \n',n)
                 % -  Upper Sediment  (s1; layer 3)    
               [s1_out(n), s1_output_sums,s1_geoResponse(n,:),s1_flux(n),s1_flux_sums,P,s1_distCount]...
                   = Huang13_cluster(n,iter,'s1',s1,simple2,cor.s1,Litho1.r(n),det,temp_P);
               s1_oc_sums = s1_oc_sums + s1_output_sums;
               s1_oc_flux_sums = s1_oc_flux_sums + s1_flux_sums; 
               temp_P = temp_P + P; % add up pressure
   
         % -  Middle Sediment  (s2; layer 4)    
               [s2_out(n), s2_output_sums,s2_geoResponse(n,:),s2_flux(n),s2_flux_sums,P,s2_distCount]...
                   = Huang13_cluster(n,iter,'s2',s2,simple2,cor.s2,Litho1.r(n),det,temp_P);
               s2_oc_sums = s2_oc_sums + s2_output_sums;               
               s2_oc_flux_sums = s2_oc_flux_sums + s2_flux_sums;    
               temp_P = temp_P + P; % add up pressure
               
         % -  Lower Sediment  (s3; layer 5)    
               [s3_out(n), s3_output_sums,s3_geoResponse(n,:),s3_flux(n),s3_flux_sums,P,s3_distCount]...
                   = Huang13_cluster(n,iter,'s3',s3,simple2,cor.s3,Litho1.r(n),det,temp_P);
               s3_oc_sums = s3_oc_sums + s3_output_sums;  
               s3_oc_flux_sums = s3_oc_flux_sums + s3_flux_sums; 
               temp_P = temp_P + P; % add up pressure
              
         % - Upper Crust (UC; layer 6)      
               [UC_out(n), UC_output_sums,UC_geoResponse(n,:),UC_flux(n),UC_flux_sums,P,UC_distCount]...
                   = Huang13_cluster(n,iter,'UC',UC,simple2,cor.UC,Litho1.r(n),det,temp_P);
              UC_oc_sums = UC_oc_sums + UC_output_sums; 
              UC_oc_flux_sums = UC_oc_flux_sums + UC_flux_sums; 
               temp_P = temp_P + P; % add up pressure   
                 
         % - Middle Crust (MC; layer 7)  (put layer "MC_oc" so it doesnt do
            % abundance calculation) "MC_oc_output_stat" is different size than "MC_output_stat"
               [MC_out(n), MC_output_sums,MC_geoResponse(n,:),MC_flux(n),MC_flux_sums,P,MC_distCount]...
                   = Huang13_cluster(n,iter,'MC_oc',MC,simple2,cor.MC,Litho1.r(n),det,temp_P);
               MC_oc_sums = MC_oc_sums + MC_output_sums; 
               MC_oc_flux_sums = MC_oc_flux_sums + MC_flux_sums; 
               temp_P = temp_P + P; % add up pressure             
                  
         % - Lower Crust (LC; layer 8)      (put layer "LC_oc" so it doesnt do
            % abundance calculation)"LC_oc_output_stat" is different size than "LC_output_stat"
               [LC_out(n), LC_output_sums,LC_geoResponse(n,:),LC_flux(n),LC_flux_sums,P,LC_distCount]...
                   = Huang13_cluster(n,iter,'LC_oc',LC,simple2,cor.LC,Litho1.r(n),det,temp_P);
               LC_oc_sums = LC_oc_sums + LC_output_sums; 
               LC_oc_flux_sums = LC_oc_flux_sums + LC_flux_sums;
               temp_P = temp_P + P; 
               %moho(n).pressure = temp_P + P; % Pressure at base of LC (i.e. moho)
               
        % - Lithospheric Mantle (LM; layer 9) (need layer = 'LM_oc' so it
        %            fills it in blank
               [LM_out(n), LM_output_sums,LM_geoResponse(n,:),LM_flux(n),LM_flux_sums,P,LM_distCount]...
                   = Huang13_cluster(n,iter,'LM_oc',LM,simple2,cor.LM,Litho1.r(n),det,temp_P);
               LM_oc_sums = LM_oc_sums + LM_output_sums; 
               LM_oc_flux_sums = LM_oc_flux_sums + LM_flux_sums;  
               %LAB(n).pressure = moho(n).pressure + P; % pressure at base of lithosphere (i.e. LAB)


    end % end of if else statement
        % -- Record Flux vs Distance information --
            s1_distCount_sums = s1_distCount_sums + s1_distCount; 
            s2_distCount_sums = s2_distCount_sums + s2_distCount; 
            s3_distCount_sums = s3_distCount_sums + s3_distCount; 
            UC_distCount_sums = UC_distCount_sums + UC_distCount;
            MC_distCount_sums = MC_distCount_sums + MC_distCount;
            LC_distCount_sums = LC_distCount_sums + LC_distCount;
            LM_distCount_sums = LM_distCount_sums + LM_distCount;

        
        % -- Record Near Field Data --
         if nf_temp(n) == true %cell is part of near field
            % Abundance/Mass values
            s1_nf_sums = s1_nf_sums + s1_output_sums; 
            s2_nf_sums = s2_nf_sums + s2_output_sums; 
            s3_nf_sums = s3_nf_sums + s3_output_sums; 
            UC_nf_sums = UC_nf_sums + UC_output_sums;
            MC_nf_sums = MC_nf_sums + MC_output_sums;
            LC_nf_sums = LC_nf_sums + LC_output_sums;
            LM_nf_sums = LM_nf_sums + LM_output_sums;
           
            % Flux values
            s1_nf_flux_sums = s1_nf_flux_sums + s1_flux_sums; 
            s2_nf_flux_sums = s2_nf_flux_sums + s2_flux_sums; 
            s3_nf_flux_sums = s3_nf_flux_sums + s3_flux_sums; 
            UC_nf_flux_sums = UC_nf_flux_sums + UC_flux_sums;
            MC_nf_flux_sums = MC_nf_flux_sums + MC_flux_sums;
            LC_nf_flux_sums = LC_nf_flux_sums + LC_flux_sums;
            LM_nf_flux_sums = LM_nf_flux_sums + LM_flux_sums;
         end
     
end % End of Monte Carlo loop



MASTER.MCRunTime = sprintf('%.1f minutes',toc/60);
%MASTER.memory = memory; %not available on a cluster
MASTER.numCells = length(Litho1.latlon); 
fprintf('Done (Time Elapsed: %.1f min) \n',toc/60)



save testing.mat UC_cc_sums
disp('worked')


return
%% 7) ---- Mantle Monte Carlo ----

if calcMantle == true
fprintf('Mantle Monte Carlo...')


% Combine data for lithosphere
   %{
     Col 
     1      mass (kg)
     2      mass of U (kg)
     3      mass of Th (kg) 
     4      mass of K40 (kg)
     5      heat flow (W/m^2)
     6      total heat production (W)
     7      heat production from U (W)
     8      heat production from Th (W)
     9      heat production from K (W)
    %}
lith = s1_cc_sums + s2_cc_sums +  s3_cc_sums + UC_cc_sums + MC_cc_sums...
    + LC_cc_sums + LM_cc_sums + s1_oc_sums + s2_oc_sums +  s3_oc_sums...
    + UC_oc_sums + MC_oc_sums + LC_oc_sums + LM_oc_sums; 

% -- Calculate Mass of BSE --
earth.mass = randist(5.97218*10^24,0.000060*10^24,iter); %(kg) Chambat et al. 2010 Tbl 1
    x = 1.835*10^24 + 9.675*10^22; %(kg) mass of inner + outer core
core.mass   = randist(x,x*0.03,iter); %(kg) Yoder 1995 Tbl 2 (arbitrary uncertainty, but probably around percent level)
bse.mass = earth.mass - core.mass; 
man.mass = bse.mass - lith(:,1); %We only need approximate for abundance calculation

if calcMantleLayered == true
    % Set Enriched Mantle (EM) as 750km thick (19% of mantle by mass from
    % PREM). 19% is calculated in Arevalo et al. 2013
    man.propEM = 0.19; %750km ~= 19% of mantle by mass from PREM

    % Define Depleted Mantle Abundance (Arevalo et al. 2013 (tbl 5))
    man.aU.dm = randist(8.3,0.25,iter,cor.bse)*10^-9;
    man.aTh.dm = randist(24,0.5,iter,cor.bse)*10^-9;
    man.aK40.dm = randist(110,5,iter,cor.bse)*10^-6*K.b; 

    % EM BSE - Lithosphere - DM 
    man.aU.em = ((bse.aU .* bse.mass) - lith(:,2) - (man.aU.dm .* man.mass*(1-man.propEM)))...
        ./(man.mass*man.propEM);
    man.aTh.em = ((bse.aTh .* bse.mass) - lith(:,3) - (man.aTh.dm .* man.mass*(1-man.propEM)))...
        ./(man.mass*man.propEM);        
    man.aK40.em = ((bse.aK40 .* bse.mass) - lith(:,4) - (man.aK40.dm .* man.mass*(1-man.propEM)))...
        ./(man.mass*man.propEM);

        % Don't allow negative abundances (negative when crust has more HPE than BSE
        man.aU.dm(man.aU.dm<0) = 0; 
        man.aTh.dm(man.aTh.dm<0) = 0 ;
        man.aK40.dm(man.aK40.dm<0) = 0;

        man.aU.em(man.aU.em<0) = 0; 
        man.aTh.em(man.aTh.em<0) = 0 ;
        man.aK40.em(man.aK40.em<0) = 0;

        

else

     % -- Calculate abundances in convecting Mantle -- (Mantle = BSE - lithosphere)
    man.aU.dm = (bse.aU .* bse.mass - lith(:,2))./man.mass; % kg/kg U
    man.aTh.dm = (bse.aTh .* bse.mass - lith(:,3))./man.mass; % kg/kg Th
    man.aK40.dm = (bse.aK40 .* bse.mass - lith(:,4))./man.mass; % kg/kg K40

    man.aU.em = 0; man.aTh.em = 0; man.aK40.em = 0;
            % Don't allow negative abundances (negative when crust has more HPE than BSE
            man.aU.dm(man.aU.dm<0) = 0; 
            man.aTh.dm(man.aTh.dm<0) = 0 ;
            man.aK40.dm(man.aK40.dm<0) = 0;
   
end

% -- Monte Carlo for Mantle (heat production, geoneutrino flux)
man_sums_dm = zeros(1,1); 
man_sums_em = zeros(1,1);
man_flux_sums_dm = zeros(1,1); 
man_flux_sums_em = zeros(1,1);
man_distCount_sums_dm = zeros(length(simple2.centers),length(simple2.energy)*2); 
man_distCount_sums_em = zeros(length(simple2.centers),length(simple2.energy)*2); 

tic
    parfor n = 1:length(Litho1.latlon)

           [man_out_dm(n), man_sums_out_dm,man_flux_dm(n),man_flux_out_dm,man_distCount_dm,...
             man_out_em(n), man_sums_out_em,man_flux_em(n),man_flux_out_em,man_distCount_em]...
           = mantleGeo(Litho1.LAB(n)*1000,Litho1.latlon(n,2),Litho1.latlon(n,1),...
           man.aU,man.aTh,man.aK40,det,Litho1.r(n),PREM,simple2);

               %fprintf('made it here')
               % Depleted Mantle
               man_sums_dm = man_sums_dm + man_sums_out_dm;  
               man_flux_sums_dm = man_flux_sums_dm + man_flux_out_dm; 
               man_distCount_sums_dm = man_distCount_sums_dm + man_distCount_dm;
               
               % Enriched Mantle
               man_sums_em = man_sums_em + man_sums_out_em;  
               man_flux_sums_em = man_flux_sums_em + man_flux_out_em; 
               man_distCount_sums_em = man_distCount_sums_em + man_distCount_em;

    end

end %end of "if calcMantleLayered == true"
MASTER.MCRunTimeMantle = sprintf('%.1f minutes',toc/60);
MASTER.EndTime = datestr(now,'mmmm dd, yyyy HH:MM AM');
fprintf('Done (Time Elapsed: %.1f min) \n',toc/60)


%% 
%{

DO NOT DELETE!!! This section defines variables used in the Huang13
function and allows for easy testing.  The n = 17114 is the cell containing
Borexino (i.e. it is a good place to test things, including the gridding). 

%%%%%%%%%%%%%
n = 17114; 
%n = 20355; 
s1 = MC;
%s1.radius = simple1.radius; 
s2 = simple2; 
cor = cor.MC; 
SurfRadius = Litho1.r(n);
detector = det; 
P = zeros(iter,1); 
layer = 'MC';

%%%%%%%%%%%%
for "geoFlux" function

geo = 1; 
a = 1;
isoParam = simple2.iso;
OscilParam = simple2.oscillate; 
energy = simple2.energy;
radius = SurfRadius - depth3; 
lat_int = lat_int3; 
lon_int = lon_int3;
lon = lon3; 
lat = lat3; 
d_int = d_int3;

%%%%%%%%%%%%%%%
For 'mantleGeo.m" function

simple2.calcMantleLayered = true;

n = 17114; %borexino

s2 = simple2; 
LAB = Litho1.LAB(n)*1000; %m (delete later)
iter = 10000; 
aU = man.aU;
aTh = man.aTh;
aK40 = man.aK40;
lat = Litho1.latlon(n,2);
lon = Litho1.latlon(n,1); 
SurfRadius = 6.3780e6; %m
detector = det; 
n = 1;







scatter(1:8,time_taken(:,1)/60,80,'filled','r')
xlabel('Number of cores used')
ylabel('Run time (min)')
title('Run time vs Pool size for 1000 iterations')
set(gca,'fontsize',20)
 


    [lon2,lat2,depth2,lon_int2,lat_int2,d_int2]...
        = miniVox(s1.lon,s1.lat,depth',1,1,thick',num,SurfRadius); 



%}

%% 8) ---- Calculate Stats of each layer ---- (hp, hf, mass)
 %determines what statistics to output (type "help stat")
method = 4;    % 4 = geometric mean with +- std
meth_norm = 5; % 5 = normal distribution maximum liklihood fit

% -- Calculate statistics (simplified by 'RM18_allocate' function)
s1.total = RM18_allocate(s1_cc_sums,meth_norm,cc,'s1'); 
s2.total = RM18_allocate(s2_cc_sums,meth_norm,cc,'s2'); 
s3.total = RM18_allocate(s3_cc_sums,meth_norm,cc,'s3'); 
UC.total = RM18_allocate(UC_cc_sums,meth_norm,cc,'UC'); 
MC.total = RM18_allocate(MC_cc_sums,method,cc,'MC'); 
LC.total = RM18_allocate(MC_cc_sums,method,cc,'LC'); 
LM.total = RM18_allocate(LM_cc_sums,method,cc,'LM'); 
mantle.dm.total = RM18_allocate(man_sums_dm,meth_norm,cc,'man'); 
mantle.em.total = RM18_allocate(man_sums_em,meth_norm,cc,'man');

% -- Move the abundance information to match other layer structures --
MC.aU = MC.total.aU; MC.aTh = MC.total.aTh; MC.aK = MC.total.aK; 
LC.aU = LC.total.aU; LC.aTh = LC.total.aTh; LC.aK = LC.total.aK; 
mantle.dm.aU = mantle.dm.total.aU; mantle.dm.aTh = mantle.dm.total.aTh; mantle.dm.aK = mantle.dm.total.aK; 
mantle.em.aU = mantle.em.total.aU; mantle.em.aTh = mantle.em.total.aTh; mantle.em.aK = mantle.em.total.aK; 
MC.total = rmfield(MC.total,{'aU','aTh','aK'});
LC.total = rmfield(LC.total,{'aU','aTh','aK'});
mantle.dm.total = rmfield(mantle.dm.total,{'aU','aTh','aK'});
mantle.em.total = rmfield(mantle.em.total,{'aU','aTh','aK'});


%% 9) ---- Reallocate individual cell data into respective structures -- 
  %{
Note: This would preferrably be done directly in loop, but parfor is
    % not friendly
    
    "ouput_stat"
    Col
    1-3     mass (kg; median, + uncertainty, - uncertainty)
    4-6     mass of U (kg)
    7-9     mass of Th (kg) 
    10-12   mass of K40 (kg)
    13-15   heat production (W)
    16-18   heat flow (W/m^2)
    19-21   abundance of U (kg/kg)
    22-24   abundance of Th (kg/kg)
    25-27   abundance of K40 (kg/kg)
    28-30   "f" fraction of felsic (only exists for MC and LC)
    31-33   temperature in center of layer (only exists for MC and LC)
    34-36   # of times repeated (only exists for MC and LC)
%}

fprintf('Re-organizing Data...')
tic

% Mass of U, Th, K, and cell ('mass')(kg)
    n = {'U';'Th';'K';'mass'};% these need to be the same as "UC_out"
    % Note: putting 'mass.U' in "n" doesn't work.
    for j = 1:length(Litho1.latlon)
        for i = 1:length(n)    
            s1.mass(:).(n{i})(j,:) = s1_out(j).mass.(n{i});
            s2.mass(:).(n{i})(j,:) = s2_out(j).mass.(n{i});
            s3.mass(:).(n{i})(j,:) = s3_out(j).mass.(n{i});
            UC.mass(:).(n{i})(j,:) = UC_out(j).mass.(n{i});
            MC.mass(:).(n{i})(j,:) = MC_out(j).mass.(n{i});
            LC.mass(:).(n{i})(j,:) = LC_out(j).mass.(n{i});
            LM.mass(:).(n{i})(j,:) = LM_out(j).mass.(n{i});
            mantle.dm.mass(:).(n{i})(j,:) = man_out_dm(j).mass.(n{i});
            mantle.em.mass(:).(n{i})(j,:) = man_out_em(j).mass.(n{i});
       end
    end
    
% Heat production (W)
    n = {'U';'Th';'K';'total'};% these need to be the same as "UC_out"
    % Note: putting 'mass.U' in "n" doesn't work.
    for j = 1:length(Litho1.latlon)
        for i = 1:length(n)    
            s1.hp(:).(n{i})(j,:) = s1_out(j).hp.(n{i});
            s2.hp(:).(n{i})(j,:) = s2_out(j).hp.(n{i});
            s3.hp(:).(n{i})(j,:) = s3_out(j).hp.(n{i});
            UC.hp(:).(n{i})(j,:) = UC_out(j).hp.(n{i});
            MC.hp(:).(n{i})(j,:) = MC_out(j).hp.(n{i});
            LC.hp(:).(n{i})(j,:) = LC_out(j).hp.(n{i});
            LM.hp(:).(n{i})(j,:) = LM_out(j).hp.(n{i});
            mantle.dm.hp(:).(n{i})(j,:) = man_out_dm(j).hp.(n{i});
            mantle.em.hp(:).(n{i})(j,:) = man_out_em(j).hp.(n{i});
       end
    end    
    
% Heat flow (W/m2)
    n = {'hf'};% these need to be the same as "UC_out"
    for j = 1:length(Litho1.latlon)
        for i = 1:length(n)    
            s1(:).(n{i})(j,:) = s1_out(j).(n{i});
            s2(:).(n{i})(j,:) = s2_out(j).(n{i});
            s3(:).(n{i})(j,:) = s3_out(j).(n{i});
            UC(:).(n{i})(j,:) = UC_out(j).(n{i});
            MC(:).(n{i})(j,:) = MC_out(j).(n{i});
            LC(:).(n{i})(j,:) = LC_out(j).(n{i});
            LM(:).(n{i})(j,:) = LM_out(j).(n{i});
            mantle.dm(:).(n{i})(j,:) = man_out_dm(j).(n{i});
            mantle.em(:).(n{i})(j,:) = man_out_em(j).(n{i});
         end
    end   
    
% Abundance of U (ppm), Th (ppm), K (wt%), temperature (C), fraction felsic
        n = {'aU';'aTh';'aK';'temp';'f';'pressure'};% these need to be the same as "UC_out"
        x = 1:length(Litho1.latlon); x = x(cc); 
    for j = 1:length(x)
        for i = 1:length(n)    
            MC(:).(n{i})(x(j),:) = MC_out(x(j)).(n{i});
            LC(:).(n{i})(x(j),:) = LC_out(x(j)).(n{i});
       end
    end  
     
fprintf('Done (Time Elapsed: %.1f min) \n',toc/60)
    
%% ---- Calculate Final Geoneutrino Flux ----

U238 = 1:length(simple2.energy); %1:19
Th232 = length(simple2.energy)+1:length(simple2.energy)*2; %20:38

% Sediment (TNU)
    sed_crust = s1_cc_flux_sums + s2_cc_flux_sums + s3_cc_flux_sums; 
    flux.sed.sums.cc.U238 = sum(bsxfun(@times,sed_crust(:,U238),TNU.U238),2); 
    flux.sed.sums.cc.Th232 = sum(bsxfun(@times,sed_crust(:,Th232),TNU.Th232),2); 

    sed_ocean = s1_oc_flux_sums + s2_oc_flux_sums + s3_oc_flux_sums; 
    flux.sed.sums.oc.U238 = sum(bsxfun(@times,sed_ocean(:,U238),TNU.U238),2); 
    flux.sed.sums.oc.Th232 = sum(bsxfun(@times,sed_ocean(:,Th232),TNU.Th232),2); 

% UC (TNU)
    flux.UC.cc.U238 = sum(bsxfun(@times,UC_cc_flux_sums(:,U238),TNU.U238),2); 
    flux.UC.cc.Th232 = sum(bsxfun(@times,UC_cc_flux_sums(:,Th232),TNU.Th232),2); 

    flux.UC.oc.U238 = sum(bsxfun(@times,UC_oc_flux_sums(:,U238),TNU.U238),2); 
    flux.UC.oc.Th232 = sum(bsxfun(@times,UC_oc_flux_sums(:,Th232),TNU.Th232),2); 

% MC (TNU)
    flux.MC.cc.U238 = sum(bsxfun(@times,MC_cc_flux_sums(:,U238),TNU.U238),2); 
    flux.MC.cc.Th232 = sum(bsxfun(@times,MC_cc_flux_sums(:,Th232),TNU.Th232),2); 

    flux.MC.oc.U238 = sum(bsxfun(@times,MC_oc_flux_sums(:,U238),TNU.U238),2); 
    flux.MC.oc.Th232 = sum(bsxfun(@times,MC_oc_flux_sums(:,Th232),TNU.Th232),2); 

% LC (TNU)
    flux.LC.cc.U238 = sum(bsxfun(@times,LC_cc_flux_sums(:,U238),TNU.U238),2); 
    flux.LC.cc.Th232 = sum(bsxfun(@times,LC_cc_flux_sums(:,Th232),TNU.Th232),2); 

    flux.LC.oc.U238 = sum(bsxfun(@times,LC_oc_flux_sums(:,U238),TNU.U238),2); 
    flux.LC.oc.Th232 = sum(bsxfun(@times,LC_oc_flux_sums(:,Th232),TNU.Th232),2); 
    
% LM (TNU)
    flux.LM.cc.U238 = sum(bsxfun(@times,LM_cc_flux_sums(:,U238),TNU.U238),2); 
    flux.LM.cc.Th232 = sum(bsxfun(@times,LM_cc_flux_sums(:,Th232),TNU.Th232),2); 

    flux.LM.oc.U238 = sum(bsxfun(@times,LM_oc_flux_sums(:,U238),TNU.U238),2); 
    flux.LM.oc.Th232 = sum(bsxfun(@times,LM_oc_flux_sums(:,Th232),TNU.Th232),2); 
   
% Near Field Crust (oc and cc combined)
    flux.s1.nf.U238 = sum(bsxfun(@times,s1_nf_flux_sums(:,U238),TNU.U238),2); 
    flux.s2.nf.U238 = sum(bsxfun(@times,s2_nf_flux_sums(:,U238),TNU.U238),2); 
    flux.s3.nf.U238 = sum(bsxfun(@times,s3_nf_flux_sums(:,U238),TNU.U238),2); 
    flux.UC.nf.U238 = sum(bsxfun(@times,UC_nf_flux_sums(:,U238),TNU.U238),2); 
    flux.MC.nf.U238 = sum(bsxfun(@times,MC_nf_flux_sums(:,U238),TNU.U238),2); 
    flux.LC.nf.U238 = sum(bsxfun(@times,LC_nf_flux_sums(:,U238),TNU.U238),2); 
    flux.LM.nf.U238 = sum(bsxfun(@times,LM_nf_flux_sums(:,U238),TNU.U238),2); 

    flux.s1.nf.Th232 = sum(bsxfun(@times,s1_nf_flux_sums(:,Th232),TNU.Th232),2); 
    flux.s2.nf.Th232 = sum(bsxfun(@times,s2_nf_flux_sums(:,Th232),TNU.Th232),2); 
    flux.s3.nf.Th232 = sum(bsxfun(@times,s3_nf_flux_sums(:,Th232),TNU.Th232),2); 
    flux.UC.nf.Th232 = sum(bsxfun(@times,UC_nf_flux_sums(:,Th232),TNU.Th232),2); 
    flux.MC.nf.Th232 = sum(bsxfun(@times,MC_nf_flux_sums(:,Th232),TNU.Th232),2); 
    flux.LC.nf.Th232 = sum(bsxfun(@times,LC_nf_flux_sums(:,Th232),TNU.Th232),2); 
    flux.LM.nf.Th232 = sum(bsxfun(@times,LM_nf_flux_sums(:,Th232),TNU.Th232),2); 
    
% Flux vs Distance Counter
    flux.s1.count.U238 = sum(bsxfun(@times,s1_distCount_sums(:,U238),TNU.U238),2); 
    flux.s2.count.U238 = sum(bsxfun(@times,s2_distCount_sums(:,U238),TNU.U238),2); 
    flux.s3.count.U238 = sum(bsxfun(@times,s3_distCount_sums(:,U238),TNU.U238),2); 
    flux.UC.count.U238 = sum(bsxfun(@times,UC_distCount_sums(:,U238),TNU.U238),2); 
    flux.MC.count.U238 = sum(bsxfun(@times,MC_distCount_sums(:,U238),TNU.U238),2); 
    flux.LC.count.U238 = sum(bsxfun(@times,LC_distCount_sums(:,U238),TNU.U238),2); 
    flux.LM.count.U238 = sum(bsxfun(@times,LM_distCount_sums(:,U238),TNU.U238),2); 
    flux.man.dm.count.U238 = sum(bsxfun(@times,man_distCount_sums_dm(:,U238),TNU.U238),2);
    flux.man.em.count.U238 = sum(bsxfun(@times,man_distCount_sums_em(:,U238),TNU.U238),2);

    flux.s1.count.Th232 = sum(bsxfun(@times,s1_distCount_sums(:,Th232),TNU.Th232),2); 
    flux.s2.count.Th232 = sum(bsxfun(@times,s2_distCount_sums(:,Th232),TNU.Th232),2); 
    flux.s3.count.Th232 = sum(bsxfun(@times,s3_distCount_sums(:,Th232),TNU.Th232),2); 
    flux.UC.count.Th232 = sum(bsxfun(@times,UC_distCount_sums(:,Th232),TNU.Th232),2); 
    flux.MC.count.Th232 = sum(bsxfun(@times,MC_distCount_sums(:,Th232),TNU.Th232),2); 
    flux.LC.count.Th232 = sum(bsxfun(@times,LC_distCount_sums(:,Th232),TNU.Th232),2); 
    flux.LM.count.Th232 = sum(bsxfun(@times,LM_distCount_sums(:,Th232),TNU.Th232),2); 
    flux.man.dm.count.Th232 = sum(bsxfun(@times,man_distCount_sums_dm(:,Th232),TNU.Th232),2);
    flux.man.em.count.Th232 = sum(bsxfun(@times,man_distCount_sums_em(:,Th232),TNU.Th232),2);

    
%% 11) ---- Reproduce Huang et al. 2013 Table 2 (flux) ---
%{
The rows of the table are reservoirs.  The columns are average with
uncertainty for each reservoir. 

%}

method = 4; % 4 = median +- 68% c.l., 6 = geometric mean +- std
method_norm = 5; %normal distribution maximum liklihood fit

U238 = 1:length(simple2.energy); %1:16
Th232 = length(simple2.energy)+1:length(simple2.energy)*2; %1:32

% Row Names
    huang.tab2.rows = {'Sed_cc';'UC';'MC';'LC';'LM';'Sed_oc';'OC';'Bulk CC';...
        'Bulk Crust';'NearField Crust';'FarField Crust';'Total LS';...
        'DM';'EM';'Grand Total'};


% Sediment (TNU) (lognormal or normal, both fit well and give ~same unc.)
    huang.tab2.U238(1,:) = stat(flux.sed.sums.cc.U238,method); 
    huang.tab2.Th232(1,:) = stat(flux.sed.sums.cc.Th232,method); 
    huang.tab2.total(1,:) = stat(flux.sed.sums.cc.U238+flux.sed.sums.cc.Th232,method); 


    huang.tab2.U238(6,:) = stat(flux.sed.sums.oc.U238,method); 
    huang.tab2.Th232(6,:) = stat(flux.sed.sums.oc.Th232,method); 
    huang.tab2.total(6,:) = stat(flux.sed.sums.oc.U238+flux.sed.sums.oc.Th232,method); 


% UC (TNU) (lognormal distribution)
    huang.tab2.U238(2,:) = stat(flux.UC.cc.U238,method); 
    huang.tab2.Th232(2,:) = stat(flux.UC.cc.Th232,method); 
    huang.tab2.total(2,:) = stat(flux.UC.cc.U238+flux.UC.cc.Th232,method); 


% MC (TNU) (lognormal distribution)
    huang.tab2.U238(3,:) = stat(flux.MC.cc.U238,method); 
    huang.tab2.Th232(3,:) = stat(flux.MC.cc.Th232,method); 
    huang.tab2.total(3,:) = stat(flux.MC.cc.U238+flux.MC.cc.Th232,method); 


% LC (TNU) (lognormal distribution)
    huang.tab2.U238(4,:) = stat(flux.LC.cc.U238,method); 
    huang.tab2.Th232(4,:) = stat(flux.LC.cc.Th232,method); 
    huang.tab2.total(4,:) = stat(flux.LC.cc.U238+flux.LC.cc.Th232,method); 


% LM (TNU) (lognormal distribution)
    huang.tab2.U238(5,:) = stat(flux.LM.cc.U238,method); 
    huang.tab2.Th232(5,:) = stat(flux.LM.cc.Th232,method); 
    huang.tab2.total(5,:) = stat(flux.LM.cc.U238+flux.LM.cc.Th232,method); 


% OC Sed (TNU)
    huang.tab2.U238(6,:) = stat(flux.sed.sums.oc.U238,method); 
    huang.tab2.Th232(6,:) = stat(flux.sed.sums.oc.Th232,method); 
    huang.tab2.total(6,:) = stat(flux.sed.sums.oc.U238+flux.sed.sums.oc.Th232,method); 


% OC (TNU) (lognormal distribution)
    U238_oc = flux.UC.oc.U238 + flux.MC.oc.U238 + flux.LC.oc.U238; 
    Th232_oc= flux.UC.oc.Th232 + flux.MC.oc.Th232 + flux.LC.oc.Th232; 

    huang.tab2.U238(7,:) = stat(U238_oc,method); 
    huang.tab2.Th232(7,:) = stat(Th232_oc,method); 
    huang.tab2.total(7,:) = stat(U238_oc + Th232_oc,method); 

% Bulk CC (Sed + UC + MC + LC) (TNU) (lognormal distribution)
    flux.bulkCC.sums.cc.U238 = flux.UC.cc.U238 + flux.MC.cc.U238 + flux.LC.cc.U238...
            + flux.sed.sums.cc.U238; 
    flux.bulkCC.sums.cc.Th232 = flux.UC.cc.Th232 + flux.MC.cc.Th232 + flux.LC.cc.Th232...
            + flux.sed.sums.cc.Th232; 

    huang.tab2.U238(8,:)  = stat(flux.bulkCC.sums.cc.U238,method); 
    huang.tab2.Th232(8,:) = stat(flux.bulkCC.sums.cc.Th232,method); 
    huang.tab2.total(8,:) = stat(flux.bulkCC.sums.cc.U238 + flux.bulkCC.sums.cc.Th232,method); 


% Bulk Crust (TNU) (lognormal distribution)
    U238_crust = flux.bulkCC.sums.cc.U238 + flux.UC.oc.U238 + flux.MC.oc.U238...
        +flux.LC.oc.U238 + flux.sed.sums.oc.U238;
    Th232_crust = flux.bulkCC.sums.cc.Th232 + flux.UC.oc.Th232 + flux.MC.oc.Th232...
        +flux.LC.oc.Th232 + flux.sed.sums.oc.Th232; 
    
    huang.tab2.U238(9,:)  = stat(U238_crust,method); 
    huang.tab2.Th232(9,:) = stat(Th232_crust,method); 
    huang.tab2.total(9,:) = stat(U238_crust + Th232_crust,method); 

    
% Near Field Crust (NFC) (TNU) (lognormal distribution)
    nf_crust.U238 = flux.s1.nf.U238 + flux.s2.nf.U238 + flux.s3.nf.U238...
        + flux.UC.nf.U238 + flux.MC.nf.U238 + flux.LC.nf.U238; 
    nf_crust.Th232 = flux.s1.nf.Th232 + flux.s2.nf.Th232 + flux.s3.nf.Th232...
        + flux.UC.nf.Th232 + flux.MC.nf.Th232 + flux.LC.nf.Th232; 
    
    huang.tab2.U238(10,:)  = stat(nf_crust.U238,method); 
    huang.tab2.Th232(10,:) = stat(nf_crust.Th232,method); 
    huang.tab2.total(10,:) = stat(nf_crust.U238 + nf_crust.Th232,method);     
    
    
% Far Field Crust (FFC) (TNU) (lognormal distribution)
    ff_crust.U238 = U238_crust - nf_crust.U238; %total signal - nearField = Far Field
    ff_crust.Th232 = Th232_crust - nf_crust.Th232; 

    huang.tab2.U238(11,:)  = stat(ff_crust.U238,method); 
    huang.tab2.Th232(11,:) = stat(ff_crust.Th232,method); 
    huang.tab2.total(11,:) = stat(ff_crust.U238 + ff_crust.Th232,method);
    
    
% Total Lithosphere  (TNU) (lognormal distribution)
    U238_litho = U238_crust + flux.LM.cc.U238; 
    Th232_litho = Th232_crust + flux.LM.cc.Th232; 

    huang.tab2.U238(12,:) = stat(U238_litho,method); 
    huang.tab2.Th232(12,:) = stat(Th232_litho,method); 
    huang.tab2.total(12,:) = stat(U238_litho + Th232_litho,method); 

    
% Depleted Mantle (DM; TNU) (normal distribution (can't do lognormal cause 0 values)    
    flux.man.dm.U238 = sum(bsxfun(@times,man_flux_sums_dm(:,U238),TNU.U238),2); 
    flux.man.dm.Th232 = sum(bsxfun(@times,man_flux_sums_dm(:,Th232),TNU.Th232),2); 

    huang.tab2.U238(13,:) = stat(flux.man.dm.U238,method_norm); 
    huang.tab2.Th232(13,:) = stat(flux.man.dm.Th232,method_norm); 
    huang.tab2.total(13,:) = stat(flux.man.dm.U238+flux.man.dm.Th232,method_norm); 

    
% Enriched Mantle (EM; TNU)
    flux.man.em.U238 = sum(bsxfun(@times,man_flux_sums_em(:,U238),TNU.U238),2); 
    flux.man.em.Th232 = sum(bsxfun(@times,man_flux_sums_em(:,Th232),TNU.Th232),2); 

    huang.tab2.U238(14,:) = stat(flux.man.em.U238,method_norm); 
    huang.tab2.Th232(14,:) = stat(flux.man.em.Th232,method_norm); 
    huang.tab2.total(14,:) = stat(flux.man.em.U238+flux.man.em.Th232,method_norm); 

    
% Total Flux at Detector (TNU) (lognormal distribution)
    flux.total.U238 = U238_litho + flux.man.dm.U238 + flux.man.em.U238; 
    flux.total.Th232 = Th232_litho + flux.man.dm.Th232 + flux.man.em.U238;
    
    huang.tab2.U238(15,:) = stat(flux.total.U238,method); 
    huang.tab2.Th232(15,:) = stat(flux.total.Th232,method); 
    huang.tab2.total(15,:) = stat(flux.total.U238 + flux.total.Th232,method); 

clear ocean U238_oc U238_crust U238_litho Th232_oc Th232_crust Th232_litho

huang.tab2.table = horzcat(huang.tab2.U238,huang.tab2.Th232,huang.tab2.total); 

%% 12) ---- Reproduce Huang et al. 2013 Table 3 (abund) ---
%{
The rows of the table are reservoirs.  The columns are average with
uncertainty for each reservoir. 

%}

m1= 4; % see help "stat"
m2 = 3; 
% note: UC, sed, and LM will always fit using normal dist. maximum
% Liklihood fit because they are always normal (see Section 4)

% Crustal Sediments
weight = s1.mass.mass(cc,1) + s2.mass.mass(cc,1) + s3.mass.mass(cc,1); 
w1 = s1.mass.mass(cc,1)./weight; % weigh each layer by mass so we can combine
w2 = s2.mass.mass(cc,1)./weight; 
w3 = s3.mass.mass(cc,1)./weight; 
sed.thick = s1.thick(cc,:) + s2.thick(cc,:) + s3.thick(cc); 
sed.rho = s1.rho(cc,:).*w1 + s2.rho(cc).*w2 + s3.rho(cc).*w3; 
sed.mass.mass = s1.mass.mass(cc) + s2.mass.mass(cc) + s3.mass.mass(cc); 
sed.mass.U = s1.mass.U(cc) + s2.mass.U(cc) + s3.mass.U(cc); 
sed.mass.Th = s1.mass.Th(cc) + s2.mass.Th(cc) + s3.mass.Th(cc);
sed.mass.K = s1.mass.K(cc) + s2.mass.K(cc) + s3.mass.K(cc);
sed_sums = s1_cc_sums + s2_cc_sums + s3_cc_sums; 

% Oceanic Sediments
weight = s1.mass.mass(oc,1) + s2.mass.mass(oc,1) + s3.mass.mass(oc,1); 
w1 = s1.mass.mass(oc,1)./weight; % weigh each layer by mass so we can combine
w2 = s2.mass.mass(oc,1)./weight; 
w3 = s3.mass.mass(oc,1)./weight; 
sed.oc.thick = s1.thick(oc) + s2.thick(oc,:) + s3.thick(oc,:); 
sed.oc.rho = s1.rho(oc,1).*w1 + s2.rho(oc,1).*w2 + s3.rho(oc,1).*w3; 
sed_oc_sums = s1_oc_sums + s2_oc_sums + s3_oc_sums; 

% Oceanic Crust
weight = UC.mass.mass(oc,1) + MC.mass.mass(oc,1) + LC.mass.mass(oc,1); 
w1 = UC.mass.mass(oc,1)./weight; % weigh each layer by mass so we can combine
w2 = MC.mass.mass(oc,1)./weight; 
w3 = LC.mass.mass(oc,1)./weight; 
ocean.thick = UC.thick(oc) + MC.thick(oc) + LC.thick(oc);
ocean.rho = UC.rho(oc).*w1 + MC.rho(oc).*w2 + LC.rho(oc).*w3; 

ocean_sums = UC_oc_sums + MC_oc_sums + LC_oc_sums; 


% Mantle Stats
if simple2.calcMantleLayered == true
    mantle.dm.thick = 2891-750-Litho1.LAB; 
    mantle.em.thick = 750; 
    mantle.dm.rho = PREM(25:2891-750,3); 
    mantle.em.rho = PREM(2891-750:2891,3);     
else
    mantle.dm.thick = 2891 - Litho1.LAB;
    mantle.em.thick = 0;
    mantle.dm.rho = PREM(25:2891,3); 
    mantle.em.rho = 0;       
end

% BSE Stats
bse_sums = s1_cc_sums + s2_cc_sums + s3_cc_sums + UC_cc_sums + MC_cc_sums...
    + LC_cc_sums + LM_cc_sums + s1_oc_sums + s2_oc_sums + s3_oc_sums + UC_oc_sums...
    + MC_oc_sums + LC_oc_sums + man_sums_dm + man_sums_em; 
bse.thick = 2891; %km
bse.rho = bse_sums(:,1)./((4/3*pi*6378137^3)-(4/3*pi*(6378137-2891000)^3));

% Continental Crust Stats
crust_sums = s1_cc_sums + s2_cc_sums + s3_cc_sums + UC_cc_sums + MC_cc_sums...
    + LC_cc_sums; 
crust.thick = s1.thick(cc) + s2.thick(cc) + s3.thick(cc) + UC.thick(cc)...
    + MC.thick(cc) + LC.thick(cc); 
x = s1.mass.mass(cc,1) + s2.mass.mass(cc,1) + s3.mass.mass(cc,1) + UC.mass.mass(cc,1)...
    + MC.mass.mass(cc,1) + LC.mass.mass(cc,1); 
crust.rho = s1.rho(cc).*(s1.mass.mass(cc,1)./x) + s2.rho(cc).*(s2.mass.mass(cc,1)./x)...
    + s3.rho(cc).*(s3.mass.mass(cc,1)./x) + UC.rho(cc).*(UC.mass.mass(cc,1)./x)...
    + MC.rho(cc).*(MC.mass.mass(cc,1)./x) + LC.rho(cc).*(LC.mass.mass(cc,1)./x); 






% Row Names
huang.tab3.rows = {'Sed';'UC';'MC';'LC';'Bulk CC';'LM';'Sed';'C';'DM';'EM';'BSE'};


% - Density (rho; g/cm3) (use 'find' to only get non-zero values as to not mess up stats)
    x = [stat(sed.rho(find(sed.rho)),m2);...
         stat(UC.rho(cc),m2);...
         stat(MC.rho(cc),m2);...
         stat(LC.rho(cc),m2);...
         stat(crust.rho,m2);...
         stat(LM.rho(cc),m2);...
         stat(sed.oc.rho(find(sed.oc.rho)),m2);...
         stat(ocean.rho,m2);...
         stat(mantle.dm.rho,m2);...
         stat(mantle.em.rho,m2);...
         stat(bse.rho,m2)]/1000; %(g/cm3) 
    huang.tab3.rho = x; 

% - Thickness (km) (use 'find' to only get non-zero values as to not mess up stats)
    x = [stat(sed.thick(find(sed.thick)),m2);...
        stat(UC.thick(cc),m2);...
        stat(MC.thick(cc),m2);...
        stat(LC.thick(cc),m2);...
        stat(crust.thick,m2);...
        stat(LM.thick(cc),m2);...
        stat(sed.oc.rho,m2);...
        stat(ocean.thick,m2);...
        stat(mantle.dm.thick(cc),m2);...
        stat(mantle.em.thick,m2);...
        stat(bse.thick,m2)]/1000; %(km)
    huang.tab3.thick = x; 

% - Mass (10^21 kg)
    x = [stat(sed_sums(:,1),m2);...
        stat(UC_cc_sums(:,1),m2);...
        stat(MC_cc_sums(:,1),m2);...
        stat(LC_cc_sums(:,1),m2);...
        stat(crust_sums(:,1),m2);...
        stat(LM_cc_sums(:,1),m2);...
        stat(sed_oc_sums(:,1),m2);...
        stat(ocean_sums(:,1),m2);...
        stat(man_sums_dm(:,1),m2);...
        stat(man_sums_em(:,1),m2);...
        stat(bse_sums(:,1),m2)]/10^21; %(10^21 kg)
    huang.tab3.mass.mass = x; 

% - Abundance of U (ppm; ug/g)
    x = [stat(sed_sums(:,2)./sed_sums(:,1),m2);...
        stat(UC_cc_sums(:,2)./UC_cc_sums(:,1),m2);...
        stat(MC_cc_sums(:,2)./MC_cc_sums(:,1),m1);...
        stat(LC_cc_sums(:,2)./LC_cc_sums(:,1),m1);...
        stat(crust_sums(:,2)./crust_sums(:,1),m1);...
        stat(LM_cc_sums(:,2)./LM_cc_sums(:,1),m1);...
        stat(sed_oc_sums(:,2)./sed_oc_sums(:,1),m2);...
        stat(ocean_sums(:,2)./ocean_sums(:,1),m2);...
        stat(man_sums_dm(:,2)./man_sums_dm(:,1),m2);...
        stat(man_sums_em(:,2)./man_sums_em(:,1),m2);...
        stat(bse_sums(:,2)./bse_sums(:,1),m1)]*10^6; %(ppm)
    huang.tab3.abund.U = x;    

% - Abundance of Th (ppm; ug/g)
    x = [stat(sed_sums(:,3)./sed_sums(:,1),m2);...
        stat(UC_cc_sums(:,3)./UC_cc_sums(:,1),m2);...
        stat(MC_cc_sums(:,3)./MC_cc_sums(:,1),m1);...
        stat(LC_cc_sums(:,3)./LC_cc_sums(:,1),m1);...
        stat(crust_sums(:,3)./crust_sums(:,1),m1);...
        stat(LM_cc_sums(:,3)./LM_cc_sums(:,1),m1);...
        stat(sed_oc_sums(:,3)./sed_oc_sums(:,1),m2);...
        stat(ocean_sums(:,3)./ocean_sums(:,1),m2);...
        stat(man_sums_dm(:,3)./man_sums_dm(:,1),m2);...
        stat(man_sums_em(:,3)./man_sums_em(:,1),m2);...
        stat(bse_sums(:,3)./bse_sums(:,1),m1)]*10^6; %(ppm)
    huang.tab3.abund.Th = x;   

% - Abundance of K (wt%)
    x = [stat(sed_sums(:,4)./sed_sums(:,1),m2);...
        stat(UC_cc_sums(:,4)./UC_cc_sums(:,1),m2);...
        stat(MC_cc_sums(:,4)./MC_cc_sums(:,1),m1);...
        stat(LC_cc_sums(:,4)./LC_cc_sums(:,1),m1);...
        stat(crust_sums(:,4)./crust_sums(:,1),m1);...
        stat(LM_cc_sums(:,4)./LM_cc_sums(:,1),m1);...
        stat(sed_oc_sums(:,4)./sed_oc_sums(:,1),m2);...
        stat(ocean_sums(:,4)./ocean_sums(:,1),m2);...
        stat(man_sums_dm(:,4)./man_sums_dm(:,1),m2);...
        stat(man_sums_em(:,4)./man_sums_em(:,1),m2);...
        stat(bse_sums(:,4)./bse_sums(:,1),m1)]*10^6; %(ppm)
    huang.tab3.abund.K = x;   
    
% - Mass of U (10^15 kg)
    x = [stat(sed_sums(:,2),m2);...
        stat(UC_cc_sums(:,2),m2);...
        stat(MC_cc_sums(:,2),m1);...
        stat(LC_cc_sums(:,2),m1);...
        stat(crust_sums(:,2),m1);...
        stat(LM_cc_sums(:,2),m1);...
        stat(sed_oc_sums(:,2),m2);...
        stat(ocean_sums(:,2),m2);...
        stat(man_sums_dm(:,2),m2);...
        stat(man_sums_em(:,2),m2);...
        stat(bse_sums(:,2),m1)]/10^15; %(10^15 kg)
    huang.tab3.mass.U = x; 
    
% - Mass of Th (10^15 kg)
    x = [stat(sed_sums(:,3),m2);...
        stat(UC_cc_sums(:,3),m2);...
        stat(MC_cc_sums(:,3),m1);...
        stat(LC_cc_sums(:,3),m1);...
        stat(crust_sums(:,3),m1);...
        stat(LM_cc_sums(:,3),m1);...
        stat(sed_oc_sums(:,3),m2);...
        stat(ocean_sums(:,3),m2);...
        stat(man_sums_dm(:,3),m2);...
        stat(man_sums_em(:,3),m2);...
        stat(bse_sums(:,3),m1)]/10^15; %(10^15 kg)
    huang.tab3.mass.Th = x; 

% - Mass of K (10^19 kg)
    x = [stat(sed_sums(:,4),m2);...
        stat(UC_cc_sums(:,4),m2);...
        stat(MC_cc_sums(:,4),m1);...
        stat(LC_cc_sums(:,4),m1);...
        stat(crust_sums(:,4),m1);...
        stat(LM_cc_sums(:,4),m1);...
        stat(sed_oc_sums(:,4),m2);...
        stat(ocean_sums(:,4),m2);...
        stat(man_sums_dm(:,4),m2);...
        stat(man_sums_em(:,4),m2);...
        stat(bse_sums(:,4),m1)]/10^19/K.b; %(10^19 kg)
    huang.tab3.mass.K = x; 

% - Heat production (TW; 10^12 W)
      x = [stat(sed_sums(:,6),m2);...
        stat(UC_cc_sums(:,6),m2);...
        stat(MC_cc_sums(:,6),m1);...
        stat(LC_cc_sums(:,6),m1);...
        stat(crust_sums(:,6),m1);...
        stat(LM_cc_sums(:,6),m1);...
        stat(sed_oc_sums(:,6),m2);...
        stat(ocean_sums(:,6),m2);...
        stat(man_sums_dm(:,6),m2);...
        stat(man_sums_em(:,6),m2);...
        stat(bse_sums(:,6),m2)]/10^12; %(10^12 W (TW))
    huang.tab3.hp = x;
    

% Horizontal concatenate and print 

    huang.tab3.table = horzcat(huang.tab3.rho, huang.tab3.thick, huang.tab3.mass.mass, huang.tab3.abund.U,...
        huang.tab3.abund.Th, huang.tab3.abund.K, huang.tab3.mass.U, huang.tab3.mass.Th,...
        huang.tab3.mass.K, huang.tab3.hp); 


   
%% Write Model information to String
% Re-order MASTER structure alphebetically
MASTER = orderfields(MASTER);

% Record Method used
if strcmp(MASTER.method, 'Huang et al. 2013') == 1
    m = 'H13Meth'; 
else
    m = 'BivarMeth';
end

% Record if calculated flux or not (and what detector)
if size(det,2)>1
    d = sprintf('flux%s',det.Properties.RowNames{1}); 
else
    d = 'noFlux';
end

% Record type of mantle
if MASTER.mantle.calc == 'true'
    if MASTER.mantle.Style == 'Layered'
        e = 'CalcMantleLayered';
    else
        e = 'CalcMantleHomo';
    end
else
    e = 'noCalcMantle';
end


sprintf('Time Elapsed: %.1f min |  Det: = %s | Meth = %s \n',toc/60, det.Properties.RowNames{1},MASTER.method)
str1 = sprintf('Results_%1.1eIter_%s_%s_%s_%s_%s',iter,MASTER.model,m,d,e,datestr(date,'ddmmmyyyy'));

MASTER.save.name = str1; MASTER.save.locate = pwd; 


%% %% 13) ---- Save data ---
    str = strcat(str1,'.mat');
    save(str,'MASTER','huang','Litho1','s1','s2','s3','UC','MC','LC','LM','mantle','flux','cc','oc','str1','iter','simple2')
   % clearvars('-except','MASTER','huang','Litho1','s1','s2','s3','UC','MC','LC','LM','flux','cc','oc'); 
    datestr(now)


fprintf('Saved data: \n %s',str)


    
%% (latex_Table 1) Format 'huang.tab3' into Latex format and write to .txt ----
clear x r h y
h = huang.tab3; 



% Always put a '&' after so it goes to the next value
for i = 1:size(h.rows,1)
    % Row Name
if i == 1 || i == 2 || i == 3 || i == 4 || i == 5 || i == 6 || i == 7 || i == 8 %push it over 1 column
    x(i).row = sprintf(' & %s &',h.rows{i}); 
else 
    x(i).row = ''; % for DM,EM, and BSE
end

    % Density
x(i).rho = sprintf('& %0.2f &',h.rho(i,1));

    % Thickness
x(i).thick = sprintf(' %0.2f $%s$%0.1f &',h.thick(i,1),'\pm',h.thick(i,2));
     
    % Mass
x(i).mass.mass = sprintf(' %0.1f $%s$%0.1f &',h.mass.mass(i,1),'\pm',h.mass.mass(i,2));    

    % U abundance
x(i).abund.U = latexUncSymAsym(h.abund.U(i,:),2);
    
    % Th abundance
x(i).abund.Th = latexUncSymAsym(h.abund.Th(i,:),2);    

    % K abundance
x(i).abund.K = latexUncSymAsym(h.abund.K(i,:),2);
 
    % U Mass
x(i).mass.U = latexUncSymAsym(h.mass.U(i,:),1);

    % Th Mass
x(i).mass.Th = latexUncSymAsym(h.mass.Th(i,:),1);

    % K Mass
x(i).mass.K = latexUncSymAsym(h.mass.K(i,:),1);

    % Heat Production
x(i).hp = latexUncSymAsym(h.hp(i,:),1,1);

% Ammend latex format stuff onto front of strings
if i == 1
   x(i).front = '\multicolumn{1}{|c|}{\multirow{5}{*}{CC}}';
elseif i == 2 || i == 3 || i == 4 || i == 5 || i == 6 || i == 8
   x(i).front = '\multicolumn{1}{|c|}{} ';
elseif i == 7
   x(i).front = '\multicolumn{1}{|c|}{\multirow{2}{*}{OC}}';
elseif i == 9
   x(i).front = '\multicolumn{2}{|c}{DM} & '; 
elseif i == 10
   x(i).front = '\multicolumn{2}{|c}{EM} & ';   
elseif i == 11
   x(i).front = '\multicolumn{2}{|c}{BSE} & ';  
end

% Ammend latex format stuff onto end of strings
if i == 1 || i == 2 || i == 3 || i == 4|| i == 5 || i == 7
    x(i).end = '\\ \cline{2-13}'; 
else 
    x(i).end = '\\ \hline';
end

        
   y{i,1} = strcat(x(i).front,x(i).row,x(i).rho,x(i).thick,x(i).mass.mass,x(i).abund.U,...
       x(i).abund.Th,x(i).abund.K,x(i).mass.U,x(i).mass.Th,x(i).mass.K,...
       x(i).hp,x(i).end); 

end

% Create Name of file
str = strcat('t1_PhysProp_',str1,'.txt');

% Ammend final lines with caption of information (so we can identify the table)
final{1,1} = '\end{tabular}}';
x = strrep(str,'_',' '); x = strrep(x,'+','$+$'); 
final{2,1} = sprintf('\\caption{%s}',x);


% Write to text file
fid = fopen(str,'wt');
fprintf(fid, '%s\n', y{1,1});
fclose(fid);

fid = fopen(str,'at');
fprintf(fid, '%s\n', y{2:end,1});
%fclose(fid);

fid = fopen(str,'at');
fprintf(fid,'%s\n',final{:,:});

fclose(fid);    
fclose('all'); % Fully releases the file (otherwise Matlab won't for some reason)
    

disp('Saved Table 1 (geophysical information)')

%% (latex_Table 2) Format 'huang.tab2' into Latex format and write to .txt ----
clear x r h y
h = huang.tab2; 



% Always put a '&' after so it goes to the next value
for i = 1:size(h.rows,1)

    if i == 1 || i == 6 
        x(i).row = sprintf('%s &',strrep(h.rows{i},'_',' '));
    else
        x(i).row = sprintf('%s &',h.rows{i}); 
    end

    if i == 1 || i == 6 || i == 7
    x(i).sU = latexUncSymAsym(h.U238(i,:),2);
    x(i).sTh = latexUncSymAsym(h.Th232(i,:),2);
    x(i).sTotal = latexUncSymAsym(h.total(i,:),2,1);
    else
        x(i).sU = latexUncSymAsym(h.U238(i,:),1);
        x(i).sTh = latexUncSymAsym(h.Th232(i,:),1);
        x(i).sTotal = latexUncSymAsym(h.total(i,:),1,1);
    end


if i == size(h.rows,1)
    x(i).end = '\\ \hline';
else
    x(i).end = '\\';
end
        
   y{i,1} = strcat(x(i).row,x(i).sU,x(i).sTh,x(i).sTotal,x(i).end); 

end


% Create Name of file
str = strcat('t2_flux_',str1,'.txt');

% Ammend Initial lines with detector
z = MASTER.detector;
start{1,1} = sprintf('\\multirow{2}{*}{} & \\multicolumn{3}{c}{%s} \\\\',z.Properties.RowNames{1});
start{2,1} = sprintf('& \\multicolumn{3}{c}{%0.2f\\textdegree, %0.2f\\textdegree} \\\\ \\hline',z{1,2},z{1,1});
start{3,1} = ' & S(U) & S(Th) & S(U+Th) \\'; 

% Ammend final lines with caption of information (so we can identify the table)
final{1,1} = '\end{tabular}}';
z = strrep(str,'_',' '); z = strrep(z,'+','$+$'); 
final{2,1} = sprintf('\\caption{%s}',z);



% Write to text file
fid = fopen(str,'wt');
fprintf(fid,'%s\n',start{:,:}); 
fclose(fid);

fid = fopen(str,'at');
fprintf(fid, '%s\n', y{:,1});
%fclose(fid);

fid = fopen(str,'at');
fprintf(fid,'%s\n',final{:,:});

fclose(fid);    
fclose('all'); % Fully releases the file (otherwise Matlab won't for some reason)
disp('Saved Table 2 (flux)')


%% (figure) Flux vs Distance

% Sediment
count.sed.U238 = flux.s1.count.U238 + flux.s2.count.U238 + flux.s3.count.U238; 
count.sed.Th232 = flux.s1.count.Th232 + flux.s2.count.Th232 + flux.s3.count.Th232; 

% Crust
count.crust.U238 = flux.UC.count.U238 + flux.MC.count.U238 + flux.LC.count.U238;
count.crust.Th232 = flux.UC.count.Th232 + flux.MC.count.Th232 + flux.LC.count.Th232;

% Lithosphere
count.litho.U238 = count.sed.U238 + count.crust.U238 + flux.LM.count.U238; 
count.litho.Th232 = count.sed.Th232 + count.crust.Th232 + flux.LM.count.Th232; 

% Convecting Mantle 
count.man.U238 = flux.man.dm.count.U238+flux.man.em.count.U238; 
count.man.Th232 = flux.man.dm.count.Th232+flux.man.em.count.Th232; 

% BSE
count.bse.U238 = count.man.U238 + count.litho.U238;
count.bse.Th232 = count.man.Th232 + count.litho.Th232; 

% Totals
count.total.sed = count.sed.U238 + count.sed.Th232; 
count.total.crust = count.crust.U238 + count.crust.Th232; 
count.total.litho = count.litho.U238 + count.litho.Th232; 
count.total.man = count.man.U238 + count.man.Th232; 
count.total.bse = count.bse.U238 + count.bse.Th232; 


% Scale values to precisely match Monte Carlo results (values of "count"
% slightly off due to how its summed)
wt = huang.tab2.total(end,1)/sum(count.total.bse);




% Plot Figure
x = simple2.centers/1000; %(km)


fig = figure('visible','off');
set(gca, 'FontName', 'Times New Roman')
title(sprintf('%s Distance vs Cumulative Flux',MASTER.detector.Properties.RowNames{1}))
set(gca,'XScale','log') % x-axis = log scale


% TNU Signal
    yyaxis left
    plot(x,cumsum(count.total.sed),':g'); hold on
    plot(x,cumsum(count.total.crust),'--','Color',[0.82031 0.41016 0.11719]);
    plot(x,cumsum(count.total.man),'-.r');
    plot(x,cumsum(count.total.bse),'-k');
    xlabel('Distance From Detector (km)')
    ylabel('Cumulative Flux (TNU)')
    set(findall(gca, 'Type', 'Line'),'LineWidth',4); %line width
    set(gca,'FontSize',35) % axes font size
    set(gca,'XMinorTick','on','YMinorTick','on') %set ticks on
    axis([0 max(x) 0 inf])
    set(gca,'YColor','k') %set axis colors to black
    grid on
    h = legend('Sediment','Crust','Mantle','Total','Location','northwest');
    h.FontSize = 40;

% Cumulative percent total signal
    yyaxis right
    plot(x,cumsum(count.total.bse)./sum(count.total.bse),'-k');
    set(gca,'YTickLabel',0:10:100) % set right y tick labels (0 - 100 percent)
    set(gca,'YColor','k') %set axis colors to black

 
 
% Save Figure
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperPosition', [0 0 20 13]); % size of image (0 0 width length)  
    str_png = strcat('f1_DistVsFlux_',str1,'.png');
    str_eps = strcat('f1_DistVsFlux_',str1,'.eps');
    print( str_eps,'-depsc');
    print(str_png,'-dpng','-r300')

disp('Saved Figure 1 (Distance vs Flux)')


   end
end

