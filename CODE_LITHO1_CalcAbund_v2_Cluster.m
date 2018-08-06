%% Information
%{
 Calculate the abundance of Heat Producing Elements in crust from LITHO1.0
% Scott A. Wipperfurth
% updated March 2018

% see http://igppweb.ucsd.edu/~gabi/crust1.html for more info


 File Litho1.types

crust1.bnds = boundary topography
Litho1.Vp = vp
crust.vs = vs
crust1.Litho1.rho = Litho1.rho 

Output File structure: 
xyz-th = thickness of layer
xyz-bd = Litho1.elevation of each layer (bottom of layer and top of water)
xyz-Litho1.rho = density of layer
xyz-vp = vp
xyz-vs = vs

The 8 crustal layers:
====================
1) water             
2) ice               
3) upper sediments   (VP, VS, Litho1.rho not defined in all cells) 
4) middle sediments  "
5) lower sediments   "
6) upper crystalline crusts
7) middle crystalline crust
8) lower crystalline crust

a ninth layer gives V_Pn, V_Sn and Litho1.rho below the Moho. The values
are associated with LLNL model G3Cv3 on continents and a thermal
model in the oceans.

The model is defined from 89.5 to -89.5 deg latitude and -179.5 to 179.5 deg
longitude.


--- REFERENCES CITED IN CODE: ---

Huang, Y., Chubakov, V., Mantovani, F., Rudnick, R.L., McDonough, W.F., 2013.
    A reference Earth model for the heat-producing elements and associated 
    geoneutrino flux. Geochem. Geophys. Geosystems 14, 2003–2029.
    https://doi.org/10.1002/ggge.20129

Laske, G., Masters, G., Ma, Z., Pasyanos, M.E., 2013. Update on CRUST1.0
    - A 1-degree Global Model of Earth’s Crust, in: Geophysical Research 
    Abstracts, EGU2013-2658. Presented at the EGU.

Olugboji, T.M., Lekic, V., McDonough, W., 2017. A statistical assessment of
     seismic models of the U.S. continental crust using Bayesian inversion of 
    ambient noise surface wave dispersion data. Tectonics 2017TC004468. 
    https://doi.org/10.1002/2017TC004468

Pasyanos, M.E., Masters, T.G., Laske, G., Ma, Z., 2014. LITHO1.0: An 
    updated crust and lithospheric model of the Earth: LITHO1.0. 
    Journal of Geophysical Research: Solid Earth 119, 2153–2173. 
    https://doi.org/10.1002/2013JB010626

Plank, T., 2014. The Chemical Composition of Subducting Sediments, 
    in: Treatise on Geochemistry. Elsevier, pp. 607–629.

Rudnick, R.L., Gao, S., 2014. Composition of the Continental Crust, 
    in: Treatise on Geochemistry. Elsevier, pp. 1–51. 

Ruedas, T., 2017. Radioactive heat production of six geologically important
    nuclides. Geochem. Geophys. Geosystems 18, 3530–3541. 
    https://doi.org/10.1002/2017GC006997

White, W.M., Klein, E.M., 2014. Composition of the Oceanic Crust,
    in: Treatise on Geochemistry. Elsevier, pp. 457–496. 
    https://doi.org/10.1016/B978-0-08-095975-7.00315-6



OTHER: 
NASA Earth Fact Sheet: https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html

LITHO1.0 Download Page: https://igppweb.ucsd.edu/~gabi/crust1.html


Stats: 

Abundance: 

10000 iterations = ~25 minutes



Abundance + GeoFlux: 
1000 iterations = 28 minutes (with coarse gridding)
50 iterations = 5.3 minutes (with coarse gridding)
1 iteration = 3.5 minutes (with coarse gridding)
1 iteration = 3.5 minutes (with coarse gridding)


30,000 iterations = > 12 hours (with coarse gridding...)

%}

%% 1) ---- LOAD DATA FROM LITHO 1.0 ----
for all_det = 1:1
fprintf('    Starting time = %s \n',datestr(now,'mmmm dd, yyyy HH:MM AM')) %print time of loop start


%close all; 
tic; clearvars -except all_det
%cd(fileparts(matlab.desktop.editor.getActiveFilename));% moves to current folder


   MASTER.StartTime = datestr(now,'mmmm dd, yyyy HH:MM AM');


%{
% -- Create  or check for parallel Pool --
    p = gcp('nocreate'); % Get info about paralel pool. If pool already, do not create new one.
    poolsize = 8; % 8 = fastest (~~40% faster than 4 cores w/ 1000 iterations)
    MASTER.poolsize = poolsize; %record data into "MASTER" file with all run information
    if isempty(p)
        parpool('local',poolsize);
        MASTER.pool = gcp; 
    elseif p.NumWorkers ~= poolsize
        delete(gcp('nocreate'))
        parpool('local',poolsize);
        MASTER.pool = gcp; 
    else
        poolsize = p.NumWorkers; 
        MASTER.pool = gcp; 
    end
    fprintf('Size of parallel pool: %d workers \n',poolsize); clear poolsize p; 
%}
    
    
% -- Prompt input for # of iterations --
%{
        Int: 
            1   |   No Monte-Carlo, only central value
           >1   |   Monte-Carlo
%}

    detectors = {137.31, 36.43, 1000, 0.7, 5.98*10^31 ;... %KamLAND - Araki etal. (2005) 
                 13.57, 42.45, 1400, 0.842, 9.76E+30;...   %Borexino - Wikipedia LNGS
                 -81.201, 46.475, 2092, 0.8, 5.76E+31;...  %SNO+ - Chen (2006), Andringa etal. (2016)
                 112.518, 22.118, 700, 0.8, 1.29E+33;...   %JUNO - An etal. (2016) JUNO Yellowbook
                 101.71, 28.15, 2400, 0.8, 2.16E+32;...    %Jinping - Beacom etal. (2016) letter of intent
                 -156.32, 19.72, 0, 0.8, 0;         %Hawaii
                 }; 
    detectors = cell2table(detectors);
    detectors.Properties.RowNames = {'KamLAND','Borexino','SNO','JUNO','Jinping','Hawaii'};
    detectors.Properties.VariableNames = {'Longitude','Latitude','Depth','Efficiency','Protons'};




%%%%%%%%%%%%%%%%%%%%%%% CHANGE THESE FOR CLUSTER ANALYSIS %%%%%%%%%%%%%%%%%














% loop through detectors by detectors(i,:), which pulls entire row




    iter = 1000;  
    simple2.meth = 1; 
    %det = detectors('Borexino',:); 
    det = detectors(all_det,:); 
    %det = 0;
    %abundance = 1; % 1 = don't calculate abundance again
    
    MASTER.detector = det; 
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    MASTER.iter = iter; 
    

% -- Load Bivariate Data --
    simple2.bivar = load('BivarData_04042018.mat'); 

     
    clear x detectors selection p poolsize
    
                    
                   
    if simple2.meth ==1
        MASTER.method = 'Huang et al. 2013'; %record method to variable with run information
    elseif simple2.meth ==2
        MASTER.method = 'Bivariate Analysis ';
    end
                   
% - Define energies we will calculate flux at -
    % this has to be here so we can pre-define variables
    E_int = 100; %(keV)
    simple2.energy(:,1) = 1806+E_int/2:E_int:3276-E_int/2; %(keV) 3272 = max but this is easier
    MASTER.energy.bin = simple2.energy; 
    MASTER.energy.binsize = '100 KeV';


% -- Load Davies 2013 data (in 1x1 degree cells) --
    %Davies = load('Davies2013_HeatFlux_Data_Table1_Eq_lonlat_64800cells.csv');
    %Davies = Davies(:,3);


% -- Load LITHO1.0 Data -- 
    load('LITHO1_BaseData_OndrejReFormat.mat')
    MASTER.model = 'LITHO1.0';
    

% ----------  Logically index cc, oc, stable, and Archean crust  --------
    %{ 
    Assign a value of 1 if the cell has a crust type that we designate as cc,
    oc, stable, or Archean crust.  Designation as part of one of these 4 broad 
    categories is performed by us based on the names given by LITHO1.0. It is
    not directly from LITHO1.0. All other values are 0. 

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
    oc_types = [26,27,28,31,36];  
    oc = any(bsxfun(@eq,oc_types,Litho1.type(:,3)),2); clear oc_types
    cc = ~oc; 

    stable_types = [3,4,5,6,7,8,9];
    stable = any(bsxfun(@eq,stable_types,Litho1.type(:,3)),2); clear stable_types

    archean_types = [3,4,5];
    archean = any(bsxfun(@eq,archean_types,Litho1.type(:,3)),2); clear archean_types
   
    
    
 clear data_load

%% 2) ---- Organize data and Pre-allocate Variables ----
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
    [s1_cc_sums,s2_cc_sums,s3_cc_sums,UC_cc_sums,MC_cc_sums,LC_cc_sums,LM_cc_sums] ...
        = deal(zeros(iter,9));   
    
    [s1_oc_sums,s2_oc_sums,s3_oc_sums,UC_oc_sums,MC_oc_sums,LC_oc_sums,LM_oc_sums] ...
        = deal(zeros(iter,9));     
    
    [s1_cc_flux_sums,s2_cc_flux_sums,s3_cc_flux_sums,UC_cc_flux_sums,MC_cc_flux_sums,LC_cc_flux_sums,LM_cc_flux_sums] ...
        = deal(zeros(iter,length(simple2.energy(:,1))*2));   
    
    [s1_oc_flux_sums,s2_oc_flux_sums,s3_oc_flux_sums,UC_oc_flux_sums,MC_oc_flux_sums,LC_oc_flux_sums,LM_oc_flux_sums] ...
        = deal(zeros(iter,length(simple2.energy(:,1))*2));     
    

    %{
    %- stats variables
    [s1_output_stat,s2_output_stat,s3_output_stat,UC_output_stat,LM_output_stat]...
        = deal(zeros(length(Litho1.latlon),18));  
    
    [s1_oc_output_stat,s2_oc_output_stat,s3_oc_output_stat,UC_oc_output_stat,...
        MC_oc_output_stat,LC_oc_output_statLM_oc_output_stat]...
        = deal(zeros(length(Litho1.latlon),18));      
    
    [MC_output_stat,LC_output_stat] = deal(zeros(length(Litho1.latlon),36)); % MC LC cc
 %}
clear n x i

%% 3) ---- Define correlation ----
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
    
%% 4) ---- Assign and define constant parameters ----

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
    
    % -- Mantle abundances -- (
   % man.dm.aU = 
    
    
    % -- BSE abundances -- (Arevalo et al. 2013 (U + Th unc. from Wipperfurth et al. 2018)
    BSE.aU = repmat([20.3,2.03,2.03]*10^-9,length(UC.lat),1); %(kg/kg) 
    BSE.aTh = repmat([79.5,7.95,7.95]*10^-9,length(UC.lat),1); %(kg/kg)
    BSE.aK = repmat([240,48,48]*10^-9,length(UC.lat),1); %(kg/kg)
    
    
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
        
    
% -- Define heat Production from kg of isotope -- (Ruedas 2017)
    simple2.hp.U238 = 94.946*10^-6;  %(Watt/kg)
    simple2.hp.U235 = 568.402*10^-6; %(Watt/kg)
    simple2.hp.U    = 98.314*10^-6;  %(Watt/kg)
    simple2.hp.Th   = 26.368*10^-6;  %(Watt/kg)
    simple2.hp.K40  = 28.761*10^-6;  %(Watt/kg) K40, not K 
    
    
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


% - Load and bin geonuetrino spectrum - (Enomoto et al. 2007)
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
        
    
    
    
    
clear prot Tmax x E enomoto top bot E_int

%% 5) ---- Perform non-loop calculations ----

% -- Find radius at detector location -- (m) (radius Earth at point 
if size(det,2)>1
    x = knnsearch([Litho1.latlon(:,1),Litho1.latlon(:,2)],[det{1,1},det{1,2}]);
    det{1,3} = Litho1.r(x) - det{1,3}; clear x 
    det.Properties.VariableNames = {'Longitude','Latitude','Radius','Efficiency','Proton'};
end
    
%{
% -- Calculate initial Gridding of cells (for geoflux calc) --
x = dis(lon, lat, det, rad, 8389.2);

grid.num = [4000,15000,25000,40000,60000]; %[15000,20000,50000,100000,150000]; %(m) ~size of grid cells  [3000,5000,10000,50000,90000]
grid.lim = [100000,160000,280000,600000,1000000]; %(m) distance from detector for each grid size
    
        for i = 1:length(x)
        if x(i) <= grid.lim(1)
            geo(2(i) = grid.num(1); % num = desired cell size
        elseif x(i) > grid.lim(1) && x(i) <= grid.lim(2)
            size(i) = grid.num(2);
        elseif x(i) > grid.lim(2) && x(i) <= grid.lim(3)
            size(i) = grid.num(3);
        elseif x(i) > grid.lim(3) && x(i) <= grid.lim(4)
            size(i) = grid.num(4);
        elseif x(i) > grid.lim(4) && x(i) <= grid.lim(5)
            size(i) = grid.num(5); 
    elseif eqn(i) >grid.lim(5)
        size(i) = 100000;
    end  
        end
%}

%% 7) ---- MONTE CARLO  ----
%{
The monte carlo is used to calculate the mass distribution and heat production
for each layer. For the middle and lower crust, the abundance of U, Th, and
K is calculated from comparison of seismic velocity and petrologic
constraints, then the heat production is calculated. 
%}
%disp('Starting Monte carlo')

rng shuffle %reseed random number generator based on time (DO NOT PUT INTO LOOP, time intensive and maybe logically not correct)

tic

parfor n = 1:length(Litho1.latlon) % n = a specific cell (out of 64,800)
%percent_done(n,length(Litho1.latlon),5)

temp_P = zeros(iter,1); %temporary pressure, reset every new cell (otherwise it will continue summing between cells


    if cc(n) == 1 % cc(n) = 1 is continental crust, 0 is oceanic crust
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
               [s1_out(n), s1_output_sums,s1_geoResponse(n,:),s1_flux(n),s1_flux_sums,P]...
                   = Huang13_cluster(n,iter,'s1',s1,simple2,cor.s1,Litho1.r(n),det,temp_P);
              s1_cc_sums = s1_cc_sums + s1_output_sums;
              s1_cc_flux_sums = s1_cc_flux_sums + s1_flux_sums;
              temp_P = temp_P + P; % add up pressure
            
         % -  Middle Sediment  (s2; layer 4)    
               [s2_out(n), s2_output_sums,s2_geoResponse(n,:),s2_flux(n),s2_flux_sums,P]...
                   = Huang13_cluster(n,iter,'s2',s2,simple2,cor.s2,Litho1.r(n),det,temp_P);
              s2_cc_sums = s2_cc_sums + s2_output_sums;  
              s2_cc_flux_sums = s2_cc_flux_sums + s2_flux_sums;
              temp_P = temp_P + P; % add up pressure
      
         % -  Lower Sediment  (s3; layer 5)    
               [s3_out(n), s3_output_sums,s3_geoResponse(n,:),s3_flux(n),s3_flux_sums,P]...
                   = Huang13_cluster(n,iter,'s3',s3,simple2,cor.s3,Litho1.r(n),det,temp_P);
              s3_cc_sums = s3_cc_sums + s3_output_sums;               
              s3_cc_flux_sums = s3_cc_flux_sums + s3_flux_sums;
              temp_P = temp_P + P; % add up pressure
           
         % - Upper Crust (UC; layer 6)      
               [UC_out(n), UC_output_sums,UC_geoResponse(n,:),UC_flux(n),UC_flux_sums,P]...
                   = Huang13_cluster(n,iter,'UC',UC,simple2,cor.UC,Litho1.r(n),det,temp_P);
              UC_cc_sums = UC_cc_sums + UC_output_sums; 
              UC_cc_flux_sums = UC_cc_flux_sums + UC_flux_sums;
              temp_P = temp_P + P; % add up pressure
              
         % - Middle Crust (MC; layer 7)      
               [MC_out(n), MC_output_sums,MC_geoResponse(n,:),MC_flux(n),MC_flux_sums,P]...
                   = Huang13_cluster(n,iter,'MC',MC,simple2,cor.MC,Litho1.r(n),det,temp_P);
              MC_cc_sums = MC_cc_sums + MC_output_sums; 
              MC_cc_flux_sums = MC_cc_flux_sums + MC_flux_sums;
              temp_P = temp_P + P; % add up pressure
              
         % - Lower Crust (LC; layer 8)      
               [LC_out(n), LC_output_sums,LC_geoResponse(n,:),LC_flux(n),LC_flux_sums,P]...
                   = Huang13_cluster(n,iter,'LC',LC,simple2,cor.LC,Litho1.r(n),det,temp_P);
              LC_cc_sums = LC_cc_sums + LC_output_sums;   
              LC_cc_flux_sums = LC_cc_flux_sums + LC_flux_sums;
              temp_P = temp_P + P; %moho(n).pressure = stat(P); % Pressure at base of LC (i.e. moho)
  
         % - Lithospheric Mantle (LM; layer 9)      
               [LM_out(n), LM_output_sums,LM_geoResponse(n,:),LM_flux(n),LM_flux_sums,P]...
                   = Huang13_cluster(n,iter,'LM',LM,simple2,cor.LM,Litho1.r(n),det,temp_P);
               LM_cc_sums = LM_cc_sums + LM_output_sums;  
               LM_cc_flux_sums = LM_cc_flux_sums + LM_flux_sums; 
               %LAB(n).pressure = stat(temp_P + P); % pressure at base of lithosphere (i.e. LAB)
 
	
    else  % ----  CALCULATE ATTRIBUTES FOR OCEANIC CRUST --------
   
                 % -  Upper Sediment  (s1; layer 3)    
               [s1_out(n), s1_output_sums,s1_geoResponse(n,:),s1_flux(n),s1_flux_sums,P]...
                   = Huang13_cluster(n,iter,'s1',s1,simple2,cor.s1,Litho1.r(n),det,temp_P);
               s1_oc_sums = s1_oc_sums + s1_output_sums;
               s1_oc_flux_sums = s1_oc_flux_sums + s1_flux_sums; 
               temp_P = temp_P + P; % add up pressure
   
         % -  Middle Sediment  (s2; layer 4)    
               [s2_out(n), s2_output_sums,s2_geoResponse(n,:),s2_flux(n),s2_flux_sums,P]...
                   = Huang13_cluster(n,iter,'s2',s2,simple2,cor.s2,Litho1.r(n),det,temp_P);
               s2_oc_sums = s2_oc_sums + s2_output_sums;               
               s2_oc_flux_sums = s2_oc_flux_sums + s2_flux_sums;    
               temp_P = temp_P + P; % add up pressure
               
         % -  Lower Sediment  (s3; layer 5)    
               [s3_out(n), s3_output_sums,s3_geoResponse(n,:),s3_flux(n),s3_flux_sums,P]...
                   = Huang13_cluster(n,iter,'s3',s3,simple2,cor.s3,Litho1.r(n),det,temp_P);
               s3_oc_sums = s3_oc_sums + s3_output_sums;  
               s3_oc_flux_sums = s3_oc_flux_sums + s3_flux_sums; 
               temp_P = temp_P + P; % add up pressure
              
         % - Upper Crust (UC; layer 6)      
               [UC_out(n), UC_output_sums,UC_geoResponse(n,:),UC_flux(n),UC_flux_sums,P]...
                   = Huang13_cluster(n,iter,'UC',UC,simple2,cor.UC,Litho1.r(n),det,temp_P);
              UC_oc_sums = UC_oc_sums + UC_output_sums; 
              UC_oc_flux_sums = UC_oc_flux_sums + UC_flux_sums; 
               temp_P = temp_P + P; % add up pressure   
                 
         % - Middle Crust (MC; layer 7)  (put layer "MC_oc" so it doesnt do
            % abundance calculation) "MC_oc_output_stat" is different size than "MC_output_stat"
               [MC_out(n), MC_output_sums,MC_geoResponse(n,:),MC_flux(n),MC_flux_sums,P]...
                   = Huang13_cluster(n,iter,'MC_oc',MC,simple2,cor.MC,Litho1.r(n),det,temp_P);
               MC_oc_sums = MC_oc_sums + MC_output_sums; 
               MC_oc_flux_sums = MC_oc_flux_sums + MC_flux_sums; 
               temp_P = temp_P + P; % add up pressure             
                  
         % - Lower Crust (LC; layer 8)      (put layer "LC_oc" so it doesnt do
            % abundance calculation)"LC_oc_output_stat" is different size than "LC_output_stat"
               [LC_out(n), LC_output_sums,LC_geoResponse(n,:),LC_flux(n),LC_flux_sums,P]...
                   = Huang13_cluster(n,iter,'LC_oc',LC,simple2,cor.LC,Litho1.r(n),det,temp_P);
               LC_oc_sums = LC_oc_sums + LC_output_sums; 
               LC_oc_flux_sums = LC_oc_flux_sums + LC_flux_sums;
               temp_P = temp_P + P; 
               %moho(n).pressure = temp_P + P; % Pressure at base of LC (i.e. moho)
               
        % - Lithospheric Mantle (LM; layer 9) (need layer = 'LM_oc' so it
        %            fills it in blank
               [LM_out(n), LM_output_sums,LM_geoResponse(n,:),LM_flux(n),LM_flux_sums,P]...
                   = Huang13_cluster(n,iter,'LM_oc',LM,simple2,cor.LM,Litho1.r(n),det,temp_P);
               LM_oc_sums = LM_oc_sums + LM_output_sums; 
               LM_oc_flux_sums = LM_oc_flux_sums + LM_flux_sums;  
               %LAB(n).pressure = moho(n).pressure + P; % pressure at base of lithosphere (i.e. LAB)


    end % end of if else statement
       
end % 1:64800



MASTER.EndTime = datestr(now,'mmmm dd, yyyy HH:MM AM');
MASTER.MCRunTime = sprintf('%.1f minutes',toc/60);
MASTER.memory = memory; 
MASTER.numCells = length(Litho1.latlon); 

%% 
%{

DO NOT DELETE!!! This section defines variables used in the Huang13
function and allows for easy testing.  The n = 17114 is the cell containing
Borexino (i.e. it is a good place to test things, including the gridding). 


n = 17114; 
s1 = UC;
%s1.radius = simple1.radius; 
s2 = simple2; 
cor = cor.UC; 
SurfRadius = Litho1.r(n);
detector = det; 
P = zeros(iter,1); 
layer = 'UC';


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



for "miniVox" function
radius = SurfRadius - depth'; 
lat = s1.lat; 
lon = s1.lon; 
lon_int = 1; 
lat_int = 1;
d_int = thick';




scatter(1:8,time_taken(:,1)/60,80,'filled','r')
xlabel('Number of cores used')
ylabel('Run time (min)')
title('Run time vs Pool size for 1000 iterations')
set(gca,'fontsize',20)



    [lon2,lat2,depth2,lon_int2,lat_int2,d_int2]...
        = miniVox(s1.lon,s1.lat,depth',1,1,thick',num,SurfRadius); 



%}

%% 8) ---- Calculate Abundances of each layer ----
method = 4; %determines what statistics to output

s1.total.hf = stat(s1_cc_sums(:,5),method); 
s1.total.hp.total = stat(s1_cc_sums(:,6),method); 
s1.total.hp.U = stat(s1_cc_sums(:,7),method); 
s1.total.hp.Th = stat(s1_cc_sums(:,8),method); 
s1.total.hp.K = stat(s1_cc_sums(:,9),method); 
s1.total.mass = stat(s1_cc_sums(:,1),method); 

s2.total.hf = stat(s2_cc_sums(:,5),method); 
s2.total.hp.total = stat(s2_cc_sums(:,6),method); 
s2.total.hp.U = stat(s2_cc_sums(:,7),method); 
s2.total.hp.Th = stat(s2_cc_sums(:,8),method); 
s2.total.hp.K = stat(s2_cc_sums(:,9),method); 
s2.total.mass = stat(s2_cc_sums(:,1),method); 

s3.total.hf = stat(s3_cc_sums(:,5),method); 
s3.total.hp.total = stat(s3_cc_sums(:,6),method); 
s3.total.hp.U = stat(s3_cc_sums(:,7),method); 
s3.total.hp.Th = stat(s3_cc_sums(:,8),method); 
s3.total.hp.K = stat(s3_cc_sums(:,9),method); 
s3.total.mass = stat(s3_cc_sums(:,1),method); 

UC.total.hf = stat(UC_cc_sums(:,5),method); 
UC.total.hp.total = stat(UC_cc_sums(:,6),method); 
UC.total.hp.U = stat(UC_cc_sums(:,7),method); 
UC.total.hp.Th = stat(UC_cc_sums(:,8),method); 
UC.total.hp.K = stat(UC_cc_sums(:,9),method); 
UC.total.mass = stat(UC_cc_sums(:,1),method); 

MC.total.aU = stat(MC_cc_sums(:,2)./(MC_cc_sums(:,1)),method)*10^6;  % aU (ppm)
MC.total.aTh = stat(MC_cc_sums(:,3)./(MC_cc_sums(:,1)),method)*10^6; % aTh (ppm)
MC.total.aK = stat(MC_cc_sums(:,4)./(MC_cc_sums(:,1)),method);  % aK (wt%) 
MC.total.hf = stat(MC_cc_sums(:,5),method); 
MC.total.hp.total = stat(MC_cc_sums(:,6),method); 
MC.total.hp.U = stat(MC_cc_sums(:,7),method); 
MC.total.hp.Th = stat(MC_cc_sums(:,8),method); 
MC.total.hp.K = stat(MC_cc_sums(:,9),method); 
MC.total.mass = stat(MC_cc_sums(:,1),method); 

LC.total.aU = stat(LC_cc_sums(:,2)./(LC_cc_sums(:,1)),method)*10^6;  % aU (ppm)
LC.total.aTh = stat(LC_cc_sums(:,3)./(LC_cc_sums(:,1)),method)*10^6; % aTh (ppm)
LC.total.aK = stat(LC_cc_sums(:,4)./(LC_cc_sums(:,1)),method);  % aK (wt%) 
LC.total.hf = stat(LC_cc_sums(:,5),method); 
LC.total.hp.total = stat(LC_cc_sums(:,6),method); 
LC.total.hp.U = stat(LC_cc_sums(:,7),method); 
LC.total.hp.Th = stat(LC_cc_sums(:,8),method); 
LC.total.hp.K = stat(LC_cc_sums(:,9),method); 
LC.total.mass = stat(LC_cc_sums(:,1),method); 

LM.total.hf = stat(LM_cc_sums(:,5),method); 
LM.total.hp.total = stat(LM_cc_sums(:,6),method); 
LM.total.hp.U = stat(LM_cc_sums(:,7),method); 
LM.total.hp.Th = stat(LM_cc_sums(:,8),method); 
LM.total.hp.K = stat(LM_cc_sums(:,9),method); 
LM.total.mass = stat(LM_cc_sums(:,1),method); 

%% 9) ---- Finish Geo-flux calculation to get TNU ---- 

iso = simple2.iso; 
eno = simple2.eno;

TNU.U238 = 1/(iso.amu.U238*iso.amu.amuKg) .*eno.U238' /sum(eno.U238)...
    *6 * iso.dc.U238 /10^4 /(7.67*10^4); % nu/cm2/s

TNU.Th232 = 1/(iso.amu.Th232*iso.amu.amuKg) .* eno.Th232' /sum(eno.Th232)...
    *4 * iso.dc.Th232 /10^4 /(2.48*10^5); % nu/cm2/s




 %clear test
 %{
TNU.U = 7.67*10^4; %conversion of nu/s/cm2 to TNU
TNU.Th = 2.48*10^5;

iso = simple2.iso; 
eno = simple2.eno;

y = UC_cc_flux_sums(:,1:length(simple2.energy)); %1 row sums to ~680 (sums to 18 after multiplication of eno)
x = UC_cc_flux_sums(:,length(simple2.energy)+1:end); 

t.fluxU = y * 1/(simple2.iso.amu.U238*1.660539*10^-27) .* eno.U238' /sum(eno.U238) * 6 * iso.dc.U238 /10^4; % nu/cm2/s
t.fluxTh = x * 1/(simple2.iso.amu.Th232*1.660539*10^-27) .* eno.Th232' /sum(eno.Th232)* 6 * iso.dc.Th232 /10^4; % nu/cm2/s


t.tnuU = t.fluxU./TNU.U; 
t.tnuTh = t.fluxTh./TNU.Th; 

z(1) = median(sum(t.tnuU,2));
z(2) = median(sum(t.tnuTh,2));
z;
clear z 
%{
% testing 
y = sort(y); 
sum(y(50,:).*eno.U238');
x = sort(t.fluxU); 
sum(x(50,:));


figure
histogram(sum(t.tnuU)); 
title(sprintf('Bin Size %d KeV',bin_size(test)))
%}
 %}
 



%% 10) ---- Reallocate individual cell data into respective structures -- 
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

%fprintf('Re-organizing Data...')

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
       end
    end
%{
% Flux of U238 and Th232 (TNU)
    n = {'U238';'Th232'};% these need to be the same as "UC_out"
    for j = 1:length(Litho1.latlon)
        for i = 1:length(n)    
            s1.flux(:).(n{i})(j,:) = s1_flux(j).(n{i}) .* TNU(:).(n{i});
            s2.flux(:).(n{i})(j,:) = s2_flux(j).(n{i}) .* TNU(:).(n{i});
            s3.flux(:).(n{i})(j,:) = s3_flux(j).(n{i}) .* TNU(:).(n{i});
            UC.flux(:).(n{i})(j,:) = UC_flux(j).(n{i}) .* TNU(:).(n{i});
            MC.flux(:).(n{i})(j,:) = MC_flux(j).(n{i}) .* TNU(:).(n{i});
            LC.flux(:).(n{i})(j,:) = LC_flux(j).(n{i}) .* TNU(:).(n{i});
            LM.flux(:).(n{i})(j,:) = LM_flux(j).(n{i}) .* TNU(:).(n{i});
       end
    end    
   %}       
    
    
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
    
    


%% 11) ---- Reproduce Huang et al. 2013 Table 2 (flux) ---
%{
The rows of the table are reservoirs.  The columns are average with
uncertainty for each reservoir. 

%}

method = 4; % 4 = median +- 68% c.l., 6 = geometric mean +- std

U238 = 1:length(simple2.energy); %1:16
Th232 = length(simple2.energy)+1:length(simple2.energy)*2; %1:32

% Row Names
huang.tab2.rows = {'Sed_cc';'UC';'MC';'LC';'LM';'Sed_oc';'OC';'Bulk CC';'Bulk Crust';'Bulk Litho'};


% Sediment (TNU)
crust = s1_cc_flux_sums + s2_cc_flux_sums + s3_cc_flux_sums; 
flux.sed.sums.cc.U238 = sum(crust(:,U238).*TNU.U238,2); 
flux.sed.sums.cc.Th232 = sum(crust(:,Th232).*TNU.Th232,2); 

ocean = s1_oc_flux_sums + s2_oc_flux_sums + s3_oc_flux_sums; 
flux.sed.sums.oc.U238 = sum(ocean(:,U238).*TNU.U238,2); 
flux.sed.sums.oc.Th232 = sum(ocean(:,Th232).*TNU.Th232,2); 

huang.tab2.U238(1,:) = stat(flux.sed.sums.cc.U238,method); 
huang.tab2.Th232(1,:) = stat(flux.sed.sums.cc.Th232,method); 
huang.tab2.total(1,:) = stat(flux.sed.sums.cc.U238+flux.sed.sums.cc.Th232,method); 


huang.tab2.U238(6,:) = stat(flux.sed.sums.oc.U238,method); 
huang.tab2.Th232(6,:) = stat(flux.sed.sums.oc.Th232,method); 
huang.tab2.total(6,:) = stat(flux.sed.sums.oc.U238+flux.sed.sums.oc.Th232,method); 


% UC (TNU)
flux.UC.sums.cc.U238 = sum(UC_cc_flux_sums(:,U238).*TNU.U238,2); 
flux.UC.sums.cc.Th232 = sum(UC_cc_flux_sums(:,Th232).*TNU.Th232,2); 

flux.UC.sums.oc.U238 = sum(UC_oc_flux_sums(:,U238).*TNU.U238,2); 
flux.UC.sums.oc.Th232 = sum(UC_oc_flux_sums(:,Th232).*TNU.Th232,2); 

huang.tab2.U238(2,:) = stat(flux.UC.sums.cc.U238,method); 
huang.tab2.Th232(2,:) = stat(flux.UC.sums.cc.Th232,method); 
huang.tab2.total(2,:) = stat(flux.UC.sums.cc.U238+flux.UC.sums.cc.Th232,method); 


% MC (TNU)
flux.MC.sums.cc.U238 = sum(MC_cc_flux_sums(:,U238).*TNU.U238,2); 
flux.MC.sums.cc.Th232 = sum(MC_cc_flux_sums(:,Th232).*TNU.Th232,2); 

flux.MC.sums.oc.U238 = sum(MC_oc_flux_sums(:,U238).*TNU.U238,2); 
flux.MC.sums.oc.Th232 = sum(MC_oc_flux_sums(:,Th232).*TNU.Th232,2); 

huang.tab2.U238(3,:) = stat(flux.MC.sums.cc.U238,method); 
huang.tab2.Th232(3,:) = stat(flux.MC.sums.cc.Th232,method); 
huang.tab2.total(3,:) = stat(flux.MC.sums.cc.U238+flux.MC.sums.cc.Th232,method); 


% LC (TNU)
flux.LC.sums.cc.U238 = sum(LC_cc_flux_sums(:,U238).*TNU.U238,2); 
flux.LC.sums.cc.Th232 = sum(LC_cc_flux_sums(:,Th232).*TNU.Th232,2); 

flux.LC.sums.oc.U238 = sum(LC_oc_flux_sums(:,U238).*TNU.U238,2); 
flux.LC.sums.oc.Th232 = sum(LC_oc_flux_sums(:,Th232).*TNU.Th232,2); 

huang.tab2.U238(4,:) = stat(flux.LC.sums.cc.U238,method); 
huang.tab2.Th232(4,:) = stat(flux.LC.sums.cc.Th232,method); 
huang.tab2.total(4,:) = stat(flux.LC.sums.cc.U238+flux.LC.sums.cc.Th232,method); 


% LM (TNU)
flux.LM.sums.cc.U238 = sum(LM_cc_flux_sums(:,U238).*TNU.U238,2); 
flux.LM.sums.cc.Th232 = sum(LM_cc_flux_sums(:,Th232).*TNU.Th232,2); 

flux.LM.sums.oc.U238 = sum(LM_oc_flux_sums(:,U238).*TNU.U238,2); 
flux.LM.sums.oc.Th232 = sum(LM_oc_flux_sums(:,Th232).*TNU.Th232,2); 

huang.tab2.U238(5,:) = stat(flux.LM.sums.cc.U238,method); 
huang.tab2.Th232(5,:) = stat(flux.LM.sums.cc.Th232,method); 
huang.tab2.total(5,:) = stat(flux.LM.sums.cc.U238+flux.LM.sums.cc.Th232,method); 


% OC Sed (TNU)
huang.tab2.U238(6,:) = stat(flux.sed.sums.oc.U238,method); 
huang.tab2.Th232(6,:) = stat(flux.sed.sums.oc.Th232,method); 
huang.tab2.total(6,:) = stat(flux.sed.sums.oc.U238+flux.sed.sums.oc.Th232,method); 


% OC (TNU)
U238_oc = flux.UC.sums.oc.U238 + flux.MC.sums.oc.U238 + flux.LC.sums.oc.U238; 
Th232_oc= flux.UC.sums.oc.Th232 + flux.MC.sums.oc.Th232 + flux.LC.sums.oc.Th232; 

huang.tab2.U238(7,:) = stat(U238_oc,method); 
huang.tab2.Th232(7,:) = stat(Th232_oc,method); 
huang.tab2.total(7,:) = stat(U238_oc + Th232_oc,method); 

% Bulk CC (Sed + UC + MC + LC) (TNU)
U238_cc = flux.UC.sums.cc.U238 + flux.MC.sums.cc.U238 + flux.LC.sums.cc.U238...
        + flux.sed.sums.cc.U238; 
Th232_cc = flux.UC.sums.cc.Th232 + flux.MC.sums.cc.Th232 + flux.LC.sums.cc.Th232...
        + flux.sed.sums.cc.Th232; 
    
huang.tab2.U238(8,:) = stat(U238_cc,method); 
huang.tab2.Th232(8,:) = stat(Th232_cc,method); 
huang.tab2.total(8,:) = stat(U238_cc + Th232_cc,method); 


% Bulk Crust (TNU)
U238_crust = U238_cc + flux.UC.sums.oc.U238 + flux.MC.sums.oc.U238...
    +flux.LC.sums.oc.U238 + flux.sed.sums.oc.U238;
Th232_crust = Th232_cc + flux.UC.sums.oc.Th232 + flux.MC.sums.oc.Th232...
    +flux.LC.sums.oc.Th232 + flux.sed.sums.oc.Th232;

huang.tab2.U238(9,:) = stat(U238_crust,method); 
huang.tab2.Th232(9,:) = stat(Th232_crust,method); 
huang.tab2.total(9,:) = stat(U238_crust + Th232_crust,method); 


% Total Lithosphere  (TNU)
U238_litho = U238_crust + flux.LM.sums.cc.U238; 
Th232_litho = Th232_crust + flux.LM.sums.cc.Th232; 

huang.tab2.U238(10,:) = stat(U238_litho,method); 
huang.tab2.Th232(10,:) = stat(Th232_litho,method); 
huang.tab2.total(10,:) = stat(U238_litho + Th232_litho,method); 




clear ocean U238_oc U238_crust U238_litho Th232_oc Th232_crust Th232_litho

huang.tab2.table = horzcat(huang.tab2.U238,huang.tab2.Th232,huang.tab2.total); 

%% 12) ---- Reproduce Huang et al. 2013 Table 3 (abund) ---
%{
The rows of the table are reservoirs.  The columns are average with
uncertainty for each reservoir. 

%}

method = 4; % see help "stat"
meth_norm = 5; 
% note: UC, sed, and LM will always fit using normal dist. maximum
% Liklihood fit because they are always normal (see Section 4)


weight = s1.mass.mass(:,1) + s2.mass.mass(:,1) + s3.mass.mass(:,1); 
w1 = s1.mass.mass(:,1)./weight; % weigh each layer by mass so we can combine
w2 = s2.mass.mass(:,1)./weight; 
w3 = s3.mass.mass(:,1)./weight; 
sed.thick = s1.thick + s2.thick + s3.thick; 
sed.rho = s1.rho.*w1 + s2.rho.*w2 + s3.rho.*w3; 
sed.mass.mass = s1.mass.mass + s2.mass.mass + s3.mass.mass; 
sed.mass.U = s1.mass.U + s2.mass.U + s3.mass.U; 
sed.mass.Th = s1.mass.Th + s2.mass.Th + s3.mass.Th;
sed.mass.K = s1.mass.K + s2.mass.K + s3.mass.K;
sed_sums = s1_cc_sums + s2_cc_sums + s3_cc_sums; 
%sed.aU = (


% Row Names
huang.tab3.rows = {'Sed';'UC';'MC';'LC';'LM'};



% - Density (rho; g/cm3)
x = [stat(sed.rho(find(sed.rho(cc))),meth_norm);...
                stat(UC.rho(cc),meth_norm); stat(MC.rho(cc),meth_norm); stat(LC.rho(cc),meth_norm);...
                stat(LM.rho(cc),meth_norm)]/1000; %(g/cm3) 
huang.tab3.rho = x; 

% - Thickness (km)
x = [stat(s1.thick(find(sed.thick(cc))),meth_norm); stat(UC.thick(cc),meth_norm); stat(MC.thick(cc),meth_norm);...
        stat(LC.thick(cc),meth_norm); stat(LM.thick(cc),meth_norm)]/1000; %(km)
huang.tab3.thick = x; 

% - Mass (10^21 kg)
x = [stat(sed_sums(:,1),meth_norm); stat(UC_cc_sums(:,1),meth_norm);stat(MC_cc_sums(:,1),meth_norm);...
        stat(LC_cc_sums(:,1),meth_norm); stat(LM_cc_sums(:,1),meth_norm)]/10^21; %(10^21 kg)
huang.tab3.mass.mass = x; 

% - Abundance of U (ppm; ug/g)
x = [stat(sed_sums(:,2)./sed_sums(:,1),meth_norm);stat(UC_cc_sums(:,2)./UC_cc_sums(:,1),meth_norm);...
        stat(MC_cc_sums(:,2)./MC_cc_sums(:,1),method );stat(LC_cc_sums(:,2)./LC_cc_sums(:,1),method );...
        stat(LM_cc_sums(:,2)./LM_cc_sums(:,1),meth_norm)]*10^6; %(ppm)
huang.tab3.abund.U = x;    

% - Abundance of Th (ppm; ug/g)
x = [stat(sed_sums(:,3)./sed_sums(:,1),meth_norm);stat(UC_cc_sums(:,3)./UC_cc_sums(:,1),meth_norm);...
        stat(MC_cc_sums(:,3)./MC_cc_sums(:,1),method );stat(LC_cc_sums(:,3)./LC_cc_sums(:,1),method );...
        stat(LM_cc_sums(:,3)./LM_cc_sums(:,1),meth_norm)]*10^6; %(ppm)
huang.tab3.abund.Th = x;   

% - Abundance of K (wt%)
x = [stat(sed_sums(:,4)./sed_sums(:,1),meth_norm);stat(UC_cc_sums(:,4)./UC_cc_sums(:,1),meth_norm);...
        stat(MC_cc_sums(:,4)./MC_cc_sums(:,1),method );stat(LC_cc_sums(:,4)./LC_cc_sums(:,1),method );...
        stat(LM_cc_sums(:,4)./LM_cc_sums(:,1),meth_norm)]/K.b*10^2; %(wt%)
huang.tab3.abund.K = x;  

% - Mass of U (10^15 kg)
x = [stat(sed_sums(:,2),meth_norm); stat(UC_cc_sums(:,2),meth_norm);stat(MC_cc_sums(:,2),method);...
        stat(LC_cc_sums(:,2),method); stat(LM_cc_sums(:,2),meth_norm)]/10^15;%(10^15 kg)
huang.tab3.mass.U = x; 

% - Mass of Th (10^15 kg)
x = [stat(sed_sums(:,3),meth_norm); stat(UC_cc_sums(:,3),meth_norm);stat(MC_cc_sums(:,3),method);...
        stat(LC_cc_sums(:,3),method); stat(LM_cc_sums(:,3),meth_norm)]/10^15; %(10^15 kg)
huang.tab3.mass.Th = x; 

% - Mass of K (10^19 kg)
x = [stat(sed_sums(:,4),meth_norm); stat(UC_cc_sums(:,4),meth_norm);stat(MC_cc_sums(:,4),method);...
        stat(LC_cc_sums(:,4),method); stat(LM_cc_sums(:,method),meth_norm)]/10^19/K.b; %(10^19 kg)
huang.tab3.mass.K = x; 

% - Heat production (TW; 10^12 W)
x = [stat(sed_sums(:,5),meth_norm); stat(UC_cc_sums(:,5),meth_norm);stat(MC_cc_sums(:,5),method);...
        stat(LC_cc_sums(:,5),method); stat(LM_cc_sums(:,5),meth_norm)]/10^12; %(10^12 W)
huang.tab3.hp = x;
    

% Horizontal concatenate and print 

huang.tab3.table = horzcat(huang.tab3.rho, huang.tab3.thick, huang.tab3.mass.mass, huang.tab3.abund.U,...
    huang.tab3.abund.Th, huang.tab3.abund.K, huang.tab3.mass.U, huang.tab3.mass.Th,...
    huang.tab3.mass.K, huang.tab3.hp); 


%% Save Testing Data
%{
load Results_TestingEnergyBins.mat 
results(test,1) = bin_size(test); 
results(test,2:4) = huang.tab2.total(end,:); 
save Results_TestingEnergyBins.mat results
fprintf('Energy Bin Size: %.1f  ||  Time Elapsed: %.1f min \n',bin_size(test),toc/60)
%}





%% 13) ---- Save data ---

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


str = sprintf('Results_%1.1eIter_%s_%s_%s_%s.mat',iter,MASTER.model,m,d,datestr(date,'ddmmmyyyy'));
save(str,'MASTER','huang','Litho1','s1','s2','s3','UC','MC','LC','LM','flux')

end

%% (function) "rand_n"

% This function will output "0" if iter = 1 so that we will draw the
% central value of a distribution (instead of randomly along the
% distribution shape). If iter > 1, then we can create random numbers along
% a distribution (normal).  These numbers are single point floating numbers
% because they use less memory. 
function [output] = rand_n(iterations)

    if iterations == 1 
        output = 0;% this will cause the creation of a distribution at the central value
    elseif iterations >1 
        output = randn(iterations,1); 
    end
end

%% (function) "stat"
function [output] = stat(PDF,whatKind)
% STAT finds the median, + 68% confidence limit, or - 68% confidence limit
%   of a distribution. 
%
%   SYNTAX: 
%   output(PDF) inds the median and 68% confidence limits of the distribution 
%               "PDF" 
%
%   output(PDF, whatKind) returns the median, + sigma, or - sigma 
%       (rounded to 5 decimal), of the distribution "PDF" correlating to
%       whatKind = 1-6 respectively (see below). Size of PDF 
%       is determined within the function.
%
%   output(PDF, whatKind, int) returns the median, + sigma, or - sigma 
%       (rounded to 5 decimal), of the distribution "PDF" correlating to
%       whatKind = 1, 2, or 3 respectively and size of distribution int. 
%
%
%   whatKind: 
%       1   |   Median +- 68.24% confidence limit (~1 sigma, but not)
%       2   |   Median +- 95.45% confidence limit (~2 sigma, but not)
%       3   |   Mean with standard deviation
%       4   |   Geometric mean with +- std of log-distribution 
%       5   |   normal dist. maximum Liklihood fit normal
%       6   |   log-normal dist. maximum Likelihood fit log-normal
%       7   |   gamma dist. maximum liklihood fit
%
%   -----            Written by Scott A. Wipperfurth             ----- 
%   -----      University of Maryland-College Park, Geology      ----- 
%   -----                        June, 2016                      -----
%   -----                last modified May, 2018             -----
%
%   See also median, mean, std, nanmedian, nanmean



% -- If only 1 function argument, i.e find median value --
if nargin == 1
    whatKind = 1;
end


% -- If 2 function argument, i.e find median value, but specify whatKind = [1-5] --

int = length(PDF); %size of PDF


% find Median, + sigma, - sigma 
if whatKind == 1 % 68.24% of data (~1 sigma)
   s_PDF = sort(PDF);
   PDF_median = nanmedian(PDF);

   output(1,1) = PDF_median;
   output(1,2) = s_PDF(round((int)/2) + round(0.3413*int)) - PDF_median; % plus sigma
   output(1,3) = PDF_median - s_PDF(round((int)/2) - round(0.3413*int)); % minus sigma

elseif whatKind == 2 % 95.45% of data (~2 sigma)
   s_PDF = sort(PDF);
   PDF_median = nanmedian(PDF);

   output(1,1) = PDF_median;
   output(1,2) = s_PDF(round((int)/2) + round(0.4773*int)) - PDF_median; % plus sigma
   output(1,3) = PDF_median - s_PDF(round((int)/2) - round(0.4773*int)); % minus sigma

elseif whatKind == 3 % mean +- std
   output(1,1) = nanmean(PDF); 
   output(1,2) = nanstd(PDF); 


elseif whatKind == 4 % mean +- from std of ln(pdf)
   output = zeros(1,3);
   s_PDF = log(PDF);
   x = std(s_PDF); 
   PDF_mean = exp(nanmean(s_PDF));

   output(1,1) = PDF_mean;
   output(1,2) = exp(nanmean(s_PDF) + x) - PDF_mean; % plus sigma
   output(1,3) = PDF_mean - exp(nanmean(s_PDF) - x); % minus sigma     

elseif whatKind == 5 % normal maximum liklihood fit
   if length(PDF) == 1
       output= exp(mean(log(PDF))); 
   else
   output = mle(PDF,'distribution','normal'); 
   output(:,3) = 0;
   end
   
elseif whatKind == 6 % log-normal maximum liklihood fit
       if length(PDF) == 1
       output= exp(mean(log(PDF))); 
   else
   x = mle(PDF,'distribution','lognormal'); 

   output(1,1) = exp(x(1));
   output(1,2) = exp(x(1) + x(2)) - exp(x(1)); % plus sigma
   output(1,3) = exp(x(1)) - exp(x(1) - x(2)); % minus sigma    
       end

elseif whatKind == 7 % gamma maximum liklihood fit
       if length(PDF) == 1
       output= exp(mean(log(PDF))); 
   else
  output = mle(PDF,'distribution','gamma');   
       end
else
     % Provide warning and stop program when wrong inputs
 error('Wrong input of "whatKind" or "int"') 
 end
end

%% (function) "logdist"
function [dist, log_mu, log_sigma] = logdist(mu,plusError,minusError,int,locate)
% LOGDIST Creates a log-normal distribution (non-gaussian) from input
% variables or selects a single value from a distribution. It is also 
% possible to select the location of the random value on the distribution 
% using the 'locate' parameter (where 0 = mean). 
%
% logdist(mu,plusError,minusError, int, locate)
%
%   mu          = mean 
%   plusError   = positive uncertainty
%   minusError  = negative uncertainty
%   int         = size of distribution ("0" if want single value)
%   locate      = location of value on PDF (can be larger than 1)
%   
%   If int = 1, dist will be the central value always. If int = 0, then 
%   output will be a single random value within the distribution, unless 
%   "locate" is larger than 1. In that case the number of random values
%   will be determined by the length of "locate".
%   
%   [dist, log_mu, log_sigma] = logdist(mu,plusError,minusError,int,locate)
%   provides the distribution (dist), the log(mean) of the distribution,
%   and log(sigma) of the distribution.  log(sigma) is a combination of
%   the positive and negative uncertainty and is calculated as: 
%
%   log_plus = log(mu + plusError) - log_mu;
%   log_minus = log_mu - log(mu - minusError);
%   log_sigma = (log_plus+log_minus)/2; 
%   
% See also randist, lognrnd, lognfit, lognpdf, randn.
%
%   -----            Written by Scott A. Wipperfurth             ----- 
%   -----      University of Maryland-College Park, Geology      ----- 
%   -----                     Created June, 2016                 ----- 
%   -----                     Updated June, 2018                 ----- 




if nargin == 4 

    if int == 1 % CENTRAL Value
              dist = mu;



    elseif int >1% MONTE CARLO 
            log_mu = log(mu);
            log_plus = log(mu + plusError) - log_mu;
            log_minus = log_mu - log(mu - minusError);
            log_sigma = (log_plus+log_minus)/2;

            dist = exp(randn(int,1).*log_sigma(:)+log_mu(:));




    elseif int ==0 %(i.e. you want 1 random value from the distribution)

            log_mu = log(mu);
            log_plus = log(mu + plusError) - log_mu;
            log_minus = log_mu - log(mu - minusError);
            log_sigma = (log_plus+log_minus)/2;

            dist = exp(randn(1,1).*log_sigma(:)+log_mu(:));


    end
    
elseif nargin == 5
    % Find random numbers when specified random distribution
    if int == 1 % CENTRAL Value
              dist = mu;



    elseif int >1% MONTE CARLO 
            log_mu = log(mu);
            log_plus = log(mu + plusError) - log_mu;
            log_minus = log_mu - log(mu - minusError);
            log_sigma = (log_plus+log_minus)/2;

            dist = exp(locate.*log_sigma(:)+log_mu(:)); 




    elseif int ==0 %(i.e. you want 1 random value from the distribution)

            log_mu = log(mu);
            log_plus = log(mu + plusError) - log_mu;
            log_minus = log_mu - log(mu - minusError);
            log_sigma = (log_plus+log_minus)/2;

            dist = exp(locate.*log_sigma(:)+log_mu(:));

    else
        warning('You entered %d input(s). You need 4 or 5.',nargin);
    end  
end
end

%% (function) "randist"

function [dist] = randist(mu,error,int,locate)
%RANDIST normal distribution.
%  RANDIST(mu,error,int,locate) creates a normal distribution (gaussian) 
%  from input variables.  mu = mean, error = +- uncertainty, int = size/type
%  of distribution output, and locate = location on PDF.  
%
%   int = 1, dist will be the central value always.
%   int = 0, returned value is single random value within distribution
%           (unless "locate" >1 in length)
%   int > 0, returns distribution comprised of "int" values
%
%   locate = "single value" will return the random number at "locate"
%       location on the distribution. This is used to provide correlated 
%       values betwen creation of different distributions. If locate is a
%       distribution it should be a uniform distribution (i.e. built using
%       the RANDN function, not RAND). 'locate' = 0 will return mu, while 
%       'locate' = 1 will return (5+3) = 8. 
%
%   
%
%   Example usage: 
%   [dist] = randist(3,2,0) returns a single value from the distribution
%       with mu = 3 and uncertainty +- = 2. 
%
%   [dist] = randist(3,2,0,0.9) returns a single value from the higher end
%       of the distribution with mu = 3 and uncertainty +- = 2. 
%
% See also randn, logdist. 
%
%   -----            Written by Scott A. Wipperfurth             ----- 
%   -----      University of Maryland-College Park, Geology      ----- 
%   -----                    Created June, 2016                  ----- 
%   -----                     Update June, 2018                  ----- 
%



if nargin == 3 
    if int == 1   % returns central value

                dist = mu;


    elseif int > 1% returns distribution of "int" length 

                dist = randn(int,1).*error + mu;


    elseif int == 0 %returns 1 random value from the distribution 

                dist = randn(1,1).*error + mu;


    end
    
elseif nargin == 4 % Find random numbers when specified random distribution
     if int == 1   % returns central value

                dist = mu;


    elseif int > 1% returns distribution of "int" length 

                dist = locate.*error + mu;


    elseif int == 0 %returns 1 random value from the distribution 

                dist = locate.*error + mu;

     end   
else
    warning('You entered %d input(s). You need 3 or 4.',nargin);
end
end