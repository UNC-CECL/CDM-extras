import os
import numpy as np
import shutil

#### model parameterization ####
## directory specification ##
parent_dir = os.path.dirname(os.path.realpath(__file__)) + '/'

overwrite=True; ## If True, this script overwrites any existing model_iter* files/folders. If False, this script continues at the next model_iter number.
if overwrite:
    itermax=0;
else:
    D=[d for d in os.listdir(parent_dir) if os.path.isdir(os.path.join(parent_dir, d))]
    itermax = len(D)

## Varied model parameters
wind_fraction_list=np.linspace(0.2, 1.0, 5)
rep=np.arange(1);

# combinations of vectors
x1, x2 = np.meshgrid(wind_fraction_list, rep); ### MODIFY THIS LINE WITH PARAMETERS BEING MANIPULATED
x1 = x1.ravel()
x2 = x2.ravel()
pairs = np.column_stack((x1, x2))
n_iter = x1.shape[0];
VarNames = ['wind_fraction','replicate']; ### MODIFY THIS LINE WITH NAME OF PARAMETERS BEING MANIPULATED

#n_iter = 1;
CDM_filename = 'CDM_v2_0'; ### MODIFY THIS LINE WITH THE NAME OF THE CDM PROGRAM FILE

### BE SURE TO ADD X1[ii] AND X1[ii] AS THE VALUE ASSOCIATED WITH EACH MANIPULATED PARAMETER BELOW.
### E.G., wind_fraction = x1[ii]

for ii in range(n_iter):
    itr=ii+itermax;
    # write metadata
    if(itr==0):
        f= open(parent_dir + "iter_metadata.txt","w")
        f.write('model_iter\t%s\t%s\n' % (VarNames[0], VarNames[1]));
        f.close()
    
    f= open(parent_dir + "iter_metadata.txt","a")
    f.write('%i\t%f\t%f\n' % (itr+1, pairs[ii,0], pairs[ii,1]))
    f.close()
    
    sub_dir='model_iter' + str(itr+1) + '/'
    data_sub_dir = parent_dir + sub_dir + 'DATA'
    
    iter_dir = parent_dir + sub_dir
    if not os.path.exists(iter_dir):
        os.mkdir(iter_dir)             # make sub-directory for iteration
    
    if not os.path.exists(data_sub_dir):
        os.mkdir(data_sub_dir);    # make data sub-sub-directory for iteration
    
    shutil.copy2(parent_dir + CDM_filename, parent_dir+sub_dir+'DUNE') # copy dune program
    
    f = open(parent_dir + sub_dir + 'DUNE_program_version_info.txt', 'w' )
    f.write('DUNE program: %s \n' % CDM_filename)
    f.close()

    ## specify grid size
    NX = 200; 					# number of grid columns
    NY = 8; #128 #64 	 		# number of grid lines (must be power of 2; also 4 is for 2D, higher for 3D) 
    dx = 1;						# (m)	such that actual lenght is NX*dx X NY*dx

    Nt = 100000; 				# number of iterations (such that total time is Nt*dt_max)

    ## save info
    save_every  = 1000;         # number of iterations between saved files (such that the number of saved files are Nt/save.every)
    save_dir    = data_sub_dir;     # name of directory where output data is stored

    ## wind model parameters
    constwind_u = 0.35;         # shear velocity (m/s) (usually from 0.2 [~transport threshold] to 0.5 [strong winds])
    wind_fraction = x1[ii];        # fraction of the year wind is above threshold (only used for time scales) such that real time is Nt*dt_max / wind.fraction
  
    ## storm model:
    calc_storm = 'true';        # initiates storms sub-routine
    storm_gamma_shapeparameter = 4.0; # gamma distribution shape parameter for TWL (above MHW)
    storm_gamma_scaleparameter = 0.18; # gamma distribution scale parameter for TWL (above MHW)
    storm_totaliter = 50; # scalar for degree of erosiveness of storms. Defaults to 50

    
    ## Vegetation model parameters
    veget_calc 	= 1;            # 0: no veget cal. 1: calc. veg.
    veget_type 	= 2;            # vegetation type: 	0 - generic vegetation in mature ecosystem (no feedback with accretion rate and no lateral propagation)
                                # 					1 - generic vegetation on new dunes (feedback with accretion rate and lateral propagation)
                                #                   2 - Biel logistic growth with lateral expansion

    ## general parameters for all veg types
    veget_xmin 	= 0;                       # vegetation limit: L_veg (m) 
    veget_zmin 	= 0.2;                      # threshold elevation for veg. growth relative to MSL: Z_veg (m) 
    veget_erosion_sensitivity  = 1; 		# plant sensitivity to erosion (m^-1) ~ inverse of the ratio of root system volume to cover area [0 - 4]

    ## vegetation type 0 (modified from Duran & Moore, PNAS 2013)
    veget_Tveg 	= 10; 						# (days) characteristic time of vegetation cover growth [3-30]
    veget_0_init = 1e-2; 					# cover fraction for initial colonization: lower values -> more difficult/lenghty colonization 

    ## vegetation type 1 (Moore, Duran & Ruggiero, Geology 2016)
    veget_Hveg 	= 0.25;                     # charateristic vegetation height (m), controls vertical growth rate (lower values faster growth)
    veget_Vlateral_factor = 100;			# lateral growth factor (dimensionless)
    veget_max_slope = 15;					# (degrees) slope limit for lateral propagation (lateral propagation only for lower angles) 
    veget_1_init = 1e-5; 					# cover fraction for initial colonization: lower values -> more difficult/lenghty colonization 

    
    ## beach parameters
    calc_shore = 'true';
    # shore dynamics
    shore_sealevelrise = 0.0;               # rate of sea level rise (m/yr) in the range 0 - 0.1
    shore_alongshore_grad = 0.0;            # shoreline erosion(+)/accretion(-) rate (m/yr) due to alongshore transport gradients in real time

    # shore geometry (initial condition, from Duran & Moore 2013)
    beach_angle = 1.0;                        # shoreface angle (degrees)
    shore_MHWL = 0.5;                       # MHWL relative to the watertable
    shore_sealevel = 0;                     # watertable elevation

    ##########################################
    # INITIAL CONDITIONS						
    ##########################################
    # initial sand surface
    Init_Surf = 'beach'; # either plain, beach, or init_h

    # ---- flat surface, initialise with constant: Init-Surf = plain
    plain_Height = 0.0;

    # ---- beach: Init-Surf = beach ----
    beach_h = 0.3; 		# MHWL relative to watertable
    
    # ----  initialisation from file: Init-Surf = init_h ----
    init_h_file = 'h.100000.dat'; 	# from h.#####.dat
    if (Init_Surf == 'init_h'):
        shutil.copy2(parent_dir + init_h_file, parent_dir + sub_dir); # copy init_h_file program
    
    ##########################################
    # initial vegetation
    veget_Init_Surf = 'plain'; # either: alea, plain or init_h (look below for details)

    # ---- flat surface, initialise with constant: Init-Surf = plain
    veget_plain_Height = 0.0;

    # ---- random seeding: Init-Surf = alea
    veget_alea_nodes_b = 10;

    # ----  initialisation from file: Init-Surf = init_h ----
    veget_init_h_file = 'veget_x.100000.dat'; 	# from veget_x.#####.dat
    veget_init_h_file_aux = 'veget_y.100000.dat'; # from veget_y.#####.dat
    veget_init_h_x_line = 0; 				# keep it the same as save.x-line
    if (veget_Init_Surf == 'init_h'):
        shutil.copy2(parent_dir + veget_init_h_file, parent_dir + sub_dir); # copy dune program
        shutil.copy2(parent_dir + veget_init_h_file_aux, parent_dir + sub_dir); # copy dune program
    

    ##########################################
    ## SAVING FIELDS
    ##########################################
    ## suppress saving of some variables:
    ## dontsave.<truncated file name> = 1
    ## for example:
    dontsave_veget = 0;
    dontsave_u = 1;
    dontsave_flux= 1;
    dontsave_flux_s= 1;
    dontsave_shear= 1;
    dontsave_shear_pert= 1;
    dontsave_stall= 1;
    dontsave_rho = 1;
    dontsave_h_deposit= 1;
    dontsave_h_nonerod= 1;
    dontsave_h_sep= 1;
    dontsave_dhdt= 1;

    ##########################################
    save_x_line = 0;     # changing to 1 will reverse the reading of rows vs. columns; for Gnu plot set to 0, to reverse set to 1
    dt_max = 1000;	# time step (constant, max is not important) (sec)

    Nt0 = 0; 					# to continue a previous simulation Nt0 = Nt of prev.

    #################################
    ## influx:  const or outflux
    influx = 'const';
    q_in = 0.0; #fraction of maximum flux (from 0 to 1)



    ## print parameter file ##
    f = open(parent_dir + sub_dir + 'params.par', 'w' )
    f.write( '# \n');
    f.write( '#  $Id: default.par,v 1.11 2004/09/22 15:13:47 schatz Exp $ \n');
    f.write( '# \n');
    f.write( '#  parameter file ( iteration %i ) \n' % (itr+1));
    f.write( '# \n');
    f.write( '################################ \n');
    f.write( '\n');
    f.write( '# grid \n');
    f.write( '# \n');
    f.write( 'NX = %i                  # number of grid columns \n' % NX);
    f.write( 'NY = %i                  # number of grid lines (must be power of 2; also 4 is for 2D, higher for 3D) \n' % NY);
    f.write( 'dx = %i                  # (m)	such that actual lenght is NX*dx X NY*dx \n' % dx);
    f.write( 'Nt = %i                  # number of iterations (such that total time is Nt*dt_max) \n' % Nt);
    f.write( '\n');
    f.write( 'save.every  = %i         # number of iterations between saved files (such that the number of saved files are Nt/save.every) \n' % save_every);
    f.write( 'save.dir    = %s         # name of directory where output data is stored \n' % save_dir);
    f.write( '\n');
    f.write( '\n');
    f.write( '############################################################################### \n');
    f.write( '# Wind model: const \n');
    f.write( '\n');
    f.write( 'constwind.u = %f         # shear velocity (m/s) (usually from 0.2 [~transport threshold] to 0.5 [strong winds]) \n' % constwind_u);
    f.write( '\n');
    f.write( 'wind.fraction = %f       # fraction of the year wind is above threshold (only used for time scales) such that real time is Nt*dt_max / wind.fraction \n' % wind_fraction);
    f.write( '\n');
    f.write( '\n');
    f.write( '############################################################################### \n');
    f.write( '# storm model: const \n');
    f.write( 'calc.storm = %s # true/false to indicate whether to run storms routine \n' % calc_storm);
    f.write( 'storm.gamma.shapeparameter = %f # gamma distribution shape parameter for TWL (above MHW) \n' % storm_gamma_shapeparameter);
    f.write( 'storm.gamma.scaleparameter = %f # gamma distribution scale parameter for TWL (above MHW) \n' % storm_gamma_scaleparameter);
    f.write( 'storm.totaliter = %i    #  scalar for degree of erosiveness of storms. Defaults to 50 \n' % storm_totaliter);
    f.write( '\n');
    f.write( '\n');
    f.write( '############################################################################### \n');
    f.write( '# Vegetation model: \n');
    f.write( '# \n');
    f.write( 'veget.calc 	= %i  	# 0: no veget cal. 1: calc. veg. \n' % veget_calc);
    f.write( '\n');
    f.write( 'veget.type 	= %i		# vegetation type: 	0 - generic vegetation in mature ecosystem (no feedback with accretion rate and no lateral propagation) \n' % veget_type);
    f.write( '					# 					1 - generic vegetation on new dunes (feedback with accretion rate and lateral propagation) \n');
    f.write( '					# 					2 - Biel logistic growth with lateral expansion \n');
    f.write( '\n');
    f.write( '## general for all veg types \n');
    f.write( 'veget.xmin 	= %f  					# vegetation limit: L_veg (m)  \n' % veget_xmin);
    f.write( 'veget.zmin 	= %f	  				# threshold elevation for veg. growth relative to MSL: Z_veg (m)  \n' % veget_zmin);
    f.write( 'veget.erosion.sensitivity  = %f 		# plant sensitivity to erosion (m^-1) ~ inverse of the ratio of root system volume to cover area [0 - 4] \n' % veget_erosion_sensitivity);
    f.write( '\n');
    f.write( '### vegetation type 0 (modified from Duran & Moore, PNAS 2013) \n');
    f.write( 'veget.Tveg 	= %i 						# (days) characteristic time of vegetation cover growth [3-30] \n' % veget_Tveg);
    f.write( 'veget.0.init = %f 					# cover fraction for initial colonization: lower values -> more difficult/lenghty colonization  \n' % veget_0_init);
    f.write( '\n');
    f.write( '### vegetation type 1 (Moore, Duran & Ruggiero, Geology 2016) \n');
    f.write( 'veget.Hveg 	= %f 		 				# charateristic vegetation height (m), controls vertical growth rate (lower values faster growth) \n' % veget_Hveg);
    f.write( 'veget.Vlateral.factor = %f			# lateral growth factor (dimensionless) \n' % veget_Vlateral_factor);
    f.write( 'veget.max.slope = %f					# (degrees) slope limit for lateral propagation (lateral propagation only for lower angles)  \n' % veget_max_slope);
    f.write( 'veget.1.init = %f 					# cover fraction for initial colonization: lower values -> more difficult/lenghty colonization  \n' % veget_1_init);
    f.write( '\n');
    f.write( '\n');
    f.write( '###############################################################################\n');
    f.write( '### beach param:\n');
    f.write( 'calc.shore = %s                          # calculate shoreline \n' % calc_shore);
    f.write( '\n');
    f.write( '## shore dynamics\n');
    f.write( 'shore.sealevelrise = %f  		# rate of sea level rise (m/yr) in the range 0 - 0.1\n' % shore_sealevelrise);
    f.write( 'shore.alongshore_grad = %f 		# shoreline erosion(+)/accretion(-) rate (m/yr) due to alongshore transport gradients in real time\n' % shore_alongshore_grad);
    f.write( '\n');
    f.write( '## shore geometry (initial condition, from Duran & Moore 2013)\n');
    f.write( 'beach.angle = %f 		# shoreface angle (degrees)\n' % beach_angle);
    f.write( 'shore.MHWL = %f 		# MHWL relative to the watertable\n' % shore_MHWL);
    f.write( 'shore.sealevel = %f 		# watertable elevation\n' % shore_sealevel);
    f.write( '\n');
    f.write( '#########################################\n');
    f.write( '# INITIAL CONDITIONS						\n');
    f.write( '#########################################\n');
    f.write( '# initial sand surface\n');
    f.write( '#\n');
    f.write( 'Init-Surf = %s # either plain, beach or init_h \n' % Init_Surf);
    f.write( '\n');
    f.write( '# ---- flat surface, initialise with constant: Init-Surf = plain\n');
    f.write( 'plain.Height = %f\n' % plain_Height);
    f.write( '\n');
    f.write( '# ---- beach: Init-Surf = beach ----\n');
    f.write( 'beach.h = %f 		# MHWL relative to watertable\n' % beach_h);
    f.write( '\n');
    f.write( '# initialisation from file: Init-Surf = init_h ----\n');
    f.write( 'init_h.file = %s 	# from h.#####.dat' % init_h_file);
    f.write( '\n');
    f.write( '\n');
    f.write( '#########################################\n');
    f.write( '# initial vegetation\n');
    f.write( '#\n');
    f.write( 'veget.Init-Surf = %s # either: alea, plain or init_h (look below for details)\n' % veget_Init_Surf);
    f.write( '\n');
    f.write( '# ---- flat surface, initialise with constant: Init-Surf = plain\n');
    f.write( 'veget.plain.Height = %f\n' % veget_plain_Height);
    f.write( '\n');
    f.write( '# ---- random seeding: Init-Surf = alea\n');
    f.write( 'veget.alea.nodes.b = %f\n' % veget_alea_nodes_b);
    f.write( '\n');
    f.write( '# ----  initialisation from file: Init-Surf = init_h ---- \n');
    f.write( 'veget.init_h.file = %s # from veget_x.#####.dat \n' % veget_init_h_file);
    f.write( 'veget.init_h.file_aux = %s # from veget_y.#####.dat \n' % veget_init_h_file_aux);
    f.write( 'veget.init_h.x-line = %i # keep it the same as save.x-line \n' % veget_init_h_x_line);
    f.write( '\n');
    f.write( '\n');
    f.write( '#########################################\n');
    f.write( '# SAVING FIELDS\n');
    f.write( '#########################################\n');
    f.write( '# suppress saving of some variables:\n');
    f.write( '# dontsave.<truncated file name> = 1\n');
    f.write( '# for example:\n');
    f.write( 'dontsave.veget = %i\n' % dontsave_veget);
    f.write( 'dontsave.u = %i\n' % dontsave_u);
    f.write( 'dontsave.flux= %i\n' % dontsave_flux);
    f.write( 'dontsave.flux_s= %i\n' % dontsave_flux_s);
    f.write( 'dontsave.shear= %i\n' % dontsave_shear);
    f.write( 'dontsave.shear_pert= %i\n' % dontsave_shear_pert);
    f.write( 'dontsave.stall= %i\n' % dontsave_stall);
    f.write( 'dontsave.rho = %i\n' % dontsave_rho);
    f.write( 'dontsave.h_deposit= %i\n' % dontsave_h_deposit);
    f.write( 'dontsave.h_nonerod= %i\n' % dontsave_h_nonerod);
    f.write( 'dontsave.h_sep= %i\n' % dontsave_h_sep);
    f.write( 'dontsave.dhdt= %i\n' % dontsave_dhdt);
    f.write( '\n');
    f.write( '#########################################\n');
    f.write( 'save.x-line = %i     # changing to 1 will reverse the reading of rows vs. columns; for Gnu plot set to 0, to reverse set to 1\n' % save_x_line);
    f.write( 'dt_max = %i	# time step (constant, max is not important) (sec)\n' % dt_max);
    f.write( '\n');
    f.write( 'Nt0 = %i 					# to continue a previous simulation Nt0 = Nt of prev.\n' % Nt0);
    f.write( '\n');
    f.write( '################################\n');
    f.write( '# influx:  const or outflux\n');
    f.write( 'influx = %s\n' % influx);
    f.write( 'q_in = %f #fraction of maximum flux (from 0 to 1)' % q_in);
    f.write( '\n');
    f.close()

