
;	Calculation method switches

 	dielectric_freeSpace = 0
 	dielectric_noPoloidal = 0
 	dispersion_freeSpace = 0
 	dispersion_jaegerPRL = 0
	dispersion_generalised = 1
 	
 	bandStorage = 0
	noIMSL = 1
	linearDevice = 0
	kjInput = 0
	kj_jP_fileName = 'kj_jP.nc'

;	Plotting switches

 	plotRunData	= 1
 	plotDispersionGeneral 	= 1
 	plotDispersionJaeger	= 0
 	plotDispersionNoPol		= 0
	plotESolution = 1
	plotHSolution = 0
	plotKdotB = 0
	plotEdotB = 0
	plotJp = 1
	plotFFTSolution = 0

;	File write switches

	writeDispersionTxt = 0

;		----------------------------

;		Variables

		r0	= 2.25d0
		aWall	= 0.7d0

		rMin	= r0 + aWall - 0.1
		rMax	= r0 + aWall
   		b0	= 3.8d0
		bR_frac	= 0.0
		bz_frac	= 0.0;15.0/90.0
		ionSpecZ	= [ 1 ]
		ionSpecAmu	= [ 2 ]
		nMax		= [ 1.0 ] * 9d16
		damping = 0.0
		freq	= 3.7e9
		nPhi = 400.0*3d0 ; kPar between 90 and 170
		kz = 100.0
		nR	= 512L
		antLoc	= 2.92
		AntennaJ_r = 1
		AntennaJ_t = 0
		AntennaJ_z = 0

		useEqdsk = 0
		useProfiles = 0
		poloidalScale = 1.0
		zSlice	= 0.0
		profile1 = 0
		profile2 = 1

;		-----------------------------


