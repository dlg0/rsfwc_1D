
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

 	plotRunData	= 0
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

		r0	= 6.6d0
		aWall	= 2.5d0

		rMin	= 4.0;r0 - aWall
		rMax	= 8.2;r0 + aWall
   		b0	= 5.3d0
		bR_frac	= 0.0
		bz_frac	= 0.0	
		ionSpecZ	= [ 1 ]
		ionSpecAmu	= [ 3 ]
		nMax		= [ 3.2 ] * 1d19
		damping = 0.02
		freq	= 53e6
		nPhi = -27.0
		kz = 0.0
		nR	= 4096L
		antLoc	= 8.0

		useEqdsk = 1
		eqdskFName = 'Scen4_bn2.57_129x129'
		useProfiles = 0
		poloidalScale = 1.0
		zSlice	= 0.0
		profile1 = 0

;		-----------------------------


