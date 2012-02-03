
;	Calculation method switches

 	dielectric_freeSpace = 0
 	dielectric_noPoloidal = 0
 	dispersion_freeSpace = 0
 	dispersion_jaegerPRL = 0
	dispersion_generalised = 1
 	
 	bandStorage = 0
	noIMSL = 1
	linearDevice = 0
	kjInput = 1
	kj_jP_fileName = 'kj_jP.nc'

;	Plotting switches

 	plotRunData	= 0
 	plotDispersionGeneral 	= 0
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

		r0	= 100d0
		aWall	= 5d0

		rMin	= r0 - aWall
		rMax	= r0 + aWall
   		b0	= 0.0001d0
		bR_frac	= 100.0
		bz_frac	= 0.0	
		ionSpecZ	= [ 1 ]
		ionSpecAmu	= [ 1 ]
		nMax		= [ 1.0 ] * 1.1d14
		damping = 0.00
		freq	= 1.15e8
		nPhi = 0.0
		kz = 0.0
		nR	= 4096L
		antLoc	= 100.0

		useEqdsk = 0
		useProfiles = 0
		poloidalScale = 1.0
		zSlice	= 0.0
		profile1 = 0

;		-----------------------------


