
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
    ar2Input = 1
    ar2InFileName = 'ar2Input.nc'
    ar2RunDataFileName = 'runData_001_+013_000001.nc'

;	Plotting switches

 	plotRunData	= 0
 	plotDispersionGeneral 	= 0
 	plotDispersionJaeger	= 0
 	plotDispersionNoPol		= 0
	plotESolution = 0
	plotHSolution = 0
	plotKdotB = 0
	plotEdotB = 0
	plotJp = 0
	plotFFTSolution = 0
	plotJdotE = 0

;	File write switches

	writeDispersionTxt = 0

;		----------------------------

;		Variables

		;r0	= 100d0
		;aWall	= 5d0

		rMin	= 0.92+100
		rMax	= 1.72+100
   		b0	= 0.0
		bR_frac	= 0.0
		bz_frac	= 0.0	
		ionSpecZ	= [ 2 ]
		ionSpecAmu	= [ 4 ]
		nMax		= [ 1.0 ] * 0
		damping = 0.00
		freq	= 42.0e6
		nPhi = 13.0
		kz = 0.0
		nR	= 2048L
		antLoc	= 1.6
		AntennaJ_r = 0
		AntennaJ_t = 0
		AntennaJ_z = 1

		antSig_r = 0.02
		antSig_t = 0.02
		antSig_z = 0.02

		jAmp = -1.0 ; This seems to be a sign difference to AORSA jA

		useEqdsk = 0
		useProfiles = 0
		poloidalScale = 1.0
		zSlice	= 0.0
		profile1 = 0

;		-----------------------------


