
;	Calculation method switches

 	dielectric_freeSpace = 0
 	dielectric_noPoloidal = 0
 	dispersion_freeSpace = 0
 	dispersion_jaegerPRL = 1
	dispersion_generalised = 1
 	
 	bandStorage = 1
	noIMSL = 1
	linearDevice = 1

;	Plotting switches

 	plotRunData	= 1
 	plotDispersionGeneral 	= 1
 	plotDispersionJaeger	= 1
 	plotDispersionNoPol		= 0
	plotSolution = 1

;	File write switches

	writeDispersionTxt = 0

;		----------------------------

;		Variables

		r0	= 100d0
		aWall	= 0.2 

		rMin	= r0 - aWall
		rMax	= r0 + aWall
   		b0	= 1.1d0
		bR_frac	= 0.0
		bz_frac	= 0.0	
		ionSpecZ	= [ 1 ]
		ionSpecAmu	= [ 1 ]
		nMax		= [ 1.0 ] * 1d8
		damping = 0.05
		freq	= 28.0e9
		nPhi = 0.0
		kz = 0.0
		nR	= 256L
		antLoc	= 0.0 

		useEqdsk = 0
		useProfiles = 0
		poloidalScale = 1.0
		zSlice	= 0.0
		profile1 = 1

;		-----------------------------


