
;	Calculation method switches

 	dielectric_freeSpace = 0
 	dielectric_noPoloidal = 0
 	dispersion_freeSpace = 0
 	dispersion_jaegerPRL = 0
	dispersion_generalised = 1
 	
 	bandStorage = 1

;	Plotting switches

 	plotRunData	= 0
 	plotDispersionGeneral 	= 0
 	plotDispersionJaeger	= 0
 	plotDispersionNoPol		= 0
	plotSolution = 1

;	File write switches

	writeDispersionTxt = 0

;		----------------------------

;		Variables

		r0	= 0.67d0
		aWall	= 0.22 

		rMin	= r0 - aWall*0.99
		rMax	= r0 + aWall*0.99
   		b0	= 5.85d0
		bR_frac	= 0.15
		bz_frac	= 0.0	
		ionSpecZ	= [ 1, 2, 1 ]
		ionSpecAmu	= [ 2, 3, 1 ]
		nMax		= [ 0.88, 0.23, 0.66 ] * 1d20
		damping = 0.00
		freq	= 80.5e6
		nPhi = 10.0
		kz = 00.0
		nR	= 4096L
		antLoc	= 0.0 

		useEqdsk = 0
		useProfiles = 0
		poloidalScale = 1.0
		zSlice	= 0.0
		profile1 = 1

;		-----------------------------


