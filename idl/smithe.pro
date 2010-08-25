
;	Calculation method switches

 	dielectric_freeSpace = 0
 	dielectric_noPoloidal = 0
 	dispersion_freeSpace = 0
 	dispersion_jaegerPRL = 1
	dispersion_generalised = 1
 	
 	bandStorage = 0

;	Plotting switches

 	plotRunData	= 1
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

		rMin	= 0.455;r0 - aWall*0.99
		rMax	= 0.885;r0 + aWall*0.99
   		b0	= 5.85d0
		bR_frac	= 0.15
		bz_frac	= 0.0	
		ionSpecZ	= [ 1, 2, 1 ]
		ionSpecAmu	= [ 2, 3, 1 ]
		nMax		= [ 0.88, 0.23, 0.66 ] * 1d20
		nLim		= [ 0.88, 0.23, 0.66 ] * 1d19
		damping = 0.00
		freq	= 80.5e6
		nPhi = 10.0
		kz = 00.0
		nR	= 1024L
		antLoc	= 0.455

		useEqdsk = 0
		useProfiles = 0
		poloidalScale = 1.0
		zSlice	= 0.0
		profile1 = 1

;		-----------------------------


