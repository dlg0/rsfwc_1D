	common switches, $
		dielectric_freeSpace, $
		dielectric_noPoloidal, $
		dispersion_generalised, $
		dispersion_freeSpace, $
		dispersion_jaegerPRL, $
        dispersion_noPoloidal, $
		bandStorage, $
		noIMSL, $
		linearDevice, $
		kjInput, $
        kjDeltaFileName, $
        ar2Input, $
        ar2InFileName, $
        ar2RunDataFileName, $
        ar2EField

	common plotSwitches, $
		plotDispersionGeneral, $
		plotDispersionJaeger, $
		plotDispersionNoPol, $
		plotESolution, $
		plotHSolution, $
		plotMovie, $
		plotFrequencies, $
		plotKdotE, $
		plotKdotB, $
		plotEdotB, $
		plotJp, $
		plotFFTSolution, $
		plotJdotE, $
        plotRHS

	common writeSwitches, $
		writeDispersionTxt

    ; Set default value for switches

		dielectric_freeSpace = 0 
		dielectric_noPoloidal = 0
		dispersion_generalised = 0
		dispersion_freeSpace = 0
		dispersion_jaegerPRL = 0
        dispersion_noPoloidal = 0
		bandStorage	= 0
		noIMSL = 1
		linearDevice = 0

		plotDispersionGeneral = 0
		plotDispersionJaeger = 0
		plotDispersionNoPol = 0
		plotSolution = 0
        plotHSolution = 0
		plotMovie = 0
		plotFrequencies = 0
		plotKdotE = 0
		plotKdotB = 0
		plotEdotB = 0
		plotJp = 0
        plotJDotE = 0
		plotFFTSolution = 0
        plotRHS = 0

		writeDispersionTxt = 0

		kjInput = 0
		kj_jP_fileName = 'kj_jP.nc'
