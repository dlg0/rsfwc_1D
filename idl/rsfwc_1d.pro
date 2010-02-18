;	Real Space Full Wave Code in 1D
;	DLG

pro rsfwc_1d, $
    eR = eR, $
	ePhi = ePhi, $
	ez = ez, $
	divD = divD, $
	rFull = rFull, rHalf = rHalf, $
	jR = jR, jPhi = jPhi, jz = jz, $
	in_kz = in_kz, $
	nMax = nMax, $
	nFac = nFac, $
	iiAnt = iiAnt, $
	gradSize = gradSize, $
	dispersionOnly = dispersionOnly, $
	kR = kR, $
	kPhi = kPhi, $
	out_kz = out_kz, $
	freq = freq, $
	rOut = rOut, $
	specData = specData, $
	stixVars = stixVars, $
	runData = runData, $
	poloidalScale = poloidalScale, $
	nPhi = nPhi

;	Parameters

    @constants
    @switches
    @load_colors

	run_setup, $
		runData = runData, $
		specData = specData, $
		in_kz = in_kz, $
		nMax = nMax, $
		nFac = nFac, $
		gradSize = gradSize, $
		freq = freq, $
		nPhi = nPhi, $
		poloidalScale = poloidalScale

	rFull 	= runData.r
	rHalf	= runData.r_
	rOut	= rFull

 	wReal	= runData.freq * 2d0 * !dpi 
	w	= dcomplexArr ( runData.nR ) + wReal
	;iiDamped	= where ( r gt 1.4 )
	w[*]	= complex ( wReal, wReal * runData.damping )

	stixVariables, wReal, specData, stixVars = stixVars

	dielectric, runData, stixVars, $
		epsilonFull = epsilon, $
		epsilonHalf = epsilon_

	dispersion, wReal, epsilon, stixVars, runData, specData, $
		kR = kR, kPhi = kPhi, kz = kz

	;out_kz = kz

	if keyword_set ( dispersionOnly ) then return

	matFill, runData.nR, runData.nPhi, runData.kz, $
		runData.r, runData.r_, epsilon, epsilon_, w, runData.dR, $
		aMat = aMat, nAll = nAll, nuca = nuca, nlca = nlca

	if runData.antLoc gt runData.rMax or runData.antLoc lt runData.rMin then begin
		print, 'ADJUSTING: Antenna location was ', runData.antLoc
		if runData.antLoc lt min ( runData.r_ ) then runData.antLoc = min ( runData.r_ )
		if runData.antLoc gt max ( runData.r_ ) then runData.antLoc = max ( runData.r_ )
		print, 'ADJUSTING: Antenna location is now ', runData.antLoc
	endif

    iiAnt   = where ( abs ( runData.r_ - runData.antLoc ) $
						eq min ( abs ( runData.r_ - runData.antLoc ) ) )
	rhs		= complexArr ( nAll )

	if not keyword_set ( jR ) then jR = 0
	if not keyword_set ( jPhi ) then jPhi = 0
	if not keyword_set ( jz ) then jz = 1

	rhs[iiAnt*3+2]	= -II * wReal * u0 * runData.dr * jz
	rhs[iiAnt*3+1]	= -II * wReal * u0 * runData.dr * jPhi 
	rhs[iiAnt*3]	= -II * wReal * u0 * runData.dr * jR 

;	Solve matrix

    print, '*** solving linear system'


    if bandStorage then begin

        eField = imsl_sp_bdSol ( rhs, nlca, nuca, aMat[*], /double )

    endif else begin

	    eField	= la_linear_equation ( aMat, rhs, status = stat )
	    print, 'lapack status: ', stat

    endelse
   
	ii_eR	= lIndGen(runData.nR)*3
	eR	= eField[ii_eR]
	ii_ePhi	= temporary(ii_eR[0:runData.nR-2]+1)
	ePhi	= eField[ii_ePhi]
	ii_ez	= temporary(ii_ePhi+1)
	ez	= eField[ii_ez]

	if keyword_set ( divD ) then begin

	;  	Calculate the Div of D @ the z,phi grid pts
	;  	but first we need to invert epsilon to get D & D_

    	divD_   = dcomplexArr ( runData.nR - 1 )

    	print, '*** calculating div D'
    	for i=1UL,runData.nR-3L do begin 
    	    
    	    Dri     = epsilon[0,0,i] * eR[i] $
    	            + ( epsilon_[1,0,i-1] * ePhi[i-1] + epsilon_[1,0,i] * ePhi[i] )/2 $
    	            + ( epsilon_[2,0,i-1] * ez[i-1] + epsilon_[2,0,i] * ez[i] )/2
    	    Dri1    =  epsilon[0,0,i+1] * eR[i+1] $
    	            + ( epsilon_[1,0,i] * ePhi[i] + epsilon_[1,0,i+1] * ePhi[i+1] )/2 $
    	            + ( epsilon_[2,0,i] * ez[i] + epsilon_[2,0,i+1] * ez[i+1] )/2
    	    Dth     = ( epsilon[0,1,i] * eR[i] + epsilon[0,1,i+1] * eR[i+1] )/2 $
    	            + epsilon[1,1,i] * ePhi[i] $
    	            + epsilon[2,1,i] * ez[i] 
    	    Dz      = ( epsilon[0,2,i] * eR[i] + epsilon[0,2,i+1] * eR[i+1] )/2 $
    	            + epsilon[1,2,i] * ePhi[i] $
    	            + epsilon[2,2,i] * ez[i] 

    	    divD_[i]    = ( r[i+1] * e0 * Dri1 - r[i] * e0 * Dri ) / ( r_[i] * dr ) $
    	                    + II * nPhi / r_[i] * e0 * Dth $
    	                    + II * kz * e0 * Dz
    	endfor

    	divE_   = dcomplexArr ( runData.nR - 1 )
    	for i=0UL,runData.nR-2L do begin
    	    divE_[i]    = II * nPhi * ePhi[i] / r_[i] $
    	                    + II * kz * ez[i] $
    	                    + ( r[i+1] * eR[i+1] - r[i] * eR[i] ) / ( r_[i] * dr ) 
    	endfor

    	rho_ = e0 * divE_

	endif

	;	Calculate the magnetic wave field

	hR		= complexArr ( runData.nR-1 )
	hPhi	= complexArr ( runData.nR )
	hz		= complexArr ( runData.nR )

	for i=0,runData.nR-1 do begin

		if i lt runData.nR-1 then $
			hR[i]	= -II * runData.kz * ePhi[i] + II * runData.nPhi * ez[i] / runData.r_[i]

		if i gt 0 and i lt runData.nR-1 then begin
		
			hPhi[i]	= II * runData.kz * eR[i] - ( ez[i] - ez[i-1] ) / runData.dR 
			hz[i]	= ( -II * runData.nPhi * eR[i] $
						+ ( runData.r_[i]*ePhi[i] $
							- runData.r_[i-1]*ePhi[i-1] ) / runData.dR ) $
					 	/ runData.r[i]
	
		endif

	endfor

	hR		= hR / ( II * wReal * u0 )
	hPhi	= hPhi / ( II * wReal * u0 )
	hz		= hz / ( II * wReal * u0 )


	if plotSolution then $
	plot_solution, runData.antLoc, runData.dR, runData.nR, $
		eR, ePhi, ez, $
		kR = kR, r_kR = runData.r, $
		r1 = runData.r, r2 = runData.r_, r3 = runData.r_

	;if plotSolution then $
	;plot_solution, runData.antLoc, runData.dR, runData.nR, $
	;	hR, hPhi, hz, $
	;	kR = kR, r_kR = runData.r, $
	;	r1 = runData.r_, r2 = runData.r, r3 = runData.r

	;	Determine the longitudinal / transvers nature of the solution

	if dispersion_generalised then begin

		ePhiFull	= complexArr ( runData.nR )
		ezFull		= complexArr ( runData.nR )
	
		for i=1,runData.nR-2 do begin
	
			ePhiFull[i]	=  ( ePhi[i] + ePhi[i-1] ) / 2.0
			ezFull[i]	=  ( ez[i] + ez[i-1] ) / 2.0
	
		endfor
	
		kMag1	= sqrt ( real_part ( kR[*,0])^2 + kPhi^2 + kz^2 )
		eMag	= sqrt ( real_part(eR)^2 + real_part(ePhiFull)^2 + real_part(ezFull)^2 )
		kDotE1	= real_part ( kR[*,0] ) * real_part ( eR ) $
					+ kPhi * real_part ( ePhiFull ) $
					+ kz * real_part ( ezFull ) 
		theta1	= aCos ( kDotE1 / ( kMag1 * eMag ) )
	
		kMag2	= sqrt ( real_part ( kR[*,1])^2 + kPhi^2 + kz^2 )
		kDotE2	= real_part ( kR[*,1] ) * real_part ( eR ) $
					+ kPhi * real_part ( ePhiFull ) $
					+ kz * real_part ( ezFull ) 
		theta2	= aCos ( kDotE2 / ( kMag2 * eMag ) )
	
		kMag3	= sqrt ( real_part ( kR[*,2])^2 + kPhi^2 + kz^2 )
		kDotE3	= real_part ( kR[*,2] ) * real_part ( eR ) $
					+ kPhi * real_part ( ePhiFull ) $
					+ kz * real_part ( ezFull ) 
		theta3	= aCos ( kDotE3 / ( kMag3 * eMag ) )
	
		kMag4	= sqrt ( real_part ( kR[*,3])^2 + kPhi^2 + kz^2 )
		kDotE4	= real_part ( kR[*,3] ) * real_part ( eR ) $
					+ kPhi * real_part ( ePhiFull ) $
					+ kz * real_part ( ezFull ) 
		theta4	= aCos ( kDotE4 / ( kMag4 * eMag ) )
	
		iPlot, runData.r, theta1 * !radeg, $
			sym_index = 4, lineStyle = 6, color = blue, $
			yRange = [0,180], window_title = 'k dot E angle'
		iPlot, runData.r, theta2 * !radeg, $
			/over, sym_index = 4, lineStyle = 6, color = green
		iPlot, runData.r, theta3 * !radeg, $
			/over, sym_index = 4, lineStyle = 6, color = red
		iPlot, runData.r, theta4 * !radeg, $
			/over, sym_index = 4, lineStyle = 6, color = purple

		eMag	= sqrt ( real_part(eR)^2 + real_part(ePhiFull)^2 + real_part(ezFull)^2 )
		eDotB	= real_part(eR) * runData.bField[*,0] $
					+ real_part(ePhiFull) * runData.bField[*,1] $
					+ real_part(ezFull) * runData.bField[*,2]
		thetaA	= aCos ( eDotB / ( eMag * runData.bMag ) )

		eMag	= sqrt ( imaginary(eR)^2 + imaginary(ePhiFull)^2 + imaginary(ezFull)^2 )
		eDotB	= imaginary(eR) * runData.bField[*,0] $
					+ imaginary(ePhiFull) * runData.bField[*,1] $
					+ imaginary(ezFull) * runData.bField[*,2]
		thetaB	= aCos ( eDotB / ( eMag * runData.bMag ) )


		iPlot, runData.r, thetaA * !radeg, $
			yRange = [0,180], $
			window_title = 'e dot B angle'		
		iPlot, runData.r, thetaB * !radeg, $
			/over	


	endif

	;	Write text file for comparison with GCC

	outData	= replicate ( { r : 0.0 , $
				eR : complex ( 0.0, 0.0 ), $
				ePhi : complex ( 0.0, 0.0 ), $
				ez : complex ( 0.0, 0.0 ) }, n_elements ( runData.r ) )

	outData.r 		= runData.r
	outData.eR		= eR
	outData[0:n_elements(ePhi)-1].ePhi	= ePhi
	outData[0:n_elements(ePhi)-1].ez		= ez


	openw, lun, 'dlg_solution.txt', /get_lun
	for i=0,n_elements(runData.r)-1 do $
		printf, lun, outData[i], $
			format = '(f8.5,2x,6(e12.4,2x))'
	close, lun
	
	close, /all


	if plotMovie then begin

		dt	= 1.0/wReal* 0.05
		rng	= max ( abs ( [eR,ePhi,ez] ) )
		for i=0,10000L do begin
			t	= i * dt
			plot, runData.r, real_part ( eR ) * cos ( wReal * t ) $
					+ imaginary ( eR ) * sin ( wReal * t ), $
				yRange = [-rng*1.5,rng*1.5]
			oplot, runData.r_, real_part ( ePhi ) * cos ( wReal * t ) $
					+ imaginary ( ePhi ) * sin ( wReal * t ), $
					lineStyle = 1
			oplot, runData.r_, real_part ( ez ) * cos ( wReal * t ) $
					+ imaginary ( ez ) * sin ( wReal * t ), $
					lineStyle = 2
			oPlot, [runData.antLoc, runData.antLoc ], $
					[-rng,rng], thick = 3
			oPlot, runData.r, specData[0].n * rng / max ( specData[0].n )
			wait, 0.025

		endfor

	endif
end
