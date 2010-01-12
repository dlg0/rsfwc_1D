;	Real Space Full Wave Code in 1D
;	DLG

pro rsfwc_1d, $
    eR = eR, $
	ePhi = ePhi, $
	ez = ez, $
    plot = plot, $
    bandStorage = bandStorage, $
	divD = divD, $
	rFull = r, rHalf = r_, $
	kR1 = kR1, kR2 = kR2

	if not keyword_set ( band ) then band = 1

;	Parameters

	common constants, e, me, mi, c, e0, u0, II
	e	= 1.60217646d-19
	me	= 9.10938188d-31
   	mi	= 1.67262158d-27
	c	= 299792458d0	
	e0	=	8.854187817d-12 
	u0	=	4d0 * !dpi * 10d0^(-7) 
	II	= dcomplex ( 0, 1 )

;	Run setup

	run_setup, runData = runData

 	wReal	= runData.freq * 2 * !pi 
	w	= dcomplexArr ( runData.nR ) + wReal
	;iiDamped	= where ( r gt 1.4 )
	w[*]	= complex ( wReal, wReal * runData.damping )

;	Get stix variables

	stixVariables, wReal, runData, stixVars = stixVars

;	Calculate dielectric

	dielectric, runData, stixVars, $
		epsilonFull = epsilon, $
		epsilonHalf = epsilon_

;	Dispersion analysis

	dispersion, wReal, epsilon, stixVars, runData

;	Fill matrix

	matFill, runData.nR, runData.nPhi, runData.kz, $
		runData.r, runData.r_, epsilon, epsilon_, w, runData.dR, $
		aMat = aMat, nAll = nAll, nuca = nuca, nlca = nlca

;	Anntenna 
 
	antLoc  = min ( runData.r_ )

	if antLoc gt runData.rMax or antLoc lt runData.rMin then begin
		print, 'ERROR: antenna outside domain, please correct.'
		stop
	endif

    iiAnt   = where ( abs ( runData.r_ - antLoc ) eq min ( abs ( runData.r_ - antLoc ) ) )
	rhs		= complexArr ( nAll )
	rhs[iiAnt*3+2]	= -II * wReal * u0 ( runData.r_[0] * runData.dr ) * 20d0
	;rhs[iiAnt*3]	= II * w * u0 / ( r_[0] * dr ) * 1d0/2
	;rhs[iiAnt*3+3]	= II * w * u0 / ( r_[0] * dr ) * 1d0/2


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

	plot_solution, runData.r, runData.r_, antLoc, runData.dR, runData.nR, $
		eR, ePhi, ez

	close, /all
end
