;	Real Space Full Wave Code in 1D
;	DLG

pro rsfwc_1d, $
    eR = eR, $
	ePhi = ePhi, $
	ez = ez, $
	divD = divD, $
	rFull = rFull, rHalf = rHalf, $
	;jA_r = jA_r, jA_t = jA_t, jA_z = jA_z, $
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
		poloidalScale = poloidalScale, $
		AntennaJ_r = AntennaJ_r, $
		AntennaJ_t = AntennaJ_t, $
		AntennaJ_z = AntennaJ_z, $
		antSig_r = antSig_r, $
		antSig_t = antSig_t, $
		antSig_z = antSig_z, $
		jAmp = jAmp


	rFull 	= runData.r
	rHalf	= runData.r_
	rOut	= rFull

 	wReal	= runData.freq * 2d0 * !dpi 
	wc	    = dcomplexArr ( runData.nR, runData.nSpec+1 ) 
	wc_	    = dcomplexArr ( runData.nR_, runData.nSpec+1 ) 

    for s=0,runData.nSpec do begin
	    wc[*,s]	= wReal * complex ( 1+fltArr(runData.nR), specData[s].nuOmg )
	    wc_[*,s]	= wReal * complex ( 1+fltArr(runData.nR_), specData[s].nuOmg_)
    endfor
    w = real_part(wc)
    w_ = real_part(wc_)

	stixVariables, w, w_, wc, wc_, specData, stixVars = stixVars, electronsOnly = 0

	dielectric, runData, stixVars, $
		epsilonFull_ = epsilon, $
		epsilonHalf_ = epsilon_, $
        epsilonFullSpec = epsilonS, $
        epsilonHalfSpec = epsilonS_, w, specData, $
        sigmaFull = sigma, $
        sigmaHalf = sigma_

	dispersionAll, wReal, epsilon, stixVars, runData, specData, $
		kR = kR, kPhi = kPhi, kz = kz

	if keyword_set ( dispersionOnly ) then return

	if kjInput or ar2EField then begin

		cdfId = ncdf_open(kj_jP_fileName)

            nCdf_varGet, cdfId, 'replace', replace
            nCdf_varGet, cdfId, 'replace_', replace_

			nCdf_varGet, cdfId, 'jP_r_re', kj_jpR_re
			nCdf_varGet, cdfId, 'jP_r_im', kj_jpR_im
			nCdf_varGet, cdfId, 'jP_p_re', kj_jpT_re
			nCdf_varGet, cdfId, 'jP_p_im', kj_jpT_im
			nCdf_varGet, cdfId, 'jP_z_re', kj_jpZ_re
			nCdf_varGet, cdfId, 'jP_z_im', kj_jpZ_im

			nCdf_varGet, cdfId, 'jP_r_re_', kj_jpR_re_
			nCdf_varGet, cdfId, 'jP_r_im_', kj_jpR_im_
			nCdf_varGet, cdfId, 'jP_p_re_', kj_jpT_re_
			nCdf_varGet, cdfId, 'jP_p_im_', kj_jpT_im_
			nCdf_varGet, cdfId, 'jP_z_re_', kj_jpZ_re_
			nCdf_varGet, cdfId, 'jP_z_im_', kj_jpZ_im_

			nCdf_varGet, cdfId, 'E_r_re', kj_ER_re
			nCdf_varGet, cdfId, 'E_r_im', kj_ER_im
			nCdf_varGet, cdfId, 'E_p_re', kj_ET_re
			nCdf_varGet, cdfId, 'E_p_im', kj_ET_im
			nCdf_varGet, cdfId, 'E_z_re', kj_EZ_re
			nCdf_varGet, cdfId, 'E_z_im', kj_EZ_im

			nCdf_varGet, cdfId, 'E_r_re_', kj_ER_re_
			nCdf_varGet, cdfId, 'E_r_im_', kj_ER_im_
			nCdf_varGet, cdfId, 'E_p_re_', kj_ET_re_
			nCdf_varGet, cdfId, 'E_p_im_', kj_ET_im_
			nCdf_varGet, cdfId, 'E_z_re_', kj_EZ_re_
			nCdf_varGet, cdfId, 'E_z_im_', kj_EZ_im_

		ncdf_close, cdfId

		kj_jpR = complex(kj_jpR_re,kj_jpR_im)
		kj_jpT = complex(kj_jpT_re,kj_jpT_im)
		kj_jpZ = complex(kj_jpZ_re,kj_jpZ_im)

		kj_jpR_ = complex(kj_jpR_re_,kj_jpR_im_)
		kj_jpT_ = complex(kj_jpT_re_,kj_jpT_im_)
		kj_jpZ_ = complex(kj_jpZ_re_,kj_jpZ_im_)

		kj_ER = complex(kj_ER_re,kj_ER_im)
		kj_ET = complex(kj_ET_re,kj_ET_im)
		kj_EZ = complex(kj_EZ_re,kj_EZ_im)

		kj_ER_ = complex(kj_ER_re_,kj_ER_im_)
		kj_ET_ = complex(kj_ET_re_,kj_ET_im_)
		kj_EZ_ = complex(kj_EZ_re_,kj_EZ_im_)

        ; Try applying a windowing function to the
        ; "replace" section to avoid the discontinuity
        ; that is causing the problem ...

        iiReplace = where(replace gt 0, iiReplaceCnt)
        HanningReplace = hanning(iiReplaceCnt)
        xGrid = fIndGen(n_elements(replace))
        xGrid_ = xGrid[0:-2]+0.5

        HanningFull = fltArr(n_elements(replace))
        HanningFull[iiReplace] = HanningReplace
        HanningHalf = interpol(HanningFull,xGrid,xGrid_,/spline)

        HanningFull = HanningFull*0+1
        HanningHalf = HanningHalf*0+1

        kjIn = { $
            jpR : kj_jpR*HanningFull, $
            jpT : kj_jpT*HanningFull, $
            jpZ : kj_jpZ*HanningFull, $
            jpR_ : kj_jpR_*HanningHalf, $
            jpT_ : kj_jpT_*HanningHalf, $
            jpZ_ : kj_jpZ_*HanningHalf, $
            ER : kj_ER*HanningFull, $
            ET : kj_ET*HanningFull, $
            EZ : kj_EZ*HanningFull, $
            ER_ : kj_ER_*HanningHalf, $
            ET_ : kj_ET_*HanningHalf, $
            EZ_ : kj_EZ_*HanningHalf, $
            replace : replace, $
            replace_ : replace_ }

	endif else begin
        replace = intArr(runData.nR)
        replace_ = intArr(runData.nR-1)
    endelse

	if runData.antLoc gt runData.rMax or runData.antLoc lt runData.rMin then begin
		print, 'ADJUSTING: Antenna location was ', runData.antLoc
		if runData.antLoc lt min ( runData.r_ ) then runData.antLoc = min ( runData.r_ )
		if runData.antLoc gt max ( runData.r_ ) then runData.antLoc = max ( runData.r_ )
		print, 'ADJUSTING: Antenna location is now ', runData.antLoc
	endif

    iiAnt   = where ( abs ( runData.r_ - runData.antLoc ) $
						eq min ( abs ( runData.r_ - runData.antLoc ) ) )

	nAll	= runData.nR + 2L * ( runData.nR - 1L )
	rhs		= complexArr ( nAll )

	;if not keyword_set ( jA_r ) then jA_R = 0
	;if not keyword_set ( jA_t ) then jA_t = 0
	;if not keyword_set ( jA_z ) then jA_z = 1

	;antSigX = (runData.rMax-runData.rMin)/n_elements(runData.r)*15
	;if jAmp gt 0 then jAmp = jAmp else jAmp 

	;antk = -2000
	;AntFluct = cos(runData.r*antk)+ii*sin(runData.r*antk)
	;AntFluct_ = cos(runData.r_*antk)+ii*sin(runData.r_*antk)

	TmpAntJr = jAmp*exp( -( (runData.r-runData.antLoc)^2/antSig_r^2 ) )
	TmpAntJr_ = jAmp*exp( -( (runData.r_-runData.antLoc)^2/antSig_r^2 ) )
	TmpAntJt = jAmp*exp( -( (runData.r-runData.antLoc)^2/antSig_t^2 ) )
	TmpAntJt_ = jAmp*exp( -( (runData.r_-runData.antLoc)^2/antSig_t^2 ) )
	TmpAntJz = jAmp*exp( -( (runData.r-runData.antLoc)^2/antSig_z^2 ) )
	TmpAntJz_ = jAmp*exp( -( (runData.r_-runData.antLoc)^2/antSig_z^2 ) )


    if AntennaJ_r then jA_r	= TmpAntJr else jA_r = TmpAntJr*0
    if AntennaJ_t then jA_t	= TmpAntJt else jA_t = TmpAntJt*0
    if AntennaJ_z then jA_z	= TmpAntJz else jA_z = TmpAntJz*0

	if AntennaJ_r then jA_r_ = TmpAntJr_ else jA_r_ = TmpAntJr_*0
    if AntennaJ_t then jA_t_ = TmpAntJt_ else jA_t_ = TmpAntJt_*0
    if AntennaJ_z then jA_z_ = TmpAntJz_ else jA_z_ = TmpAntJz_*0


   	; NOTE: the equations were derived as +curlcurlE - k0edotE = +iwuJa
	; and as such, the sign in the below is +ve, NOT -ve as in the aorsa
	; formulation ;)	
   	for i=0,runData.nR-2 do begin

		rhs[i*3+2]	= II * wReal * u0 * jA_z_[i]
		rhs[i*3+1]	= II * wReal * u0 * jA_t_[i]
		rhs[i*3]	= II * wReal * u0 * jA_r[i]

		if kjInput then begin
			rhs[i*3+2]	+= II * wReal * u0 * kjIn.jpZ_[i]
			rhs[i*3+1]	+= II * wReal * u0 * kjIn.jpT_[i]
			rhs[i*3]	+= II * wReal * u0 * kjIn.jpR[i]
		endif

		;if ar2EField then begin
        ;    if replace[i] then begin
        ;       ;print, 'Adding AORSA driving term ...'
		;       rhs[i*3]	+= II * wReal * u0 * kjIn.jpR[i]
		;    endif
        ;    if replace[i] and replace[i+1] then begin
		;       rhs[i*3+2]	+= II * wReal * u0 * kjIn.jpZ_[i]
		;       rhs[i*3+1]	+= II * wReal * u0 * kjIn.jpT_[i]
        ;    endif	
        ;endif

	endfor

    ;ii3 = IndGen(RunData.nR)*3
    ;rhs_R = rhs[ii3]
    ;rhs_T = rhs[ii3[0:-2]+1]
    ;rhs_Z = rhs[ii3[0:-2]+2]

    ;p=plot(RunData.r,rhs_R)
    ;p=plot(RunData.r,imaginary(rhs_R),/over,color='r')
    ;p=plot(RunData.r_,rhs_T,/over)
    ;p=plot(RunData.r_,imaginary(rhs_T),color='r',/over)
    ;p=plot(RunData.r_,rhs_Z,/over)
    ;p=plot(RunData.r_,imaginary(rhs_z),color='r',/over)

	matFill, runData.nR, runData.nPhi, runData.kz, $
		runData.r, runData.r_, epsilon, epsilon_, w, runData.dR, $
		aMat = aMat, nAll = nAll, nuca = nuca, nlca = nlca, $
        replaceFull = replace, replaceHalf = replace_, rhs=rhs, kjIn = kjIn

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

	hR		= hR / ( II * w[*,0] * u0 )
	hPhi	= hPhi / ( II * w[*,0] * u0 )
	hz		= hz / ( II * w[*,0] * u0 )

	ePhiFull	= complexArr ( runData.nR )
	ezFull		= complexArr ( runData.nR )
	hRFull  	= complexArr ( runData.nR )

    ;; Remove the Inident Field from the solution
    ;if ar2EField then begin
    ;    eR = eR - kjIn.eR*kjIn.replace
    ;    ePhi = ePhi - kjIn.eT_*kjIn.replace_
    ;    eZ = eZ - kjIn.eZ_*kjIn.replace_
    ;endif

	for i=1,runData.nR-2 do begin
	
		ePhiFull[i]	=  ( ePhi[i] + ePhi[i-1] ) / 2.0
		ezFull[i]	=  ( ez[i] + ez[i-1] ) / 2.0
		hRFull[i]	=  ( hR[i] + hR[i-1] ) / 2.0

	endfor

	if plotESolution then $
	rs_plot_solution, runData.antLoc, runData.dR, runData.nR, $
		eR, ePhi, ez, $
		kR = kR, r_kR = runData.r, $
		r1 = runData.r, r2 = runData.r_, r3 = runData.r_, kjIn=kjIn

	if plotHSolution then begin
	    ;rs_plot_solution, runData.antLoc, runData.dR, runData.nR, $
	    ;	hR, hPhi, hz, $
	    ;	kR = kR, r_kR = runData.r, $
	    ;	r1 = runData.r_, r2 = runData.r, r3 = runData.r
        p=plot(runData.r,hr,layout=[1,3,1],title='h_r',window_title='rsfwc_1d')
        p=plot(runData.r,imaginary(hr),/over,color='r')
        p=plot(runData.r,hphi,layout=[1,3,2],/current,title='h_t')
        p=plot(runData.r,imaginary(hphi),/over,color='r')
        p=plot(runData.r,hz,layout=[1,3,3],/current,title='h_z')
        p=plot(runData.r,imaginary(hz),/over,color='r')
       
        e_r = eR
        e_t = ePhiFull
        e_z = eZFull
        h_r = hRFull
        h_t = hPhi
        h_z = hZ

        S_r = e_t*conj(h_z)-e_z*conj(h_t)
        S_t = -(e_r*conj(h_z)-e_z*conj(h_r))
        S_z = e_r*conj(h_t)-e_t*conj(h_r)

        S_r = 0.5*real_part(S_r)
        S_t = 0.5*real_part(S_t)
        S_z = 0.5*real_part(S_z)
    
        SRange = max(abs([S_r,S_t,S_z])) 
        p=plot(runData.r,S_r,layout=[1,3,1],title='Time Averaged Poynting Vector S_r',$
                window_title='rsfwc_1d', yRange=[-SRange,SRange])
        p=plot(runData.r,S_t,layout=[1,3,2],/current,title='S_t', yRange=[-SRange,SRange])
        p=plot(runData.r,S_z,layout=[1,3,3],/current,title='S_z', yRange=[-SRange,SRange])

    endif

	;	Determine the longitudinal / transvers nature of the solution

	if dispersion_generalised then begin

		if plotKdotE then begin	
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
		endif

		if plotEdotB then begin
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


			p=plot(runData.r, thetaA * !radeg, $
					lineStyle=6, symbol="D",yRange = [0,180], title = 'e dot B angle')
			p=plot(runData.r, thetaB * !radeg, $
					/over, lineStyle=6, symbol="S")
		endif

	endif

	e_r	= eR
	e_t	= ePhiFull
	e_z	= eZFull

	; Calculate plasma current
	; ------------------------

    jP_r_S = complexArr(runData.nR,runData.nSpec+1)
    jP_t_S = complexArr(runData.nR,runData.nSpec+1)
    jP_z_S = complexArr(runData.nR,runData.nSpec+1)
    
    jP_r = complexArr(runData.nR)
    jP_t = complexArr(runData.nR)
    jP_z = complexArr(runData.nR)

    for i=1,runData.nR-2 do begin
    
    	if kjInput and replace[i] then begin
    
    		jP_r[i]	= kj_jpR[i]
    		jP_t[i]	= kj_jpT[i]
    		jP_z[i]	= kj_jpZ[i]
            print, 'REPLACING THINGS ...'    
    	endif else begin
    
            for s=0,runData.nSpec do begin
                    ;thisSigma = (epsilonS[*,*,s,i]-identity(3))*wReal*e0/II
                    thisSigma = sigma[*,*,s,i]
                    thisE = [[e_r[i]],[e_t[i]],[e_z[i]]]
    
                    ; Changing this now matches the AORSA
                    ; jP, BUT is it the right one, or are they BOTH 
                    ; wrong now?
    
                    ;this_jP = thisSigma # transpose(thisE)
                    this_jP = thisSigma ## thisE
  
                    jP_r_S[i,s] = this_jP[0]
                    jP_t_S[i,s] = this_jP[1]
                    jP_z_S[i,s] = this_jP[2]
    
            endfor

            jP_r[i] = total(jP_r_S[i,*],2)
            jP_t[i] = total(jP_t_S[i,*],2)
            jP_z[i] = total(jP_z_S[i,*],2)

    	endelse
    endfor


	; jDotE
	; -----

    jPDotE_S = fltArr(runData.nR,runData.nSpec+1)
    for s=0,runData.nSpec do begin
        jPDotE_S[*,s] = 0.5 * real_part ( jP_r_S[*,s] * conj(e_r) $
				+ jP_t_S[*,s] * conj(e_t) $
				+ jP_z_S[*,s] * conj(e_z) )
    endfor

	jPDotE	= 0.5 * real_part ( jP_r * conj(e_r) $
				+ jP_t * conj(e_t) $
				+ jP_z * conj(e_z) )

	if plotJdotE then begin

		p = plot (runData.r,jPDotE,color='b',thick=1,$
				title='J dot E',name='jDote_0',font_size=10,$
				layout=[1,3,1],window_title='rsfwc_1d')

        colorArr = ['blue','red','green','orange']
        for s=0,n_elements(jPDotE_S[0,*])-1 do begin
            p = plot(runData.r,jPDotE_S[*,s],color=colorArr[s],/over,thick=3,linestyle='--',transparency=50)
        endfor

        p=plot(runData.r,specData[0].nuOmg,title='nuOmg',/current,layout=[1,3,2])

		; jAnt
		; ---

		pr = plot (runData.r,jA_r,color='b',thick=3,transparency=50,$
				title='jAnt',name='jAnt_r',font_size=10,$
				layout=[1,3,3],/current)
		pt = plot (runData.r,jA_t,color='r',thick=3,transparency=50,$
				name='jAnt_t',/over)
		pz = plot (runData.r,jA_z,color='g',thick=3,transparency=50,$
				name='jAnt_z',/over)

        if ar2EField then begin
            p = plot(runData.r,kjIn.jPr*kjIn.replace,linestyle='--',/over)
            p = plot(runData.r,imaginary(kjIn.jPr)*kjIn.replace,linestyle='-',/over)
            p = plot(runData.r,kjIn.jPt*kjIn.replace,linestyle='--',/over,color='r')
            p = plot(runData.r,imaginary(kjIn.jPt)*kjIn.replace,linestyle='-',/over,color='r')
            p = plot(runData.r,kjIn.jPz*kjIn.replace,linestyle='--',/over,color='g')
            p = plot(runData.r,imaginary(kjIn.jPz)*kjIn.replace,linestyle='-',/over,color='g')
        endif



		;l = legend(target=[pr,pt,pz],position=[0.8,0.4],/norm,font_size=10)

	endif

	;save, $
	;	jP_r, jP_t, jP_z, $
	;	e_r, e_t, e_z, $
	;	fileName = 'solutionVals.sav'

	; Write netCDF file

	nc_id = nCdf_create ('rsfwc_1d_'+runIdent+'.nc', /clobber )

	nCdf_control, nc_id, /fill
	
	nr_id = nCdf_dimDef ( nc_id, 'nR', runData.nR )
	nrH_id = nCdf_dimDef ( nc_id, 'nR_', runData.nR-1 )
	scalar_id = nCdf_dimDef ( nc_id, 'scalar', 1 )

	freq_id = nCdf_varDef ( nc_id, 'freq', scalar_id, /float )
	r_id = nCdf_varDef ( nc_id, 'r', nr_id, /float )
	rH_id = nCdf_varDef ( nc_id, 'r_', nrH_id, /float )
	z_id = nCdf_varDef ( nc_id, 'z', nr_id, /float )
	zH_id = nCdf_varDef ( nc_id, 'z_', nrH_id, /float )

	B0_r_id = nCdf_varDef ( nc_id, 'B0_r', nr_id, /float )
	B0_p_id = nCdf_varDef ( nc_id, 'B0_p', nr_id, /float )
	B0_z_id = nCdf_varDef ( nc_id, 'B0_z', nr_id, /float )

	e_r_re_id = nCdf_varDef ( nc_id, 'e_r_re', nr_id, /float )
	e_r_im_id = nCdf_varDef ( nc_id, 'e_r_im', nr_id, /float )
	e_p_re_id = nCdf_varDef ( nc_id, 'e_p_re', nr_id, /float )
	e_p_im_id = nCdf_varDef ( nc_id, 'e_p_im', nr_id, /float )
	e_z_re_id = nCdf_varDef ( nc_id, 'e_z_re', nr_id, /float )
	e_z_im_id = nCdf_varDef ( nc_id, 'e_z_im', nr_id, /float )

	jP_r_re_id = nCdf_varDef ( nc_id, 'jP_r_re', nr_id, /float )
	jP_r_im_id = nCdf_varDef ( nc_id, 'jP_r_im', nr_id, /float )
	jP_p_re_id = nCdf_varDef ( nc_id, 'jP_p_re', nr_id, /float )
	jP_p_im_id = nCdf_varDef ( nc_id, 'jP_p_im', nr_id, /float )
	jP_z_re_id = nCdf_varDef ( nc_id, 'jP_z_re', nr_id, /float )
	jP_z_im_id = nCdf_varDef ( nc_id, 'jP_z_im', nr_id, /float )

	jA_r_re_id = nCdf_varDef ( nc_id, 'jA_r_re', nr_id, /float )
	jA_r_im_id = nCdf_varDef ( nc_id, 'jA_r_im', nr_id, /float )
	jA_p_re_id = nCdf_varDef ( nc_id, 'jA_p_re', nr_id, /float )
	jA_p_im_id = nCdf_varDef ( nc_id, 'jA_p_im', nr_id, /float )
	jA_z_re_id = nCdf_varDef ( nc_id, 'jA_z_re', nr_id, /float )
	jA_z_im_id = nCdf_varDef ( nc_id, 'jA_z_im', nr_id, /float )

	nCdf_control, nc_id, /enDef

	nCdf_varPut, nc_id, freq_id, runData.freq

	nCdf_varPut, nc_id, r_id, runData.r
	nCdf_varPut, nc_id, rH_id, runData.r_

	nCdf_varPut, nc_id, z_id, runData.z
	nCdf_varPut, nc_id, zH_id, runData.z_

	nCdf_varPut, nc_id, B0_r_id, runData.BField[*,0]
	nCdf_varPut, nc_id, B0_p_id, runData.BField[*,1]
	nCdf_varPut, nc_id, B0_z_id, runData.BField[*,2]

	nCdf_varPut, nc_id, e_r_re_id, real_part(e_r) 
	nCdf_varPut, nc_id, e_r_im_id, imaginary(e_r) 
	nCdf_varPut, nc_id, e_p_re_id, real_part(e_t) 
	nCdf_varPut, nc_id, e_p_im_id, imaginary(e_t) 
	nCdf_varPut, nc_id, e_z_re_id, real_part(e_z) 
	nCdf_varPut, nc_id, e_z_im_id, imaginary(e_z) 

	nCdf_varPut, nc_id, jP_r_re_id, real_part(jP_r)
	nCdf_varPut, nc_id, jP_r_im_id, imaginary(jP_r) 
	nCdf_varPut, nc_id, jP_p_re_id, real_part(jP_t) 
	nCdf_varPut, nc_id, jP_p_im_id, imaginary(jP_t) 
	nCdf_varPut, nc_id, jP_z_re_id, real_part(jP_z) 
	nCdf_varPut, nc_id, jP_z_im_id, imaginary(jP_z) 

	nCdf_varPut, nc_id, jA_r_re_id, jA_r 
	nCdf_varPut, nc_id, jA_r_im_id, jA_r*0 
	nCdf_varPut, nc_id, jA_p_re_id, jA_t
	nCdf_varPut, nc_id, jA_p_im_id, jA_t*0
	nCdf_varPut, nc_id, jA_z_re_id, jA_z
	nCdf_varPut, nc_id, jA_z_im_id, jA_z*0

	nCdf_close, nc_id


	if plotJp then begin

		jpRange = max(abs([abs(jp_r),abs(jp_t),abs(jp_z)]))
		yRange = [-jPRange,jPRange]*0.5

		p_re = plot (rFull,jP_r,thick=2,$
				title='jP_r',name='Jp_re',font_size=10,$
				layout=[1,3,1],yRange=yRange,transparency=50,window_title='rsfwc_1d')
		p_im = plot (rFull,imaginary(jP_r),color='r',thick=2,transparency=50,$
				name='Jp_re',font_size=10,/over)
	   	;l = legend(target=[p_re,p_im],position=[0.99,0.95],/norm,font_size=10,horizontal_alignment='RIGHT')

		p_re = plot (rFull,jP_t,thick=2,$
				title='jP_t',name='Jp_re',font_size=10,$
				layout=[1,3,2],/current,yRange=yRange,transparency=50)
		p_im = plot (rFull,imaginary(jP_t),color='r',thick=2,transparency=50,$
				name='Jp_re',font_size=10,/over)
	   	;l = legend(target=[p_re,p_im],position=[0.99,0.63],/norm,font_size=10,horizontal_alignment='RIGHT')

		p_re = plot (rFull,jP_z,thick=2,$
				title='jP_z',name='Jp_re',font_size=10,$
				layout=[1,3,3],/current,yRange=yRange,transparency=50)
		p_im = plot (rFull,imaginary(jP_z),color='r',thick=2,transparency=50,$
				name='Jp_re',font_size=10,/over)
	   	;l = legend(target=[p_re,p_im],position=[0.99,0.25],/norm,font_size=10,horizontal_alignment='RIGHT')


	endif

	;;	Write text file for comparison with GCC

	;outData	= replicate ( { r : 0.0 , $
	;			eR : complex ( 0.0, 0.0 ), $
	;			ePhi : complex ( 0.0, 0.0 ), $
	;			ez : complex ( 0.0, 0.0 ) }, n_elements ( runData.r ) )

	;outData.r 		= runData.r
	;outData.eR		= eR
	;outData[0:n_elements(ePhi)-1].ePhi	= ePhi
	;outData[0:n_elements(ePhi)-1].ez		= ez

	;openw, lun, 'dlg_solution.txt', /get_lun
	;for i=0,n_elements(runData.r)-1 do $
	;	printf, lun, outData[i], $
	;		format = '(f8.5,2x,6(e12.4,2x))'
	;close, lun
	;
	;close, /all


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
