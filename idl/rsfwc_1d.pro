;	Real Space Full Wave Code in 1D
;	DLG

pro rsfwc_1d, $
    eR = eR, $
	eP = eP, $
	ez = ez, $
	divD = divD, $
	r = r, r_ = r_, $
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
	nPhi = nPhi, $
    plot_Jp_spec = plot_Jp_spec, $
    cartesian_offset = _cartesian_offset
  
    if keyword_set(_cartesian_offset) then cartesian_offset = _cartesian_offset else cartesian_offset = 0

;	Parameters

    @dlg_constants
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
		jAmp = jAmp, ar2 = ar2


	r 	= runData.r
	r_	= runData.r_
    nR = n_elements(r)
    nR_ = n_elements(r_)
	rOut	= r

 	wReal	= runData.freq * 2d0 * !dpi 
	wc	    = dcomplexArr ( runData.nR, runData.nSpec ) 
	wc_	    = dcomplexArr ( runData.nR_, runData.nSpec ) 

    for s=0,runData.nSpec-1 do begin
	    wc[*,s]	    = wReal 
	    wc_[*,s]	= wReal 
    endfor
    w = real_part(wc)
    w_ = real_part(wc_)

	stixVariables, w, w_, wc, wc_, specData, stixVars = stixVars

	dielectric, runData, stixVars, w, specData, $
		epsilonFull_ = epsilon, $
		epsilonHalf_ = epsilon_, $
        epsilonFullSpec = epsilonS, $
        epsilonHalfSpec = epsilonS_, $
        sigmaFull_ = sigma, $
        sigmaHalf_ = sigma_, $
		sigmaFullSpec_abp = sigma_abp, $
		sigmaHalfSpec_abp = sigma_abp_

	dispersionAll, wReal, epsilon, stixVars, runData, specData, $
		kR = kR, kPhi = kPhi, kz = kz

	if keyword_set ( dispersionOnly ) then return

	if kjInput then begin
        kjIn = kj_read_kj_jP(kjDeltaFileName)
    endif else if ar2EField then begin
        replace = intArr(runData.nR)*0+1
        replace_ = intArr(runData.nR-1)*0+1
	endif else begin
        replace = intArr(runData.nR)
        replace_ = intArr(runData.nR-1)
    endelse


	nAll	= runData.nR + 2L * ( runData.nR - 1L )
	rhs		= dcomplexArr ( nAll )


    jA_r = ar2('jAnt_r')
    jA_t = ar2('jAnt_t')
    jA_z = ar2('jAnt_z')

    jA_r_ = interpol(jA_r,r,r_,/spline) 
    jA_t_ = interpol(jA_t,r,r_,/spline) 
    jA_z_ = interpol(jA_z,r,r_,/spline) 

   	; NOTE: the equations were derived as +curlcurlE - k0edotE = +iwuJa
	; and as such, the sign in the below is +ve, NOT -ve as in the aorsa
	; formulation ;)	

   	for i=0,runData.nR-1 do begin

		rhs[i*3]	= _II * wReal * _u0 * jA_r[i]
        if i lt nR-1 then begin
		rhs[i*3+1]	= _II * wReal * _u0 * jA_t_[i]
		rhs[i*3+2]	= _II * wReal * _u0 * jA_z_[i]
        endif

		if kjInput then begin
			rhs[i*3]	+= _II * wReal * _u0 * total(kjIn.jpR[i,0,*],3)
            if i lt nR-1 then begin
			rhs[i*3+1]	+= _II * wReal * _u0 * total(kjIn.jpT_[i,0,*],3)
			rhs[i*3+2]	+= _II * wReal * _u0 * total(kjIn.jpZ_[i,0,*],3)
            endif
		endif

	endfor

    if plotRHS then begin

        p=plot(r,jA_r,layout=[1,3,1], title='RHS (jA and delta from kJ)')
        p=plot(r,imaginary(jA_r),color='r',/over)
        if kjInput then begin
            p=plot(r,kjIn.jPr,/over,thick=2)
            p=plot(r,imaginary(kjIn.jPr),/over,color='r',thick=2)
        endif

        p=plot(r,jA_t,layout=[1,3,2], /current)
        p=plot(r,imaginary(jA_t),color='r',/over)
        if kjInput then begin
            p=plot(r,kjIn.jPt,/over,thick=2)
            p=plot(r,imaginary(kjIn.jPt),/over,color='r',thick=2)
        endif

        p=plot(r,jA_z,layout=[1,3,3], /current)
        p=plot(r,imaginary(jA_z),color='r',/over)
        if kjInput then begin
            p=plot(r,kjIn.jPz,/over,thick=2)
            p=plot(r,imaginary(kjIn.jPz),/over,color='r',thick=2)
        endif

    endif

;   Write the run-data to file

    runDataHash = HASH( $
            "freq",runData.freq, $
            "br",runData.bField[*,0], $
            "bt",runData.bField[*,1], $
            "bz",runData.bField[*,2], $
            "r",runData.r, $
            "nPhi",runData.nPhi,$
            "densitySpec",specData.n,$
            "nuOmg",specData.nuOmg,$
            "kz", runData.kz )

    runDataHash = runDataHash + HASH('jA_r_re',real_part(jA_r));
    runDataHash = runDataHash + HASH('jA_r_im',imaginary(jA_r));
    runDataHash = runDataHash + HASH('jA_t_re',real_part(jA_t));
    runDataHash = runDataHash + HASH('jA_t_im',imaginary(jA_t));
    runDataHash = runDataHash + HASH('jA_z_re',real_part(jA_z));
    runDataHash = runDataHash + HASH('jA_z_im',imaginary(jA_z));

    NCDF_PUT, 'output/rs-rundata.nc', /NEW, VARIABLES=runDataHash

	matFill, runData.nR, runData.nPhi, runData.kz, $
		runData.r, runData.r_, epsilon, epsilon_, wReal, runData.dR, $
		aMat = aMat, nAll = nAll, nuca = nuca, nlca = nlca, $
        replaceFull = replace, replaceHalf = replace_, rhs=rhs

;	Solve matrix

    print, '*** solving linear system'

    if bandStorage then begin

        eField = imsl_sp_bdSol ( rhs, nlca, nuca, aMat[*], /double )

    endif else begin

	    eField	= la_linear_equation ( aMat, rhs, status = stat, /double )
	    print, 'lapack status: ', stat

    endelse

;   Calculate residual

    lhs = aMat ## eField
    res = lhs - rhs
    resNorm = norm(res)

    print, 'Residual: ', resNorm


; Store A in sparse storage (re and im parts seperately - stupid IDL)

    a_re = real_part(aMat)
    a_im = imaginary(aMat)

    aSparse_re = sprsIn(a_re,/double)
    aSparse_im = sprsIn(a_im,/double)

    aFull_re = fulstr(aSparse_re)
    aFull_im = fulstr(aSparse_im)

    save, aSparse_re, aSparse_im, r, r_, eField, rhs, fileName = 'rs-amat.sav'


    eP_LeftBoundary_ = eField[0]
    eZ_LeftBoundary_ = eField[1]
    eP_RightBoundary_ = eField[-2]
    eZ_RightBoundary_ = eField[-1]

    ; Remove BC layers
    rhs = rhs[2:-3]
    eField = eField[2:-3]
  
	ii_eR	= lIndGen(runData.nR)*3
	eR	= eField[ii_eR]
	ii_eP	= ii_eR[0:runData.nR-2]+1
	eP_	= eField[ii_eP]
	ii_ez	= ii_eP+1
	ez_	= eField[ii_ez]


	if keyword_set ( divD ) then begin

	;  	Calculate the Div of D @ the z,phi grid pts
	;  	but first we need to invert epsilon to get D & D_

    	divD_   = dcomplexArr ( runData.nR - 1 )

    	print, '*** calculating div D'
    	for i=1UL,runData.nR-3L do begin 
    	    
    	    Dri     = epsilon[0,0,i] * eR[i] $
    	            + ( epsilon_[1,0,i-1] * eP_[i-1] + epsilon_[1,0,i] * eP_[i] )/2 $
    	            + ( epsilon_[2,0,i-1] * ez_[i-1] + epsilon_[2,0,i] * ez_[i] )/2
    	    Dri1    =  epsilon[0,0,i+1] * eR[i+1] $
    	            + ( epsilon_[1,0,i] * eP_[i] + epsilon_[1,0,i+1] * eP_[i+1] )/2 $
    	            + ( epsilon_[2,0,i] * ez_[i] + epsilon_[2,0,i+1] * ez_[i+1] )/2
    	    Dth     = ( epsilon[0,1,i] * eR[i] + epsilon[0,1,i+1] * eR[i+1] )/2 $
    	            + epsilon[1,1,i] * eP_[i] $
    	            + epsilon[2,1,i] * ez_[i] 
    	    Dz      = ( epsilon[0,2,i] * eR[i] + epsilon[0,2,i+1] * eR[i+1] )/2 $
    	            + epsilon[1,2,i] * eP_[i] $
    	            + epsilon[2,2,i] * ez_[i] 

    	    divD_[i]    = ( r[i+1] * e0 * Dri1 - r[i] * e0 * Dri ) / ( r_[i] * dr ) $
    	                    + II * nPhi / r_[i] * e0 * Dth $
    	                    + II * kz * e0 * Dz
    	endfor

    	divE_   = dcomplexArr ( runData.nR - 1 )
    	for i=0UL,runData.nR-2L do begin
    	    divE_[i]    = II * nPhi * eP_[i] / r_[i] $
    	                    + II * kz * ez_[i] $
    	                    + ( r[i+1] * eR[i+1] - r[i] * eR[i] ) / ( r_[i] * dr ) 
    	endfor

    	rho_ = e0 * divE_

	endif

	;	Calculate the magnetic wave field

	hR_	= complexArr ( runData.nR-1 )
	hP	= complexArr ( runData.nR )
	hz	= complexArr ( runData.nR )

	for i=0,runData.nR-1 do begin

		if i lt runData.nR-1 then $
			hR_[i]	= -_II * runData.kz * eP_[i] + _II * runData.nPhi * ez_[i] / runData.r_[i]

		if i gt 0 and i lt runData.nR-1 then begin
		
			hP[i]	= _II * runData.kz * eR[i] - ( ez_[i] - ez_[i-1] ) / runData.dR 
			hz[i]	= ( -_II * runData.nPhi * eR[i] $
						+ ( runData.r_[i]*eP_[i] $
							- runData.r_[i-1]*eP_[i-1] ) / runData.dR ) $
					 	/ runData.r[i]
	
		endif

	endfor

	hR_	= hR_ / ( _II * w[*,0] * _u0 )
	hP	= hP / ( _II * w[*,0] * _u0 )
	hz	= hz / ( _II * w[*,0] * _u0 )

	eP	= complexArr ( runData.nR )
	eZ	= complexArr ( runData.nR )
	hR 	= complexArr ( runData.nR )

	for i=1,runData.nR-2 do begin
		eP[i]	=  ( eP_[i] + eP_[i-1] ) / 2.0
		eZ[i]	=  ( ez_[i] + ez_[i-1] ) / 2.0
		hR[i]	=  ( hR_[i] + hR_[i-1] ) / 2.0
	endfor
    eP[0] = (eP_LeftBoundary_+eP_[0])/2.0
    eZ[0] = (eZ_LeftBoundary_+eZ_[0])/2.0
    eP[-1] = (eP_RightBoundary_+eP_[-1])/2.0
    eZ[-1] = (eZ_RightBoundary_+eZ_[-1])/2.0


	if plotESolution then $
	rs_plot_solution, runData.antLoc, runData.dR, runData.nR, $
		eR, eP_, ez_, $
		kR = kR, r_kR = runData.r, $
		r1 = runData.r, r2 = runData.r_, r3 = runData.r_

	if plotHSolution then begin
        p=plot(runData.r,hr,layout=[1,3,1],title='h_r',window_title='rsfwc_1d')
        p=plot(runData.r,imaginary(hr),/over,color='r')
        p=plot(runData.r,hP,layout=[1,3,2],/current,title='h_t')
        p=plot(runData.r,imaginary(hP),/over,color='r')
        p=plot(runData.r,hz,layout=[1,3,3],/current,title='h_z')
        p=plot(runData.r,imaginary(hz),/over,color='r')
       
        e_r = eR
        e_t = eP
        e_z = eZ
        h_r = hR
        h_t = hP
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
			eMag	= sqrt ( real_part(eR)^2 + real_part(eP)^2 + real_part(eZ)^2 )
			kDotE1	= real_part ( kR[*,0] ) * real_part ( eR ) $
						+ kPhi * real_part ( eP ) $
						+ kz * real_part ( eZ ) 
			theta1	= aCos ( kDotE1 / ( kMag1 * eMag ) )
	
			kMag2	= sqrt ( real_part ( kR[*,1])^2 + kPhi^2 + kz^2 )
			kDotE2	= real_part ( kR[*,1] ) * real_part ( eR ) $
						+ kPhi * real_part ( eP ) $
						+ kz * real_part ( eZ ) 
			theta2	= aCos ( kDotE2 / ( kMag2 * eMag ) )
	
			kMag3	= sqrt ( real_part ( kR[*,2])^2 + kPhi^2 + kz^2 )
			kDotE3	= real_part ( kR[*,2] ) * real_part ( eR ) $
						+ kPhi * real_part ( eP ) $
						+ kz * real_part ( eZ ) 
			theta3	= aCos ( kDotE3 / ( kMag3 * eMag ) )
	
			kMag4	= sqrt ( real_part ( kR[*,3])^2 + kPhi^2 + kz^2 )
			kDotE4	= real_part ( kR[*,3] ) * real_part ( eR ) $
						+ kPhi * real_part ( eP ) $
						+ kz * real_part ( eZ ) 
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
			eMag	= sqrt ( real_part(eR)^2 + real_part(eP)^2 + real_part(eZ)^2 )
			eDotB	= real_part(eR) * runData.bField[*,0] $
						+ real_part(eP) * runData.bField[*,1] $
						+ real_part(eZ) * runData.bField[*,2]
			thetaA	= aCos ( eDotB / ( eMag * runData.bMag ) )

			eMag	= sqrt ( imaginary(eR)^2 + imaginary(eP)^2 + imaginary(eZ)^2 )
			eDotB	= imaginary(eR) * runData.bField[*,0] $
						+ imaginary(eP) * runData.bField[*,1] $
						+ imaginary(eZ) * runData.bField[*,2]
			thetaB	= aCos ( eDotB / ( eMag * runData.bMag ) )


			p=plot(runData.r, thetaA * !radeg, $
					lineStyle=6, symbol="D",yRange = [0,180], title = 'e dot B angle')
			p=plot(runData.r, thetaB * !radeg, $
					/over, lineStyle=6, symbol="S")
		endif

	endif

	e_r	= eR
	e_t	= eP
	e_z	= eZ

    h_r = hR
    h_t = hP
    h_z = hZ

    b_r = _u0 * h_r 
    b_t = _u0 * h_t 
    b_z = _u0 * h_z 

	; Calculate plasma current
	; ------------------------

    jP_r_S = complexArr(runData.nR,runData.nSpec)
    jP_t_S = complexArr(runData.nR,runData.nSpec)
    jP_z_S = complexArr(runData.nR,runData.nSpec)
    
    jP_a_S = complexArr(runData.nR,runData.nSpec)
    jP_b_S = complexArr(runData.nR,runData.nSpec)
    jP_p_S = complexArr(runData.nR,runData.nSpec)
    
    for i=0,runData.nR-1 do begin
    
            thisE = [[e_r[i]],[e_t[i]],[e_z[i]]]
			B_rtz = RunData.bField[i,*]
            Bu_rtz = B_rtz / sqrt(B_rtz[0]^2+B_rtz[1]^2+B_rtz[2]^2)
			Rot_abp_to_rtz = Get_RotMat_abp_to_rtz(Bu_rtz)
			thisE_abp = transpose(Rot_abp_to_rtz)##ThisE

            if total(transpose(Rot_abp_to_rtz)) ne total(transpose(Rot_abp_to_rtz)) then stop

            for s=0,runData.nSpec-1 do begin

                    ThisSigma = sigma[*,*,s,i]
                    ThisSigma_abp = sigma_abp[*,*,s,i]
 
                    this_jP = thisSigma ## thisE
                    this_jP_abp = thisSigma_abp ## thisE_abp
          
    				jP_r_S[i,s] = this_jP[0]
                    jP_t_S[i,s] = this_jP[1]
                    jP_z_S[i,s] = this_jP[2]

    				jP_a_S[i,s] = this_jP_abp[0]
                    jP_b_S[i,s] = this_jP_abp[1]
                    jP_p_S[i,s] = this_jP_abp[2]

            endfor

    endfor

    jP_r = Total(jP_r_S,2)
    jP_t = Total(jP_t_S,2)
    jP_z = Total(jP_z_S,2)

    jP_a = Total(jP_a_S,2)
    jP_b = Total(jP_b_S,2)
    jP_p = Total(jP_p_S,2)

	; jDotE
	; -----

    jPDotE_S = fltArr(runData.nR,runData.nSpec)
    for s=0,runData.nSpec-1 do begin
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

	endif

	; Write netCDF file

    _jP_r_S = complexArr(runData.nR,1,runData.nSpec)
    _jP_t_S = complexArr(runData.nR,1,runData.nSpec)
    _jP_z_S = complexArr(runData.nR,1,runData.nSpec)

    for s=0,runData.nSpec-1 do begin
        _jP_r_S[*,0,s] = jP_r_S[*,s]
        _jP_t_S[*,0,s] = jP_t_S[*,s]
        _jP_z_S[*,0,s] = jP_z_S[*,s]
    endfor

    TmpDensity = FltArr(RunData.nR,RunData.nSpec)
    for s=0,RunData.nSpec-1 do begin
        TmpDensity[*,s] = SpecData[s].n
    endfor

    solHash = HASH() 

    solHash = solHash + HASH('r', runData.r-cartesian_offset);
    solHash = solHash + HASH('z', runData.z);
    solHash = solHash + HASH('freq',runData.freq);

    solHash = solHash + HASH('density_m3', TmpDensity);

    solHash = solHash + HASH('B0_r',runData.BField[*,0]);
    solHash = solHash + HASH('B0_t',runData.BField[*,1]);
    solHash = solHash + HASH('B0_z',runData.BField[*,2]);

    solHash = solHash + HASH('E_r_re',real_part(e_r));
    solHash = solHash + HASH('E_r_im',imaginary(e_r));
    solHash = solHash + HASH('E_t_re',real_part(e_t));
    solHash = solHash + HASH('E_t_im',imaginary(e_t));
    solHash = solHash + HASH('E_z_re',real_part(e_z));
    solHash = solHash + HASH('E_z_im',imaginary(e_z));
    
    solHash = solHash + HASH('jP_r_re',real_part(jP_r));
    solHash = solHash + HASH('jP_r_im',imaginary(jP_r));
    solHash = solHash + HASH('jP_t_re',real_part(jP_t));
    solHash = solHash + HASH('jP_t_im',imaginary(jP_t));
    solHash = solHash + HASH('jP_z_re',real_part(jP_z));
    solHash = solHash + HASH('jP_z_im',imaginary(jP_z));

    solHash = solHash + HASH('jP_r_spec_re',real_part(_jP_r_S));
    solHash = solHash + HASH('jP_r_spec_im',imaginary(_jP_r_S));
    solHash = solHash + HASH('jP_t_spec_re',real_part(_jP_t_S));
    solHash = solHash + HASH('jP_t_spec_im',imaginary(_jP_t_S));
    solHash = solHash + HASH('jP_z_spec_re',real_part(_jP_z_S));
    solHash = solHash + HASH('jP_z_spec_im',imaginary(_jP_z_S));

    solHash = solHash + HASH('B_r_re',real_part(b_r));
    solHash = solHash + HASH('B_r_im',imaginary(b_r));
    solHash = solHash + HASH('B_t_re',real_part(b_t));
    solHash = solHash + HASH('B_t_im',imaginary(b_t));
    solHash = solHash + HASH('B_z_re',real_part(b_z));
    solHash = solHash + HASH('B_z_im',imaginary(b_z));

    solHash = solHash + HASH('jA_r_re',real_part(jA_r));
    solHash = solHash + HASH('jA_r_im',imaginary(jA_r));
    solHash = solHash + HASH('jA_t_re',real_part(jA_t));
    solHash = solHash + HASH('jA_t_im',imaginary(jA_t));
    solHash = solHash + HASH('jA_z_re',real_part(jA_z));
    solHash = solHash + HASH('jA_z_im',imaginary(jA_z));

    NCDF_PUT, 'output/rs-solution.nc', /NEW, VARIABLES=solHash


	if plotJp then begin

        nS = runData.nSpec
		jpRange = max(abs([abs(jp_r),abs(jp_t),abs(jp_z)]))

        cnt = 0
        for s=0,nS-1 do begin
            This_amu_str = ', amu: '+string(round(SpecData[s].m/_mi),format='(i1.1)')
            This_Z_str = ', Z: '+string(round(SpecData[s].q/_e),format='(i+2.1)')
            if cnt eq 0 then current = 0 else current = 1
		    p_re = plot (r,jP_r_s[*,s],thick=2,$
		    		title='jP_r'+This_amu_str+This_Z_str,name='Jp_re',font_size=10,$
		    		layout=[nS,3,1+s],yRange=yRange,transparency=50,current=current )
		    p_im = plot (r,imaginary(jP_r_s[*,s]),color='r',thick=2,transparency=50,$
		    		name='Jp_re',font_size=10,/over)

		    p_re = plot (r,jP_t_s[*,s],thick=2,$
		    		title='jP_t',name='Jp_re',font_size=10,$
		    		layout=[nS,3,1+1*nS+s],/current,yRange=yRange,transparency=50)
		    p_im = plot (r,imaginary(jP_t_s[*,s]),color='r',thick=2,transparency=50,$
		    		name='Jp_re',font_size=10,/over)

		    p_re = plot (r,jP_z_s[*,s],thick=2,$
		    		title='jP_z',name='Jp_re',font_size=10,$
		    		layout=[nS,3,1+2*nS+s],/current,yRange=yRange,transparency=50)
		    p_im = plot (r,imaginary(jP_z_s[*,s]),color='r',thick=2,transparency=50,$
		    		name='Jp_re',font_size=10,/over)
            cnt++

        endfor

        p=plot(r,jP_r,layout=[1,3,1], title='jP (summed over species)')
        p=plot(r,imaginary(jP_r),color='r',/over)
        if kjInput then begin
            p=plot(r_,kjIn.jPr_,/over, thick=2)
            p=plot(r_,imaginary(kjIn.jPr_),/over,color='r',thick=2)
        endif

        p=plot(r,jP_t,layout=[1,3,2], /current)
        p=plot(r,imaginary(jP_t),color='r',/over)
        if kjInput then begin
            p=plot(r_,kjIn.jPt_,/over, thick=2)
            p=plot(r_,imaginary(kjIn.jPt_),/over,color='r',thick=2)
        endif

        p=plot(r,jP_z,layout=[1,3,3], /current)
        p=plot(r,imaginary(jP_z),color='r',/over)
        if kjInput then begin
            p=plot(r_,kjIn.jPz_,/over, thick=2)
            p=plot(r_,imaginary(kjIn.jPz_),/over,color='r',thick=2)
        endif

	endif

    plotJp_abp = 0
	if plotJp_abp then begin

        nS = runData.nSpec
		jpRange = max(abs([abs(jp_r),abs(jp_t),abs(jp_z)]))

        p=plot([0,0],[0,0],/noData, layout=[nS,3,1],window_title='rsfwc_1d - Jp_abp',dimensions=[1200,800])

        for s=0,nS-1 do begin
            This_amu_str = ', amu: '+string(SpecData[s].m/mi,format='(i1.1)')
            This_Z_str = ', Z: '+string(SpecData[s].q/e,format='(i+2.1)')

		    p_re = plot (r,jP_a_s[*,s],thick=2,$
		    		title='jP_a'+This_amu_str+This_Z_str,name='Jp_re',font_size=10,$
		    		layout=[nS,3,1+s],yRange=yRange,transparency=50,/current)
		    p_im = plot (r,imaginary(jP_a_s[*,s]),color='r',thick=2,transparency=50,$
		    		name='Ja_re',font_size=10,/over)

		    p_re = plot (r,jP_b_s[*,s],thick=2,$
		    		title='jP_b',name='Jp_re',font_size=10,$
		    		layout=[nS,3,1+1*nS+s],/current,yRange=yRange,transparency=50)
		    p_im = plot (r,imaginary(jP_b_s[*,s]),color='r',thick=2,transparency=50,$
		    		name='Jp_re',font_size=10,/over)

		    p_re = plot (r,jP_p_s[*,s],thick=2,$
		    		title='jP_p',name='Jp_re',font_size=10,$
		    		layout=[nS,3,1+2*nS+s],/current,yRange=yRange,transparency=50)
		    p_im = plot (r,imaginary(jP_p_s[*,s]),color='r',thick=2,transparency=50,$
		    		name='Jp_re',font_size=10,/over)
        endfor

        p=plot(r,jP_a,layout=[1,3,1], title='jP_abp (summed over species)')
        p=plot(r,imaginary(jP_r),color='r',/over)
        if kjInput then begin
                p=plot(r_,kjIn.jPr_,/over,thick=2)
                p=plot(r_,imaginary(kjIn.jPr_),/over,color='r',thick=2)
        endif

        p=plot(r,jP_b,layout=[1,3,2], /current)
        p=plot(r,imaginary(jP_t),color='r',/over)

        p=plot(r,jP_p,layout=[1,3,3], /current)
        p=plot(r,imaginary(jP_z),color='r',/over)

	endif



	if plotMovie then begin

		dt	= 1.0/wReal* 0.05
		rng	= max ( abs ( [eR,eP,ez] ) )
		for i=0,10000L do begin
			t	= i * dt
			plot, runData.r, real_part ( eR ) * cos ( wReal * t ) $
					+ imaginary ( eR ) * sin ( wReal * t ), $
				yRange = [-rng*1.5,rng*1.5]
			oplot, runData.r_, real_part ( eP ) * cos ( wReal * t ) $
					+ imaginary ( eP ) * sin ( wReal * t ), $
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
stop
end
