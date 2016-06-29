;	Real Space Full Wave Code in 1D
;	DLG

pro rsfwc_1d, $
    eR = eR, $
	eP = eP, $
	ez = ez, $
	divD = divD, $
	r = r, r_ = r_, $
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
	nPhi = nPhi, $
    plot_Jp_spec = plot_Jp_spec

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


	r 	= runData.r
	r_	= runData.r_
    nR = n_elements(r)
    nR_ = n_elements(r_)
	rOut	= r

 	wReal	= runData.freq * 2d0 * !dpi 
	wc	    = dcomplexArr ( runData.nR, runData.nIonSpec+1 ) 
	wc_	    = dcomplexArr ( runData.nR_, runData.nIonSpec+1 ) 

    for s=0,runData.nIonSpec do begin
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
        sigmaFull_ = sigma, $
        sigmaHalf_ = sigma_, $
		sigmaFullSpec_abp = sigma_abp, $
		sigmaHalfSpec_abp = sigma_abp_

	dispersionAll, wReal, epsilon, stixVars, runData, specData, $
		kR = kR, kPhi = kPhi, kz = kz

	if keyword_set ( dispersionOnly ) then return

	if kjInput then begin
        kjIn = kj_read_kj_jP(kj_jP_fileName, r, rH=r_)
        runData.nIonSpec = n_elements(kjIn.jPr_spec[0,*])-1
    endif else if ar2EField then begin
        replace = intArr(runData.nR)*0+1
        replace_ = intArr(runData.nR-1)*0+1
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
	rhs		= dcomplexArr ( nAll )

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

   	for i=0,runData.nR-1 do begin

		rhs[i*3]	= II * wReal * u0 * jA_r[i]
        if i lt nR-1 then begin
		rhs[i*3+1]	= II * wReal * u0 * jA_t_[i]
		rhs[i*3+2]	= II * wReal * u0 * jA_z_[i]
        endif

		if kjInput then begin
			rhs[i*3]	+= II * wReal * u0 * kjIn.jpR[i]
            if i lt nR-1 then begin
			rhs[i*3+1]	+= II * wReal * u0 * kjIn.jpT_[i]
			rhs[i*3+2]	+= II * wReal * u0 * kjIn.jpZ_[i]
            endif
		endif

	endfor

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

    eP_LeftBoundary_ = eField[0]
    eZ_LeftBoundary_ = eField[1]
    eP_RightBoundary_ = eField[-2]
    eZ_RightBoundary_ = eField[-1]

    ; Remove BC layers
    rhs = rhs[2:-3]
    eField = eField[2:-3]
  
	ii_eR	= lIndGen(runData.nR)*3
	eR	= eField[ii_eR]
	ii_eP	= temporary(ii_eR[0:runData.nR-2]+1)
	eP_	= eField[ii_eP]
	ii_ez	= temporary(ii_eP+1)
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
			hR_[i]	= -II * runData.kz * eP_[i] + II * runData.nPhi * ez_[i] / runData.r_[i]

		if i gt 0 and i lt runData.nR-1 then begin
		
			hP[i]	= II * runData.kz * eR[i] - ( ez_[i] - ez_[i-1] ) / runData.dR 
			hz[i]	= ( -II * runData.nPhi * eR[i] $
						+ ( runData.r_[i]*eP_[i] $
							- runData.r_[i-1]*eP_[i-1] ) / runData.dR ) $
					 	/ runData.r[i]
	
		endif

	endfor

	hR_	= hR_ / ( II * w[*,0] * u0 )
	hP	= hP / ( II * w[*,0] * u0 )
	hz	= hz / ( II * w[*,0] * u0 )

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

    b_r = u0 * h_r 
    b_t = u0 * h_t 
    b_z = u0 * h_z 

	; Calculate plasma current
	; ------------------------

    if kjInput then begin
    
    	jP_r	= kjIn.jpR
    	jP_t	= kjIn.jpT
    	jP_z	= kjIn.jpZ

    	jP_r_S	= kjIn.jpR_spec
    	jP_t_S	= kjIn.jpT_spec
    	jP_z_S	= kjIn.jpZ_spec

    endif else begin

        jP_r_S = complexArr(runData.nR,runData.nIonSpec+1)
        jP_t_S = complexArr(runData.nR,runData.nIonSpec+1)
        jP_z_S = complexArr(runData.nR,runData.nIonSpec+1)
        
        jP_a_S = complexArr(runData.nR,runData.nIonSpec+1)
        jP_b_S = complexArr(runData.nR,runData.nIonSpec+1)
        jP_p_S = complexArr(runData.nR,runData.nIonSpec+1)
        
        for i=0,runData.nR-1 do begin
        
                thisE = [[e_r[i]],[e_t[i]],[e_z[i]]]
				B_rtz = RunData.bField[i,*]
                Bu_rtz = B_rtz / sqrt(B_rtz[0]^2+B_rtz[1]^2+B_rtz[2]^2)
				Rot_abp_to_rtz = Get_RotMat_abp_to_rtz(Bu_rtz)
				thisE_abp = transpose(Rot_abp_to_rtz)##ThisE

                if total(transpose(Rot_abp_to_rtz)) ne total(transpose(Rot_abp_to_rtz)) then stop

                for s=0,runData.nIonSpec do begin

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

    endelse

	; jDotE
	; -----

    jPDotE_S = fltArr(runData.nR,runData.nIonSpec+1)
    for s=0,runData.nIonSpec do begin
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

	nc_id = nCdf_create ('rsfwc_1d.nc', /clobber )

	nCdf_control, nc_id, /fill
	
	nr_id = nCdf_dimDef ( nc_id, 'nR', runData.nR )
	nrH_id = nCdf_dimDef ( nc_id, 'nR_', runData.nR-1 )
	scalar_id = nCdf_dimDef ( nc_id, 'scalar', 1 )
	nIonSpec_id = nCdf_dimDef ( nc_id, 'nSpec', runData.nIonSpec+1 )

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
    
	b_r_re_id = nCdf_varDef ( nc_id, 'b_r_re', nr_id, /float )
	b_r_im_id = nCdf_varDef ( nc_id, 'b_r_im', nr_id, /float )
	b_p_re_id = nCdf_varDef ( nc_id, 'b_p_re', nr_id, /float )
	b_p_im_id = nCdf_varDef ( nc_id, 'b_p_im', nr_id, /float )
	b_z_re_id = nCdf_varDef ( nc_id, 'b_z_re', nr_id, /float )
	b_z_im_id = nCdf_varDef ( nc_id, 'b_z_im', nr_id, /float )

	jP_r_re_id = nCdf_varDef ( nc_id, 'jP_r_re', nr_id, /float )
	jP_r_im_id = nCdf_varDef ( nc_id, 'jP_r_im', nr_id, /float )
	jP_p_re_id = nCdf_varDef ( nc_id, 'jP_p_re', nr_id, /float )
	jP_p_im_id = nCdf_varDef ( nc_id, 'jP_p_im', nr_id, /float )
	jP_z_re_id = nCdf_varDef ( nc_id, 'jP_z_re', nr_id, /float )
	jP_z_im_id = nCdf_varDef ( nc_id, 'jP_z_im', nr_id, /float )

	jP_r_re_spec_id = nCdf_varDef ( nc_id, 'jP_r_re_spec', [nr_id,nIonSpec_id], /float )
	jP_r_im_spec_id = nCdf_varDef ( nc_id, 'jP_r_im_spec', [nr_id,nIonSpec_id], /float )
	jP_p_re_spec_id = nCdf_varDef ( nc_id, 'jP_p_re_spec', [nr_id,nIonSpec_id], /float )
	jP_p_im_spec_id = nCdf_varDef ( nc_id, 'jP_p_im_spec', [nr_id,nIonSpec_id], /float )
	jP_z_re_spec_id = nCdf_varDef ( nc_id, 'jP_z_re_spec', [nr_id,nIonSpec_id], /float )
	jP_z_im_spec_id = nCdf_varDef ( nc_id, 'jP_z_im_spec', [nr_id,nIonSpec_id], /float )

	jA_r_re_id = nCdf_varDef ( nc_id, 'jA_r_re', nr_id, /float )
	jA_r_im_id = nCdf_varDef ( nc_id, 'jA_r_im', nr_id, /float )
	jA_p_re_id = nCdf_varDef ( nc_id, 'jA_p_re', nr_id, /float )
	jA_p_im_id = nCdf_varDef ( nc_id, 'jA_p_im', nr_id, /float )
	jA_z_re_id = nCdf_varDef ( nc_id, 'jA_z_re', nr_id, /float )
	jA_z_im_id = nCdf_varDef ( nc_id, 'jA_z_im', nr_id, /float )

	Density_id = nCdf_varDef ( nc_id, 'density_m3', [nr_id,nIonSpec_id], /float )

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

	nCdf_varPut, nc_id, b_r_re_id, real_part(b_r) 
	nCdf_varPut, nc_id, b_r_im_id, imaginary(b_r) 
	nCdf_varPut, nc_id, b_p_re_id, real_part(b_t) 
	nCdf_varPut, nc_id, b_p_im_id, imaginary(b_t) 
	nCdf_varPut, nc_id, b_z_re_id, real_part(b_z) 
	nCdf_varPut, nc_id, b_z_im_id, imaginary(b_z) 

	nCdf_varPut, nc_id, jP_r_re_id, real_part(jP_r)
	nCdf_varPut, nc_id, jP_r_im_id, imaginary(jP_r) 
	nCdf_varPut, nc_id, jP_p_re_id, real_part(jP_t) 
	nCdf_varPut, nc_id, jP_p_im_id, imaginary(jP_t) 
	nCdf_varPut, nc_id, jP_z_re_id, real_part(jP_z) 
	nCdf_varPut, nc_id, jP_z_im_id, imaginary(jP_z) 

	nCdf_varPut, nc_id, jP_r_re_spec_id, real_part(jP_r_S)
	nCdf_varPut, nc_id, jP_r_im_spec_id, imaginary(jP_r_S) 
	nCdf_varPut, nc_id, jP_p_re_spec_id, real_part(jP_t_S)
	nCdf_varPut, nc_id, jP_p_im_spec_id, imaginary(jP_t_S) 
	nCdf_varPut, nc_id, jP_z_re_spec_id, real_part(jP_z_S)
	nCdf_varPut, nc_id, jP_z_im_spec_id, imaginary(jP_z_S) 

	nCdf_varPut, nc_id, jA_r_re_id, jA_r 
	nCdf_varPut, nc_id, jA_r_im_id, jA_r*0 
	nCdf_varPut, nc_id, jA_p_re_id, jA_t
	nCdf_varPut, nc_id, jA_p_im_id, jA_t*0
	nCdf_varPut, nc_id, jA_z_re_id, jA_z
	nCdf_varPut, nc_id, jA_z_im_id, jA_z*0

    TmpDensity = FltArr(RunData.nR,RunData.nIonSpec+1)
    for s=0,RunData.nIonSpec do begin
        TmpDensity[*,s] = SpecData[s].n
    endfor

    nCdf_varPut, nc_id, Density_id, TmpDensity

	nCdf_close, nc_id


	if plotJp then begin

        nS = runData.nIonSpec+1
		jpRange = max(abs([abs(jp_r),abs(jp_t),abs(jp_z)]))
		;yRange = [-1,1];[-jPRange,jPRange]*0.5

        p=plot([0,0],[0,0],/noData, layout=[nS,3,1],window_title='rsfwc_1d',dimensions=[1200,800])

        for s=0,nS-1 do begin
            This_amu_str = ', amu: '+string(SpecData[s].m/mi,format='(i1.1)')
            This_Z_str = ', Z: '+string(SpecData[s].q/e,format='(i+2.1)')

		    p_re = plot (r,jP_r_s[*,s],thick=2,$
		    		title='jP_r'+This_amu_str+This_Z_str,name='Jp_re',font_size=10,$
		    		layout=[nS,3,1+s],yRange=yRange,transparency=50,/current)
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
        endfor

        p=plot(r,jP_r,layout=[1,3,1], title='jP (summed over species)')
        p=plot(r,imaginary(jP_r),color='r',/over)
        if kjInput then begin
                p=plot(r_,kjIn.jPr_,/over)
                p=plot(r_,imaginary(kjIn.jPr_),/over,color='r')
        endif

        p=plot(r,jP_t,layout=[1,3,2], /current)
        p=plot(r,imaginary(jP_t),color='r',/over)

        p=plot(r,jP_z,layout=[1,3,3], /current)
        p=plot(r,imaginary(jP_z),color='r',/over)

	endif

    plotJp_abp = 1
	if plotJp_abp then begin

        nS = runData.nIonSpec+1
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
                p=plot(r_,kjIn.jPr_,/over)
                p=plot(r_,imaginary(kjIn.jPr_),/over,color='r')
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
