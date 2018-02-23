pro run_setup, $
	runData = runData, $
	specData = specData, $
	in_kz = kz, $
	nMax = nMax, $
	nFac = nfac, $
	gradSize = gradSize, $
	freq = freq, $
	poloidalScale = poloidalScale, $
	nPhi = nPhi, $
	AntennaJ_r = AntennaJ_r, $
	AntennaJ_t = AntennaJ_t, $
	AntennaJ_z = AntennaJ_z, $
	antSig_r = antSig_r, $
	antSig_t = antSig_t, $
	antSig_z = antSig_z, $
	jAmp = jAmp


	@dlg_constants
	common switches
	common plotSwitches
	common writeSwitches

	common dlg_colors, ct12, $
		red, $
		blue, $
		green, $
		purple

	loadct, 12, /sil, rgb_table = ct12
	red	= transpose ( ct12[12*16-1,*] )
	blue	= transpose ( ct12[8*16-1,*] )
	green	= transpose ( ct12[2*16-1,*] )
	purple	= transpose ( ct12[9*16-1,*] )


;	Defaults
;	--------

	r0	= 1.0d0
	aWall	= 0.22 

	rMin	= 0.7
	rMax	= 1.3
   	b0	= 0.53
	bR_frac	= 0.1
	bz_frac	= 0.0	
	specZ	= [ 2 ]
	specAmu	= [ 4 ]
	if not keyword_set ( nMax ) then nMax = [ 2.0 ] * 1d18
	if not keyword_set ( nLim ) then nLim = !NULL 
	if not keyword_set ( nFac ) then nFac = 1.0 
	damping = 0.06
	if not keyword_set ( freq ) then freq = 30.0e6
	if not keyword_set ( nPhi ) then nPhi = 5.0
	if not keyword_set ( in_kz ) then kz = 0.0 else kz = in_kz 
	nR	= 512L
	AntennaJ_r = 0
	AntennaJ_t = 0
	AntennaJ_z = 1
	antLoc	= 1.6
	jAmp = 1.0
	antSig_r = 0.02 
	antSig_t = 0.02  
	antSig_z = 0.02 
	useEqdsk = 0
	useProfiles = 0
	if not keyword_set ( poloidalScale ) then poloidalScale = 1.0
	zSlice	= 0.00 
	sliceSlope = 0.0
	profile1 = 0
	profile2 = 0
    kjInput = 0
    ar2EField = 0
    specifyE = 0
    ar2Input = 1

;	-----------------------------

;	Read input file

	@rs_input 

    if ar2Input gt 0 then begin 
        ar2 = ar2_read_ar2Input('./', fileName = ar2InFileName)
        rMin = ar2.rMin
        rMax = ar2.rMax
    endif

    print, 'rMin: ', rMin
    print, 'rMax: ', rMax
    print, 'nPhi: ', nPhi
    print, 'kz: ', kz

;		Grid

		dR	= ( rMax - rMin ) / ( nR - 1 )
		r	= dIndGen ( nR ) * dR + rMin
		;	half grid variables are _
		r_	= r[1:*] - dR / 2

		z	= r * sliceSlope + zSlice
		z_	= r_ * sliceSlope + zSlice

;		Setup plasma and field profiles

		bPhi	= dblArr(nR) + b0 / r * r0
		bPhi_	= dblArr(nR-1) + b0 / r_ * r0

		bR		= dblArr ( nR ) + bPhi * bR_frac 
		bR_		= dblArr ( nR-1 ) + bPhi * bR_frac 
		bz		= dblArr ( nR ) + bPhi * bz_frac 
		bz_		= dblArr ( nR-1 ) + bPhi * bz_frac 

        if ar2Input gt 0 then begin 

            ar2_nR = n_elements(ar2.r)
            ar2_nz = n_elements(ar2.z)

			bR  = interpolate ( ar2.br, $
					( r - ar2.rMin ) / (ar2.rMax-ar2.rMin) * (ar2_nR-1), $
    		    ( z - ar2.zMin ) / (ar2.zMax-ar2.zMin) * (ar2_nZ-1), $
				cubic = -0.5 )
    		bPhi  = interpolate ( ar2.bt, $
					( r - ar2.rMin ) / (ar2.rMax-ar2.rMin) * (ar2_nR-1), $
    		    ( z - ar2.zMin ) / (ar2.zMax-ar2.zMin) * (ar2_nZ-1), $
				cubic = -0.5 )
     		bZ  = interpolate ( ar2.bz, $
					( r - ar2.rMin ) / (ar2.rMax-ar2.rMin) * (ar2_nR-1), $
    		    ( z - ar2.zMin ) / (ar2.zMax-ar2.zMin) * (ar2_nZ-1), $
				cubic = -0.5 )

    		bR_  = interpolate ( ar2.br, $
					( r_ - ar2.rMin ) / (ar2.rMax-ar2.rMin) * (ar2_nR-1), $
    		    ( z_ - ar2.zMin ) / (ar2.zMax-ar2.zMin) * (ar2_nZ-1), $
				cubic = -0.5 )
    		bPhi_  = interpolate ( ar2.bt, $
					( r_ - ar2.rMin ) / (ar2.rMax-ar2.rMin) * (ar2_nR-1), $
    		    ( z_ - ar2.zMin ) / (ar2.zMax-ar2.zMin) * (ar2_nZ-1), $
				cubic = -0.5 )
     		bZ_  = interpolate ( ar2.bz, $
					( r_ - ar2.rMin ) / (ar2.rMax-ar2.rMin) * (ar2_nR-1), $
    		    ( z_ - ar2.zMin ) / (ar2.zMax-ar2.zMin) * (ar2_nZ-1), $
				cubic = -0.5 )
   
        endif else begin
            print, 'ERROR: You must have an ar2Input.nc file now'                    
            stop
        endelse


		bMag	= sqrt ( bR^2 + bPhi^2 + bz^2 )
		bMag_	= sqrt ( bR_^2 + bPhi_^2 + bz_^2 )

        if ar2Input then begin
            specZ = ar2.AtomicZ
            specAmu = ar2.amu
        endif

		nSpec	= n_elements ( specZ )

		specData	= replicate ( $
				{ 	q : 0d0, $
					m : dComplexArr(nR), $
                    m_ : dComplexArr(nR-1), $
                    amu : 0d0, $
                    Z : 0d0, $
					wp : dComplexarr ( nr ), $
					wp_ : dComplexArr ( nR - 1 ), $
					wc : dComplexArr ( nR ), $
					wc_ : dComplexArr ( nR - 1 ), $
					n : dblArr ( nR ), $
					n_ : dblArr ( nR - 1 ), $
					wLH : dblArr ( nR ), $
					wUH : dblArr ( nR ), $
                    nuOmg: dblArr(nR), $
                    nuOmg_ : dblArr(nR) }, nSpec )

    
        if ar2Input then begin

            for s=0,nSpec-1 do begin

			    this_n  = interpolate ( ar2.Density_m3[*,*,s], $
			    		( r - ar2.rMin ) / (ar2.rMax-ar2.rMin) * (ar2_nR-1), $
    		        ( z - ar2.zMin ) / (ar2.zMax-ar2.zMin) * (ar2_nZ-1), $
			    	cubic = -0.5 )
			    this_n_  = interpolate ( ar2.Density_m3[*,*,s], $
			    		( r_ - ar2.rMin ) / (ar2.rMax-ar2.rMin) * (ar2_nR-1), $
    		        ( z_ - ar2.zMin ) / (ar2.zMax-ar2.zMin) * (ar2_nZ-1), $
			    	cubic = -0.5 )

			    this_nuOmg  = interpolate ( ar2.nuOmg[*,*,s], $
			    		( r - ar2.rMin ) / (ar2.rMax-ar2.rMin) * (ar2_nR-1), $
    		        ( z - ar2.zMin ) / (ar2.zMax-ar2.zMin) * (ar2_nZ-1), $
			    	cubic = -0.5 )
			    this_nuOmg_  = interpolate ( ar2.nuOmg[*,*,s], $
			    		( r_ - ar2.rMin ) / (ar2.rMax-ar2.rMin) * (ar2_nR-1), $
    		        ( z_ - ar2.zMin ) / (ar2.zMax-ar2.zMin) * (ar2_nZ-1), $
			    	cubic = -0.5 )

                specData[s].n = this_n
                specData[s].n_ = this_n_

                specData[s].nuOmg = this_nuOmg
                specData[s].nuOmg_ = this_nuOmg_
 
                specData[s].m = ar2.amu[s]*_amu * (1 + _ii*this_nuOmg)
                specData[s].m_ = ar2.amu[s]*_amu * (1 + _ii*this_nuOmg_)
                specData[s].amu = ar2.amu[s]
                specData[s].q = ar2.AtomicZ[s]*_e
                specData[s].Z = ar2.AtomicZ[s]

            endfor 
        endif


		for s = 0, nSpec-1 do begin

			specData[s].wp	= sqrt ( specData[s].n * specData[s].q^2 $
    	                        / ( specData[s].m*_e0 ))
			specData[s].wc	= specData[s].q * bMag / specData[s].m
			specData[s].wp_	= sqrt ( specData[s].n_ * specData[s].q^2 $
    	                        / ( specData[s].m_*_e0 ))
			specData[s].wc_	= specData[s].q * bMag_ / specData[s].m_

		endfor

		runData = { bField : [ [bR], [bPhi], [bz] ], $
					bField_ : [ [bR_], [bPhi_], [bz_] ], $
					bMag : bMag, $
					bMag_ : bMag_, $
					nR : nR, $
                    nR_ : nR-1, $
					r : r, $
					r_ : r_, $
					z : z, $
					z_ : z_, $
					dR : dR, $
					b0 : b0, $
					r0 : r0, $
					rMin : rMin, $
					rMax : rMax, $
					aWall : aWall, $
					damping : damping, $
					freq : freq, $
					nPhi : nPhi, $
					kz : kz, $
					nSpec : nSpec, $
					antLoc : antLoc, $
                    antSig_r : antSig_r, $
                    antSig_t : antSig_t, $
                    antSig_z : antSig_z, $
                    jAmp : jAmp $
                    }

    	if plotRunData then begin 
		    p=Plot( r, runData.bField[*,1], $
    	        layout = [2,1,1], $
    	        title = 'b Field Components')
		    p=plot( r, runData.bField[*,0], $
				/over, $
				color = red)
		    p=plot( r, runData.bField[*,2], $
				/over, $
				color = blue)
		    p=plot( r, specData[0].n, $
    	        layout=[2,1,2],/current, $
    	        title = 'density 0', yRange=[1e10,1e22],/ylog)

			for s=1,nSpec-1 do p=plot( r, specData[s].n, /over )

            p=plot(r,specData[0].nuOmg, title='nuOmg')
            for s=0,nSpec-1 do p=plot(r,specData[s].nuOmg, /over)
    	endif


end
