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


	common constants
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

;	Calculation method switches

 	dielectric_freeSpace = 0
 	dielectric_noPoloidal = 0
 	dispersion_freeSpace = 0
 	dispersion_jaegerPRL = 0
	dispersion_generalised = 1
 	
 	bandStorage = 1

;	Plotting switches

 	plotRunData	= 0
 	plotDispersionGeneral 	= 1
 	plotDispersionJaeger	= 0
 	plotDispersionNoPol		= 0
	plotSolution = 0
	plotMovie	= 0
	plotFrequencies = 0
	plotKdotE = 0
	plotEdotB = 0
	plotKdotB = 0
	plotJp = 0
	plotJdotE = 0

;	File write switches

	writeDispersionTxt = 0


;	Variables
;	----------------------------

	r0	= 1.0d0
	aWall	= 0.22 

	rMin	= 0.7;0.5;1.08;r0 - aWall*0.99
	rMax	= 1.3;1.7;1.14;r0 + aWall*0.99
   	b0	= 0.53
	bR_frac	= 0.1
	bz_frac	= 0.0	
	ionSpecZ	= [ 2 ]
	ionSpecAmu	= [ 4 ]
	if not keyword_set ( nMax ) then nMax = [ 2.0 ] * 1d18
	if not keyword_set ( nLim ) then nLim = !NULL 
	if not keyword_set ( nFac ) then nFac = 1.0 
	damping = 0.06
	if not keyword_set ( freq ) then freq = 30.0e6
	if not keyword_set ( nPhi ) then nPhi = -5.0;-22.0
	if not keyword_set ( in_kz ) then kz = 0.0 else kz = in_kz;51.5;63.6 
	nR	= 512L
	AntennaJ_r = 0
	AntennaJ_t = 0
	AntennaJ_z = 1
	antLoc	= 1.6
	jAmp = 1.0
	useEqdsk = 1
	useProfiles = 1
	if not keyword_set ( poloidalScale ) then poloidalScale = 1.0
	zSlice	= 1.06 
	sliceSlope = 0.0
	profile1 = 0
	profile2 = 0
    kjInput = 0
    ar2EField = 0
    specifyE = 0

;	-----------------------------

;	Benchmarking cases

	;@smithe
	;@steffi
	;@iter
	;@langmuir
	;@klepper
    ;@ar2_vorpal_right_simple
    ;@ar2_vorpal_left_simple
    ;@ar2_vorpal_simple_full
    @colestock-kashuba
	@rsfwc_input

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

		if linearDevice then begin

			rm = 0.3
			rm0 = r0 - rm

			bPhi = b0 / (r-rm0) * (r[0]-rm0)
			bPhi_ = b0 / (r_-rm0) * (r[0]-rm0)

			bPhi = bPhi + reverse(bPhi) 
			bPhi_ = bPhi_ + reverse(bPhi_)

			bPhi = bPhi / bPhi[0] * b0
			bPhi_ = bPhi_ / bPhi[0] * b0

			bR		= dblArr ( nR ) + bPhi * bR_frac 
			bR_		= dblArr ( nR-1 ) + bPhi * bR_frac 
			bz		= dblArr ( nR ) + bPhi * bz_frac 
			bz_		= dblArr ( nR-1 ) + bPhi * bz_frac 

		endif

		if useEqdsk then begin

			eqdsk	= readgeqdsk ( eqdskFName, /noTor )
			

			bR  = interpolate ( eqdsk.bR, $
					( r - eqdsk.rleft ) / eqdsk.rdim * (eqdsk.nW-1.0), $
    		    ( z - min ( eqdsk.z ) ) / eqdsk.zdim * (eqdsk.nH-1.0), $
				cubic = -0.5 )
    		bPhi  = interpolate ( eqdsk.bPhi, $
					( r - eqdsk.rleft ) / eqdsk.rdim * (eqdsk.nW-1.0), $
    		    ( z - min ( eqdsk.z ) ) / eqdsk.zdim * (eqdsk.nH-1.0), $
				cubic = -0.5 )
    		bz  = interpolate ( eqdsk.bz, $
					( r - eqdsk.rleft ) / eqdsk.rdim * (eqdsk.nW-1.0), $
    		    ( z - min ( eqdsk.z ) ) / eqdsk.zdim * (eqdsk.nH-1.0), $
				cubic = -0.5 )

			bR_  = interpolate ( eqdsk.bR, $
					( r_ - eqdsk.rleft ) / eqdsk.rdim * (eqdsk.nW-1.0), $
    		    ( z_ - min ( eqdsk.z ) ) / eqdsk.zdim * (eqdsk.nH-1.0), $
				cubic = -0.5 )
    		bPhi_  = interpolate ( eqdsk.bPhi, $
					( r_ - eqdsk.rleft ) / eqdsk.rdim * (eqdsk.nW-1.0), $
    		    ( z_ - min ( eqdsk.z ) ) / eqdsk.zdim * (eqdsk.nH-1.0), $
				cubic = -0.5 )
    		bz_  = interpolate ( eqdsk.bz, $
					( r_ - eqdsk.rleft ) / eqdsk.rdim * (eqdsk.nW-1.0), $
    		    ( z_ - min ( eqdsk.z ) ) / eqdsk.zdim * (eqdsk.nH-1.0), $
				cubic = -0.5 )

		endif 

        if ar2Input gt 0 then begin 

            ar2_read_ar2Input, ar2InFileName, ar2=ar2

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
    
        endif

		if poloidalScale ne 1 then begin

			print, 'Scaling poloidal field components'
			bR	= bR * poloidalScale
			bR_	= bR_ * poloidalScale 
			bz	= bz * poloidalScale
			bz_	= bz_ * poloidalScale

		endif

		bMag	= sqrt ( bR^2 + bPhi^2 + bz^2 )
		bMag_	= sqrt ( bR_^2 + bPhi_^2 + bz_^2 )

        if ar2Input then begin
            ionSpecZ = ar2.AtomicZ[1:-1]
            ionSpecAmu = ar2.amu[1:-1]
        endif

		nSpec	= n_elements ( ionSpecZ )

		specData	= replicate ( $
				{ 	q : 0d0, $
					m : 0d0, $
					wp : dblArr ( nR ), $
					wp_ : dblArr ( nR - 1 ), $
					wc : dblArr ( nR ), $
					wc_ : dblArr ( nR - 1 ), $
					n : dblArr ( nR ), $
					n_ : dblArr ( nR - 1 ), $
					wLH : dblArr ( nR ), $
					wUH : dblArr ( nR ), $
                    nuOmg: dblArr(nR), $
                    nuOmg_ : dblArr(nR) }, nSpec + 1 )

        if not ar2Input then begin
		for i=0,nSpec do begin

			nProfile	= r * 0 + 1
			nProfile_	= r_* 0 + 1

			if profile1 and i lt nSpec then begin

				;nProfile	= 1.0 - ( abs ( r - r0 ) / aWall )^2 
				;nProfile_	= 1.0 - ( abs ( r_ - r0 ) / aWall )^2

				if nLim ne !NULL then scaleFac = (nMax[i]-nLim[i])/aWall^2
    			if nLim eq !NULL then scaleFac = (nMax[i])/aWall^2
                    
                nProfile  = ( nMax[i]-scaleFac*(r-r0)^2 )  / nMax[i]
                nProfile_  = ( nMax[i]-scaleFac*(r_-r0)^2 )  / nMax[i]

			endif

			if profile2 and i lt nSpec then begin

				nProfile = (exp(-70*((r-max(r)))))
				nProfile_ = (exp(-70*((r_-max(r)))))

			endif

			; ions
			if i lt nSpec then begin

				specData[i].q 		= ionSpecZ[i] * e
				specData[i].m 		= ionSpecAmu[i] * mi
				specData[i].n 		= nMax[i] * nProfile
				specData[i].n_		= nMax[i] * nProfile_

			; electrons for charge neutrality
			endif else begin

				specData[i].q 		= -e
				specData[i].m 		= me
				if nSpec gt 1 then begin
					specData[i].n 		= total ( specData[0:nSpec-1].n $
											* transpose(rebin(ionSpecZ,nSpec,nR)), 2 ) 
					specData[i].n_		= total ( specData[0:nSpec-1].n_ $
											* transpose(rebin(ionSpecZ,nSpec,nR-1)), 2 ) 
				endif else begin
					specData[i].n 		= specData[0:nSpec-1].n $
											* transpose(rebin(ionSpecZ,nSpec,nR)) 
					specData[i].n_		= specData[0:nSpec-1].n_ $
											* transpose(rebin(ionSpecZ,nSpec,nR-1)) 
				endelse

			endelse

		endFor
        endif

		;	import density profiles 

		if useProfiles then begin

			nstx_profile	= '../profiles/dlg_profiles_x1.00.nc' 
			cdfId = ncdf_open ( nstx_profile, /noWrite ) 
    	    ncdf_varGet, cdfId, 'R_binCenters', xMap_R
    	    ncdf_varGet, cdfId, 'z_binCenters', xMap_z
    	    ncdf_varGet, cdfId, 'ne', xMap_ne
    	    ncdf_varGet, cdfId, 'te', xMap_te
    	    ncdf_varGet, cdfId, 'ti', xMap_ti
    	    ncdf_varGet, cdfId, 'bMag', xMap_B
    	    ncdf_varGet, cdfId, 'r2d', xMap_R2D
    	    ncdf_varGet, cdfId, 'z2d', xMap_z2D
    	    nCdf_close, cdfId
    	         
    	    nX  = n_elements ( xMap_B[*,0] )
    	    nY  = n_elements ( xMap_B[0,*] )

			ne_  = interpolate ( xMap_ne, $
					( r - min(xMap_R2d) ) / (max(xMap_R2d)-min(xMap_R2d)) * (nX-1.0), $
    		    ( z - min ( xMap_z2d ) ) / (max(xMap_z2d)-min(xMap_z2d)) * (nY-1.0), $
				cubic = -0.5 )
    		ne__  = interpolate ( xMap_ne, $
					( r_ - min(xMap_R2d) ) / (max(xMap_R2d)-min(xMap_R2d)) * (nX-1.0), $
    		    ( z_ - min ( xMap_z2d ) ) / (max(xMap_z2d)-min(xMap_z2d)) * (nY-1.0), $
				cubic = -0.5 )

			if keyword_set ( nFac ) then begin
				ne_	= ne_ * nFac 
				ne__	= ne__ * nFac
			endif

			specData[1].n	= ne_
			specData[1].n_	= ne__
    
			specData[0].n	= ne_ / ionSpecZ[0]
			specData[0].n_	= ne__ / ionSpecZ[0]

		endif
        
        if ar2Input then begin

            for s=0,nSpec do begin

                index = s+1
                if s eq nSpec then index = 0

			    this_n  = interpolate ( ar2.Density_m3[*,*,index], $
			    		( r - ar2.rMin ) / (ar2.rMax-ar2.rMin) * (ar2_nR-1), $
    		        ( z - ar2.zMin ) / (ar2.zMax-ar2.zMin) * (ar2_nZ-1), $
			    	cubic = -0.5 )
			    this_n_  = interpolate ( ar2.Density_m3[*,*,index], $
			    		( r_ - ar2.rMin ) / (ar2.rMax-ar2.rMin) * (ar2_nR-1), $
    		        ( z_ - ar2.zMin ) / (ar2.zMax-ar2.zMin) * (ar2_nZ-1), $
			    	cubic = -0.5 )

			    this_nuOmg  = interpolate ( ar2.nuOmg[*,*,index], $
			    		( r - ar2.rMin ) / (ar2.rMax-ar2.rMin) * (ar2_nR-1), $
    		        ( z - ar2.zMin ) / (ar2.zMax-ar2.zMin) * (ar2_nZ-1), $
			    	cubic = -0.5 )
			    this_nuOmg_  = interpolate ( ar2.nuOmg[*,*,index], $
			    		( r_ - ar2.rMin ) / (ar2.rMax-ar2.rMin) * (ar2_nR-1), $
    		        ( z_ - ar2.zMin ) / (ar2.zMax-ar2.zMin) * (ar2_nZ-1), $
			    	cubic = -0.5 )

                specData[s].n = this_n
                specData[s].n_ = this_n_

                specData[s].nuOmg = this_nuOmg
                specData[s].nuOmg_ = this_nuOmg_
 
                specData[s].m = ar2.amu[index]*mi
                if s eq nSpec then specData[s].m = me
                specData[s].q = ar2.AtomicZ[index]*e
                ;if s eq nSpec then specData[s].q = -e

            endfor 
        endif


		for i = 0, nSpec do begin

			specData[i].wp	= sqrt ( specData[i].n * specData[i].q^2 $
    	                        / ( specData[i].m*e0 ))
			specData[i].wc	= specData[i].q * bMag / specData[i].m
			specData[i].wp_	= sqrt ( specData[i].n_ * specData[i].q^2 $
    	                        / ( specData[i].m*e0 ))
			specData[i].wc_	= specData[i].q * bMag_ / specData[i].m

		
		endfor

		for i = 0, nSpec - 1 do begin

			specData[i].wLH	= sqrt ( specData[i].wc^2 $
				+ specData[i].wp^2 / ( 1.0 + specData[nSpec].wp^2 / specData[nSpec].wc^2 ) )
			;specData[i].wLH	= sqrt ( specData[i].wc * abs(specData[nSpec].wc) $
			;	* ( specData[nSpec].wp^2 + specData[nSpec].wc * specData[i].wc) $
			;	/ ( specData[nSpec].wp^2 + specData[nSpec].wc^2 ) )
			;specData[i].wLH	= ( 1d0 / (specData[i].wc * abs(specData[nSpec].wc)) $
			;					+ 1d0/(specData[i].wp^2) )^(-0.5)
			specData[i].wUH	= sqrt ( specData[nSpec].wp^2 + specData[nSpec].wc^2 ) 
	
		endfor

		if plotFrequencies then begin
			iplot, abs(specData[1].wp), /ylog, thick = 3, insert_legend = 1, name = 'wpe' 
			iplot, abs(specData[1].wc), /over, thick = 3, lineStyle = 2, insert_legend = 1, name = 'wce'
			iplot, abs(specData[0].wp), /over, insert_legend = 1, name = 'wpi'
			iplot, abs(specData[0].wc), /over, lineStyle = 2, insert_legend = 1, name = 'wci'
		endif

        
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
					antLoc : antLoc }



    	if plotRunData then begin 
		    p=Plot( r, runData.bField[*,1], $
    	        layout = [2,1,1], $
    	        title = 'bPhi')
		    p=plot( r, runData.bField[*,0], $
				/over, $
				color = red)
		    p=plot( r, runData.bField[*,2], $
				/over, $
				color = blue)
		    p=plot( r, specData[nSpec].n, $
    	        layout=[2,1,2],/current, $
    	        title = 'density 0', yRange=[1e10,1e22],/ylog)
			for i=0,nSpec-1 do $
		    p=plot( r, specData[i].n, $
    	        /over )
    	endif


end
