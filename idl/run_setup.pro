pro run_setup, $
	runData = runData, $
	specData = specData, $
	kz = kz, $
	nMax = nMax

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

 	plotRunData	= 1
 	plotDispersionGeneral 	= 1
 	plotDispersionJaeger	= 0
 	plotDispersionNoPol		= 0
	plotSolution = 1

;	File write switches

	writeDispersionTxt = 0

;		----------------------------
;	Variables

	r0	= 0.67d0
	aWall	= 0.22 

	rMin	= 0.2;r0 - aWall*0.99
	rMax	= 1.7;r0 + aWall*0.99
   	b0	= 5.85d0
	bR_frac	= 0.0
	bz_frac	= 0.0	
	ionSpecZ	= [ 2 ]
	ionSpecAmu	= [ 2 ]
	if not keyword_set ( nMax ) then nMax = [ 5.5 ] * 1d18
	damping = 0.04
	freq	= 30.0e6
	nPhi = -22.0
	if not keyword_set ( kz ) then kz = 98.5
	nR	= 1024L
	antLoc	= 1.2

	useEqdsk = 1
	useProfiles = 0
	poloidalScale = 1.0
	zSlice	= 1.0 
	profile1 = 0

;	-----------------------------



;	Benchmarking cases

;	@smithe


;		Grid

		dR	= ( rMax - rMin ) / ( nR - 1 )
		r	= dIndGen ( nR ) * dR + rMin
		;	half grid variables are _
		r_	= r[1:*] - dR / 2

		z	= r * 0 + zSlice
		z_	= r_ * 0 + zSlice

;		Setup plasma and field profiles

		bPhi	= dblArr(nR) + b0 / r * r0
		bPhi_	= dblArr(nR-1) + b0 / r_ * r0

		bR		= dblArr ( nR ) + bPhi * bR_frac 
		bR_		= dblArr ( nR-1 ) + bPhi * bR_frac 
		bz		= dblArr ( nR ) + bPhi * bz_frac 
		bz_		= dblArr ( nR-1 ) + bPhi * bz_frac 

		if useEqdsk then begin

			nstx_eqdsk	= '../eqdsk/g120740.00275.EFIT02.mds.uncorrected.qscale_1.00000.dlgMod_1.67'
			eqdsk	= readgeqdsk ( nstx_eqdsk )
			
		
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

			;for i=0,5000 do begin
			;	bPhi	= smooth ( bPhi, nR / 10, /edge ) 
			;	bPhi_	= smooth ( bPhi_, nR / 10, /edge ) 
			;endfor
			;bPhi	= bPhi * 0.5
			;bPhi_	= bPhi_ * 0.5
	
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

		nSpec	= n_elements ( ionSpecZ )

		specData	= replicate ( $
				{ 	q : 0d0, $
					m : 0d0, $
					wp : dblArr ( nR ), $
					wp_ : dblArr ( nR - 1 ), $
					wc : dblArr ( nR ), $
					wc_ : dblArr ( nR - 1 ), $
					n : dblArr ( nR ), $
					n_ : dblArr ( nR - 1 ) }, nSpec + 1 )

		for i=0,nSpec do begin

			nProfile	= r * 0 + 1
			nProfile_	= r_* 0 + 1

			if profile1 then begin

				nProfile	= 1.0 - ( abs ( r - r0 ) / aWall )^2 
				nProfile_	= 1.0 - ( abs ( r_ - r0 ) / aWall )^2
	
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

		;;	sinusoidal density gradient

		;gradSize = 0.3
		;gradStart	= 0.85
		;gradFreq	= 1 / (2*gradSize)
		;ne_profile	= fltArr ( nR ) 
		;ii1	= where ( r lt gradStart )
		;ii2	= where ( r ge gradStart and r lt gradStart+gradSize )
		;ne_profile[*]	= -1
		;ne_profile[ii1]	= 1
		;ne_profile[ii2]	= cos ( 2*!pi*gradFreq*(r[ii2]-gradStart) )
		;ne_profile	= (ne_profile + 1) * 25e17 + nMax[0]

		;specData[0].n	= ne_profile
		;specData[0].n_	= ((specData[0].n)[0:nR-2] + (specData[0].n)[1:nR-1])/2
		;
		;specData[1].n 		= specData[0].n / 2
		;specData[1].n_		= specData[0].n_ / 2
    	   
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

			;for i=0,100 do begin
			;	ne_	= smooth ( ne_, nR / 10, /edge )
			;	ne__	= smooth ( ne__, nR / 10, /edge )
			;endfor
			;ne_	= ne_ * 1.5
			;ne__	= ne__ * 1.5

			specData[1].n	= ne_
			specData[1].n_	= ne__
    
			specData[0].n	= ne_ / ionSpecZ[0]
			specData[0].n_	= ne__ / ionSpecZ[0]
 	
		endif

		for i = 0, nSpec do begin

			specData[i].wp	= sqrt ( specData[i].n * specData[i].q^2 $
    	                        / ( specData[i].m*e0 ))
			specData[i].wc	= specData[i].q * bMag / specData[i].m
			specData[i].wp_	= sqrt ( specData[i].n_ * specData[i].q^2 $
    	                        / ( specData[i].m*e0 ))
			specData[i].wc_	= specData[i].q * bMag_ / specData[i].m

		endfor

		runData = { bField : [ [bR], [bPhi], [bz] ], $
					bField_ : [ [bR_], [bPhi_], [bz_] ], $
					bMag : bMag, $
					bMag_ : bMag_, $
					nR : nR, $
					r : r, $
					r_ : r_, $
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
		    iPlot, r, runData.bField[*,1], $
    	        view_grid = [2,1], $
    	        title = 'bPhi', $
    	        /zoom_on_resize, $
				id = runDataPlot 
		    iPlot, r, runData.bField[*,0], $
				/over, $
				color = red
		    iPlot, r, runData.bField[*,2], $
				/over, $
				color = blue
		    iPlot, r, specData[nSpec].n, $
    	        /view_next, $
    	        title = 'density 0'
			for i=0,nSpec-1 do $
		    iPlot, r, specData[i].n, $
    	        /over 
    	endif


end
