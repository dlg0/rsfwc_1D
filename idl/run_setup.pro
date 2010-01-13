pro run_setup, $
	runData = runData, $
	specData = specData

	common constants
	common switches
	common plotSwitches

	common dlg_colors, ct12, $
		red, $
		blue, $
		green

	loadct, 12, /sil, rgb_table = ct12
	red	= transpose ( ct12[12*16-1,*] )
	blue	= transpose ( ct12[8*16-1,*] )
	green	= transpose ( ct12[2*16-1,*] )

;	Calculation method switches

 	dielectric_freeSpace = 0
 	dielectric_noPoloidal = 0
 	dispersion_freeSpace = 0
 	dispersion_jaegerPRL = 0
 	
 	bandStorage = 1

;	Plotting switches

 	plotRunData	= 1
 	plotDispersionGeneral 	= 1
 	plotDispersionJaeger	= 1
 	plotDispersionNoPol		= 0

	runBench = 0
	if runBench then begin

		restore, 'bench/bench_smithe_icw.sav'

	endif else begin

;		----------------------------

;		Variables

		r0	= 0.67d0
		aWall	= 0.22 

		rMin	= r0 - aWall*0.95
		rMax	= r0 + aWall*0.95
   		b0	= 5.85d0
		bR_frac	= 0.15
		bz_frac	= 0.0	
		ionSpecZ	= [ 1, 2, 1]
		ionSpecAmu	= [ 2, 3, 1 ]
		nMax		= [ 0.88, 0.23, 0.66 ] * 1d20
		damping = 0.00
		freq	= 80.5e6
		nPhi = 10.0
		kz = 0.1
		nR	= 4096L
		antLoc	= 1.8

		useEqdsk = 0
		useProfiles = 0
		poloidalScale = 1
		zSlice	= 0.8
		profile1 = 1

;		-----------------------------


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

		;gradSize = 0.1
		;gradStart	= 1.2
		;gradFreq	= 1 / (2*gradSize)
		;ne_profile	= fltArr ( nR ) 
		;ii1	= where ( r lt gradStart )
		;ii2	= where ( r ge gradStart and r lt gradStart+gradSize )
		;ne_profile[*]	= -1
		;ne_profile[ii1]	= 1
		;ne_profile[ii2]	= cos ( 2*!pi*gradFreq*(r[ii2]-gradStart) )
		;ne_profile	= (ne_profile + 1) * 50e17 + nMax

		;specData[0].n	= ne_profile
		;specData[0].n_	= ((specData[0].n)[0:nR-2] + (specData[0].n)[1:nR-1])/2
		
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

			specData[1].n	= ne_
			specData[1].n_	= ne__
    
			specData[0].n	= ne_ / ionSpecAmu[0]
			specData[0].n_	= ne__ / ionSpecAmu[0]
 	
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

		;	Create a benchmark save file from this data

		;save, specData, runData, fileName =  
	endelse

end
