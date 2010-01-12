pro run_setup, $
	useEqdsk = useEqdsk, $
	useProfiles = useProfiles, $
	runData = runData

	common constants

;	Variables

	r0	= 0.67d0
	aWall	= 0.22 

	rMin	= r0 - aWall*0.95
	rMax	= r0 + aWall*0.95
   	b0	= 5.85d0
	bR_frac	= 0.1
	bz_frac	= 0.0	
	ionSpecAmu	= [ 2, 3, 1 ]
	ionSpecZ	= [ 1, 2, 1 ]
	damping = 0.0
	freq	= 80.5d6
	nPhi = 10
	kz = 0
	nR	= 256L

	common switches, $
		dielectric_freeSpace, $
		dielectric_noPoloidal, $
		dispersion_freeSpace, $
		dispersion_jaegerPRL, $
		bandStorage	

	dielectric_freeSpace = 0
	dielectric_noPoloidal = 0
	dispersion_freeSpace = 0
	dispersion_jaegerPRL = 0
	bandStorage = 1

;	Grid

	dR	= ( rMax - rMin ) / ( nR - 1 )
	r	= dIndGen ( nR ) * dR + rMin
	;	half grid variables are _
	r_	= r[1:*] - dR / 2

;	Setup plasma and field profiles

	bPhi	= dblArr(nR) + b0 / r * r0
	bPhi_	= dblArr(nR-1) + b0 / r_ * r0

	bR		= dblArr ( nR ) + bPhi * bR_frac 
	bR_		= dblArr ( nR-1 ) + bPhi * bR_frac 
	bz		= dblArr ( nR ) + bPhi * bz_frac 
	bz_		= dblArr ( nR-1 ) + bPhi * bz_frac 

	if keyword_set ( useEqdsk ) then begin

		nstx_eqdsk	= '../eqdsk/g120740.00275.EFIT02.mds.uncorrected.qscale_1.00000.dlgMod_1.67'
		eqdsk	= readgeqdsk ( nstx_eqdsk )
		
		zSlice	= 0.0
		z	= r * 0 + zSlice
		z_	= r_ * 0 + zSlice

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

	if keyword_set ( poloidalScale ) then begin

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

	nMax	= 2.0d20

	for i=0,nSpec do begin

		nProfile	= 1.0 - ( abs ( r - r0 ) / aWall )^2
		nProfile_	= 1.0 - ( abs ( r_ - r0 ) / aWall )^2

		; ions
		if i lt nSpec then begin

			specData[i].q 		= ionSpecZ[i] * e
			specData[i].m 		= ionSpecAmu[i] * mi
			specData[i].n 		= nMax* nProfile
			specData[i].n_		= nMax* nProfile_

		; electrons for charge neutrality
		endif else begin

			specData[i].q 		= -e
			specData[i].m 		= me
			specData[i].n 		= total ( specData[0:nSpec-1].n * ionSpecZ ) 
			specData[i].n_		= total ( specData[0:nSpec-1].n_ * ionSpecZ ) 

		endelse

		specData[i].wp	= sqrt ( specData[i].n * specData[i].q^2 $
                            / ( specData[i].m*e0 ))
		specData[i].wc	= specData[i].q * bMag / specData[i].m
		specData[i].wp_	= sqrt ( specData[i].n_ * specData[i].q^2 $
                            / ( specData[i].m*e0 ))
		specData[i].wc_	= specData[i].q * bMag_ / specData[i].m

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

	if keyword_set ( useProfiles ) then begin

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

		specData[0].n	= ne_
		specData[0].n_	= ne__
    
	endif

	runData = { specData : specData, $
				bField : [ [bR], [bPhi], [bz] ], $
				bField_ : [ [bR_], [bPhi_], [bz_] ], $
				bMag : bMag, $
				bMag_ : bMag, $
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
				nSpec : nSpec }

	common dlg_colors, ct12, $
		red, $
		blue, $
		green

	loadct, 12, /sil, rgb_table = ct12
	red	= transpose ( ct12[12*16-1,*] )
	blue	= transpose ( ct12[8*16-1,*] )
	green	= transpose ( ct12[2*16-1,*] )

    if keyword_set ( plot ) then begin 
	    iPlot, r, bPhi, $
            view_grid = [3,2], $
            title = 'bPhi', $
            /zoom_on_resize, $
			id = panel
	    iPlot, r, bR, $
			/over, $
			color = red
	    iPlot, r, bz, $
			/over, $
			color = blue
	    iPlot, r, specData[0].n, $
            /view_next, $
            title = 'density 0', $
			/yLog, $
			yRange = [1e17,1e20]
	    iPlot, r, specData[1].n, $
            /over
    endif


end
