pro run_setup, $
	runData = runData, $
	specData = specData, $
	in_kz = kz, $
	nMax = nMax, $
	nFac = nfac, $
	gradSize = gradSize, $
	freq = freq, $
	poloidalScale = poloidalScale, $
	nPhi = nPhi

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
 	dispersion_jaegerPRL = 1
	dispersion_generalised = 1
 	
 	bandStorage = 1

;	Plotting switches

 	plotRunData	= 1
 	plotDispersionGeneral 	= 1
 	plotDispersionJaeger	= 1
 	plotDispersionNoPol		= 0
	plotSolution = 1
	plotMovie	= 0
	plotFrequencies = 0

;	File write switches

	writeDispersionTxt = 0

;		----------------------------
;	Variables

	r0	= 1.0d0
	aWall	= 0.22 

	rMin	= 0.2;0.5;1.08;r0 - aWall*0.99
	rMax	= 1.8;1.7;1.14;r0 + aWall*0.99
   	b0	= 0.525718819
	bR_frac	= 0.1
	bz_frac	= 0.0	
	ionSpecZ	= [ 1 ]
	ionSpecAmu	= [ 2 ]
	if not keyword_set ( nMax ) then nMax = [ 2.0 ] * 1d18
	if not keyword_set ( nFac ) then nFac = 1.0 
	damping = 0.06
	if not keyword_set ( freq ) then freq = 30.0e6
	if not keyword_set ( nPhi ) then nPhi = 13.0;-22.0
	if not keyword_set ( in_kz ) then kz = 20.0 else kz = in_kz;51.5;63.6 
	nR	= 512L
	antLoc	= 1.6

	useEqdsk = 0
	useProfiles = 0
	if not keyword_set ( poloidalScale ) then poloidalScale = 1.0
	zSlice	= 1.0 
	sliceSlope = 0.0
	profile1 = 0

;	-----------------------------



;	Benchmarking cases

;	@smithe


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
			;bPhi	= bPhi * 0.6
			;bPhi_	= bPhi_ * 0.6

			;;	z towards x
			;rotTh2	=	-aTan ( max(r)-min(r), max(z)-min(z) )
			;rot_y	= [ [ cos ( rotTh2 ), 0, -sin ( rotTh2 ) ], $
			;			[ 0, 1, 0 ], $
			;			[ sin ( rotTh2 ), 0, cos ( rotTh2 ) ] ]
			;inv_rot_y	= transpose ( rot_y )

			;bVec	= rot_y ## [ [bR],[bPhi],[bz] ]
			;bVec_	= rot_y ## [ [bR_],[bPhi_],[bz_] ]

			;bR	= bVec[*,0]
			;bz	= bVec[*,2]

			;bR_	= bVec_[*,0]
			;bz_	= bVec_[*,2]

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
					n_ : dblArr ( nR - 1 ), $
					wLH : dblArr ( nR ), $
					wUH : dblArr ( nR ) }, nSpec + 1 )

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

		;if not keyword_set ( gradSize ) then gradSize = 0.21 ;0.26
		;gradEnd		= 1.2
		;gradFreq	= 1 / (2*gradSize)
		;nGrad	= fix(gradSize / ( dR ) + 1)
		;ne_profile	= fltArr ( nR ) 
		;r_grad	= fIndGen ( nGrad ) * dR
		;ii1	= where ( r lt gradEnd )
		;ii3	= where ( abs ( r - gradEnd ) eq min ( abs ( r - gradEnd ) ) ) + 1
		;ne_profile[*]	= -1
		;ne_profile[ii1]	= 1
		;ne_profile[ii3[0]-nGrad:ii3[0]-1]	= $
		;	cos ( 2*!pi*gradFreq*r_grad )
		;ne_profile	= (ne_profile + 1) * 40e17 + nMax[0]

		;specData[1].n	= ne_profile
		;specData[1].n_	= ((specData[1].n)[0:nR-2] + (specData[1].n)[1:nR-1])/2
		;
		;specData[0].n 		= specData[1].n / ionSpecZ[0]
		;specData[0].n_		= specData[1].n_ / ionSpecZ[0]
    	   
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

			if keyword_set ( nFac ) then begin
				ne_	= ne_ * nFac 
				ne__	= ne__ * nFac
			endif

			;iplot, r, specData[1].n, thick = 2, trans = 50
			;iPlot, r, specData[0].n, thick = 2, /over, trans = 50
		
			specData[1].n	= ne_
			specData[1].n_	= ne__
    
			specData[0].n	= ne_ / ionSpecZ[0]
			specData[0].n_	= ne__ / ionSpecZ[0]

			;iPlot, r, ne_, /over
			;iplot, r, ne_ / ionSPecZ[0], /over
 
			;iPlot, eqdsk.rbbbs, eqdsk.zbbbs, /iso
			;iPlot, r, z, /over
;stop	
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
