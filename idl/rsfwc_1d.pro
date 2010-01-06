;	Real Space Full Wave Code in 1D
;	DLG

;function plotLog10, y
;
;	iiPos	= where ( y gt 0, iiPosCnt )
;	if iiPosCnt gt 0 then $
;		y[iiPos] = aLog10 ( y[iiPos] )
;	iiNeg	= where ( y lt 0, iiNegCnt )
;	if iiNegCnt gt 0 then $
;		y[iiNeg] = -aLog10 ( -y[iiNeg] )
;
;	return, y
;
;end

function dlg_interpB, r, bStruct, bMag = bMag

	bR  = interpolate ( bStruct.bR, $
			( r[0] - bStruct.rleft ) / bStruct.rdim * (bStruct.nW-1.0), $
        ( r[2] - min ( bStruct.z ) ) / bStruct.zdim * (bStruct.nH-1.0), $
		cubic = -0.5 )
    bPhi  = interpolate ( bStruct.bPhi, $
			( r[0] - bStruct.rleft ) / bStruct.rdim * (bStruct.nW-1.0), $
        ( r[2] - min ( bStruct.z ) ) / bStruct.zdim * (bStruct.nH-1.0), $
		cubic = -0.5 )
    bz  = interpolate ( bStruct.bz, $
			( r[0] - bStruct.rleft ) / bStruct.rdim * (bStruct.nW-1.0), $
        ( r[2] - min ( bStruct.z ) ) / bStruct.zdim * (bStruct.nH-1.0), $
		cubic = -0.5 )

	bMag	= sqrt ( bR^2 + bPhi^2 + bz^2 )

	return, [ bR, bPhi, bz ]

end


pro rsfwc_1d, $
    eR = eR, $
	ePhi = ePhi, $
	ez = ez, $
    w = wReal, $
    nR = nR, $
    plot = plot, $
    freeSpace = freeSpace, $
    bandStorage = bandStorage, $
	divD = divD, $
	kz = kz, $
	nPhi = nPhi, $
	rFull = r, rHalf = r_, $
	kR1 = kR1, kR2 = kR2, $
	damping = damping, $
	useEqdsk = useEqdsk, $
	useProfiles = useProfiles, $
	noPoloidal = noPoloidal
	

	if not keyword_set ( band ) then band = 1
;	Parameters

	e	= 1.60217646d-19
	me	= 9.10938188d-31
   	mi	= 1.67262158d-27
	c	= 299792458d0	
	e0	=	8.854187817d-12 
	u0	=	4d0 * !dpi * 10d0^(-7) 

;	Setup Grid

	if not keyword_set ( nR ) then $
        nR = 256L

	r0	= 1.00
	rMin	= 0.25 
	rMax	= 1.7

	dR	= ( rMax - rMin ) / ( nR - 1 )
	r	= dIndGen ( nR ) * dR + rMin
	;	half grid variables are _
	r_	= r[1:*] - dR / 2

;	Setup plasma and field profiles

	if keyword_set ( useEqdsk ) then begin

		nstx_eqdsk	= '../eqdsk/g120740.00275.EFIT02.mds.uncorrected.qscale_1.00000.dlgMod_1.67'
		eqdsk	= readgeqdsk ( nstx_eqdsk )

		z	= r * 0  + 0.8
		z_	= r_ * 0 + 0.8 

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
	
	endif else begin

    	if keyword_set ( freeSpace ) then $
    	    bMax = 0.0 $
    	else $
		    bMax	= 0.375d0

		bPhi	= dblArr(nR)+bMax / r * r0
		bPhi_	= dblArr(nR-1)+bMax / r_ * r0

		bR		= dblArr ( nR ) - bPhi * 0.15;+ cos ( 15 * (r - r0) ) * bPhi * 0.15 
		bR_		= dblArr ( nR-1 ) - bPhi * 0.15;+ cos ( 15 * (r - r0) ) * bPhi_ * 0.15 
		bz		= dblArr ( nR ) + bPhi * 0.15 ;+ sin ( 15 * (r - r0) ) * bPhi * 0.15
		bz_		= dblArr ( nR-1 ) + bPhi * 0.15;+ sin ( 15 * (r_ - r0) ) * bPhi_ * 0.15 

	endelse

	if keyword_set ( noPoloidal ) then begin

		bR[*]	= 0
		bR_[*]	= 0
		bz[*]	= 0
		bz_[*]	= 0

	endif

	bMag	= sqrt ( bR^2 + bPhi^2 + bz^2 )
	bMag_	= sqrt ( bR_^2 + bPhi_^2 + bz_^2 )

	nSpec	= 2
	specData	= replicate ( $
			{ 	q : 0d0, $
				m : 0d0, $
				wp : dblArr ( nR ), $
				wp_ : dblArr ( nR - 1 ), $
				wc : dblArr ( nR ), $
				wc_ : dblArr ( nR - 1 ), $
				n : dblArr ( nR ), $
				n_ : dblArr ( nR - 1 ) }, nSpec )

	nPeakR	= 2.65d0
	nMax	= 16d18

	;	electrons

	specData[0].q 		= -1 * e
	specData[0].m 		= me
	specData[0].n 		= nMax;((-(r-nPeakR)^2+(r[0]-nPeakR)^2) / (r[0]-nPeakR)^2 * nMax)>0
	specData[0].n_		= nMax;((-(r_-nPeakR)^2+(r[0]-nPeakR)^2) / (r[0]-nPeakR)^2 * nMax)>0

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

	;;   deuterium

	;specData[1].q 		= 1 * e
	;specData[1].m 		= 2 * mi 
	;specData[1].n		= specData[0].n
	;specData[1].n_		= specData[0].n_

	;	helium 4	

	specData[1].q 		= 2 * e
	specData[1].m 		= 4 * mi 
	specData[1].n		= specData[0].n / 2
	specData[1].n_		= specData[0].n_ / 2

	;;	hydrogen

	;specData[3].q 		= 1 * e
	;specData[3].m 		= 1 * mi 
	;specData[3].n		= 0.33 * specData[0].n
	;specData[3].n_		= 0.33 * specData[0].n_


   	;	anntenna location
 
	antLoc  = 0.8;max ( r_ )

	if antLoc gt rMax or antLoc lt rMin then begin
		print, 'ERROR: antenna outside domain, please correct.'
		stop
	endif


	for i = 0, nSpec - 1 do begin

		specData[i].wp	= sqrt ( specData[i].n * specData[i].q^2 $
                            / ( specData[i].m*e0 ))
		specData[i].wc	= specData[i].q * bMag / specData[i].m
		specData[i].wp_	= sqrt ( specData[i].n_ * specData[i].q^2 $
                            / ( specData[i].m*e0 ))
		specData[i].wc_	= specData[i].q * bMag_ / specData[i].m

	endfor

	if keyword_set ( plot ) then begin
    	loadct, 12, /sil, rgb_table = ct12
		red	= transpose ( ct12[12*16-1,*] )
		blue	= transpose ( ct12[8*16-1,*] )
		green	= transpose ( ct12[2*16-1,*] )
	endif

   
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

;	Dispersion analysis

	if not keyword_set ( damping ) then $
		damping = 0.0


    if not keyword_set ( wReal ) then $
	    wReal	= 30e6 * 2 * !pi 
	w	= dcomplex ( wReal, wReal * damping )
	w	= dcomplexArr ( nR ) + wReal
	iiDamped	= where ( r gt 1.4 )
	w[iiDamped]	= complex ( wReal, wReal * damping )

    if keyword_set ( freeSpace ) then begin
        wReal   = 3000d6 * 2d0 * !dpi 
        w  = wReal
    endif

	if not keyword_set ( nPhi ) then $
		nPhi = 0 
	kPar	= nPhi / r
	if not keyword_set ( kz ) then $
		kz = 0 
	nPar	= kPar * c / wReal

    if keyword_set ( freeSpace ) then begin

        k   = wReal / c
        kPer1   = sqrt ( k^2 - kPar^2 )
        kPer2   = kPer1 * 0
        iPlot, r, kPer1, $
            title = 'free space dispersion', $
            /view_next

    endif else begin

	    stixR	= 1d0 - total ( specData.wp^2 / ( wReal * ( wReal + specData.wc ) ),2)
	    stixL	= 1d0 - total ( specData.wp^2 / ( wReal * ( wReal - specData.wc ) ),2)
	    stixS	= 0.5d0 * ( stixR + stixL )
	    stixD	= 0.5d0 * ( stixR - stixL )
	    stixP	= 1d0 - total ( specData.wp^2 / wReal^2, 2 )

	    stixR_	= 1d0 - total ( specData.wp_^2 / ( wReal * ( wReal + specData.wc_ ) ),2)
	    stixL_	= 1d0 - total ( specData.wp_^2 / ( wReal * ( wReal - specData.wc_ ) ),2)
	    stixS_	= 0.5d0 * ( stixR_ + stixL_ )
	    stixD_	= 0.5d0 * ( stixR_ - stixL_ )
	    stixP_	= 1d0 - total ( specData.wp_^2 / wReal^2, 2 )

        part1   = -complex(4d0,0d0) * stixP * ( -stixD^2 + ( nPar^2 - stixS )^2 ) * stixS $
                    + ( stixD^2 + ( nPar^2 - stixS ) * ( stixP + stixS ) )^2
        part2_1   = stixD^2 + nPar^2 * stixP + nPar^2 * stixS - stixP * stixS - stixS^2
        part2_2   = -stixD^2 - nPar^2 * stixP - nPar^2 * stixS + stixP * stixS + stixS^2

        kPerpSq1_   =  ( part2_2 + sqrt ( part1 ) ) / ( 2 * stixS ) * wReal / c 
        kPerpSq2_   =  ( -part2_1 - sqrt ( part1 ) ) / ( 2 * stixS ) * wReal / c 

	    AA	= stixS
	    BB	= -1d0 * ( stixR * stixL + stixP * stixS - nPar^2 * ( stixP + stixS ) )
	    CC	= stixP * ( nPar^2 - stixR ) * ( nPar^2 - stixL )
	    B24AC	= BB^2 - 4d0 * AA * CC
	    kPerpSq1	=  ( -BB $
            + sqrt ( complex ( B24AC, B24AC * 0 ) ) ) / ( 2d0 * AA ) * wReal / c
	    kPerpSq2	= ( -BB $
            - sqrt ( complex ( B24AC, B24AC * 0 ) ) ) / ( 2d0 * AA ) * wReal / c

		kPerp1	= sqrt ( kPerpSq1 )
		kPerp2	= sqrt ( kPerpSq2 )

        kR1  = sqrt ( kPerpSq1 - kz^2 ) * 0
        kR2  = sqrt ( kPerpSq2 - kz^2 ) * 0

		;	poloidal field dispersion analysis

		kPhi	= nPhi / r

		nPhi_	= kPhi * c / wReal	
		bPol	= sqrt ( bR^2 + bz^2 ) / bMag
		bTor	= bPhi / bMag

		c0		= -stixD^2 * stixP $
					+ nPhi_^4 * stixP / bTor^4 $
					- 2 * nPhi_^2 * stixP * stixS / bTor^2 $
					+ stixP * stixS^2
		c1		= 4 * bPol * nPhi_^3 * stixP / bTor^4 $
					- 4 * bPol * nPhi_ * stixP * stixS / bTor^2
		c2		= stixD^2 $
					+ 6 * bPol^2 * nPhi_^2 * stixP / bTor^4 $
					+ nPhi_^2 * stixP / bTor^2 $
					+ nPhi_^2 * stixS / bTor^2 $
					- stixP * stixS $
					- 2 * bPol^2 * stixP * stixS / bTor^2 $
					- stixS^2
		c3		= 4 * bPol^3 * nPhi_ * stixP / bTor^4 $
					+ 2 * bPol * nPhi_ * stixP / bTor^2 $
					+ 2 * bPol * nPhi_ * stixS / bTor^2
		c4		= bPol^4 * stixP / bTor^4 $
					+ bPol^2 * stixP / bTor^2 $
					+ stixS $
					+ bPol^2 * stixS / bTor^2

		kPerp__	= complexArr ( nR, 4 )	
		for i = 0, nR - 1 do begin	
	
			c_	= [ c0[i], c1[i], c2[i], c3[i], c4[i] ]
			kPerp__[i,*]	= imsl_zeroPoly ( c_, /double ) * wReal / c

		endfor

		sqrtArg	= kPerp__^2 - kz^2
		kR	= sqrt ( sqrtArg )
		kR2	= -sqrt ( -sqrtArg)


		if keyword_set ( plot ) then begin

			iPlot, r,kPerp__[*,0], lineStyle = 6, sym_index = 4, /zoom_on_reSize, $
				view_grid = [2,2]
			iPlot, r,imaginary ( kPerp__[*,0] ), lineStyle = 6, sym_index = 3, /over
			iPlot, r,kPerp__[*,1], lineStyle = 6, sym_index = 4, /over
			iPlot, r,imaginary ( kPerp__[*,1] ), lineStyle = 6, sym_index = 3, /over
			iPlot, r,kPerp__[*,2], lineStyle = 6, sym_index = 4, /over
			iPlot, r,imaginary ( kPerp__[*,2] ), lineStyle = 6, sym_index = 3, /over
			iPlot, r,kPerp__[*,3], lineStyle = 6, sym_index = 4, /over
			iPlot, r,imaginary ( kPerp__[*,3] ), lineStyle = 6, sym_index = 3, /over

			for i = 1, nSpec - 1 do begin

				iiRes	= where ( abs ( specData[i].wc - wReal ) $
							eq min ( abs ( specData[i].wc - wReal ) ) )

				if iiRes ne 0 and iiRes ne nR-1 then begin

					iPlot, [ r[iiRes], r[iiRes] ], [ -1e3, 1e3 ], thick = 10, trans = 50, /over
	
				endif			

			endfor
	
			iPlot, r,abs(real_part(kPerp__[*,0])), lineStyle = 6, sym_index = 4, /zoom_on_reSize, $
				/view_next, /ylog, yRange = [1,1e4]
			iPlot, r,abs(imaginary ( kPerp__[*,0] )), lineStyle = 6, sym_index = 3, /over
			iPlot, r,abs(real_part(kPerp__[*,1])), lineStyle = 6, sym_index = 4, /over
			iPlot, r,abs(imaginary ( kPerp__[*,1] )), lineStyle = 6, sym_index = 3, /over
			iPlot, r,abs(real_part(kPerp__[*,2])), lineStyle = 6, sym_index = 4, /over
			iPlot, r,abs(imaginary ( kPerp__[*,2] )), lineStyle = 6, sym_index = 3, /over
			iPlot, r,abs(real_part(kPerp__[*,3])), lineStyle = 6, sym_index = 4, /over
			iPlot, r,abs(imaginary ( kPerp__[*,3] )), lineStyle = 6, sym_index = 3, /over

			iPlot, r,real_part(kR[*,0]), lineStyle = 6, sym_index = 4, /zoom_on_reSize, $
				/view_next
			iPlot, r,imaginary ( kR[*,0] ), lineStyle = 6, sym_index = 3, /over
			iPlot, r,real_part(kR[*,1]), lineStyle = 6, sym_index = 4, /over
			iPlot, r,imaginary ( kR[*,1] ), lineStyle = 6, sym_index = 3, /over
			iPlot, r,real_part(kR[*,2]), lineStyle = 6, sym_index = 4, /over
			iPlot, r,imaginary ( kR[*,2] ), lineStyle = 6, sym_index = 3, /over
			iPlot, r,real_part(kR[*,3]), lineStyle = 6, sym_index = 4, /over
			iPlot, r,imaginary ( kR[*,3] ), lineStyle = 6, sym_index = 3, /over

			iPlot, r,abs(real_part(kR[*,0])), lineStyle = 6, sym_index = 4, /zoom_on_reSize, $
				/view_next, /ylog, yRange = [1,1e4]
			iPlot, r,abs(imaginary ( kR[*,0] )), lineStyle = 6, sym_index = 3, /over
			iPlot, r,abs(real_part(kR[*,1])), lineStyle = 6, sym_index = 4, /over
			iPlot, r,abs(imaginary ( kR[*,1] )), lineStyle = 6, sym_index = 3, /over
			iPlot, r,abs(real_part(kR[*,2])), lineStyle = 6, sym_index = 4, /over
			iPlot, r,abs(imaginary ( kR[*,2] )), lineStyle = 6, sym_index = 3, /over
			iPlot, r,abs(real_part(kR[*,3])), lineStyle = 6, sym_index = 4, /over
			iPlot, r,abs(imaginary ( kR[*,3] )), lineStyle = 6, sym_index = 3, /over

	
		endif

    endelse

;	;Seperate visualisation of the dispersion relation

	;kPerpSq1_re	= plotLog10 ( real_part ( kPerpSq1 ) )
	;kPerpSq1_im	= plotLog10 ( imaginary ( kPerpSq1 ) )
	;kPerpSq2_re	= plotLog10 ( real_part ( kPerpSq2 ) )
	;kPerpSq2_im	= plotLog10 ( imaginary ( kPerpSq2 ) )

	;
	;iPlot, r, kPerpSq1_re, $
	;	/zoom_on_resize
	;iPlot, r, kPerpSq1_im, $
	;	/over, $
	;	lineStyle = 2
	;iPlot, r, kPerpSq2_re, $
	;	/over, $
	;	color = blue
	;iPlot, r, kPerpSq2_im, $
	;	/over, $
	;	color = blue, $
	;	lineStyle = 2


;	Calculate dielectric tensor

	II	= dcomplex ( 0, 1 )

	epsilon		= dcomplexArr ( 3, 3, nR )
	epsilon_	= dcomplexArr ( 3, 3, nR-1 )

    if keyword_set ( freeSpace ) then begin

        epsilon[0,0,*]  = complex ( 1, 0 )
        epsilon[1,1,*]  = complex ( 1, 0 )
        epsilon[2,2,*]  = complex ( 1, 0 )

        epsilon_[0,0,*]  = complex ( 1, 0 )
        epsilon_[1,1,*]  = complex ( 1, 0 )
        epsilon_[2,2,*]  = complex ( 1, 0 )

    endif else begin

	    epsilon[0,0,*]	= stixS
	    epsilon[0,2,*]	= -II * stixD
	    epsilon[1,1,*]	= stixP
	    epsilon[2,0,*]	= II * stixD
	    epsilon[2,2,*]	= stixS

	    epsilon_[0,0,*]	= stixS_
	    epsilon_[0,2,*]	= -II * stixD_
	    epsilon_[1,1,*]	= stixP_
	    epsilon_[2,0,*]	= II * stixD_
	    epsilon_[2,2,*]	= stixS_

    endelse

;	Generic dielectric for arbitrary magnetic field direction

	bUnit_cyl	= [ [ bR / bMag ], $	
					[ bPhi / bMag ], $
					[ bz / bMag ] ]
	bUnit_cyl_	= [ [ bR_ / bMag_ ], $	
					[ bPhi_ / bMag_ ], $
					[ bz_ / bMag_ ] ]

	;	rotate to cartesian

	phi	= 0d0
	bUnit_car	= dblArr ( nR, 3 )
	bUnit_car_	= dblArr ( nR-1, 3 )

	for i = 0, nR - 1 do begin

		;	Here the rotation from cylindrical to cartesian and
		;	its inverse is not used since for phi = 0, the 
		;	rotation is simply the identity matrix.

		;cyl2car	= [ [ cos ( phi ), -sin ( phi ), 0 ], $
		;			[ sin ( phi ), cos ( phi ), 0 ], $
		;			[ 0, 0, 1 ] ]

		;car2cyl	= invert ( cyl2car )

		;bUnit_car[i,*]	= cyl2car ## bUnit_cyl[i,*]
		;bUnit_car_[i,*]	= cyl2car ## bUnit_cyl[i,*]

		bUnit_car[i,*]	= bUnit_cyl[i,*]	

		;	build rotation matrices for dependant on b field

		;	y towards z
		rotTh1	=	-aTan ( bUnit_car[i,1], bUnit_car[i,2] ) 
		rot_x	= [ [ 1, 0, 0], $
					[ 0, cos ( rotTh1 ), sin ( rotTh1 ) ], $
					[ 0, -sin ( rotTh1 ), cos ( rotTh1 ) ] ]
		inv_rot_x	= transpose ( rot_x )

		;	z towards x
		rotTh2	=	aTan ( bUnit_car[i,0], bUnit_car[i,2] ) 
		rot_y	= [ [ cos ( rotTh2 ), 0, -sin ( rotTh2 ) ], $
					[ 0, 1, 0 ], $
					[ sin ( rotTh2 ), 0, cos ( rotTh2 ) ] ]
		inv_rot_y	= transpose ( rot_y )

		;	build dielectric in x,y,z with z assumed along b

		epsilon_stix	= [ [ stixS[i], II * stixD[i], 0 ], $
							[ -II * stixD[i], stixS[i], 0 ], $
							[ 0, 0, stixP[i] ] ]

		;	rotate

		epsilon[*,*,i]	= rot_x ## epsilon_stix ## inv_rot_x 
		epsilon[*,*,i]	= rot_y ## epsilon[*,*,i] ## inv_rot_y
		;;	rotate to cylindrical
		;epsilon[*,*,i]	=	car2cyl ## epsilon_car ## cyl2car


		;	same for 1/2 grid

		if i lt nR - 1 then begin

			bUnit_car_[i,*]	= bUnit_cyl_[i,*]	

			rotTh1_	=	-aTan ( bUnit_car_[i,1], bUnit_car_[i,2] ) 
			rot_x_	= [ [ 1, 0, 0], $
						[ 0, cos ( rotTh1_ ), sin ( rotTh1_ ) ], $
						[ 0, -sin ( rotTh1_ ), cos ( rotTh1_ ) ] ]
			inv_rot_x_	= invert ( rot_x_ )

			rotTh2_	=	aTan ( bUnit_car_[i,0], bUnit_car_[i,2] ) 
			rot_y_	= [ [ cos ( rotTh2_ ), 0, -sin ( rotTh2_ ) ], $
						[ 0, 1, 0 ], $
						[ sin ( rotTh2_ ), 0, cos ( rotTh2_ ) ] ]
			inv_rot_y_	= invert ( rot_y_ )

			epsilon_stix_	= [ [ stixS_[i], II * stixD_[i], 0 ], $
								[ -II * stixD_[i], stixS_[i], 0 ], $
								[ 0, 0, stixP_[i] ] ]

			epsilon_[*,*,i]	= rot_x_ ## epsilon_stix_ ## inv_rot_x_
			epsilon_[*,*,i]	= rot_y_ ## epsilon_[*,*,i] ## inv_rot_y_

		endif
	
	endfor

;	Build matrix

    print, '*** filling matrix'

	nAll	= nR + 2L * ( nR - 1L )
    if keyword_set ( bandStorage ) then begin

        nuca    = 5L
        nlca    = 3L
        m   = nlca + nuca + 1L
        aMat_bandStorage    = complexArr ( nAll, m )

        print, 'aMat: ', (nAll*2.0)^2*16.0/(1024.0^2)
        print, 'aMat_bandStorage: ', (m*1d0*nAll*2.0)*16.0/(1024.0^2) 

    endif else begin

	    aMat	= complexArr ( nAll, nAll )

    endelse

	rhs		= complexArr ( nAll )

	for i = 0UL, nR - 1L do begin

		;	r component
            if keyword_set ( bandStorage ) then begin

                j_  = 3*i
                if i gt 0 then begin
                    i_  = 3*i-2
			    	aMat_bandStorage[i_,nuca+j_-i_]	= -II * nPhi * r_[i-1] / ( r[i]^2 * dr ) $
                                        - w[i]^2*r_[i-1]/(2*c^2*r[i])*epsilon_[1,0,i-1]
                    i_  = 3*i-1
			    	aMat_bandStorage[i_,nuca+j_-i_]	= -II * kz / dr $
                                        - w[i]^2*r_[i-1]/(2*c^2*r[i])*epsilon_[2,0,i-1]
                endif
                    
                    i_=3*i    
			    	aMat_bandStorage[i_,nuca+j_-i_]	= nPhi^2 / r[i]^2 + kz^2 $
                                        - w[i]^2 / c^2 * epsilon[0,0,i]
                if i lt nR-1 then begin
                    i_=3*i+1
			    	aMat_bandStorage[i_,nuca+j_-i_]	= II * nPhi * r_[i] / ( r[i]^2 * dr ) $
                                        - w[i]^2*r_[i]/(2*c^2*r[i])*epsilon_[1,0,i]
                    i_=3*i+2
			    	aMat_bandStorage[i_,nuca+j_-i_]	= II * kz / dr $
                                        - w[i]^2*r_[i]/(2*c^2*r[i])*epsilon_[2,0,i]
                endif

            endif else begin

                if i gt 0 then begin
			    	aMat[3*i-2,3*i]	= -II * nPhi * r_[i-1] / ( r[i]^2 * dr ) $
                                        - w[i]^2*r_[i-1]/(2*c^2*r[i])*epsilon_[1,0,i-1]
			    	aMat[3*i-1,3*i]	= -II * kz / dr $
                                        - w[i]^2*r_[i-1]/(2*c^2*r[i])*epsilon_[2,0,i-1]
                endif
			    	aMat[3*i,3*i]	= nPhi^2 / r[i]^2 + kz^2 $
                                        - w[i]^2 / c^2 * epsilon[0,0,i]
                if i lt nR-1 then begin
			    	aMat[3*i+1,3*i]	= II * nPhi * r_[i] / ( r[i]^2 * dr ) $
                                        - w[i]^2*r_[i]/(2*c^2*r[i])*epsilon_[1,0,i]
			    	aMat[3*i+2,3*i]	= II * kz / dr $
                                        - w[i]^2*r_[i]/(2*c^2*r[i])*epsilon_[2,0,i]
                endif

            endelse


		;	phi component
			if i lt nR-1 then begin
           
                if keyword_set ( bandStorage ) then begin
            
                    j_=3*i+1
                    i_=3*i-2
                    if i gt 0 then $
				    aMat_bandStorage[i_,nuca+j_-i_] = -r_[i-1] / ( r[i] * dr^2 ) 
                    i_=3*i
				    aMat_bandStorage[i_,nuca+j_-i_]	= -II*nPhi/(r[i]*dr) $
                                        - w[i]^2*r[i]/(2*c^2*r_[i])*epsilon[0,1,i]
                    i_=3*i+1
				    aMat_bandStorage[i_,nuca+j_-i_]	= kz^2 + r_[i]/(r[i+1]*dr^2) $
                                            + r_[i]/(r[i]*dr^2) $
                                            - w[i]^2/c^2*epsilon_[1,1,i]
                    i_=3*i+2
				    aMat_bandStorage[i_,nuca+j_-i_]	= -kz*nPhi/r_[i] - w[i]^2/(c^2)*epsilon_[2,1,i]
                    i_=3*i+3
				    aMat_bandStorage[i_,nuca+j_-i_]	= II * nPhi / (r[i+1]*dr) $
                                            - w[i]^2*r[i+1]/(2*c^2*r_[i])*epsilon[0,1,i+1]
                    i_=3*i+4
				    if i lt nR-2 then $
				    aMat_bandStorage[i_,nuca+j_-i_]	= -r_[i+1]/(r[i+1]*dr^2)

                endif else begin

                    if i gt 0 then $
				    aMat[3*i-2,3*i+1] = -r_[i-1] / ( r[i] * dr^2 ) 
				    aMat[3*i,3*i+1]	= -II*nPhi/(r[i]*dr) $
                                        - w[i]^2*r[i]/(2*c^2*r_[i])*epsilon[0,1,i]
				    aMat[3*i+1,3*i+1]	= kz^2 + r_[i]/(r[i+1]*dr^2) $
                                            + r_[i]/(r[i]*dr^2) $
                                            - w[i]^2/c^2*epsilon_[1,1,i]
				    aMat[3*i+2,3*i+1]	= -kz*nPhi/r_[i] - w[i]^2/(c^2)*epsilon_[2,1,i]
				    aMat[3*i+3,3*i+1]	= II * nPhi / (r[i+1]*dr) $
                                            - w[i]^2*r[i+1]/(2*c^2*r_[i])*epsilon[0,1,i+1]
				    if i lt nR-2 then $
				    aMat[3*i+4,3*i+1]	= -r_[i+1]/(r[i+1]*dr^2)
                
                endelse

		;	z component	


                if keyword_set ( bandStorage ) then begin
                    j_=3*i+2
                    i_=3*i-1
                    if i gt 0 then $
				    aMat_bandStorage[i_,nuca+j_-i_] = -r[i]/(r_[i]*dr^2)
                    i_=3*i
				    aMat_bandStorage[i_,nuca+j_-i_] = -II*kz*r[i]/(r_[i]*dr)$
                                        -w[i]^2*r[i]/(2*c^2*r_[i])*epsilon[0,2,i]
                    i_=3*i+1
				    aMat_bandStorage[i_,nuca+j_-i_] = -nPhi*kz/r_[i] - w[i]^2/c^2*epsilon_[1,2,i]
                    i_=3*i+2
				    aMat_bandStorage[i_,nuca+j_-i_] = (r[i+1]+r[i])/(r_[i]*dr^2) + nPhi^2/r_[i]^2 $
                                            - w[i]^2/c^2*epsilon_[2,2,i]
                    i_=3*i+3
				    aMat_bandStorage[i_,nuca+j_-i_] = II * kz * r[i+1] / (r_[i]*dr) $
                                            - w[i]^2*r[i+1]/(2*c^2*r_[i])*epsilon[0,2,i+1]
                    i_=3*i+5
				    if i lt nR-2 then $
				    aMat_bandStorage[i_,nuca+j_-i_] = -r[i+1] / ( r_[i] * dr^2 ) 

                endif else begin
                
	                if i gt 0 then $
				    aMat[3*i-1,3*i+2] = -r[i]/(r_[i]*dr^2)
				    aMat[3*i,3*i+2] = -II*kz*r[i]/(r_[i]*dr)$
                                        -w[i]^2*r[i]/(2*c^2*r_[i])*epsilon[0,2,i]
				    aMat[3*i+1,3*i+2] = -nPhi*kz/r_[i] - w[i]^2/c^2*epsilon_[1,2,i]
				    aMat[3*i+2,3*i+2] = (r[i+1]+r[i])/(r_[i]*dr^2) + nPhi^2/r_[i]^2 $
                                            - w[i]^2/c^2*epsilon_[2,2,i]
				    aMat[3*i+3,3*i+2] = II * kz * r[i+1] / (r_[i]*dr) $
                                            - w[i]^2*r[i+1]/(2*c^2*r_[i])*epsilon[0,2,i+1]
				    if i lt nR-2 then $
				    aMat[3*i+5,3*i+2] = -r[i+1] / ( r_[i] * dr^2 ) 

                endelse
	
			endif

	endfor

    iiAnt   = where ( abs ( r_ - antLoc ) eq min ( abs ( r_ - antLoc ) ) )
    
	rhs[iiAnt*3+2]	= -II * wReal * u0 ( r_[0] * dr ) * 20d0
	;rhs[iiAnt*3]	= II * w * u0 / ( r_[0] * dr ) * 1d0/2
	;rhs[iiAnt*3+3]	= II * w * u0 / ( r_[0] * dr ) * 1d0/2


;	Solve matrix

    print, '*** solving linear system'


    if keyword_set ( bandStorage ) then begin

        eField = imsl_sp_bdSol ( rhs, nlca, nuca, aMat_bandStorage[*], /double )

    endif else begin

	    eField	= la_linear_equation ( aMat, rhs, status = stat )
	    print, 'lapack status: ', stat

    endelse
   
	ii_eR	= lIndGen(nR)*3
	eR	= eField[ii_eR]
	ii_ePhi	= temporary(ii_eR[0:nR-2]+1)
	ePhi	= eField[ii_ePhi]
	ii_ez	= temporary(ii_ePhi+1)
	ez	= eField[ii_ez]

	if keyword_set ( divD ) then begin

	;  	Calculate the Div of D @ the z,phi grid pts
	;  	but first we need to invert epsilon to get D & D_

    	divD_   = dcomplexArr ( nR - 1 )

    	print, '*** calculating div D'
    	for i=1UL,nR-3L do begin 
    	    
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

    	divE_   = dcomplexArr ( nR - 1 )
    	for i=0UL,nR-2L do begin
    	    divE_[i]    = II * nPhi * ePhi[i] / r_[i] $
    	                    + II * kz * ez[i] $
    	                    + ( r[i+1] * eR[i+1] - r[i] * eR[i] ) / ( r_[i] * dr ) 
    	endfor

    	rho_ = e0 * divE_

	endif


    if keyword_set ( plot ) then begin

   
    ;	Visualise solution
    
        loadct, 12, /sil
		eMag	= sqrt ( abs(er)^2 + abs(ePhi)^2 + abs(ez)^2 )
    
        iPlot, r, eR, $
            view_number = 4, $
            title = 'eR', $
			over = panel
        iPlot, r, imaginary ( eR ), $
			view_number = 4, $
            over = panel, $
            color = transpose ( ct12[8*16-1,*] )
		iPlot, [antLoc,antLoc], [-max(abs(eR)),max(abs(eR))], $
			transparency = 80, $
			thick = 10, $
			over = panel, $
			yRange = [ -mean ( eMag ), mean ( eMag ) ] * 2
	
    	iPlot, [r_[0]-dr,r_,r_[nR-2]+dr], [0,ePhi,0], $
            psym = -4, $
            view_number = 5, $
            title = 'ePhi', $
			over = panel
    	iPlot, r_, imaginary(ePhi), $
            psym = -4, $
            over = panel, $
			view_number = 5, $
            color = transpose ( ct12[8*16-1,*] )
		iPlot, [antLoc,antLoc], [-max(abs(ePhi)),max(abs(ePhi))], $
			transparency = 80, $
			thick = 10, $
			over = panel, $
			yRange = [ -mean ( eMag ), mean ( eMag ) ] * 2

    	iPlot, [r_[0]-dr,r_,r_[nR-2]+dr], [0,ez,0], $
            psym = -4, $
            view_number = 6, $
            title = 'ez', $
			over = panel
    	iPlot, [r_[0]-dr,r_,r_[nR-2]+dr], [0,imaginary(ez),0], $
            psym = -4, $
            over = panel, $
            color = transpose ( ct12[8*16-1,*] )
		iPlot, [antLoc,antLoc], [-max(abs(ez)),max(abs(ez))], $
			transparency = 80, $
			thick = 10, $
			over = panel, $
			yRange = [ -mean ( eMag ), mean ( eMag ) ] * 2

		print, 'FFTing solution ...'
		nFFT	= 128 
		fftData_eR	= complexArr ( nFFT )
		fftData_ePhi	= complexArr ( nFFT )
		fftData_ez	= complexArr ( nFFT )
		rFFTData	= fltArr ( nFFT )
		fftStep	= nR / nFFT 
		fftWidth	= nR / 10 
		fftWidth	+= (fftWidth mod 2)
		osf	= 10

		fftData2d	= complexArr ( nFFT, fftWidth*osf )
		hannWin	= hanning ( fftWidth )
		kRFFT	= findgen ( fftWidth * osf ) / ( fftWidth * osf * dr ) * ( 2 * !pi )
		eTmp	= complexArr ( fftWidth * osf )

		for i=0,nFFT-1 do begin
		
			iiLow	= i*fftStep-fftWidth/2
			iiHig	= i*fftStep+fftWidth/2-1

			rFFTData[i]	= r[i*fftStep]

			if iiLow gt 0 and iiHig lt nR-2 then begin

				ezTmp	= ( ez[iiLow:iiHig] + ez[iiLow+1:iiHig+1] )	/ 2
				ePhiTmp	= ( ePhi[iiLow:iiHig] + ePhi[iiLow+1:iiHig+1] )	/ 2
				eRTmp	= eR[iiLow:iiHig] 
				eTmp[fftWidth*osf/2-fftWidth/2:fftWidth*osf/2+fftWidth/2-1]	= $
					eRTmp * hannWin
				fftTmp	= abs ( fft ( eTmp ) )^2
				fftData2D[i,*]	= fftTmp

			endif

		endfor
		print, 'DONE'
	
		nLevs	= 21	
		levels	= 10.0^fIndGen ( nLevs )/1e10
		colors	= 255 - ( bytScl ( fIndGen(nLevs), top = 253 ) + 1 )
		kMax	= 300 
		iiPlotk	= where ( kRFFT gt 0 and kRFFT le kMax )
		iContour, fftData2D[*,iiPlotk]/max(fftData2D), rFFTData, kRFFT[iiPlotk], $
			yRange = [1,kMax], $
			view_number = 3, $
			xRange = [min(r),max(r)], $
			rgb_indices = colors, $
			rgb_table = 1, $
			/fill, $
			c_value = levels, $
			/zoom_on_resize, $
			over = panel


		iPlot, [antLoc,antLoc], [1,1e3], $
			transparency = 80, $
			thick = 10, $
			over = panel
	

		if keyword_set ( divD ) then begin

    		;   plot the numerical and real charge density
 
        	iPlot, r_, rho_, $
        	    color = transpose ( ct12[8*16-1,*] ), $
        	    /over, $
        	    title = 'rho', $
        	    /view_next
        	iPlot, r_, divD_, $
        	    /over
        	iPlot, r_, imaginary ( divD_ ), $
        	    /over, $
        	    lineStyle = 1

        	iPlot, r_, divD_, $
        	    title = 'div D', $
        	    /over, $
        	    /view_next
        	iPlot, r_, imaginary ( divD_ ), $
        	    /over, $
        	    lineStyle = 1

		endif

    endif

	close, /all
end
