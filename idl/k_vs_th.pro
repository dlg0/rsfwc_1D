pro k_vs_th

    @constants

    ionSpecZ    = [ 1.0 ]
    ionSpecAmu  = [ 2.0 ]

    nIonSpec	= n_elements ( ionSpecZ )

	r   	= 1.4
    z       = -0.00
	
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

    r0  = 1.0
    b0  = 0.53
    bPhi    = b0 / r * r0
    bR  = 0.1 * bPhi
    bz  = 0.0 * bPhi
 
    bMag    = sqrt ( bR^2 + bPhi^2 + bz^2 )

    nstx_profile	= '../profiles/dlg_profiles_x1.00.nc' 
	cdfId = ncdf_open ( nstx_profile, /noWrite ) 
    ncdf_varGet, cdfId, 'ne', xMap_ne
    ncdf_varGet, cdfId, 'r2d', xMap_R2D
    ncdf_varGet, cdfId, 'z2d', xMap_z2D
    nCdf_close, cdfId
         
    nX  = n_elements ( xMap_R2D[*,0] )
    nY  = n_elements ( xMap_R2D[0,*] )

	ne_  = interpolate ( xMap_ne, $
			( r - min(xMap_R2d) ) / (max(xMap_R2d)-min(xMap_R2d)) * (nX-1.0), $
        ( z - min ( xMap_z2d ) ) / (max(xMap_z2d)-min(xMap_z2d)) * (nY-1.0), $
		cubic = -0.5 )

    nMax    = [ ne_ / ionSpecZ[0] ]
    nMax    = 2d19

    create_specData, ionSpecZ, ionSpecAmu, nMax, bMag, $
        specData = specData 

    nFreq   = 1 
	;freq	= 10d0^(fIndGen(nFreq)/(nFreq-1)*5+6)
    freq    = 30d6
    fNorm   = specData[0].wc / ( 2 * !pi )

    iiNSTX  = where ( abs ( freq - 30d6 ) eq min ( abs ( freq - 30d6 ) ) )

    sTh = 0.0
    eTh = 360.0 
    nTheta  = 3600 
    dTh = (eTh-sTh)/(nTheta-1)
    thetaArr   = sTh + fIndGen ( nTheta ) * dTh

    kSq1 = complexArr ( nFreq, nTheta )
    kSq2 = complexArr ( nFreq, nTheta )
    uSq1 = complexArr ( nFreq, nTheta )
    uSq2 = complexArr ( nFreq, nTheta )

    for f=0,nFreq-1 do begin

        wReal   = 2d0 * !dpi * freq[f]

        stixVariables, wReal, specData, $
	        stixVars = stixVars, /noHalf

        for i = 0, nTheta-1 do begin

            AA   = stixVars.stixS * sin ( thetaArr[i] * !dtor )^2 $
                    + stixVars.stixP * cos ( thetaArr[i] * !dtor )^2
            BB   = stixVars.stixR * stixVars.stixL * sin ( thetaArr[i] * !dtor )^2 $
                    + stixVars.stixP * stixVars.stixS * ( 1d0 + cos ( thetaArr[i]*!dtor )^2 )
            CC   = stixVars.stixP * stixVars.stixR * stixVars.stixL

            c_	= [ CC, -BB, AA ]
	        ksqTmp	= imsl_zeroPoly ( c_, /double ) * wReal^2 / c^2

            kSq1[f,i]    = ksqTmp[0]
            kSq2[f,i]    = ksqTmp[1]

            c_	= [ AA, -BB, CC ]
	        uSqTmp	= imsl_zeroPoly ( c_, /double ) 

            uSq1[f,i]    = uSqTmp[0]
            uSq2[f,i]    = uSqTmp[1]

        endfor

    endfor
stop
    u1  = sqrt(uSq1);wReal / ( sqrt(kSq1) * c )
    u2  = sqrt(uSq2);wReal / ( sqrt(kSq2) * c )
    u1  = 1.0 / sqrt(kSq1); * wReal / c
    u2  = 1.0 / sqrt(kSq2); * wReal / c

    ;loadct, 13, file = 'davect.tbl', /sil, rgb_table = ct1
    ;cSF = 253.0/(nTheta/4)    
    ;device, decomposed = 0

    ;for i=0,nTheta/4-1 do begin

    ;    kPlot	= real_part ( kSq1[*,i] )
	;    iiNeg	= where ( kPlot lt 0 )
	;    kPlot[iiNeg]	= -1.0 * kPlot[iiNeg] 
	;    kPlot	= ( kPlot )^(1d0/2d0)
	;    kPlot[iiNeg]	= -1.0 * kPlot[iiNeg]

    ;    iplot, freq / fNorm, kSq1[*,i],$;kPlot, $
    ;        /xLog, $
    ;        sym_index = 4, lineStyle = 6, $
    ;        /zoom_on_resize, /stretch_to_fit, $
    ;        /over, yRange = [-1d4,1d4], sym_size = 0.4, $
    ;        color = transpose ( ct1[((i*cSF-1)>1)<254,*] )

    ;    kPlot	= real_part ( kSq2[*,i] )
	;    iiNeg	= where ( kPlot lt 0 )
	;    kPlot[iiNeg]	= -1.0 * kPlot[iiNeg] 
	;    kPlot	= ( kPlot )^(1d0/2d0)
	;    kPlot[iiNeg]	= -1.0 * kPlot[iiNeg]

    ;    iplot, freq / fNorm, kSq2[*,i],$;kPlot, $
    ;        /xLog, $
    ;        sym_index = 4, lineStyle = 6, $
    ;        /zoom_on_resize, /stretch_to_fit, $
    ;        /over, yRange = [-1d4,1d4], sym_size = 0.4, $
    ;        color = transpose ( ct1[((i*cSF-1)>1)<254,*] )

    ;endfor

    ;iPlot, [specData[0].wLH/(2*!pi)/fNorm,specData[0].wLH/(2*!pi)/fNorm], $
    ;        [-1d4,1d4], /over, thick = 6, trans = 80
    ;iPlot, [specData[0].wc/(2*!pi)/fNorm,specData[0].wc/(2*!pi)/fNorm], $
    ;        [-1d4,1d4], /over, thick = 6, trans = 80
    ;iPlot, [30d6/fNorm,30d6/fNorm], $
    ;        [-1d4,1d4], /over, thick = 6, trans = 80


    @load_colors

    iPlot, real_part(u1[iiNSTX,*]), thetaArr*!dtor, /polar, color = red
    iPlot, imaginary(u1[iiNSTX,*]), thetaArr*!dtor, /polar, /over, color = red, lineStyle = 2
    iPlot, real_part(u2[iiNSTX,*]), thetaArr*!dtor, /polar, /over, color = blue
    iPlot, imaginary(u2[iiNSTX,*]), thetaArr*!dtor, /polar, /over, color = blue, lineStyle = 2

    save, u1, u2, iinstx, thetaarr, file = 'u1u2.sav'
  stop 
    ;iPlot, [0.05,0.05], [0,0], /polar, /over, color = black, sym_index = 4, sym_thick = 4

    ;   CMA diagram

    wReal   = 2d0 * !dpi * 30d6 

    wpe_nstx    = sqrt ( specData[nIonSpec].n * specData[nIonSpec].q^2 $
                                    / ( specData[nIonSpec].m*e0 )) 
    wpi_nstx    = sqrt ( specData[0].n * specData[0].q^2 $
                                    / ( specData[0].m*e0 )) 

    wce_nstx	= specData[nIonSpec].q * bMag / specData[nIonSpec].m

    myX = (wpe_nstx^2+wpi_nstx^2)/wReal^2
    myY = abs ( wce_nstx ) / wReal

    stop
    nDens   = 200
	nMax	= 10d0^(fIndGen(nDens)/(nDens-1)*15+7)
    nB  = 300
    bMag    = 10d0^(fIndGen(nB)/(nB-1)*8-4)

    S2D = fltArr ( n_elements(nMax), n_elements(bMag) )
    L2D = fltArr ( n_elements(nMax), n_elements(bMag) )
    P2D = fltArr ( n_elements(nMax), n_elements(bMag) )
    RLPS2D = fltArr ( n_elements(nMax), n_elements(bMag) )

    bAxis   = fltArr ( n_elements ( bMag ) )
    nAxis   = fltArr ( n_elements ( nMax ) )

    for ii=0,n_elements(nMax)-1 do begin
            print, ii
        for jj=0,n_elements(bMag)-1 do begin


	        for i=0,nIonSpec do begin

	        	; ions
	        	if i lt nIonSpec then begin

	        		specData[i].q 		= ionSpecZ[i] * e
	        		specData[i].m 		= ionSpecAmu[i] * mi
	        		specData[i].n 		= nMax[ii] 

	        	; electrons for charge neutrality
	        	endif else begin

	        		specData[i].q 		= -e
	        		specData[i].m 		= me
	        		if nIonSpec gt 1 then begin
	        			specData[i].n 		= total ( specData[0:nIonSpec-1].n $
	        									* ionSpecZ ) 
	        		endif else begin
	        			specData[i].n 		= specData[0:nIonSpec-1].n $
	        									* ionSpecZ 
	        		endelse

	        	endelse

	        endFor


            for i = 0, nIonSpec do begin

	        	specData[i].wp	= sqrt ( specData[i].n * specData[i].q^2 $
                                    / ( specData[i].m*e0 ))
	        	specData[i].wc	= specData[i].q * bMag[jj] / specData[i].m
	        
	        endfor


	        stixR	= 1d0 - total ( specData.wp^2 / ( wReal * ( wReal + specData.wc ) ))
	        stixL	= 1d0 - total ( specData.wp^2 / ( wReal * ( wReal - specData.wc ) ))
	        stixS	= 0.5d0 * ( stixR + stixL )
	        stixD	= 0.5d0 * ( stixR - stixL )
	        stixP	= 1d0 - total ( specData.wp^2 / wReal^2 )

            S2D[ii,jj]    = stixS
            L2D[ii,jj]    = stixL
            P2D[ii,jj]    = stixP
            RLPS2D[ii,jj]    = (stixR*stixL)-(stixP*stixS) 

            bAxis[jj]    = abs ( specData[nIonSpec].wc ) / wReal
            nAxis[ii]    = ( specData[nIonSpec].wp^2 + specData[0].wp^2 ) / wReal^2

        endfor
    endfor
        
contour, l2d, naxis, baxis, levels = [0.0], /xlog, /ylog, yRange=[1d-1,1d4],xRange=[1d-1,1d6]         
contour, s2d, naxis, baxis, levels = [0.0], /over 
contour, p2d, naxis, baxis, levels = [0.0], /over 
contour, rlps2d, naxis, baxis, levels = [0.0], /over 
oPlot, [1d-2,1d6], [1,1]
tmp = specdata[0].q*bmag/specdata[0].m/wreal 
iiWci   = where ( abs ( 1.0-tmp) eq min ( abs ( 1.0 - tmp ) ) )
oPlot, [1d-2,1d6],[bAxis[iiWci],bAxis[iiWci]]
oPlot, [myX,myX], [myY,myY], psym = 4, symSize = 4
stop
end
