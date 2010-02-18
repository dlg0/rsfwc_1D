pro k_vs_th

    @constants

    ionSpecZ    = [ 1 ]
    ionSpecAmu  = [ 2 ]

    nIonSpec	= n_elements ( ionSpecZ )

	specData	= replicate ( $
			{ 	q : 0d0, $
				m : 0d0, $
				wp : 0.0, $
				wc : 0.0, $
				n : 0.0, $
				wLH : 0.0, $
				wUH : 0.0 }, nIonSpec + 1 )


    nMax     = [ 2.7e18 ]
    bMag   = 0.53

	for i=0,nIonSpec do begin

		; ions
		if i lt nIonSpec then begin

			specData[i].q 		= ionSpecZ[i] * e
			specData[i].m 		= ionSpecAmu[i] * mi
			specData[i].n 		= nMax[i] 

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
		specData[i].wc	= specData[i].q * bMag / specData[i].m
	
	endfor

	for i = 0, nIonSpec - 1 do begin

		specData[i].wLH	= sqrt ( specData[i].wc^2 $
			+ specData[i].wp^2 / ( 1.0 + specData[nIonSpec].wp^2 / specData[nIonSpec].wc^2 ) )
		specData[i].wUH	= sqrt ( specData[nIonSpec].wp^2 + specData[nIonSpec].wc^2 ) 
	
	endfor


    nFreq   = 500
	freq	= 10d0^(fIndGen(nFreq)/(nFreq-1)*5+6)
    fNorm   = specData[0].wc / ( 2 * !pi )

    iiNSTX  = where ( abs ( freq - 30d6 ) eq min ( abs ( freq - 30d6 ) ) )

    sTh = 0.0
    eTh = 360.0 
    nTheta  = 360 
    dTh = (eTh-sTh)/(nTheta-1)
    thetaArr   = sTh + fIndGen ( nTheta ) * dTh

    kSq1 = complexArr ( nFreq, nTheta )
    kSq2 = complexArr ( nFreq, nTheta )
    uSq1 = complexArr ( nFreq, nTheta )
    uSq2 = complexArr ( nFreq, nTheta )

    for f=0,nFreq-1 do begin

        wReal   = 2d0 * !dpi * freq[f]

	    stixR	= 1d0 - total ( specData.wp^2 / ( wReal * ( wReal + specData.wc ) ))
	    stixL	= 1d0 - total ( specData.wp^2 / ( wReal * ( wReal - specData.wc ) ))
	    stixS	= 0.5d0 * ( stixR + stixL )
	    stixD	= 0.5d0 * ( stixR - stixL )
	    stixP	= 1d0 - total ( specData.wp^2 / wReal^2 )
    
        for i = 0, nTheta-1 do begin

            AA   = stixS * sin ( thetaArr[i] * !dtor )^2 + stixP * cos ( thetaArr[i] * !dtor )^2
            BB   = stixR * stixL * sin ( thetaArr[i] * !dtor )^2 + stixP * stixS * ( 1d0 + cos ( thetaArr[i]*!dtor )^2 )
            CC   = stixP * stixR * stixL

            c_	= [ CC, -BB, AA ]
	        ksqTmp	= imsl_zeroPoly ( c_, /double ) * wReal / c

            kSq1[f,i]    = ksqTmp[0]
            kSq2[f,i]    = ksqTmp[1]

            c_	= [ AA, -BB, CC ]
	        uSqTmp	= imsl_zeroPoly ( c_, /double ) 

            uSq1[f,i]    = uSqTmp[0]
            uSq2[f,i]    = uSqTmp[1]

        endfor

    endfor

    u1  = sqrt(uSq1);wReal / ( sqrt(kSq1) * c )
    u2  = sqrt(uSq2);wReal / ( sqrt(kSq2) * c )

   
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

    iPlot, real_part(u1[iiNSTX,*]), thetaArr*!dtor, /polar, color = red, /iso
    iPlot, imaginary(u1[iiNSTX,*]), thetaArr*!dtor, /polar, /over, color = red, lineStyle = 2
    iPlot, (u2[iiNSTX,*]), thetaArr*!dtor, /polar, /over, color = blue
    iPlot, imaginary(u2[iiNSTX,*]), thetaArr*!dtor, /polar, /over, color = blue, lineStyle = 2
    iPlot, [0.05,0.05], [0,0], /polar, /over, color = black, sym_index = 4, sym_thick = 4

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
