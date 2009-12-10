;	Real Space Full Wave Code in 1D
;	DLG

pro rsfwc_1d, $
    eR = eR, $
    w = wReal, $
    nR = nR, $
    noPlot = noPlot, $
    freeSpace = freeSpace, $
    bandStorage = bandStorage

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
	rMin	= 1.0d0
	rMax	= 2.0d0
	dR	= ( rMax - rMin ) / ( nR - 1 )
	r	= dIndGen ( nR ) * dR + rMin
	;	half grid variables are _
	r_	= r[1:*] - dR / 2

;	Setup plasma and field profiles

    if keyword_set ( freeSpace ) then $
        bMax = 0.0 $
    else $
	    bMax	= 0.5d0
	bPhi	= dblArr(nR)+bMax / r 
	bPhi_	= dblArr(nR-1)+bMax / r_

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

	nPeakR	= 1.5d0
	nMax	= 1.0d19

	;	electrons

	specData[0].q 		= -1 * e
	specData[0].m 		= me
	specData.n 		= ((-(r-nPeakR)^2+(r[0]-nPeakR)^2) / (r[0]-nPeakR)^2 * nMax)
	specData.n_		= ((-(r_-nPeakR)^2+(r[0]-nPeakR)^2) / (r[0]-nPeakR)^2 * nMax)

	;	dueterium

	specData[1].q 		= 1 * e
	specData[1].m 		= 2 * mi 

	for i = 0, nSpec - 1 do begin

		specData[i].wp	= sqrt ( specData[i].n * specData[i].q^2 $
                            / ( specData[i].m*e0 ))
		specData[i].wc	= specData[i].q * bPhi / specData[i].m
		specData[i].wp_	= sqrt ( specData[i].n_ * specData[i].q^2 $
                            / ( specData[i].m*e0 ))
		specData[i].wc_	= specData[i].q * bPhi_ / specData[i].m

	endfor

   
    if not keyword_set ( noPlot ) then begin 
	    iPlot, r, bPhi, $
            view_grid = [3,3], $
            title = 'bPhi', $
            /zoom_on_resize
	    iPlot, r, specData[0].n, $
            /view_next, $
            title = 'density 0'
	    iPlot, r, specData[1].n, $
            /over
    endif

;	Dispersion analysis

    if not keyword_set ( wReal ) then $
	    wReal	= 30d6 * 2d0 * !dpi
	w	= dcomplex ( wReal, wReal * 0.00 )

    if keyword_set ( freeSpace ) then begin
        wReal   = 3000d6 * 2d0 * !dpi 
        w  = wReal
    endif

	nPhi	= 15.0 
	kPar	= nPhi / r
	kz		= 0.0 
	nPar	= kPar * c / wReal

    loadct, 12, /sil, rgb_table = ct12
    device, decomposed = 0

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

        part1   = -4 * stixP * ( -stixD^2 + ( nPar^2 - stixS )^2 ) * stixS $
                    + ( stixD^2 + ( nPar^2 - stixS ) * ( stixP + stixS ) )^2
        part2_1   = stixD^2 + nPar^2 * stixP + nPar^2 * stixS - stixP * stixS - stixS^2
        part2_2   = -stixD^2 - nPar^2 * stixP - nPar^2 * stixS + stixP * stixS + stixS^2

        kPer1   =  sqrt ( complex ( ( -part2_1 - sqrt ( part1 ) ) / ( 2 * stixS ), fltArr ( nR ) ) * wReal / c )
        kPer2   =  sqrt ( complex ( ( part2_2 + sqrt ( part1 ) ) / ( 2 * stixS ), fltArr ( nR ) ) * wReal / c )

	    ;AA	= stixS
	    ;BB	= -1d0 * ( stixR * stixL + stixP * stixS - nPar^2 * ( stixP + stixS ) )
	    ;CC	= stixP * ( nPar^2 - stixR ) * ( nPar^2 - stixL )
	    ;B24AC	= BB^2 - 4d0 * AA * CC
	    ;kPerp1	= sqrt ( ( -BB $
        ;    + sqrt ( complex ( B24AC, B24AC * 0 ) ) ) / ( 2d0 * AA ) ) * wReal / c
	    ;kPerp2	= sqrt ( ( -BB $
        ;    - sqrt ( complex ( B24AC, B24AC * 0 ) ) ) / ( 2d0 * AA ) ) * wReal / c

        kR1  = sqrt ( kPer1^2 - kz^2 )
        kR2  = sqrt ( kPer2^2 - kz^2 )

        if not keyword_set ( noPlot ) then begin

	        iPlot, r, real_part ( kPer1 ), $
                /yLog, $
                yRange = [1, 1e3], $
                /view_next, $
                title = 'dispersion root 0'
	        iPlot, r, imaginary ( kPer1 ), $
	        	lineStyle = 1, $
                /over
	        iPlot, r, real_part ( kR1 ), $
                /over, $
                color = transpose ( ct12[8*16-1,*] )
	        iPlot, r, imaginary ( kR1 ), $
	        	lineStyle = 1, $
                /over, $
                color = transpose ( ct12[8*16-1,*] )


	        iPlot, r, real_part ( kPer2 ), $
                /over
	        iPlot, r, imaginary ( kPer2 ), $
	        	lineStyle = 1, $
                /over
	        iPlot, r, real_part ( kR2 ), $
                /over, $
                color = transpose ( ct12[8*16-1,*] )
	        iPlot, r, imaginary ( kR2 ), $
	        	lineStyle = 1, $
                /over, $
                color = transpose ( ct12[8*16-1,*] )

        endif

	    !p.multi = 0

    endelse

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

	for i = 0L, nR - 1L do begin

		;	r component
            if keyword_set ( bandStorage ) then begin

                j_  = 3*i
                if i gt 0 then begin
                    i_  = 3*i-2
			    	aMat_bandStorage[i_,nuca+j_-i_]	= -II * nPhi * r_[i-1] / ( r[i]^2 * dr ) $
                                        - w^2*r_[i-1]/(2*c^2*r[i])*epsilon_[1,0,i-1]
                    i_  = 3*i-1
			    	aMat_bandStorage[i_,nuca+j_-i_]	= -II * kz / dr $
                                        - w^2*r_[i-1]/(2*c^2*r[i])*epsilon_[2,0,i-1]
                endif
                    
                    i_=3*i    
			    	aMat_bandStorage[i_,nuca+j_-i_]	= nPhi^2 / r[i]^2 + kz^2 $
                                        - w^2 / c^2 * epsilon[0,0,i]
                if i lt nR-1 then begin
                    i_=3*i+1
			    	aMat_bandStorage[i_,nuca+j_-i_]	= II * nPhi * r_[i] / ( r[i]^2 * dr ) $
                                        - w^2*r_[i]/(2*c^2*r[i])*epsilon_[1,0,i]
                    i_=3*i+2
			    	aMat_bandStorage[i_,nuca+j_-i_]	= II * kz / dr $
                                        - w^2*r_[i]/(2*c^2*r[i])*epsilon_[2,0,i]
                endif

            endif else begin

                if i gt 0 then begin
			    	aMat[3*i-2,3*i]	= -II * nPhi * r_[i-1] / ( r[i]^2 * dr ) $
                                        - w^2*r_[i-1]/(2*c^2*r[i])*epsilon_[1,0,i-1]
			    	aMat[3*i-1,3*i]	= -II * kz / dr $
                                        - w^2*r_[i-1]/(2*c^2*r[i])*epsilon_[2,0,i-1]
                endif
			    	aMat[3*i,3*i]	= nPhi^2 / r[i]^2 + kz^2 $
                                        - w^2 / c^2 * epsilon[0,0,i]
                if i lt nR-1 then begin
			    	aMat[3*i+1,3*i]	= II * nPhi * r_[i] / ( r[i]^2 * dr ) $
                                        - w^2*r_[i]/(2*c^2*r[i])*epsilon_[1,0,i]
			    	aMat[3*i+2,3*i]	= II * kz / dr $
                                        - w^2*r_[i]/(2*c^2*r[i])*epsilon_[2,0,i]
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
                                        - w^2*r[i]/(2*c^2*r_[i])*epsilon[0,1,i]
                    i_=3*i+1
				    aMat_bandStorage[i_,nuca+j_-i_]	= kz^2 + r_[i]/(r[i+1]*dr^2) $
                                            + r_[i]/(r[i]*dr^2) $
                                            - w^2/c^2*epsilon_[1,1,i]
                    i_=3*i+2
				    aMat_bandStorage[i_,nuca+j_-i_]	= -kz*nPhi/r_[i] - w^2/(c^2)*epsilon_[2,1,i]
                    i_=3*i+3
				    aMat_bandStorage[i_,nuca+j_-i_]	= II * nPhi / (r[i+1]*dr) $
                                            - w^2*r[i+1]/(2*c^2*r_[i])*epsilon[0,1,i+1]
                    i_=3*i+4
				    if i lt nR-2 then $
				    aMat_bandStorage[i_,nuca+j_-i_]	= -r_[i+1]/(r[i+1]*dr^2)

                endif else begin

                    if i gt 0 then $
				    aMat[3*i-2,3*i+1] = -r_[i-1] / ( r[i] * dr^2 ) 
				    aMat[3*i,3*i+1]	= -II*nPhi/(r[i]*dr) $
                                        - w^2*r[i]/(2*c^2*r_[i])*epsilon[0,1,i]
				    aMat[3*i+1,3*i+1]	= kz^2 + r_[i]/(r[i+1]*dr^2) $
                                            + r_[i]/(r[i]*dr^2) $
                                            - w^2/c^2*epsilon_[1,1,i]
				    aMat[3*i+2,3*i+1]	= -kz*nPhi/r_[i] - w^2/(c^2)*epsilon_[2,1,i]
				    aMat[3*i+3,3*i+1]	= II * nPhi / (r[i+1]*dr) $
                                            - w^2*r[i+1]/(2*c^2*r_[i])*epsilon[0,1,i+1]
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
                                        -w^2*r[i]/(2*c^2*r_[i])*epsilon[0,2,i]
                    i_=3*i+1
				    aMat_bandStorage[i_,nuca+j_-i_] = -nPhi*kz/r_[i] - w^2/c^2*epsilon_[1,2,i]
                    i_=3*i+2
				    aMat_bandStorage[i_,nuca+j_-i_] = (r[i+1]+r[i])/(r_[i]*dr^2) + nPhi^2/r_[i]^2 $
                                            - w^2/c^2*epsilon_[2,2,i]
                    i_=3*i+3
				    aMat_bandStorage[i_,nuca+j_-i_] = II * kz * r[i+1] / (r_[i]*dr) $
                                            - w^2*r[i+1]/(2*c^2*r_[i])*epsilon[0,2,i+1]
                    i_=3*i+5
				    if i lt nR-2 then $
				    aMat_bandStorage[i_,nuca+j_-i_] = -r[i+1] / ( r_[i] * dr^2 ) 

                endif else begin
                
	                if i gt 0 then $
				    aMat[3*i-1,3*i+2] = -r[i]/(r_[i]*dr^2)
				    aMat[3*i,3*i+2] = -II*kz*r[i]/(r_[i]*dr)$
                                        -w^2*r[i]/(2*c^2*r_[i])*epsilon[0,2,i]
				    aMat[3*i+1,3*i+2] = -nPhi*kz/r_[i] - w^2/c^2*epsilon_[1,2,i]
				    aMat[3*i+2,3*i+2] = (r[i+1]+r[i])/(r_[i]*dr^2) + nPhi^2/r_[i]^2 $
                                            - w^2/c^2*epsilon_[2,2,i]
				    aMat[3*i+3,3*i+2] = II * kz * r[i+1] / (r_[i]*dr) $
                                            - w^2*r[i+1]/(2*c^2*r_[i])*epsilon[0,2,i+1]
				    if i lt nR-2 then $
				    aMat[3*i+5,3*i+2] = -r[i+1] / ( r_[i] * dr^2 ) 

                endelse
	
			endif

	endfor

;    if keyword_set ( bandStorage ) then begin
;
;        test_bandStorage    = complexArr ( nAll, m )
;        for i=0,nAll-1 do begin
;            for j=0,nAll-1 do begin
;
;                if abs ( aMat[i,j] ) ne 0 then begin
;                    test_bandStorage[i,nuca+j-i] = aMat[i,j] 
;                    ;if test_bandStorage[(i-j+nuca+1)*i+j] ne aMat_bandStorage[(i-j+nuca+1)*i+j] then begin
;                    ;    print, aMat[i,j], test_bandStorage[(i-j+nuca+1)*i+j], aMat_bandStorage[(i-j+nuca+1)*i+j], i, j
;                    ;endif
;                endif
;
;            endfor
;        endfor
;    stop
;    endif

    ;endelse

	rhs[(nR-2)*3+2]	= II * w * u0 / ( r_[0] * dr ) * 1d0

;	Solve matrix

    print, '*** solving linear system'


    if keyword_set ( bandStorage ) then begin

        eField = imsl_sp_bdSol ( rhs, nlca, nuca, aMat_bandStorage[*], /double )

    endif else begin

	    eField	= la_linear_equation ( aMat, rhs, status = stat )
	    print, 'lapack status: ', stat

    endelse
    
	;eFieldTmp	= dcomplexArr ( nAll )
	;eFieldTmp	= eField
	;eField	= temporary ( eFieldTmp )

	ii_eR	= indGen(nR, /long)*3
	ii_ePhi	= indGen(nR-1, /long)*3+1
	ii_ez	= indGen(nR-1, /long)*3+2
	eR	= eField[ii_eR]
	ePhi	= eField[ii_ePhi]
	ez	= eField[ii_ez]

;   Calculate the Div of D @ the z,phi grid pts
;   but first we need to invert epsilon to get D & D_

    divD_   = dcomplexArr ( nR - 1 )

    print, '*** calculating div D'
    for i=1,nR-3 do begin 
        
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
    for i=0,nR-2 do begin
        divE_[i]    = II * nPhi * ePhi[i] / r_[i] $
                        + II * kz * ez[i] $
                        + ( r[i+1] * eR[i+1] - r[i] * eR[i] ) / ( r_[i] * dr ) 
    endfor

    rho_ = e0 * divE_


    if not keyword_set ( noPlot ) then begin

   
    ;	Visualise solution
    
        loadct, 12, /sil
    	!p.charSize = 3
    
        iPlot, r, eR, $
            /view_next, $
            yTickFont_size = 30, $
            title = 'eR'
        iPlot, r, imaginary ( eR ), $
            /over, $
            color = transpose ( ct12[8*16-1,*] )
    	iPlot, [r_[0]-dr,r_,r_[nR-2]+dr], [0,ePhi,0], $
            psym = -4, $
            /view_next, $
            title = 'ePhi'
    	iPlot, r_, imaginary(ePhi), $
            psym = -4, $
            /over, $
            color = transpose ( ct12[8*16-1,*] )
    	iPlot, [r_[0]-dr,r_,r_[nR-2]+dr], [0,ez,0], $
            psym = -4, $
            /view_next, $
            title = 'ez'
    	iPlot, [r_[0]-dr,r_,r_[nR-2]+dr], [0,imaginary(ez),0], $
            psym = -4, $
            /over, $
            color = transpose ( ct12[8*16-1,*] )

        iPlot, findgen(nR)*dr*nR*2*!pi,(abs(fft(eR)))>0.001,$
            xrange=[0,900], $
            /ylog, $
            yrange=[0.001,1e3], $
            /over, $
            /view_next, $
            title = 'Soln spectrum'

        iPlot, findgen(nR-1)*dr*nr*2*!pi,(abs(fft(ePhi)))>0.001, /over
        iPlot, findgen(nR-1)*dr*nr*2*!pi,(abs(fft(ez)))>0.001, /over

        iPlot, [median(abs(kPer1)),median(abs(kPer1))], [0.001,1e3], /over
        iPlot, [mean(abs(kPer2)),mean(abs(kPer2))], [0.001,1e3], /over
        iPlot, [median(abs(kPer1)),median(abs(kPer1))], [0.001,1e3], /over
        iPlot, [mean(abs(kPer2)),mean(abs(kPer2))], [0.001,1e3], /over

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
stop
end
