;	Real Space Full Wave Code in 1D
;	DLG

pro rsfwc_1d, $
    eR = eR, $
    w = wReal, $
    nR = nR, $
    noPlot = noPlot

;	Parameters

	e	= 1.60217646d-19
	me	= 9.10938188d-31
   	mi	= 1.67262158d-27
	c	= 299792458d0	
	e0	=	8.854187817d-12 
	u0	=	4d0 * !dpi * 10d0^(-7) 

;	Setup Grid

	if not keyword_set ( nR ) then $
        nR = 1024 
	rMin	= 1d0
	rMax	= 1.8d0
	dR	= ( rMax - rMin ) / ( nR - 1 )
	r	= dIndGen ( nR ) * dR + rMin
	;	half grid variables are _
	r_	= r[1:*] - dR / 2

;	Setup plasma and field profiles

	bMax	= 0.6d0
	bPhi	= dblArr(nR)+bMax; / r 
	bPhi_	= dblArr(nR-1)+bMax; / r_

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

	nPeakR	= 1.1d0
	nMax	= 0.05d19

	;	electrons

	specData[0].q 		= -1 * e
	specData[0].m 		= me
	specData.n 		= nMax;(-(r-nPeakR)^2+(r[0]-nPeakR)^2) / (r[0]-1.1)^2 * nMax
	specData.n_		= nMax;(-(r_-nPeakR)^2+(r[0]-nPeakR)^2) / (r[0]-1.1)^2 * nMax

	;	dueterium

	specData[1].q 		= 1 * e
	specData[1].m 		= 2 * mi 

	for i = 0, nSpec - 1 do begin

		specData[i].wp		= sqrt ( specData[i].n * specData[i].q^2 / ( specData[i].m * e0 ) )
		specData[i].wc		= specData[i].q * bPhi / specData[i].m
		specData[i].wp_		= sqrt ( specData[i].n_ * specData[i].q^2 / ( specData[i].m * e0 ) )
		specData[i].wc_		= specData[i].q * bPhi_ / specData[i].m

	endfor

   
    if not keyword_set ( noPlot ) then begin 
	    iPlot, r, bPhi, $
            view_grid = [3,1]
	    iPlot, r, specData[0].n, $
            /view_next
	    iPlot, r, specData[1].n, $
            /view_next
    endif

;	Dispersion analysis

    if not keyword_set ( wReal ) then $
	    wReal	= 30d6 * 2d0 * !dpi
	w	= complex ( wReal, wReal * 0.00 )

	nPhi	= 5 
	kPar	= nPhi / r
	kz		= 0 
	nPar	= kPar * c / wReal

	stixR	= 1d0 - total ( specData.wp^2 / ( wReal * ( wReal + specData.wc ) ), 2 )
	stixL	= 1d0 - total ( specData.wp^2 / ( wReal * ( wReal - specData.wc ) ), 2 )
	stixS	= 0.5d0 * ( stixR + stixL )
	stixD	= 0.5d0 * ( stixR - stixL )
	stixP	= 1d0 - total ( specData.wp^2 / wReal^2, 2 )

	stixR_	= 1d0 - total ( specData.wp_^2 / ( wReal * ( wReal + specData.wc_ ) ), 2 )
	stixL_	= 1d0 - total ( specData.wp_^2 / ( wReal * ( wReal - specData.wc_ ) ), 2 )
	stixS_	= 0.5d0 * ( stixR_ + stixL_ )
	stixD_	= 0.5d0 * ( stixR_ - stixL_ )
	stixP_	= 1d0 - total ( specData.wp_^2 / wReal^2, 2 )


	AA	= stixS
	BB	= -1d0 * ( stixR * stixL + stixP * stixS - nPar^2 * ( stixP + stixS ) )
	CC	= stixP * ( nPar^2 - stixR ) * ( nPar^2 - stixL )
	B24AC	= BB^2 - 4d0 * AA * CC
	kPerp1	= sqrt ( ( -BB + sqrt ( complex ( B24AC, B24AC * 0 ) ) ) / ( 2d0 * AA ) ) * wReal / c
	kPerp2	= sqrt ( ( -BB - sqrt ( complex ( B24AC, B24AC * 0 ) ) ) / ( 2d0 * AA ) ) * wReal / c

    kR1  = sqrt ( kPerp1^2 - kz^2 )
    kR2  = sqrt ( kPerp2^2 - kz^2 )


    if not keyword_set ( noPlot ) then begin
    loadct, 12, /sil, rgb_table = ct12
    device, decomposed = 0

	    iPlot, r, real_part ( kPerp1 ), $
            /yLog, $
            yRange = [1, 1e3], $
            view_grid = [2,1]
	    iPlot, r, imaginary ( kPerp1 ), $
	    	lineStyle = 1, $
            /over
	    iPlot, r, real_part ( kR1 ), $
            /over, $
            color = transpose ( ct12[8*16-1,*] )
	    iPlot, r, imaginary ( kR1 ), $
	    	lineStyle = 1, $
            /over, $
            color = transpose ( ct12[8*16-1,*] )


	    iPlot, r, real_part ( kPerp2 ), $
            /yLog, $
            yRange = [1, 1e3], $
            /view_next
	    iPlot, r, imaginary ( kPerp2 ), $
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

;	Calculate dielectric tensor

	II	= dcomplex ( 0, 1 )

	epsilon		= dcomplexArr ( 3, 3, nR )
	epsilon_	= dcomplexArr ( 3, 3, nR-1 )

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

;	Build matrix

	nAll	= nR + 2 * ( nR - 1 )
	aMat	= complexArr ( nAll, nAll )
	rhs		= complexArr ( nAll )

	;	matrix fill

    print, '*** filling matrix'
	for i = 0, nR - 1 do begin

		;	r component
            if i gt 0 then begin
				aMat[3*i-2,3*i]	= II * nPhi / (2 * r[i]^2) - II * nPhi / ( r[i] * dr ) $
                                    - w^2/(2*c^2)*epsilon_[1,0,i-1]
				aMat[3*i-1,3*i]	= -II * kz / dr - w^2 / (2 * c^2) * epsilon_[2,0,i-1]
            endif
				aMat[3*i,3*i]	= nPhi / r[i]^2 + kz^2 - w^2 / c^2 * epsilon[0,0,i]
            if i lt nR-1 then begin
				aMat[3*i+1,3*i]	= II * nPhi / ( r[i] * dr ) + II * nPhi / ( 2 * r[i]^2 ) $
                                    - w^2/(2*c^2)*epsilon_[1,0,i]
				aMat[3*i+2,3*i]	= II * kz / dr - w^2 / (2*c^2) * epsilon_[2,0,i]
            endif

		;	phi component
			if i lt nR-1 then begin
                if i gt 0 then $
				aMat[3*i-2,3*i+1] = 1 / ( 2 * r_[i] * dr ) - 1/ dr^2
				aMat[3*i,3*i+1]	= -II*nPhi/(r_[i]*dr) $
                                    - II*nPhi/(2*r_[i]^2)-w^2/(2*c^2)*epsilon[0,1,i]
				aMat[3*i+1,3*i+1]	= 2 / dr^2 + kz^2 + 1 / r_[i]^2 -w^2/c^2*epsilon_[1,1,i]
				aMat[3*i+2,3*i+1]	= -kz*nPhi/r_[i] - w^2/(c^2)*epsilon_[2,1,i]
				aMat[3*i+3,3*i+1]	= II * nPhi / (r_[i]*dr) - II*nPhi/(2*r_[i]^2) $
                                        - w^2/(2*c^2)*epsilon[0,1,i+1]
				if i lt nR-2 then $
				aMat[3*i+4,3*i+1]	= -1/dr^2 - 1/(2*r_[i]*dr)
			endif

		;	z component	
			if i lt nR-1 then begin
                if i gt 0 then $
				aMat[3*i-1,3*i+2] = 1/(2*r_[i]*dr)-1/dr^2
				aMat[3*i,3*i+2] = II*kz/(2*r_[i])-II*kz/dr-w^2/(2*c^2)*epsilon[0,2,i]
				aMat[3*i+1,3*i+2] = -nPhi*kz/r_[i] - w^2/c^2*epsilon_[1,2,i]
				aMat[3*i+2,3*i+2] = 2/dr^2 + nPhi^2/r_[i]^2 - w^2/c^2*epsilon_[2,2,i]
				aMat[3*i+3,3*i+2] = II * kz / dr + II * kz / (2*r_[i]) -w^2/(2*c^2)*epsilon[0,2,i+1]
				if i lt nR-2 then $
				aMat[3*i+5,3*i+2] = -1/dr^2 -1/(2*r_[i]*dr)
			endif

	endfor

	rhs[3*(nR-20)+2]	= II * w * u0 * r_[nR-20] * dr * 1d0

;	Solve matrix

    print, '*** solving linear system', (nAll*2.0)^2*16.0/(1024.0^2)
	eField	= la_linear_equation ( aMat, rhs, status = stat, /double )
	print, 'lapack status: ', stat

	eFieldTmp	= dcomplexArr ( nAll )
	eFieldTmp	= eField
	eField	= temporary ( eFieldTmp )

	ii_eR	= indGen(nR)*3
	ii_ePhi	= indGen(nR-1)*3+1
	ii_ez	= indGen(nR-1)*3+2
	eR	= eField[ii_eR]
	ePhi	= eField[ii_ePhi]
	ez	= eField[ii_ez]

;   Calculate the Div of D @ the z,phi grid pts
;   but first we need to invert epsilon to get D & D_

    D   = dcomplexArr ( nR, 3 )
    D_  = dcomplexArr ( nR-1, 3 )

    print, '*** calculating D'
    for i=0,nR-1 do begin
       
        if i eq 0 then begin

            D[i,0]   = epsilon[0,0,i] * eR[i] $
                        + ( epsilon_[1,0,i] * ePhi[i] ) / 2 $
                        + ( epsilon_[2,0,i] * ez[i] ) / 2
            D[i,1]   = epsilon[0,1,i] * eR[i] $
                        + ( epsilon_[1,1,i] * ePhi[i] ) / 2 $
                        + ( epsilon_[2,1,i] * ez[i] ) / 2
            D[i,2]   = epsilon[0,2,i] * eR[i] $
                        + ( epsilon_[1,2,i] * ePhi[i] ) / 2 $
                        + ( epsilon_[2,2,i] * ez[i] ) / 2

        endif else if i eq nR-1 then begin

            D[i,0]   = epsilon[0,0,i] * eR[i] $
                        + ( epsilon_[1,0,i-1] * ePhi[i-1] ) / 2 $
                        + ( epsilon_[2,0,i-1] * ez[i-1] ) / 2
            D[i,1]   = epsilon[0,1,i] * eR[i] $
                        + ( epsilon_[1,1,i-1] * ePhi[i-1] ) / 2 $
                        + ( epsilon_[2,1,i-1] * ez[i-1] ) / 2
            D[i,2]   = epsilon[0,2,i] * eR[i] $
                        + ( epsilon_[1,2,i-1] * ePhi[i-1] ) / 2 $
                        + ( epsilon_[2,2,i-1] * ez[i-1] ) / 2

        endif else begin
    
            D[i,0]   = epsilon[0,0,i] * eR[i] $
                        + ( epsilon_[1,0,i-1] * ePhi[i-1] + epsilon_[1,0,i] * ePhi[i] ) / 2 $
                        + ( epsilon_[2,0,i-1] * ez[i-1] + epsilon_[2,0,i] * ez[i] ) / 2
            D[i,1]   = epsilon[0,1,i] * eR[i] $
                        + ( epsilon_[1,1,i-1] * ePhi[i-1] + epsilon_[1,1,i] * ePhi[i] ) / 2 $
                        + ( epsilon_[2,1,i-1] * ez[i-1] + epsilon_[2,1,i] * ez[i] ) / 2
            D[i,2]   = epsilon[0,2,i] * eR[i] $
                        + ( epsilon_[1,2,i-1] * ePhi[i-1] + epsilon_[1,2,i] * ePhi[i] ) / 2 $
                        + ( epsilon_[2,2,i-1] * ez[i-1] + epsilon_[2,2,i] * ez[i] ) / 2

        endelse

        if i lt nR-1 then begin
 
            D_[i,0]   = ( epsilon[0,0,i] * eR[i] + epsilon[0,0,i+1] * eR[i+1] ) / 2 $
                        + epsilon_[1,0,i] * ePhi[i] $
                        + epsilon_[2,0,i] * ez[i]  
            D_[i,1]   = ( epsilon[0,1,i] * eR[i] + epsilon[0,1,i+1] * eR[i+1] ) / 2 $
                        + epsilon_[1,1,i] * ePhi[i] $
                        + epsilon_[2,1,i] * ez[i]  
            D_[i,2]   = ( epsilon[0,2,i] * eR[i] + epsilon[0,2,i+1] * eR[i+1] ) / 2 $
                        + epsilon_[1,2,i] * ePhi[i] $
                        + epsilon_[2,2,i] * ez[i]  

        endif
 
    endfor

    D   = e0 * D
    D_  = e0 * D_

    divD_   = dcomplexArr ( nR - 1 )

    print, '*** calculating div D'
    for i=0,nR-2 do begin 
        divD_[i]    = D_[i,0] / r_[i] $
                        + II * nPhi * D_[i,1] / r_[i] $
                        + II * kz * D_[i,2] $
                        + D[i+1,0] / dr $
                        - D[i,0] / dr
    endfor

    ;   and @ the r grid pts

    divD    = dcomplexArr ( nR )

    for i=1,nR-2 do begin 
        divD[i]    = D[i,0] / r[i] $
                        + II * nPhi * D[i,1] / r[i] $
                        + II * kz * D[i,2] $
                        - D[i-1,0] / ( 2 * dr ) $
                        + D[i,0] / ( 2 * dr )
    endfor

    divE_   = dcomplexArr ( nR - 1 )
    for i=0,nR-2 do begin
        divE_[i]    = II * nPhi * ePhi[i] / r_[i] $
                        + II * kz * ez[i] $
                        + eR[i+1] * ( 0.5 + r_[i] / dr ) / r_[i] $
                        + eR[i] * ( 0.5 - r_[i] / dr ) / r_[i] 
    endfor

    rho_ = e0 * divE_


    if not keyword_set ( noPlot ) then begin

    ;   plot the numerical and real charge density
    
        iPlot, r_, divD_, $
            /zoom_on_resize
        iPlot, r, divD, $
            /over, $
            color = transpose ( ct12[12*16-1,*] )
        iPlot, r_, imaginary ( divD_ ), $
            /over, $
            lineStyle = 1
        iPlot, r, imaginary ( divD ), $
            /over, $
            color = transpose ( ct12[12*16-1,*] ), $
            lineStyle = 1
        iPlot, r_, rho_, $
            color = transpose ( ct12[8*16-1,*] ), $
            /over
    
    ;	Visualise solution
    
        loadct, 12, /sil
    	!p.charSize = 3
    
        iPlot, r, eR, $
            view_grid = [3,1], $
            yTickFont_size = 30, $
            /zoom_on_resize
        iPlot, r, imaginary ( eR ), $
            /over, $
            color = transpose ( ct12[8*16-1,*] )
    	iPlot, [r_[0]-dr,r_,r_[nR-2]+dr], [0,ePhi,0], $
            psym = -4, /view_next
    	iPlot, r_, imaginary(ePhi), $
            psym = -4, $
            /over, $
            color = transpose ( ct12[8*16-1,*] )
    	iPlot, [r_[0]-dr,r_,r_[nR-2]+dr], [0,ez,0], $
            psym = -4, /view_next
    	iPlot, r_, imaginary(ez), $
            psym = -4, $
            /over, $
            color = transpose ( ct12[8*16-1,*] )
    endif

end
