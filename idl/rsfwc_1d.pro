;	Real Space Full Wave Code in 1D
;	DLG

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
	damping = damping
	

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
	rMin	= 0.2d0
	rMax	= 2.7d0
	dR	= ( rMax - rMin ) / ( nR - 1 )
	r	= dIndGen ( nR ) * dR + rMin
	;	half grid variables are _
	r_	= r[1:*] - dR / 2

;	Setup plasma and field profiles

    if keyword_set ( freeSpace ) then $
        bMax = 0.0 $
    else $
	    bMax	= 0.45d0

	r0	= 1

	bPhi	= dblArr(nR)+bMax / r^0.1
	bPhi_	= dblArr(nR-1)+bMax / r_^0.1

	bR		= dblArr ( nR ) 
	bR_		= dblArr ( nR-1 ) 
	bz		= dblArr ( nR ) + sin ( r - r0 ) * 0.3
	bz_		= dblArr ( nR-1 ) + sin ( r_ - r0 ) * 0.3

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
	nMax	= 1.0d17

	;	electrons

	specData[0].q 		= -1 * e
	specData[0].m 		= me
	specData[0].n 		= nMax;((-(r-nPeakR)^2+(r[0]-nPeakR)^2) / (r[0]-nPeakR)^2 * nMax)>0
	specData[0].n_		= nMax;((-(r_-nPeakR)^2+(r[0]-nPeakR)^2) / (r[0]-nPeakR)^2 * nMax)>0

	;	sinusoidal density gradient

	gradSize = 0.15
	gradStart	= 1.4
	gradFreq	= 1 / (2*gradSize)
	ne_profile	= fltArr ( nR ) 
	ii1	= where ( r lt gradStart )
	ii2	= where ( r ge gradStart and r lt gradStart+gradSize )
	ne_profile[*]	= -1
	ne_profile[ii1]	= 1
	ne_profile[ii2]	= cos ( 2*!pi*gradFreq*(r[ii2]-gradStart) )
	ne_profile	= (ne_profile + 1) * 200e17 + nMax

	specData[0].n	= ne_profile
	specData[0].n_	= ((specData[0].n)[0:nR-2] + (specData[0].n)[1:nR-1])/2
	
	specData[1].n 		= specData[0].n / 2
	specData[1].n_		= specData[0].n_ / 2

	;   helium	

	specData[1].q 		= 2 * e
	specData[1].m 		= 2 * mi 

	for i = 0, nSpec - 1 do begin

		specData[i].wp	= sqrt ( specData[i].n * specData[i].q^2 $
                            / ( specData[i].m*e0 ))
		specData[i].wc	= specData[i].q * bPhi / specData[i].m
		specData[i].wp_	= sqrt ( specData[i].n_ * specData[i].q^2 $
                            / ( specData[i].m*e0 ))
		specData[i].wc_	= specData[i].q * bPhi_ / specData[i].m

	endfor

	if keyword_set ( plot ) then begin
    	loadct, 12, /sil, rgb_table = ct12
    	device, decomposed = 0
		red	= transpose ( ct12[12*16-1,*] )
		blue	= transpose ( ct12[8*16-1,*] )
		green	= transpose ( ct12[2*16-1,*] )
	endif

   
    if keyword_set ( plot ) then begin 
	    iPlot, r, bPhi, $
            view_grid = [3,2], $
            title = 'bPhi', $
            /zoom_on_resize, $
			yRange = [-1,1]
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
	    wReal	= 1.118468d11 
	w	= dcomplex ( wReal, wReal * damping )

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

        ;part1   = -4 * stixP * ( -stixD^2 + ( nPar^2 - stixS )^2 ) * stixS $
        ;            + ( stixD^2 + ( nPar^2 - stixS ) * ( stixP + stixS ) )^2
        ;part2_1   = stixD^2 + nPar^2 * stixP + nPar^2 * stixS - stixP * stixS - stixS^2
        ;part2_2   = -stixD^2 - nPar^2 * stixP - nPar^2 * stixS + stixP * stixS + stixS^2

        ;kPer1   =  sqrt ( complex ( ( -part2_1 - sqrt ( part1 ) ) / ( 2 * stixS ), fltArr ( nR ) ) * wReal / c )
        ;kPer2   =  sqrt ( complex ( ( part2_2 + sqrt ( part1 ) ) / ( 2 * stixS ), fltArr ( nR ) ) * wReal / c )

	    AA	= stixS
	    BB	= -1d0 * ( stixR * stixL + stixP * stixS - nPar^2 * ( stixP + stixS ) )
	    CC	= stixP * ( nPar^2 - stixR ) * ( nPar^2 - stixL )
	    B24AC	= BB^2 - 4d0 * AA * CC
	    kPerp1	= sqrt ( ( -BB $
            + sqrt ( complex ( B24AC, B24AC * 0 ) ) ) / ( 2d0 * AA ) ) * wReal / c
	    kPerp2	= sqrt ( ( -BB $
            - sqrt ( complex ( B24AC, B24AC * 0 ) ) ) / ( 2d0 * AA ) ) * wReal / c

        kR1  = sqrt ( kPerp1^2 - kz^2 )
        kR2  = sqrt ( kPerp2^2 - kz^2 )

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

;	Generic dielectric for arbitrary magnetic field direction


	bMag	= sqrt ( bR^2 + bPhi^2 + bz^2 )
	bMag_	= sqrt ( bR_^2 + bPhi_^2 + bz_^2 )

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
		inv_rot_x	= invert ( rot_x )

		;	z towards x
		rotTh2	=	aTan ( bUnit_car[i,0], bUnit_car[i,2] ) 
		rot_y	= [ [ cos ( rotTh2 ), 0, -sin ( rotTh2 ) ], $
					[ 0, 1, 0 ], $
					[ sin ( rotTh2 ), 0, cos ( rotTh2 ) ] ]
		inv_rot_y	= invert ( rot_y )

		;	build dielectric in x,y,z with z assumed along b

		epsilon_stix	= [ [ stixS[i], II * stixD[i], 0 ], $
							[ -II * stixD[i], stixS[i], 0 ], $
							[ 0, 0, stixP[i] ] ]

		;	rotate

		epsilon[*,*,i]	= rot_x ## rot_y ## epsilon_stix ## inv_rot_y ## inv_rot_x
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

			epsilon_[*,*,i]	= rot_x_ ## rot_y_ ## epsilon_stix_ ## inv_rot_y_ ## inv_rot_x_

		endif
	
	endfor
stop

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

    antLoc  = 1.7d0
    iiAnt   = where ( abs ( r_ - antLoc ) eq min ( abs ( r_ - antLoc ) ) )
    
	rhs[iiAnt*3+2]	= II * w * u0 / ( r_[0] * dr ) * 1d0
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
    
        iPlot, r, eR, $
            view_number = 4, $
            title = 'eR'
        iPlot, r, imaginary ( eR ), $
            /over, $
            color = transpose ( ct12[8*16-1,*] )
		iPlot, [antLoc,antLoc], [-max(abs(eR)),max(abs(eR))], $
			transparency = 80, $
			thick = 10, $
			/over
	
    	iPlot, [r_[0]-dr,r_,r_[nR-2]+dr], [0,ePhi,0], $
            psym = -4, $
            view_number = 5, $
            title = 'ePhi'
    	iPlot, r_, imaginary(ePhi), $
            psym = -4, $
            /over, $
            color = transpose ( ct12[8*16-1,*] )
		iPlot, [antLoc,antLoc], [-max(abs(ePhi)),max(abs(ePhi))], $
			transparency = 80, $
			thick = 10, $
			/over
	
    	iPlot, [r_[0]-dr,r_,r_[nR-2]+dr], [0,ez,0], $
            psym = -4, $
            view_number = 6, $
            title = 'ez'
    	iPlot, [r_[0]-dr,r_,r_[nR-2]+dr], [0,imaginary(ez),0], $
            psym = -4, $
            /over, $
            color = transpose ( ct12[8*16-1,*] )
		iPlot, [antLoc,antLoc], [-max(abs(ez)),max(abs(ez))], $
			transparency = 80, $
			thick = 10, $
			/over
	
		print, 'FFTing solution ...'
		nFFT	= 128 
		fftData_eR	= complexArr ( nFFT )
		fftData_ePhi	= complexArr ( nFFT )
		fftData_ez	= complexArr ( nFFT )
		rFFTData	= fltArr ( nFFT )
		fftStep	= nR / nFFT 
		fftWidth	= nR / 5 
		fftWidth	+= (fftWidth mod 2)
		osf	= 10

		fftData2d_eR	= complexArr ( nFFT, fftWidth*osf )
		hannWin	= hanning ( fftWidth )
		kRFFT	= findgen ( fftWidth * osf ) / ( fftWidth * osf * dr ) * ( 2 * !pi )
		eRTmp2	= complexArr ( fftWidth * osf )

		for i=0,nFFT-1 do begin
		
			iiLow	= i*fftStep-fftWidth/2
			iiHig	= i*fftStep+fftWidth/2-1

			rFFTData[i]	= r[i*fftStep]

			if iiLow gt 0 and iiHig lt nR-1 then begin

				eRTmp	= eR[iiLow:iiHig] * hannWin
				eRTmp2[fftWidth*osf/2-fftWidth/2:fftWidth*osf/2+fftWidth/2-1]	= eRTmp
				fftTmp	= abs ( fft ( eRTmp2 ) )^2
				fftData2D_eR[i,*]	= fftTmp

			endif

		endfor
		print, 'DONE'
	
		nLevs	= 11	
		levels	= 10.0^fIndGen ( nLevs )/1e3
		colors	= 255 - ( bytScl ( fIndGen(nLevs), top = 253 ) + 1 )
		iiPlotk	= where ( kRFFT gt 0 and kRFFT le 1e3 )
		iContour, fftData2D_eR[*,iiPlotk], rFFTData, kRFFT[iiPlotk], $
			/yLog, $
			yRange = [1,1e3], $
			view_number = 3, $
			xRange = [min(r),max(r)], $
			rgb_indices = colors, $
			rgb_table = 1, $
			/fill, $
			c_value = levels

		iPlot, r, real_part ( kR1 ), $
            /yLog, $
            /over, $
            title = 'dispersion', $
			xTitle = 'R[m]', $
			yTitle = 'kR', $
			xRange = [min(r),max(r)], $
			thick = 2, $
			trans	= 50

	    iPlot, r, imaginary ( kR1 ), $
	    	lineStyle = 1, $
            /over, $
			thick = 2, $
			trans = 50

	    iPlot, r, real_part ( kR2 ), $
            /over, $
            color = transpose ( ct12[12*16-1,*] ), $
			thick = 2, $
			trans = 50
	    iPlot, r, imaginary ( kR2 ), $
	    	lineStyle = 1, $
            /over, $
            color = transpose ( ct12[12*16-1,*] ), $
			thick = 2, $
			trans = 50

		iPlot, [antLoc,antLoc], [1,1e3], $
			transparency = 80, $
			thick = 10, $
			/over
	

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

end
