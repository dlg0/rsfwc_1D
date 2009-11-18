;	Real Space Full Wave Code in 1D
;	DLG

pro rsfwc_1d

;	Parameters

	e	= 1.60217646d-19
	me	= 9.10938188d-31
   	mi	= 1.67262158d-27
	c	= 299792458d0	
	e0	=	8.854187817d-12 
	u0	=	4d0 * !dpi * 10d0^(-7) 

;	Setup Grid

	nR		= 256 
	rMin	= 1.0
	rMax	= 1.2
	dR	= ( rMax - rMin ) / ( nR - 1 )
	r	= fIndGen ( nR ) * dR + rMin
	;	half grid variables are _
	r_	= r[1:*] - dR / 2

;	Setup plasma and field profiles

	bMax	= 0.6
	bPhi	= fltArr(nR)+bMax; / r 
	bPhi_	= fltArr(nR-1)+bMax; / r_

	nSpec	= 2
	specData	= replicate ( $
			{ 	q : 0.0, $
				m : 0.0, $
				wp : fltArr ( nR ), $
				wp_ : fltArr ( nR - 1 ), $
				wc : fltArr ( nR ), $
				wc_ : fltArr ( nR - 1 ), $
				n : fltArr ( nR ), $
				n_ : fltArr ( nR - 1 ) }, nSpec )

	nPeakR	= 1.1
	nMax	= 0.6e19

	;	electrons

	specData[0].q 		= -1 * e
	specData[0].m 		= me
	specData.n 		= (-(r-nPeakR)^2+(r[0]-nPeakR)^2) / (r[0]-1.1)^2 * nMax
	specData.n_		= (-(r_-nPeakR)^2+(r[0]-nPeakR)^2) / (r[0]-1.1)^2 * nMax

	;	dueterium

	specData[1].q 		= 1 * e
	specData[1].m 		= 2 * mi 

	for i = 0, nSpec - 1 do begin

		specData[i].wp		= sqrt ( specData[i].n * specData[i].q^2 / ( specData[i].m * e0 ) )
		specData[i].wc		= specData[i].q * bPhi / specData[i].m
		specData[i].wp_		= sqrt ( specData[i].n_ * specData[i].q^2 / ( specData[i].m * e0 ) )
		specData[i].wc_		= specData[i].q * bPhi_ / specData[i].m

	endfor

	winNo	= 0
	;window, winNo
	;!p.multi = [0,2,2]
    
	iPlot, r, bPhi, $
        view_grid = [3,1]
	iPlot, r, specData[0].n, $
        /view_next
	iPlot, r, specData[1].n, $
        /view_next

;	Dispersion analysis

	wReal	= 30d6 * 2d0 * !dpi
	w	= complex ( wReal, wReal * 0.1 )

	nPhi	= 15 
	kPar	= nPhi; / r
	kz		= 15 
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

	;++winNo
	;window, winNo
	;!p.multi = [0,2,2]
	iPlot, r, real_part ( kPerp1 ), $
        /yLog, $
        yRange = [1, 1e3], $
        view_grid = [2,1]
	iPlot, r, imaginary ( kPerp1 ), $
		lineStyle = 1, $
        /over
	iPlot, r, real_part ( kPerp2 ), $
        /yLog, $
        yRange = [1, 1e3], $
        /view_next
	iPlot, r, imaginary ( kPerp2 ), $
		lineStyle = 1, $
        /over

	!p.multi = 0

;	Calculate dielectric tensor

	II	= complex ( 0, 1 )

	epsilon		= fltArr ( 3, 3, nR )
	epsilon_	= fltArr ( 3, 3, nR-1 )

	epsilon[0,0,*]	= stixS
	epsilon[0,2,*]	= II * stixD
	epsilon[1,1,*]	= stixP
	epsilon[0,2,*]	= -II * stixD
	epsilon[2,2,*]	= stixS

	epsilon_[0,0,*]	= stixS_
	epsilon_[0,2,*]	= II * stixD_
	epsilon_[1,1,*]	= stixP_
	epsilon_[0,2,*]	= -II * stixD_
	epsilon_[2,2,*]	= stixS_

;	Build matrix

	nAll	= nR + 2 * ( nR - 1 )
	aMat	= complexArr ( nAll, nAll )
	rhs		= complexArr ( nAll )

	;	matrix fill

	for i = 0, nR - 1 do begin

		;	forward difference pts

			if i eq 0 then begin

				;	r component

				;	aMat[3*i,3*i]	= kz^2 + nPhi^2/r[i]^2 - w^2/c^2*epsilon[0,0,i] 
				;	aMat[3*i+1,3*i]	= 7 * II * nPhi / ( 2 * r[i] * dr )
				;	aMat[3*i+2,3*i]	= 7 * II * kz / ( 2 * dr )
				;	aMat[3*i+3,3*i]	= 0
				;	aMat[3*i+4,3*i]	= -II * nPHi / ( 2 * r[i] * dr )
				;	aMat[3*i+5,3*i]	= -II * kz / ( 2 * dr )

					;aMat[3*i,3*i]	= kz^2 + nPhi^2/r[i]^2 - w^2/c^2*epsilon[0,0,i] 
					;aMat[3*i+1,3*i]	= II * nPhi * ( 2 * r[i] + dr ) / ( 2 * r[i]^2 + dr ) - w^2 * epsilon_[1,0,i] / ( 2 * c^2 )
					;aMat[3*i+2,3*i]	= II * kz / dr - w^2 * epsilon_[2,0,i] / ( 2 * c^2 )
					;aMat[3*i+3,3*i]	= 0
					;aMat[3*i+4,3*i]	= 0 
					;aMat[3*i+5,3*i]	= 0 


			;	;	phi component

			;		aMat[3*i,3*i+1] 	= -II * nPhi * ( 2 * r_[i] + dr ) / ( 2 * r_[i]^2 * dr ) - w^2 * epsilon[0,1,i] / ( 2 * c^2 ) 
			;		aMat[3*i+1,3*i+1]	= kz^2 + 1/r_[i]^2 - 2 / dr^2 + 3 / ( 2 * r_[i] * dr ) - w^2 * epsilon_[1,1,i] / c^2
			;		aMat[3*i+2,3*i+1]	= -kz * nPhi / r_[i] - w^2 * epsilon_[2,1,i] / c^2
			;		aMat[3*i+3,3*i+1]	= -II * nPhi / ( 2 * r_[i]^2 ) + II * nPhi / ( r_[i] * dr ) - w^2 * epsilon[0,1,i+1] / ( 2 * c^2 )
			;		aMat[3*i+4,3*i+1]	= ( 5 * r_[i] - 2 * dr ) / ( r_[i] * dr^2 )
			;		aMat[3*i+5,3*i+1]	= 0
			;	   	aMat[3*i+6,3*i+1]	= 0
			;		aMat[3*i+7,3*i+1]	= ( -8 * r_[i] + dr ) / ( 2 * r_[i] * dr^2 )
			;		aMat[3*i+8,3*i+1]	= 0
			;	   	aMat[3*i+9,3*i+1]	= 0
			;		aMat[3*i+10,3*i+1]	= 1 / dr^2 

					;;aMat[3*i+1-3,3*i+1] = 1 / ( 2 * r_[i] * dr ) - 1/ dr^2
					;aMat[3*i,3*i+1]	= -II*nPhi/(r_[i]*dr) -II*nPhi/(2*r_[i]^2)-w^2/(2*c^2)*epsilon[0,1,i]
					;aMat[3*i+1,3*i+1]	= 2 / dr^2 + kz^2 + 1 / r_[i]^2 -w^2/c^2*epsilon_[1,1,i]
					;aMat[3*i+2,3*i+1]	= -kz*nPhi/r_[i] - w^2/(c^2)*epsilon[2,1,i]
					;aMat[3*i+3,3*i+1]	= II * nPhi / (r_[i]*dr) - II*nPhi/(2*r_[i]^2) - w^2/(2*c^2)*epsilon[0,1,i]
					;aMat[3*i+4,3*i+1]	= -1/dr^2 - 1/(2*r_[i]*dr)
	
			;	;	z component

			;		aMat[3*i,3*i+2] 	= II * kz / 2 * ( 1 / r_[i] - 2 / dr ) - w^2 * epsilon[0,2,i] / ( 2 * c^2 )  
			;		aMat[3*i+1,3*i+2]	= -kz * nPhi / r_[i] - w^2 * epsilon_[1,2,i] / c^2 
			;		aMat[3*i+2,3*i+2]	= nPhi^2 / r_[i]^2 - 2 / dr^2 + 3 / ( 2 * r_[i] * dr ) - w^2 * epsilon_[2,2,i] / c^2 
			;		aMat[3*i+3,3*i+2]	= II * kz / 2 * ( 1 / r_[i] + 2 / dr ) - w^2 * epsilon[0,2,i] / ( 2 * c^2 ) 
			;		aMat[3*i+4,3*i+2]	= 0 
			;		aMat[3*i+5,3*i+2]	= ( 5 * r_[i] - 2 * dr ) / ( r_[i] * dr^2 ) 
			;	   	aMat[3*i+6,3*i+2]	= 0 
			;		aMat[3*i+7,3*i+2]	= 0
			;		aMat[3*i+8,3*i+2]	= ( -8 * r_[i] + dr ) / ( 2 * r_[i] * dr^2 ) 
			;	   	aMat[3*i+9,3*i+2]	= 0
			;		aMat[3*i+10,3*i+2]	= 0 
			;		aMat[3*i+11,3*i+2]	= 1 / dr^2

					;;aMat[3*i-1,3*i+2] = 1/(2*r_[i]*dr)-1/dr^2
					;aMat[3*i,3*i+2] = II*kz/(2*r_[i])-II*kz/dr-w^2/(2*c^2)*epsilon[0,2,i]
					;aMat[3*i+1,3*i+2] = -nPhi*kz/r_[i] - w^2/c^2*epsilon[1,2,i]
					;aMat[3*i+2,3*i+2] = 2/dr^2 + nPhi^2/r_[i]^2 - w^2/c^2*epsilon[2,2,i]
					;aMat[3*i+3,3*i+2] = II * kz / dr + II * kz / (2*r_[i]) -w^2/(2*c^2)*epsilon[0,2,i+1]
					;aMat[3*i+5,3*i+2] = -1/dr^2 -1/(2*r_[i]*dr)
	
			endif

		;	r component
            if i gt 0 then begin
				aMat[3*i-2,3*i]	= II * nPhi / (2 * r[i]^2) - II * nPhi / ( r[i] * dr ) - w^2/(2*c^2)*epsilon_[1,0,i-1]
				aMat[3*i-1,3*i]	= -II * kz / dr - w^2 / (2 * c^2) * epsilon_[2,0,i-1]
            endif
				aMat[3*i,3*i]	= nPhi / r[i]^2 + kz^2 - w^2 / c^2 * epsilon[0,0,i]
            if i lt nR-1 then begin
				aMat[3*i+1,3*i]	= II * nPhi / ( r[i] * dr ) + II * nPhi / ( 2 * r[i]^2 ) - w^2/(2*c^2)*epsilon_[1,0,i]
				aMat[3*i+2,3*i]	= II * kz / dr - w^2 / (2*c^2) * epsilon_[2,0,i]
            endif

		;	phi component
			if i lt nR-1 then begin
                if i gt 0 then $
				aMat[3*i-2,3*i+1] = 1 / ( 2 * r_[i] * dr ) - 1/ dr^2
				aMat[3*i,3*i+1]	= -II*nPhi/(r_[i]*dr) -II*nPhi/(2*r_[i]^2)-w^2/(2*c^2)*epsilon[0,1,i]
				aMat[3*i+1,3*i+1]	= 2 / dr^2 + kz^2 + 1 / r_[i]^2 -w^2/c^2*epsilon_[1,1,i]
				aMat[3*i+2,3*i+1]	= -kz*nPhi/r_[i] - w^2/(c^2)*epsilon[2,1,i]
				aMat[3*i+3,3*i+1]	= II * nPhi / (r_[i]*dr) - II*nPhi/(2*r_[i]^2) - w^2/(2*c^2)*epsilon[0,1,i]
				if i lt nR-2 then $
				aMat[3*i+4,3*i+1]	= -1/dr^2 - 1/(2*r_[i]*dr)
			endif

		;	z component	
			if i lt nR-1 then begin
                if i gt 0 then $
				aMat[3*i-1,3*i+2] = 1/(2*r_[i]*dr)-1/dr^2
				aMat[3*i,3*i+2] = II*kz/(2*r_[i])-II*kz/dr-w^2/(2*c^2)*epsilon[0,2,i]
				aMat[3*i+1,3*i+2] = -nPhi*kz/r_[i] - w^2/c^2*epsilon[1,2,i]
				aMat[3*i+2,3*i+2] = 2/dr^2 + nPhi^2/r_[i]^2 - w^2/c^2*epsilon[2,2,i]
				aMat[3*i+3,3*i+2] = II * kz / dr + II * kz / (2*r_[i]) -w^2/(2*c^2)*epsilon[0,2,i+1]
				if i lt nR-2 then $
				aMat[3*i+5,3*i+2] = -1/dr^2 -1/(2*r_[i]*dr)
			endif

			if i eq nR-2 then begin

				;;	phi component
				;aMat[3*i+3,3*i+1]	= -II * nPhi / ( 2 * r_[i]^2 ) + II * nPhi / ( r_[i] * dr ) - w^2 * epsilon[0,1,i+1] / ( 2 * c^2 ) 
				;aMat[3*i+2,3*i+1]	= -kz * nPhi / r_[i] - w^2 * epsilon_[2,1,i] / c^2
				;aMat[3*i+1,3*i+1]	= kz^2 + 1 / r_[i]^2 - 2 / dr^2 - 3 / ( 2 * r_[i] * dr ) - w^2 * epsilon_[1,1,i] / c^2
				;aMat[3*i+0,3*i+1]	= -II * nPhi * ( 2 * r_[i] + dr ) / ( 2 * r_[i]^2 * dr ) - w^2 * epsilon_[0,1,i] / ( 2 * c^2 ) 
				;aMat[3*i-1,3*i+1]	= 0 
				;aMat[3*i-2,3*i+1]	= ( 5 * r_[i] + 2 * dr ) / ( r_[i] * dr^2 ) 
				;aMat[3*i-3,3*i+1]	= 0
				;aMat[3*i-4,3*i+1]	= 0 
				;aMat[3*i-5,3*i+1]	= -1 * ( 8 + dr / r_[i] ) / ( 2 * dr^2 ) 
				;aMat[3*i-6,3*i+1]	= 0 
				;aMat[3*i-7,3*i+1]	= 0 
				;aMat[3*i-8,3*i+1]	= 1 / dr^2 
				
				;;	z component
				;aMat[3*i+3,3*i+2]	= II * kz / 2 * ( 1 / r_[i] + 2 / dr ) - w^2 * epsilon[0,2,i+1] / ( 2 * c^2 )  
				;aMat[3*i+2,3*i+2]	= nPhi^2 / r_[i]^2 - 2 / dr^2 - 3 / ( 2 * r_[i] * dr ) - w^2 * epsilon_[2,2,i] / c^2
				;aMat[3*i+1,3*i+2]	= -kz * nPhi / r_[i] - w^2 * epsilon[1,2,i] / c^2 
				;aMat[3*i+0,3*i+2]	= II * kz / 2 * ( 1 / r_[i] - 2 / dr ) - w^2 * epsilon[0,2,i] / ( 2 * c^2 ) 
				;aMat[3*i-1,3*i+2]	= ( 5 * r_[i] + 2 * dr ) / ( r_[i] * dr^2 ) 
				;aMat[3*i-2,3*i+2]	= 0 
				;aMat[3*i-3,3*i+2]	= 0
				;aMat[3*i-4,3*i+2]	= -1 * ( 8 * r_[i] + dr ) / ( 2 * r_[i] * dr^2 ) 
				;aMat[3*i-5,3*i+2]	= 0
				;aMat[3*i-6,3*i+2]	= 0
				;aMat[3*i-7,3*i+2]	= 1 / dr^2 

			endif

			if i eq nR-1 then begin
					
				;	r component
				;aMat[3*i,3*i]	= kz^2 + nPhi^2 / r[i]^2 - w^2 * epsilon[0,0,i] / c^2
				;aMat[3*i-1,3*i]	= -9 * II * kz / ( 2 * dr ) 
				;aMat[3*i-2,3*i]	= -7 * II * nPhi / ( 2 * r[i] * dr ) 
				;aMat[3*i-3,3*i]	= 0 
				;aMat[3*i-4,3*i]	= -II * kz / ( 2 * dr ) 
				;aMat[3*i-5,3*i]	= II * nPhi / ( 2 * r[i] * dr ) 

				;aMat[3*i,3*i]	= kz^2 + nPhi^2 / r[i]^2 - w^2 * epsilon[0,0,i] / c^2
				;aMat[3*i-1,3*i]	= -II * kz / dr - w^2 * epsilon_[2,0,i-1] / ( 2 * c^2 ) 
				;aMat[3*i-2,3*i]	= II * nPhi * ( 1 - 2 * r[i] / dr ) / ( 2 * r[i]^2 ) - w^2 * epsilon_[1,0,i-1] / ( 2 * c^2 ) 
				;aMat[3*i-3,3*i]	= 0 
				;aMat[3*i-4,3*i]	= 0 
				;aMat[3*i-5,3*i]	= 0 


			endif

	endfor

	;	over determine the matrix by adding two extra equations to each end such that
	;	the boundary conditions are enforced on the phi and z components on the full 
	;	(not half) end grid points. This will then require SVD to get a least squares
	;	solution. 

	ePhiLeftBry	= complexArr ( nAll )
	ezLeftBry	= complexArr ( nAll )
	ePhiRightBry	= complexArr ( nAll )
	ezRightBry		= complexArr ( nAll )

	ePhiLeftBry[0]	= -II * nPhi * ( 3 * r[0] + 2 * dr ) / ( 2 * r[0]^2 * dr ) - w^2 * epsilon[0,1,0] / ( c^2 )
	ePhiLeftBry[1]	= ( 24 - 7 * dr / r[0] ) / ( 2 * dr^2 )
	ePhiLeftBry[2]	= 0
	ePhiLeftBry[3]	= 2 * II * nPhi / ( r[0] * dr )
	ePhiLeftBry[4]	= ( -8 * r[0] + dr ) / ( 2 * r[0] * dr^2 )
	ePhiLeftBry[5]	= 0
	ePhiLeftBry[6]	= -II * nPhi / ( 2 * r[0] * dr )

	ezLeftBry[0]	= II * kz / r[0] - 3 * II * kz / ( 2 * dr ) - w^2 * epsilon[0,2,0] / c^2  
	ezLeftBry[1]	= 0 
	ezLeftBry[2]	= ( 24 - 7 * dr / r[0] ) / ( 2 * dr^2 ) 
	ezLeftBry[3]	= 2 * II * kz / dr 
	ezLeftBry[4]	= 0 
	ezLeftBry[5]	= ( -8 * r[0] + dr ) / ( 2 * r[0] * dr^2 )
	ezLeftBry[6]	= -II * kz / ( 2 * dr ) 

	ePhiRightBry[nAll-1]	= -II * nPhi / r[nR-1]^2 +  3 * II * nPhi / ( 2 * r[nR-1] * dr ) - w^2 * epsilon[0,1,nR-1] / ( c^2 )
	ePhiRightBry[nAll-2]	= 0 
	ePhiRightBry[nAll-3]	= ( 24 + 7 * dr / r[nR-1] ) / ( 2 * dr^2 ) 
	ePhiRightBry[nAll-4]	= -2 * II * nPhi / ( r[nR-1] * dr )
	ePhiRightBry[nAll-5]	= 0 
	ePhiRightBry[nAll-6]	= -1 * ( 8 * r[nR-1] + dr ) / ( 2 * r[nR-1] * dr^2 ) 
	ePhiRightBry[nAll-7]	= II * nPhi / ( 2 * r[nR-1] * dr )

	ezRightBry[nAll-1]	= II * kz / r[nR-1] - 3 * II * kz / ( 2 * dr ) - w^2 * epsilon[0,2,nR-1] / c^2  
	ezRightBry[nAll-2]	= II * kz / ( 2 * dr ) 
	ezRightBry[nAll-3]	= 0 
	ezRightBry[nAll-4]	= -2 * II * kz / dr 
	ezRightBry[nAll-5]	= -1 * ( 8 * r[nR-1] + dr ) / ( 2 * r[nR-1] * dr^2 ) 
	ezRightBry[nAll-6]	= 0
	ezRightBry[nAll-7]	= II * kz / ( 2 * dr ) 

	rhs[3*nR*0.6+2]	= II * w * u0 * 1 

;	Solve matrix

	inner_aMat	= aMat;[0:nAll-5,0:nAll-5]
	inner_rhs	= rhs;[0:nAll-5]
	eField	= la_linear_equation ( inner_aMat, inner_rhs, status = stat )
	print, 'lapack status: ', stat

	eFieldTmp	= complexArr ( nAll )
	;eFieldTmp[0:nAll-5]	= eField
	eFieldTmp	= eField
	eField	= temporary ( eFieldTmp )

	ii_eR	= indGen(nR)*3
	ii_ePhi	= indGen(nR-1)*3+1
	ii_ez	= indGen(nR-1)*3+2
	eR	= eField[ii_eR]
	ePhi	= eField[ii_ePhi]
	ez	= eField[ii_ez]

	;	try using the least square over determined system

	OD_aMat	= [ [ePhiLeftBry], [ezLeftBry], [aMat], [ePhiRightBry], [ezRightBry] ]
	OD_rhs	= [ 0, 0, rhs, 0, 0 ]
	OD_eField	= la_least_squares ( OD_aMat, OD_rhs, $
			status = OD_stat, $
		   	method = 0 )
	print, 'OD lapack status: ', OD_stat

	ii_eR	= indGen(nR)*3+0
	ii_ePhi	= indGen(nR-1)*3+1
	ii_ez	= indGen(nR-1)*3+2
	OD_eR	= OD_eField[ii_eR]
	OD_ePhi	= OD_eField[ii_ePhi]
	OD_ez	= OD_eField[ii_ez]


;   Calculate the Div of D @ the z,phi grid pts
;   but first we need to invert epsilon to get D & D_

    D   = complexArr ( nR, 3 )
    D_  = complexArr ( nR-1, 3 )

    for i=1,nR-3 do begin
        
        D[i,0]   = epsilon[0,0,i] * eR[i] $
                    + ( epsilon_[1,0,i-1] * ePhi[i-1] + epsilon_[1,0,i] * ePhi[i] ) / 2 $
                    + ( epsilon_[2,0,i-1] * ez[i-1] + epsilon_[2,0,i] * ez[i] ) / 2
        D[i,1]   = epsilon[0,1,i] * eR[i] $
                    + ( epsilon_[1,1,i-1] * ePhi[i-1] + epsilon_[1,1,i] * ePhi[i] ) / 2 $
                    + ( epsilon_[2,1,i-1] * ez[i-1] + epsilon_[2,1,i] * ez[i] ) / 2
        D[i,2]   = epsilon[0,2,i] * eR[i] $
                    + ( epsilon_[1,2,i-1] * ePhi[i-1] + epsilon_[1,2,i] * ePhi[i] ) / 2 $
                    + ( epsilon_[2,2,i-1] * ez[i-1] + epsilon_[2,2,i] * ez[i] ) / 2
 
        D_[i,0]   = ( epsilon[0,0,i] * eR[i] + epsilon[0,0,i+1] * eR[i+1] ) / 2 $
                    + epsilon_[1,0,i] * ePhi[i] $
                    + epsilon_[2,0,i] * ez[i]  
        D_[i,1]   = ( epsilon[0,1,i] * eR[i] + epsilon[0,1,i+1] * eR[i+1] ) / 2 $
                    + epsilon_[1,1,i] * ePhi[i] $
                    + epsilon_[2,1,i] * ez[i]  
        D_[i,2]   = ( epsilon[0,2,i] * eR[i] + epsilon[0,2,i+1] * eR[i+1] ) / 2 $
                    + epsilon_[1,2,i] * ePhi[i] $
                    + epsilon_[2,2,i] * ez[i]  
 
    endfor

    D   = e0 * D
    D_  = e0 * D_

    divD_   = complexArr ( nR - 1 )

    for i=1,nR-3 do begin 
        divD_[i]    = D_[i,0] / r_[i] $
                        + II * nPhi * D_[i,1] / r_[i] $
                        + II * kz * D_[i,2] $
                        + D[i+1,0] / dr $
                        - D[i,0] / dr
    endfor

    ;   and @ the r grid pts

    divD    = complexArr ( nR )

    for i=1,nR-2 do begin 
        divD[i]    = D[i,0] / r[i] $
                        + II * nPhi * D[i,1] / r[i] $
                        + II * kz * D[i,2] $
                        - D[i-1,0] / ( 2 * dr ) $
                        + D[i,0] / ( 2 * dr )
    endfor

    loadct, 12, /sil, rgb_table = ct12
    device, decomposed = 0

    iPlot, r_, divD_
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

;	Visualise solution

	++winNo
    loadct, 12, /sil
	!p.charSize = 3

    iPlot, r, eR, $
        view_grid = [3,1]
    iPlot, r, imaginary ( eR ), $
        /over, $
        color = transpose ( ct12[8*16-1,*] )
	iPlot, [r_[0]-dr,r_,r_[nR-2]+dr], [0,ePhi,0], $
        psym = -4, /view_next
	iPlot, [r_[0]-dr,r_,r_[nR-2]+dr], imaginary([0,ePhi,0]), $
        psym = -4, $
        /over, $
        color = transpose ( ct12[8*16-1,*] )
	iPlot, [r_[0]-dr,r_,r_[nR-2]+dr], [0,ez,0], $
        psym = -4, /view_next
	iPlot, [r_[0]-dr,r_,r_[nR-2]+dr], imaginary([0,ez,0]), $
        psym = -4, $
        /over, $
        color = transpose ( ct12[8*16-1,*] )

	;++winNo
	;window, winNo, ySize = 800
	;!p.multi = [0,2,3]
	;!p.charSize = 3
	;plot, r, OD_eR
	;plot, r, imaginary (OD_eR)
	;plot, r_, OD_ePhi
	;plot, r_, imaginary (OD_ePhi)
	;plot, r_, OD_ez
	;plot, r_, imaginary(OD_ez)
stop
end
