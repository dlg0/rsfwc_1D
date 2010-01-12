pro dielectric, runData, stixVars, $
	epsilonFull = epsilon, $
	epsilonHalf = epsilon_

	common switches
	common constants

;	Calculate dielectric tensor


	epsilon		= dcomplexArr ( 3, 3, runData.nR )
	epsilon_	= dcomplexArr ( 3, 3, runData.nR-1 )

    if dielectric_freeSpace then begin

        epsilon[0,0,*]  = complex ( 1, 0 )
        epsilon[1,1,*]  = complex ( 1, 0 )
        epsilon[2,2,*]  = complex ( 1, 0 )

        epsilon_[0,0,*]  = complex ( 1, 0 )
        epsilon_[1,1,*]  = complex ( 1, 0 )
        epsilon_[2,2,*]  = complex ( 1, 0 )

    endif else if dielectric_noPoloidal then begin

	    epsilon[0,0,*]	= stixVars.stixS
	    epsilon[0,2,*]	= -II * stixVars.stixD
	    epsilon[1,1,*]	= stixVars.stixP
	    epsilon[2,0,*]	= II * stixVars.stixD
	    epsilon[2,2,*]	= stixVars.stixS

	    epsilon_[0,0,*]	= stixVars.stixS_
	    epsilon_[0,2,*]	= -II * stixVars.stixD_
	    epsilon_[1,1,*]	= stixVars.stixP_
	    epsilon_[2,0,*]	= II * stixVars.stixD_
	    epsilon_[2,2,*]	= stixVars.stixS_

		;epsilonNoPol	= epsilon
		;epsilonNoPol_	= epsilon_

    endif else begin

		;	Generic dielectric for arbitrary magnetic field direction
	
		bUnit_cyl	= [ [ runData.bField[0] / runData.bMag ], $	
						[ runData.bField[1] / runData.bMag ], $
						[ runData.bField[2] / runData.bMag ] ]
		bUnit_cyl_	= [ [ runData.bField_[0] / runData.bMag_ ], $	
						[ runData.bField_[1] / runData.bMag_ ], $
						[ runData.bField_[2] / runData.bMag_ ] ]
	
		;	rotate to cartesian
	
		phi	= 0d0
		bUnit_car	= dblArr ( runData.nR, 3 )
		bUnit_car_	= dblArr ( runData.nR-1, 3 )
	
		for i = 0L, runData.nR - 1L do begin
	
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
	
			epsilon_stix	= [ [ stixVars.stixS[i], II * stixVars.stixD[i], 0d0 ], $
								[ -II * stixVars.stixD[i], stixVars.stixS[i], 0d0 ], $
								[ 0d0, 0d0, stixVars.stixP[i] ] ]
	
			;;	rotate to cylindrical
			;epsilon[*,*,i]	=	car2cyl ## epsilon_car ## cyl2car
	
			epsilon[*,*,i]	= rotateEpsilon ( epsilon_stix, bUnit_car[i,*] )
	
			;	same for 1/2 grid
	
			if i lt runData.nR - 1 then begin
	
				bUnit_car_[i,*]	= bUnit_cyl_[i,*]	
		
				epsilon_stix_	= [ [ stixVars.stixS_[i], II * stixVars.stixD_[i], 0d0 ], $
									[ -II * stixVars.stixD_[i], stixVars.stixS_[i], 0d0 ], $
									[ 0d0, 0d0, stixVars.stixP_[i] ] ]
	
				epsilon_[*,*,i]	= rotateEpsilon ( epsilon_stix_, bUnit_car_[i,*] )
	
	
			endif
		
		endfor

	endelse

;	Check epsilon for NaNs

	iiNaN	= where ( epsilon ne epsilon, iiNaNCnt )
	if iiNaNCnt gt 0 then begin

		print, '______________________________'
		print, 'ERROR: NaN detected in epsilon'
		iiNeg	= where ( runData.specData.n le 0, iiNegCnt )
		if iiNegCnt gt 0 then begin
			
			print, '_____________________________'	
			print, 'NEGATIVE DENSITY DETECTED :-('

		endif
		stop	
	endif

end
