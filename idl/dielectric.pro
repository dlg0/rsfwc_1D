pro dielectric, runData, stixVars, $
	epsilonFull_ = epsilon, $
	epsilonHalf_ = epsilon_, $
    noHalfGrid = noHalfGrid, $
    epsilonFullSpec = epsilonFullSpec, $
    epsilonHalfSpec = epsilonHalfSpec, w, specData, $
    sigmaFull_ = sigma, $
    sigmaHalf_ = sigma_, $
	sigmaFullSpec_abp = sigma_abp, $
	sigmaHalfSpec_abp = sigma_abp_

	common switches
	common constants

;	Calculate dielectric tensor


    nSpec = n_elements(stixVars.StixR[0,*])
	epsilon		= dcomplexArr ( 3, 3, nSpec, runData.nR )
    if not keyword_set(noHalfGrid) then $
	epsilon_	= dcomplexArr ( 3, 3, nSpec, runData.nR-1 )

	sigma		= dcomplexArr ( 3, 3, nSpec, runData.nR )
    if not keyword_set(noHalfGrid) then $
	sigma_	= dcomplexArr ( 3, 3, nSpec, runData.nR-1 )

			sigma_abp = dComplexArr(3,3,nSpec,RunData.nR)
			sigma_abp_ = dComplexArr(3,3,nSpec,RunData.nR-1)


    if dielectric_freeSpace then begin

        epsilon[0,0,*,*]  = complex ( 1, 0 )
        epsilon[1,1,*,*]  = complex ( 1, 0 )
        epsilon[2,2,*,*]  = complex ( 1, 0 )

        if not keyword_set(noHalfGrid) then begin
        epsilon_[0,0,*,*]  = complex ( 1, 0 )
        epsilon_[1,1,*,*]  = complex ( 1, 0 )
        epsilon_[2,2,*,*]  = complex ( 1, 0 )
        endif

    endif else if dielectric_noPoloidal then begin

	    epsilon[0,0,*,*]	= stixVars.stixS
	    epsilon[0,2,*,*]	= -II * stixVars.stixD
	    epsilon[1,1,*,*]	= stixVars.stixP
	    epsilon[2,0,*,*]	= II * stixVars.stixD
	    epsilon[2,2,*,*]	= stixVars.stixS

        if not keyword_set(noHalfGrid) then begin
	    epsilon_[0,0,*,*]	= stixVars.stixS_
	    epsilon_[0,2,*,*]	= -II * stixVars.stixD_
	    epsilon_[1,1,*,*]	= stixVars.stixP_
	    epsilon_[2,0,*,*]	= II * stixVars.stixD_
	    epsilon_[2,2,*,*]	= stixVars.stixS_
        endif

		;epsilonNoPol	= epsilon
		;epsilonNoPol_	= epsilon_

    endif else begin

		;	Generic dielectric for arbitrary magnetic field direction
	
		bUnit_cyl	= [ [ runData.bField[*,0] / runData.bMag ], $	
						[ runData.bField[*,1] / runData.bMag ], $
						[ runData.bField[*,2] / runData.bMag ] ]

        if not keyword_set(noHalfGrid) then begin
		bUnit_cyl_	= [ [ runData.bField_[*,0] / runData.bMag_ ], $	
						[ runData.bField_[*,1] / runData.bMag_ ], $
						[ runData.bField_[*,2] / runData.bMag_ ] ]
        endif
	
		;	rotate to cartesian
	
		phi	= 0d0
        if not keyword_set(noHalfGrid) then $

		for i = 0L, runData.nR - 1L do begin
        for s = 0, nSpec-1 do begin

			;	Here the rotation from cylindrical to cartesian and
			;	its inverse is not used since for phi = 0, the 
			;	rotation is simply the identity matrix.
	
			;cyl2car	= [ [ cos ( phi ), -sin ( phi ), 0 ], $
			;			[ sin ( phi ), cos ( phi ), 0 ], $
			;			[ 0, 0, 1 ] ]
	
			;car2cyl	= invert ( cyl2car )
	
	
	
			epsilon_stix	= [ [ stixVars.stixS[i,s], II * stixVars.stixD[i,s], 0d0 ], $
								[ -II * stixVars.stixD[i,s], stixVars.stixS[i,s], 0d0 ], $
								[ 0d0, 0d0, stixVars.stixP[i,s] ] ]

			sigma_stix	    = [ [ stixVars.sig1[i,s], stixVars.sig2[i,s], 0d0 ], $
								[ -stixVars.sig2[i,s], stixVars.sig1[i,s], 0d0 ], $
								[ 0d0, 0d0, stixVars.sig3[i,s] ] ]


	        thisSigma = (epsilon_stix-identity(3))*w[i,s]*e0/II
            thisEpsilon = identity(3) + II / (w[i,s]*e0) * sigma_stix

            epsilon_stix = thisEpsilon

            ;print, 'r: ', runData.r[i]
            ;print, ''
            ;print, 'sigma_abp: '
            ;print, thisSigma

			epsilon[*,*,s,i]	= rotateEpsilon ( epsilon_stix, bUnit_cyl[i,*], w[i] )
    		sigma[*,*,s,i]	= rotateEpsilon ( sigma_stix, bUnit_cyl[i,*], w[i] )
			sigma_abp[*,*,s,i] = sigma_stix
        
            ;print, 'r: ', runData.r[i]
            ;print, 'bu: ', bUnit_cyl[i,*]
            ;print, 'rot: ', rot_
            ;stop    
			;	same for 1/2 grid

            if not keyword_set(noHalfGrid) then begin    
			if i lt runData.nR - 1 then begin
	
		
				epsilon_stix_	= [ [ stixVars.stixS_[i,s], II * stixVars.stixD_[i,s], 0d0 ], $
									[ -II * stixVars.stixD_[i,s], stixVars.stixS_[i,s], 0d0 ], $
									[ 0d0, 0d0, stixVars.stixP_[i,s] ] ]

				sigma_stix_	    = [ [ stixVars.sig1_[i,s], stixVars.sig2_[i,s], 0d0 ], $
								[ -stixVars.sig2_[i,s], stixVars.sig1_[i,s], 0d0 ], $
								[ 0d0, 0d0, stixVars.sig3_[i,s] ] ]

	            thisSigma_ = (epsilon_stix_-identity(3))*w[i,s]*e0/II
                thisEpsilon_ = identity(3) + II / (w[i,s]*e0) * sigma_stix_

                epsilon_stix_ = thisEpsilon_


				epsilon_[*,*,s,i]	= rotateEpsilon ( epsilon_stix_, bUnit_cyl_[i,*], w[i] )
				sigma_[*,*,s,i]	= rotateEpsilon ( sigma_stix_, bUnit_cyl_[i,*], w[i] )
				sigma_abp_[*,*,s,i] = sigma_stix_
	
			endif
            endif
		
		endfor
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

    epsilonFullSpec = epsilon
    epsilonHalfSpec = epsilon_

	SigmaFullSpec_abp = sigma_abp
	SigmaHalfSpec_abp = sigma_abp_

    epsilon = total(epsilon,3)
    epsilon_ = total(epsilon_,3)

end
