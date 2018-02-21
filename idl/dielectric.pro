pro dielectric, runData, stixVars, w, specData, $
	epsilonFull_ = epsilon, $
	epsilonHalf_ = epsilon_, $
    noHalfGrid = noHalfGrid, $
    epsilonFullSpec = epsilonFullSpec, $
    epsilonHalfSpec = epsilonHalfSpec, $
    sigmaFull_ = sigma, $
    sigmaHalf_ = sigma_, $
	sigmaFullSpec_abp = sigma_abp, $
	sigmaHalfSpec_abp = sigma_abp_

	common switches
	@dlg_constants

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
	    epsilon[0,2,*,*]	= -_II * stixVars.stixD
	    epsilon[1,1,*,*]	= stixVars.stixP
	    epsilon[2,0,*,*]	= _II * stixVars.stixD
	    epsilon[2,2,*,*]	= stixVars.stixS

        if not keyword_set(noHalfGrid) then begin
	    epsilon_[0,0,*,*]	= stixVars.stixS_
	    epsilon_[0,2,*,*]	= -_II * stixVars.stixD_
	    epsilon_[1,1,*,*]	= stixVars.stixP_
	    epsilon_[2,0,*,*]	= _II * stixVars.stixD_
	    epsilon_[2,2,*,*]	= stixVars.stixS_
        endif

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
	
		phi	= 0d0
        if not keyword_set(noHalfGrid) then $

		for i = 0L, runData.nR - 1L do begin
        for s = 0, nSpec-1 do begin

            f = runData.freq
            amu = specData[s].amu 
            atomicZ = specData[s].Z
            B = runData.bMag[i]
            density = (specData[s].n)[i]
            nu_omg = (specData[s].nuOmg)[i]

            epsilon_stix = kj_epsilon_cold(f, amu, atomicZ, B, density, nu_omg, sigma = ThisSigma)
        
            sigma_stix = ThisSigma

			epsilon[*,*,s,i]    = rotateEpsilon ( epsilon_stix, bUnit_cyl[i,*] )
    		sigma[*,*,s,i]	    = rotateEpsilon ( sigma_stix, bUnit_cyl[i,*] )
			sigma_abp[*,*,s,i]  = sigma_stix
          
            if i eq 256 then begin
                    print, s
                    print, 'abp : '
                    print, sigma_stix
                    print, 'rtz : '
                    print, sigma[*,*,s,i]
            endif

    
			;	same for 1/2 grid

            if not keyword_set(noHalfGrid) then begin    
			if i lt runData.nR - 1 then begin
	
                f = runData.freq
                amu = specData[s].amu 
                atomicZ = specData[s].Z
                B = runData.bMag_[i]
                density = (specData[s].n_)[i]
                nu_omg = (specData[s].nuOmg_)[i]

                epsilon_stix_ = kj_epsilon_cold(f, amu, atomicZ, B, density, nu_omg, sigma = ThisSigma_)
 
                sigma_stix_ = ThisSigma_
               
				epsilon_[*,*,s,i]	= rotateEpsilon ( epsilon_stix_, bUnit_cyl_[i,*] )
				sigma_[*,*,s,i]	= rotateEpsilon ( sigma_stix_, bUnit_cyl_[i,*] )
				sigma_abp_[*,*,s,i] = sigma_stix_
	
			endif
            endif
		
		endfor
        endfor

	endelse

;	Check epsilon for NaNs

	iiNaN	= where ( epsilon ne epsilon, iiNaNCnt )
	if iiNaNCnt gt 0 then begin

		print, 'ERROR: NaN detected in epsilon'
		iiNeg	= where ( runData.specData.n le 0, iiNegCnt )
		if iiNegCnt gt 0 then begin
			
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
