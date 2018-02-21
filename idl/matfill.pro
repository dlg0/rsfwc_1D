pro matfill, nR, nPhi, kz, r, r_, epsilon, epsilon_, w, dr, $
	aMat = aMat, nAll = nAll, nuca = nuca, nlca = nlca, $
    replaceFull = replace, replaceHalf = replace_, rhs=rhs

	@dlg_constants
	common switches

	;	Build matrix

    print, '*** filling matrix'

    if bandStorage then begin

        nuca    = 5L
        nlca    = 3L
        m   = nlca + nuca + 1L
        aMat_bandStorage    = complexArr ( nAll, m )

        print, 'aMat: ', (nAll*2.0)^2*16.0/(1024.0^2)
        print, 'aMat_bandStorage: ', (m*1d0*nAll*2.0)*16.0/(1024.0^2) 

    endif else begin

        ; The extra are to apply the PEC BCs
	    aMat	= dcomplexArr ( nAll+4, nAll+4 )
        aMat[*] = 0d0

    endelse

    ;if kjInput then begin
    ;
    ;        epsilonCopy = epsilon
    ;        epsilonCopy_ = epsilon_

    ;        epsilon = IntArr(3,3,nR)
    ;        epsilon_ = IntArr(3,3,nR-1)

    ;        for i=0,nR-1 do begin
    ;            epsilon[*,*,i] = identity(3)
    ;        endfor  
    ;        for i=0,nR-2 do begin
    ;            epsilon_[*,*,i] = identity(3)
    ;        endfor
    ;endif

    conservative_form = 0

    ; PEC BC eqns (Tangential components == 0)

    ; Left boundary
    ; eT_[i-1] = -eT_[i] sets eT[i] == 0, i.e., eT_[i-1] + eT_[i] == 0 (so set RHS == 0)
    aMat[0,0] = 1
    aMat[3,0] = 1
    ; eZ_[i-1] = -eZ_[i] sets eZ[i] == 0
    aMat[1,1] = 1
    aMat[4,1] = 1

    ; Right boundary
    ; eT_[i-1] = -eT_[i] sets eT[i] == 0, i.e., eT_[i-1] + eT_[i] == 0 (so set RHS == 0)
    aMat[-5,-2] = 1
    aMat[-2,-2] = 1
    ; eZ_[i-1] = -eZ_[i] sets eZ[i] == 0
    aMat[-4,-1] = 1
    aMat[-1,-1] = 1


    AOffSet = 2

    rhs = [0,0,rhs,0,0]

    ; Create an overlap layer
    rBak = r
    rBak_ = r
    epsilonBak = epsilon
    epsilonBak_ = epsilon

    r = [r[0]-dr,r,r[-1]+dr]
    r_ = [r_[0]-dr,r_,r_[-1]+dr]

    epsilon = ComplexArr(3,3,n_elements(epsilonBak[0,0,*])+2)
    epsilon[*,*,1:-2] = epsilonBak
    epsilon[*,*,0] = epsilonBak[0]
    epsilon[*,*,-1] = epsilonBak[-1]

    epsilon_ = ComplexArr(3,3,n_elements(epsilonBak_[0,0,*])+2)
    epsilon_[*,*,1:-2] = epsilonBak_
    epsilon_[*,*,0] = epsilonBak_[0]
    epsilon_[*,*,-1] = epsilonBak_[-1]

    SpatialOffSet = 1
	for i = 0 +SpatialOffSet, nR-1 +SpatialOffSet do begin

            jj = 3*(i-SpatialOffSet)+AOffSet

		    ;	r component
       
                if conservative_form then begin 
   				    ; eT i-1/2
			        t1 = -II * nPhi * r_[i-1] / ( r[i]^2 * dr ) $
                                        - w^2*r_[i-1]/(2*c^2*r[i])*epsilon_[1,0,i-1]
				    ; eZ i-1/2
			        t2 = -II * kz / dr - w^2*r_[i-1]/(2*c^2*r[i])*epsilon_[2,0,i-1]

                    aMat[jj-2,jj] = t1
                    aMat[jj-1,jj] = t2

				    ; eR i
				    t3	= nPhi^2 / r[i]^2 + kz^2 - w^2 / c^2 * epsilon[0,0,i]
			        aMat[jj,jj]	= t3

			        ; eT i+1/2
			        t4 = II * nPhi * r_[i] / ( r[i]^2 * dr ) $
                                        - w^2*r_[i]/(2*c^2*r[i])*epsilon_[1,0,i]
				    ; eZ i+1/2
			        t5 = II * kz / dr - w^2*r_[i]/(2*c^2*r[i])*epsilon_[2,0,i]

                    aMat[jj+1,jj] = t4             
                    aMat[jj+2,jj] = t5

                endif else begin

                    t1 = 0
                    t2 = 0
                    t3 = 0
                    t4 = 0
                    t5 = 0
                    t6 = 0
                    t7 = 0
                    t8 = 0
                    t9 = 0

                    t2 = _II*nPhi*(1-(2*r[i]/dr))/(2*r[i]^2) - w^2/(2*_c^2)*r_[i-1]/r[i]*epsilon_[1,0,i-1]
                    t3 = -_II*kz/dr - w^2/(2*_c^2)*r_[i-1]/r[i] * epsilon_[2,0,i-1]
                    t4 = kz^2 + nPhi^2/r[i]^2 - w^2/_c^2 * epsilon[0,0,i]
                    t5 = _II*nPhi*(2*r[i]+dr)/(2*r[i]^2*dr) - w^2/(2*_c^2)*r_[i]/r[i] * epsilon_[1,0,i]
                    t6 = _II*kz/dr - w^2/(2*_c^2)*r_[i]/r[i] * epsilon_[2,0,i]

                    aMat[jj-3,jj+0] = t1
				    aMat[jj-2,jj+0] = t2
				    aMat[jj-1,jj+0] = t3 

				    aMat[jj+0,jj+0] = t4 
				    aMat[jj+1,jj+0] = t5 
				    aMat[jj+2,jj+0] = t6 

				    ;aMat[jj+3,jj+0] = t7 
				    ;aMat[jj+4,jj+0] = t8 
				    ;aMat[jj+5,jj+0] = t9 

                endelse 


		;	phi component
			if i lt nR-1+SpatialOffSet then begin
         
                    if conservative_form then begin 
                        ;if i gt 0 then begin
				            t1 = -r_[i-1] / ( r[i] * dr^2 ) 
                            aMat[jj-2,jj+1] = t1
                        ;endif

				        t2 = -II*nPhi/(r[i]*dr) $
                                            - w^2*r[i]/(2*c^2*r_[i])*epsilon[0,1,i]
				        t3 = II * nPhi / (r[i+1]*dr) $
                                            - w^2*r[i+1]/(2*c^2*r_[i])*epsilon[0,1,i+1]
				        t4 = kz^2 + r_[i]/(r[i+1]*dr^2) $
                                                + r_[i]/(r[i]*dr^2) $
                                                - w^2/c^2*epsilon_[1,1,i]
				        t5 = -kz*nPhi/r_[i] - w^2/(c^2)*epsilon_[2,1,i]

                        aMat[jj+0,jj+1] = t2
                        aMat[jj+1,jj+1] = t4
                        aMat[jj+2,jj+1] = t5
                        aMat[jj+3,jj+1] = t3

				        ;if i lt nR-2 then begin
				            t6 = -r_[i+1]/(r[i+1]*dr^2)
                            aMat[jj+4,jj+1] = t6
                        ;endif

                    endif else begin

                        t1 = 0
                        t2 = 0
                        t3 = 0
                        t4 = 0
                        t5 = 0
                        t6 = 0
                        t7 = 0
                        t8 = 0
                        t9 = 0

                        t2 = (-2*r_[i]+dr)/(2*r_[i]*dr^2)
                        t4 = -_II*nPhi*(2*r_[i]+dr)/(2*r_[i]^2*dr) - w^2/(2*_c^2)*r[i]/r_[i]*epsilon[0,1,i]
                        t5 = kz^2 + 1/r_[i]^2 + 2/dr^2 - w^2/_c^2 * epsilon_[1,1,i]
                        t6 = -kz*nPhi/r_[i] - w^2/_c^2 * epsilon_[2,1,i]
                        t7 = -_II*nPhi/(2*r_[i]^2) + _II*nPhi/(r_[i]*dr) - w^2/(2*_c^2)*r[i+1]/r_[i]*epsilon[0,1,i+1]
                        t8 = -(2+dr/r_[i])/(2*dr^2)
                        
                        aMat[jj-3,jj+1] = t1
				        aMat[jj-2,jj+1] = t2
				        aMat[jj-1,jj+1] = t3 

				        aMat[jj+0,jj+1] = t4 
				        aMat[jj+1,jj+1] = t5 
				        aMat[jj+2,jj+1] = t6 

				        aMat[jj+3,jj+1] = t7 
				        aMat[jj+4,jj+1] = t8 
				        aMat[jj+5,jj+1] = t9 

                    endelse 
 
		;	z component	

                    if conservative_form then begin

					    t1	= -r[i]/(r_[i]*dr^2)

					    t2	= -II*kz*r[i]/(r_[i]*dr)$
                                        -w^2*r[i]/(2*c^2*r_[i])*epsilon[0,2,i]
                        t5	= II * kz * r[i+1] / (r_[i]*dr) $
                                            - w^2*r[i+1]/(2*c^2*r_[i])*epsilon[0,2,i+1]
                        
					    t3	= -nPhi*kz/r_[i] - w^2/c^2*epsilon_[1,2,i]
					    t4	= (r[i+1]+r[i])/(r_[i]*dr^2) + nPhi^2/r_[i]^2 $
                                                    - w^2/c^2*epsilon_[2,2,i]
					    t6	= -r[i+1] / ( r_[i] * dr^2 ) 



	                    ;if i gt 0 then begin 
				            aMat[jj-1,jj+2] = t1
                        ;endif
				            aMat[jj+0,jj+2] = t2
				            aMat[jj+1,jj+2] = t3
				            aMat[jj+2,jj+2] = t4
				            aMat[jj+3,jj+2] = t5
			            ;if i lt nR-2 then begin
				            aMat[jj+5,jj+2] = t6
                        ;endif

                    endif else begin

                        t1 = 0
                        t2 = 0
                        t3 = 0
                        t4 = 0
                        t5 = 0
                        t6 = 0
                        t7 = 0
                        t8 = 0
                        t9 = 0

                        t3 = (-2*r_[i]+dr)/(2*r_[i]*dr^2)
                        t9 = -(2+dr/r_[i])/(2*dr^2)
                        t7 = 0.5*_II*kz*(1/r_[i]+2/dr) - w^2*r[i+1]*epsilon[0,2,i+1] / (2*_c^2*r_[i])
                        t6 = nPhi^2/r_[i]^2 + 2/dr^2 - w^2/_c^2*epsilon_[2,2,i]
                        t5 = -kz*nPhi/r_[i] - w^2/_c^2 * epsilon_[1,2,i]
                        t4 = 0.5*_II*kz*(1/r_[i]-2/dr) - w^2*r[i]*epsilon[0,2,i]/(2*_c^2*r_[i])

                        aMat[jj-3,jj+2] = t1
				        aMat[jj-2,jj+2] = t2
				        aMat[jj-1,jj+2] = t3 

				        aMat[jj+0,jj+2] = t4 
				        aMat[jj+1,jj+2] = t5 
				        aMat[jj+2,jj+2] = t6 

				        aMat[jj+3,jj+2] = t7 
				        aMat[jj+4,jj+2] = t8 
				        aMat[jj+5,jj+2] = t9 

                    endelse 

			endif
	endfor

    r = rBak
    r_ = rBak_
    epsilon = epsilonBak 
    epsilon_ = epsilonBak_ 

end
