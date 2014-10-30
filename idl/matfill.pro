pro matfill, nR, nPhi, kz, r, r_, epsilon, epsilon_, w, dr, $
	aMat = aMat, nAll = nAll, nuca = nuca, nlca = nlca, $
    replaceFull = replace, replaceHalf = replace_, rhs=rhs

	common constants
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
	    aMat	= complexArr ( nAll+4, nAll+4 )
        aMat[*] = 0

    endelse

    if kjInput then begin
    
            epsilonCopy = epsilon
            epsilonCopy_ = epsilon_

            epsilon = IntArr(3,3,nR)
            epsilon_ = IntArr(3,3,nR-1)

            for i=0,nR-1 do begin
                epsilon[*,*,i] = identity(3)
            endfor  
            for i=0,nR-2 do begin
                epsilon_[*,*,i] = identity(3)
            endfor
    endif

    ;t1Tmp = ComplexArr(nR)
    ;t2Tmp = ComplexArr(nR)
    ;t3Tmp = ComplexArr(nR)
    ;t4Tmp = ComplexArr(nR)
    ;t5Tmp = ComplexArr(nR)
    ;t6Tmp = ComplexArr(nR)

    conservative_form = 1

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

            ;if ar2EField and i gt 0 and i lt nR-2 then begin
            ;   this_kjE = [$
            ;       kjIn.eR[i-1],kjIn.eT_[i-1],kjIn.eZ_[i-1], $
            ;       kjIn.eR[i],  kjIn.eT_[i],  kjIn.eZ_[i], $
            ;       kjIn.eR[i+1],kjIn.eT_[i+1],kjIn.eZ_[i+1]]
            ;endif
		    ;	r component
        
                ;; Set the Hperp == 0 at the boundary (PEC)
                ;if i eq 0 then begin

                ;    t1 = 0
                ;    t2 = 0
                ;    t3 = 0
                ;    t4 = 0
                ;    t5 = 0
                ;    t6 = 0
                ;    t7 = 0
                ;    t8 = 0
                ;    t9 = 0

                ;    

				;    if i gt 0 then begin
                ;    aMat[3*i-3,3*i+2] = t1
				;    aMat[3*i-2,3*i+2] = t2
				;    aMat[3*i-1,3*i+2] = t3 
                ;    endif

				;    aMat[3*i+0,3*i+2] = t4 
				;    aMat[3*i+1,3*i+2] = t5 
				;    aMat[3*i+2,3*i+2] = t6 
				;    aMat[3*i+3,3*i+2] = t7 

                ;    if i lt nR-2 then begin
				;    aMat[3*i+4,3*i+2] = t8 
				;    aMat[3*i+5,3*i+2] = t9 
                ;    endif


                ;endif

                ;if i gt 0 then begin
   				    ; eT i-1/2
			        t1 = -II * nPhi * r_[i-1] / ( r[i]^2 * dr ) $
                                        - w^2*r_[i-1]/(2*c^2*r[i])*epsilon_[1,0,i-1]
				    ; eZ i-1/2
			        t2 = -II * kz / dr - w^2*r_[i-1]/(2*c^2*r[i])*epsilon_[2,0,i-1]

                    aMat[jj-2,jj] = t1
                    aMat[jj-1,jj] = t2
                ;endif

					; eR i
				    t3	= nPhi^2 / r[i]^2 + kz^2 - w^2 / c^2 * epsilon[0,0,i]
			    	aMat[jj,jj]	= t3

                ;if i lt nR-1 then begin
			    	; eT i+1/2
			        t4 = II * nPhi * r_[i] / ( r[i]^2 * dr ) $
                                        - w^2*r_[i]/(2*c^2*r[i])*epsilon_[1,0,i]
				    ; eZ i+1/2
			        t5 = II * kz / dr - w^2*r_[i]/(2*c^2*r[i])*epsilon_[2,0,i]

                    aMat[jj+1,jj] = t4             
                    aMat[jj+2,jj] = t5

                ;endif

                    ;if kjInput then begin 
                    ;extraTerm = -II*w*u0*kjIn.jpR[i]
                    ;aMat[*,3*i] += extraTerm
                    ;endif


                ;;; TF/SF Method
                ;;; First add the incident field to the RHS
                ;;; but only for the "total" field points, not "scattered" field, 
                ;;; i.e., where replace[i]==1
                ;;if ar2EField and i gt 0 and i lt nR-1 then begin
                ;;if replace[i] then begin
                ;;    this_aMat = aMat[3*i-3:3*i+5,3*i]
                ;;    rhs[3*i] -= total(this_aMat * this_kjE)
                ;;endif
                ;;endif
                ;; At the left TF/SF boundary we have to add the incident
                ;; field to ensure the derivative terms are between TOTAL fields, 
                ;; i.e., they are just TOTAL field differences 
                ;if ar2EField and i gt 0 then begin
                ;if replace[i-1] eq 0 and replace[i] then begin
                ;    this_aMat = aMat[3*i-3:3*i+5,3*i]
                ;    this_aMat[3:-1] = 0
                ;    rhs[3*i] -= total(this_aMat * this_kjE)
                ;endif 
                ;endif
                ;; At the right TF/SF boundary we have to add the incident
                ;; field to ensure the derivative terms are between TOTAL fields, 
                ;; i.e., they are just TOTAL field differences 
                ;if ar2EField and i lt nR-1 then begin
                ;if replace[i] and replace[i+1] eq 0 then begin
                ;    this_aMat = aMat[3*i-3:3*i+5,3*i]
                ;    this_aMat[0:3] = 0
                ;    rhs[3*i] -= total(this_aMat * this_kjE)
                ;endif 
                ;endif
 
		;	phi component
			if i lt nR-1+SpatialOffSet then begin
          
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

                    ;if kjInput then begin 
                    ;extraTerm = -II*w*u0*kjIn.jpT_[i]
                    ;aMat[*,3*i+1] += extraTerm
                    ;endif

                    ;;; TF/SF Method
                    ;;if ar2EField and i gt 0 and i lt nR-1 then begin
                    ;;if replace[i-1] and replace[i] and replace[i+1] and replace[i+2] then begin
                    ;;    this_aMat = aMat[3*i-3:3*i+5,3*i+1]
                    ;;    rhs[3*i+1] -= total(this_aMat * this_kjE)
                    ;;endif
                    ;;endif
                    ;; Both the half points inside and outside of the full
                    ;; point interface layer need to be corrected for since 
                    ;; unlinke the r component, the th and z components have
                    ;; a dependence on a larger stencil. 
                    ;; 
                    ;; Left TF/SF boundary
                    ;; at the half point outside we need to SUBTRACT the incident
                    ;; field so this is an equation for SCATTERED fields
                    ;if ar2EField and i gt 0 then begin ; outside
                    ;if replace[i] eq 0 and replace[i+1] then begin
                    ;    this_aMat = -aMat[3*i-3:3*i+5,3*i+1]
                    ;    this_aMat[0:5] = 0
                    ;    rhs[3*i+1] -= total(this_aMat * this_kjE)
                    ;endif 
                    ;endif
                    ;; at the half point inside we need to ADD the incident
                    ;; field so this is an equation for TOTAL fields
                    ;if ar2EField and i gt 0 then begin ; inside
                    ;if replace[i-1] eq 0 and replace[i] and replace[i+1] then begin
                    ;    this_aMat = aMat[3*i-3:3*i+5,3*i+1]
                    ;    this_aMat[3:-1] = 0
                    ;    rhs[3*i+1] -= total(this_aMat * this_kjE)
                    ;endif 
                    ;endif
                    ;; Right TF/SF boundary
                    ;; at the half point outside we need to SUBTRACT the incident
                    ;; field so this is an equation for SCATTERED fields
                    ;if ar2EField and i lt nR-1 then begin ; outside
                    ;if replace[i] and replace[i+1] eq 0 then begin
                    ;    this_aMat = -aMat[3*i-3:3*i+5,3*i+1]
                    ;    this_aMat[4:-1] = 0
                    ;    rhs[3*i+1] -= total(this_aMat * this_kjE)
                    ;endif 
                    ;endif
                    ;; at the half point inside we need to ADD the incident
                    ;; field so this is an equation for TOTAL fields
                    ;if ar2EField and i lt nR-2 then begin ; inside
                    ;if replace[i] and replace[i+1] and replace[i+2] eq 0 then begin
                    ;    this_aMat = aMat[3*i-3:3*i+5,3*i+1]
                    ;    this_aMat[0:3] = 0
                    ;    rhs[3*i+1] -= total(this_aMat * this_kjE)
                    ;endif 
                    ;endif
 
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
                        t7 = 0.5*II*kz*(1/r_[i]+2/dr) - w^2*r[i+1]*epsilon[0,2,i+1] / (2*c^2*r_[i])
                        t6 = nPhi^2/r_[i]^2 + 2/dr^2 - w^2/c^2*epsilon_[2,2,i]
                        t5 = -kz*nPhi/r_[i] - w^2/c^2 * epsilon_[1,2,i]
                        t4 = 0.5*II*kz*(1/r_[i]-2/dr) - w^2*r[i]*epsilon[0,2,i]/(2*c^2*r_[i])

				        ;if i gt 0 then begin
                        aMat[jj-3,jj+2] = t1
				        aMat[jj-2,jj+2] = t2
				        aMat[jj-1,jj+2] = t3 
                        ;endif

				        aMat[jj+0,jj+2] = t4 
				        aMat[jj+1,jj+2] = t5 
				        aMat[jj+2,jj+2] = t6 
				        aMat[jj+3,jj+2] = t7 

                        ;if i lt nR-2 then begin
				        aMat[jj+4,jj+2] = t8 
				        aMat[jj+5,jj+2] = t9 
                        ;endif

                    endelse 

                    ;if kjInput then begin 
                    ;extraTerm = -II*w*u0*kjIn.jpZ_[i]
                    ;aMat[*,3*i+2] += extraTerm
                    ;endif

                    ;;; TF/SF Method
                    ;;if ar2EField and i gt 0 and i lt nR-1 then begin
                    ;;if replace[i-1] and replace[i] and replace[i+1] and replace[i+2] then begin
                    ;;    this_aMat = aMat[3*i-3:3*i+5,3*i+2]
                    ;;    rhs[3*i+2] -= total(this_aMat * this_kjE)
                    ;;endif
                    ;;endif
                    ;; Left TF/SF boundary
                    ;; at the half point outside we need to SUBTRACT the incident
                    ;; field so this is an equation for SCATTERED fields
                    ;if ar2EField and i gt 0 then begin
                    ;if replace[i] eq 0 and replace[i+1] then begin
                    ;    this_aMat = -aMat[3*i-3:3*i+5,3*i+2]
                    ;    this_aMat[0:5] = 0
                    ;    rhs[3*i+2] -= total(this_aMat * this_kjE)
                    ;endif 
                    ;endif
                    ;; at the half point inside we need to ADD the incident
                    ;; field so this is an equation for TOTAL fields
                    ;if ar2EField and i gt 0 then begin
                    ;if replace[i-1] eq 0 and replace[i] and replace[i+1] then begin
                    ;    this_aMat = aMat[3*i-3:3*i+5,3*i+2]
                    ;    this_aMat[3:-1] = 0
                    ;    rhs[3*i+2] -= total(this_aMat * this_kjE)
                    ;endif 
                    ;endif
                    ;; Right TF/SF boundary
                    ;; at the half point outside we need to SUBTRACT the incident
                    ;; field so this is an equation for SCATTERED fields
                    ;if ar2EField and i lt nR-1 then begin
                    ;if replace[i] and replace[i+1] eq 0 then begin
                    ;    this_aMat = -aMat[3*i-3:3*i+5,3*i+2]
                    ;    this_aMat[4:-1] = 0
                    ;    rhs[3*i+2] -= total(this_aMat * this_kjE)
                    ;endif 
                    ;endif
                    ;; at the half point inside we need to ADD the incident
                    ;; field so this is an equation for TOTAL fields
                    ;if ar2EField and i lt nR-2 then begin
                    ;if replace[i] and replace[i+1] and replace[i+2] eq 0 then begin
                    ;    this_aMat = aMat[3*i-3:3*i+5,3*i+2]
                    ;    this_aMat[0:3] = 0
                    ;    rhs[3*i+2] -= total(this_aMat * this_kjE)
                    ;endif 
                    ;endif
 
			endif
	endfor

    r = rBak
    r_ = rBak_
    epsilon = epsilonBak 
    epsilon_ = epsilonBak_ 

stop
   ; check the solve wth another routine
   ; find the zero crossing point in amat
end
