pro matfill, nR, nPhi, kz, r, r_, epsilon, epsilon_, w, dr, $
	aMat = aMat, nAll = nAll, nuca = nuca, nlca = nlca, $
    replaceFull = replace, replaceHalf = replace_, rhs=rhs, $
    kjIn = kjIn

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

	    aMat	= complexArr ( nAll, nAll )
        aMat[*] = 0

    endelse


	for i = 0UL, nR - 1L do begin

            if ar2EField and i gt 0 and i lt nR-2 then begin
               this_kjE = [$
                   kjIn.eR[i-1],kjIn.eT_[i-1],kjIn.eZ_[i-1], $
                   kjIn.eR[i],  kjIn.eT_[i],  kjIn.eZ_[i], $
                   kjIn.eR[i+1],kjIn.eT_[i+1],kjIn.eZ_[i+1]]
            endif
		    ;	r component
        

                if i gt 0 then begin
                if kjInput then begin
				    	; eT i-1/2
			        	aMat[3*i-2,3*i]	= -II * nPhi * r_[i-1] / ( r[i]^2 * dr ) 
				    	; eZ i-1/2
			        	aMat[3*i-1,3*i]	= -II * kz / dr 
                endif else begin
				    	; eT i-1/2
			        	aMat[3*i-2,3*i]	= -II * nPhi * r_[i-1] / ( r[i]^2 * dr ) $
                                            - w[i]^2*r_[i-1]/(2*c^2*r[i])*epsilon_[1,0,i-1]
				    	; eZ i-1/2
			        	aMat[3*i-1,3*i]	= -II * kz / dr $
                                            - w[i]^2*r_[i-1]/(2*c^2*r[i])*epsilon_[2,0,i-1]
                endelse
                endif

                if kjInput then begin
                    ; eR i
			    	aMat[3*i,3*i]	= nPhi^2 / r[i] + kz^2 - w[i]^2 / c^2
                endif else begin
					; eR i
				    t3	= nPhi^2 / r[i] + kz^2 - w[i]^2 / c^2 * epsilon[0,0,i]
			    	aMat[3*i,3*i]	= t3;nPhi^2 / r[i]^2 + kz^2 $
                                        ;- w[i]^2 / c^2 * epsilon[0,0,i]
                endelse

                if i lt nR-1 then begin
                if kjInput then begin
				    	; eT i+1/2
			        	aMat[3*i+1,3*i]	= II * nPhi * r_[i] / ( r[i]^2 * dr ) 
				    	; eZ i+1/2
			        	aMat[3*i+2,3*i]	= II * kz / dr 
                endif else begin
				    	; eT i+1/2
			        	aMat[3*i+1,3*i]	= II * nPhi * r_[i] / ( r[i]^2 * dr ) $
                                            - w[i]^2*r_[i]/(2*c^2*r[i])*epsilon_[1,0,i]
				    	; eZ i+1/2
			        	aMat[3*i+2,3*i]	= II * kz / dr $
                                            - w[i]^2*r_[i]/(2*c^2*r[i])*epsilon_[2,0,i]
                endelse
                endif

                ;; TF/SF Method
                ;; First add the incident field to the RHS
                ;; but only for the "total" field points, not "scattered" field, 
                ;; i.e., where replace[i]==1
                ;if ar2EField and i gt 0 and i lt nR-1 then begin
                ;if replace[i] then begin
                ;    this_aMat = aMat[3*i-3:3*i+5,3*i]
                ;    rhs[3*i] -= total(this_aMat * this_kjE)
                ;endif
                ;endif
                ; At the left TF/SF boundary we have to add the incident
                ; field to ensure the derivative terms are between TOTAL fields, 
                ; i.e., they are just TOTAL field differences 
                if ar2EField and i gt 0 then begin
                if replace[i-1] eq 0 and replace[i] then begin
                    this_aMat = aMat[3*i-3:3*i+5,3*i]
                    this_aMat[3:-1] = 0
                    rhs[3*i] -= total(this_aMat * this_kjE)
                endif 
                endif
                ; At the right TF/SF boundary we have to add the incident
                ; field to ensure the derivative terms are between TOTAL fields, 
                ; i.e., they are just TOTAL field differences 
                if ar2EField and i lt nR-1 then begin
                if replace[i] and replace[i+1] eq 0 then begin
                    this_aMat = aMat[3*i-3:3*i+5,3*i]
                    this_aMat[0:3] = 0
                    rhs[3*i] -= total(this_aMat * this_kjE)
                endif 
                endif
 
		;	phi component
			if i lt nR-1 then begin
          
                    if i gt 0 then begin
				            aMat[3*i-2,3*i+1] = -r_[i-1] / ( r[i] * dr^2 ) 
                    endif

                    if kjInput then begin
                        aMat[3*i,3*i+1]	= -II*nPhi/(r[i]*dr) 
				        aMat[3*i+3,3*i+1]	= II * nPhi / (r[i+1]*dr) 
                    endif else begin
				        aMat[3*i,3*i+1]	= -II*nPhi/(r[i]*dr) $
                                            - w[i]^2*r[i]/(2*c^2*r_[i])*epsilon[0,1,i]
				        aMat[3*i+3,3*i+1]	= II * nPhi / (r[i+1]*dr) $
                                            - w[i]^2*r[i+1]/(2*c^2*r_[i])*epsilon[0,1,i+1]
                    endelse
                    if kjInput then begin
				        aMat[3*i+1,3*i+1]	= kz^2 + r_[i]/(r[i+1]*dr^2) $
                                            + r_[i]/(r[i]*dr^2) $
                                            - w[i]^2/c^2

				        aMat[3*i+2,3*i+1]	= -kz*nPhi/r_[i] 
                    endif else begin
				        aMat[3*i+1,3*i+1]	= kz^2 + r_[i]/(r[i+1]*dr^2) $
                                                + r_[i]/(r[i]*dr^2) $
                                                - w[i]^2/c^2*epsilon_[1,1,i]

				        aMat[3*i+2,3*i+1]	= -kz*nPhi/r_[i] - w[i]^2/(c^2)*epsilon_[2,1,i]
                    endelse
				    if i lt nR-2 then begin
				            aMat[3*i+4,3*i+1]	= -r_[i+1]/(r[i+1]*dr^2)
                    endif

                    ;; TF/SF Method
                    ;if ar2EField and i gt 0 and i lt nR-1 then begin
                    ;if replace[i-1] and replace[i] and replace[i+1] and replace[i+2] then begin
                    ;    this_aMat = aMat[3*i-3:3*i+5,3*i+1]
                    ;    rhs[3*i+1] -= total(this_aMat * this_kjE)
                    ;endif
                    ;endif
                    ; Both the half points inside and outside of the full
                    ; point interface layer need to be corrected for since 
                    ; unlinke the r component, the th and z components have
                    ; a dependence on a larger stencil. 
                    ; 
                    ; Left TF/SF boundary
                    ; at the half point outside we need to SUBTRACT the incident
                    ; field so this is an equation for SCATTERED fields
                    if ar2EField and i gt 0 then begin ; outside
                    if replace[i] eq 0 and replace[i+1] then begin
                        this_aMat = -aMat[3*i-3:3*i+5,3*i+1]
                        this_aMat[0:5] = 0
                        rhs[3*i+1] -= total(this_aMat * this_kjE)
                    endif 
                    endif
                    ; at the half point inside we need to ADD the incident
                    ; field so this is an equation for TOTAL fields
                    if ar2EField and i gt 0 then begin ; inside
                    if replace[i-1] eq 0 and replace[i] and replace[i+1] then begin
                        this_aMat = aMat[3*i-3:3*i+5,3*i+1]
                        this_aMat[3:-1] = 0
                        rhs[3*i+1] -= total(this_aMat * this_kjE)
                    endif 
                    endif
                    ; Right TF/SF boundary
                    ; at the half point outside we need to SUBTRACT the incident
                    ; field so this is an equation for SCATTERED fields
                    if ar2EField and i lt nR-1 then begin ; outside
                    if replace[i] and replace[i+1] eq 0 then begin
                        this_aMat = -aMat[3*i-3:3*i+5,3*i+1]
                        this_aMat[4:-1] = 0
                        rhs[3*i+1] -= total(this_aMat * this_kjE)
                    endif 
                    endif
                    ; at the half point inside we need to ADD the incident
                    ; field so this is an equation for TOTAL fields
                    if ar2EField and i lt nR-2 then begin ; inside
                    if replace[i] and replace[i+1] and replace[i+2] eq 0 then begin
                        this_aMat = aMat[3*i-3:3*i+5,3*i+1]
                        this_aMat[0:3] = 0
                        rhs[3*i+1] -= total(this_aMat * this_kjE)
                    endif 
                    endif
 
		;	z component	

					t1	= -r[i]/(r_[i]*dr^2)
                    if kjInput then begin
					    t2	= -II*kz*r[i]/(r_[i]*dr) ; Propagate this "dr" correction
					    t5	= II * kz * r[i+1] / (r_[i]*dr) ; And this "dr" correction 
                    endif else begin
					    t2	= -II*kz*r[i]/(r_[i]*dr)$
                                        -w[i]^2*r[i]/(2*c^2*r_[i])*epsilon[0,2,i]
                        t5	= II * kz * r[i+1] / (r_[i]*dr) $
                                            - w[i]^2*r[i+1]/(2*c^2*r_[i])*epsilon[0,2,i+1]
                    endelse
                    if kjInput then begin
					    t3	= -nPhi*kz/r_[i] ; Also propagate this "r_[i]" correction
					    t4	= (r[i+1]+r[i])/(r_[i]*dr^2) + nPhi^2/r_[i]^2 - w[i]^2/c^2 ; Also propagate this "r_[i]^2" correction
                    endif else begin
					    t3	= -nPhi*kz/r_[i] - w[i]^2/c^2*epsilon_[1,2,i]
					    t4	= (r[i+1]+r[i])/(r_[i]*dr^2) + nPhi^2/r_[i]^2 $
                                                - w[i]^2/c^2*epsilon_[2,2,i]
                    endelse
					t6	= -r[i+1] / ( r_[i] * dr^2 ) 

	                if i gt 0 then begin 
				            aMat[3*i-1,3*i+2] = t1
                    endif
				        aMat[3*i,3*i+2]   = t2 
				        aMat[3*i+1,3*i+2] = t3
				        aMat[3*i+2,3*i+2] = t4
				        aMat[3*i+3,3*i+2] = t5
			        if i lt nR-2 then begin
				            aMat[3*i+5,3*i+2] = t6 
                    endif

                    ;; TF/SF Method
                    ;if ar2EField and i gt 0 and i lt nR-1 then begin
                    ;if replace[i-1] and replace[i] and replace[i+1] and replace[i+2] then begin
                    ;    this_aMat = aMat[3*i-3:3*i+5,3*i+2]
                    ;    rhs[3*i+2] -= total(this_aMat * this_kjE)
                    ;endif
                    ;endif
                    ; Left TF/SF boundary
                    ; at the half point outside we need to SUBTRACT the incident
                    ; field so this is an equation for SCATTERED fields
                    if ar2EField and i gt 0 then begin
                    if replace[i] eq 0 and replace[i+1] then begin
                        this_aMat = -aMat[3*i-3:3*i+5,3*i+2]
                        this_aMat[0:5] = 0
                        rhs[3*i+2] -= total(this_aMat * this_kjE)
                    endif 
                    endif
                    ; at the half point inside we need to ADD the incident
                    ; field so this is an equation for TOTAL fields
                    if ar2EField and i gt 0 then begin
                    if replace[i-1] eq 0 and replace[i] and replace[i+1] then begin
                        this_aMat = aMat[3*i-3:3*i+5,3*i+2]
                        this_aMat[3:-1] = 0
                        rhs[3*i+2] -= total(this_aMat * this_kjE)
                    endif 
                    endif
                    ; Right TF/SF boundary
                    ; at the half point outside we need to SUBTRACT the incident
                    ; field so this is an equation for SCATTERED fields
                    if ar2EField and i lt nR-1 then begin
                    if replace[i] and replace[i+1] eq 0 then begin
                        this_aMat = -aMat[3*i-3:3*i+5,3*i+2]
                        this_aMat[4:-1] = 0
                        rhs[3*i+2] -= total(this_aMat * this_kjE)
                    endif 
                    endif
                    ; at the half point inside we need to ADD the incident
                    ; field so this is an equation for TOTAL fields
                    if ar2EField and i lt nR-2 then begin
                    if replace[i] and replace[i+1] and replace[i+2] eq 0 then begin
                        this_aMat = aMat[3*i-3:3*i+5,3*i+2]
                        this_aMat[0:3] = 0
                        rhs[3*i+2] -= total(this_aMat * this_kjE)
                    endif 
                    endif
 
			endif
	endfor
end
