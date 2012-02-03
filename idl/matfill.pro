pro matfill, nR, nPhi, kz, r, r_, epsilon, epsilon_, w, dr, $
	aMat = aMat, nAll = nAll, nuca = nuca, nlca = nlca

	common constants
	common switches

	;	Build matrix

    print, '*** filling matrix'

	nAll	= nR + 2L * ( nR - 1L )
    if bandStorage then begin

        nuca    = 5L
        nlca    = 3L
        m   = nlca + nuca + 1L
        aMat_bandStorage    = complexArr ( nAll, m )

        print, 'aMat: ', (nAll*2.0)^2*16.0/(1024.0^2)
        print, 'aMat_bandStorage: ', (m*1d0*nAll*2.0)*16.0/(1024.0^2) 

    endif else begin

	    aMat	= complexArr ( nAll, nAll )

    endelse


	for i = 0UL, nR - 1L do begin

		;	r component
            if bandStorage then begin

                j_  = 3*i
                if i gt 0 then begin
                    i_  = 3*i-2
					; eT i-1/2
			    	aMat_bandStorage[i_,nuca+j_-i_]	= -II * nPhi * r_[i-1] / ( r[i]^2 * dr ) $
                                        - w[i]^2*r_[i-1]/(2*c^2*r[i])*epsilon_[1,0,i-1]
                    i_  = 3*i-1

					; eZ i-1/2
			    	aMat_bandStorage[i_,nuca+j_-i_]	= -II * kz / dr $
                                        - w[i]^2*r_[i-1]/(2*c^2*r[i])*epsilon_[2,0,i-1]
                endif
                    
                    i_=3*i    
					; eR i
			    	aMat_bandStorage[i_,nuca+j_-i_]	= nPhi^2 / r[i]^2 + kz^2 $
                                        - w[i]^2 / c^2 * epsilon[0,0,i]
                if i lt nR-1 then begin
                    i_=3*i+1
					; eT i+1/2
			    	aMat_bandStorage[i_,nuca+j_-i_]	= II * nPhi * r_[i] / ( r[i]^2 * dr ) $
                                        - w[i]^2*r_[i]/(2*c^2*r[i])*epsilon_[1,0,i]
                    i_=3*i+2
					; eZ i+1/2
			    	aMat_bandStorage[i_,nuca+j_-i_]	= II * kz / dr $
                                        - w[i]^2*r_[i]/(2*c^2*r[i])*epsilon_[2,0,i]
                endif

            endif else if kjInput then begin

                if i gt 0 then begin
					; eT i-1/2
			    	aMat[3*i-2,3*i]	= -II * nPhi * r_[i-1] / ( r[i]^2 * dr ) 
					; eZ i-1/2
			    	aMat[3*i-1,3*i]	= -II * kz / dr 
                endif
					; eR i
			    	aMat[3*i,3*i]	= nPhi^2 / r[i] + kz^2 - w[i]^2 / c^2

                if i lt nR-1 then begin
					; eT i+1/2
			    	aMat[3*i+1,3*i]	= II * nPhi * r_[i] / ( r[i]^2 * dr ) 
					; eZ i+1/2
			    	aMat[3*i+2,3*i]	= II * kz / dr 
                endif


			endif else begin

				t3	= nPhi^2 / r[i] + kz^2 - w[i]^2 / c^2 * epsilon[0,0,i]

                if i gt 0 then begin
					; eT i-1/2
			    	aMat[3*i-2,3*i]	= -II * nPhi * r_[i-1] / ( r[i]^2 * dr ) $
                                        - w[i]^2*r_[i-1]/(2*c^2*r[i])*epsilon_[1,0,i-1]
					; eZ i-1/2
			    	aMat[3*i-1,3*i]	= -II * kz / dr $
                                        - w[i]^2*r_[i-1]/(2*c^2*r[i])*epsilon_[2,0,i-1]
                endif
					; eR i
			    	aMat[3*i,3*i]	= t3;nPhi^2 / r[i]^2 + kz^2 $
                                        ;- w[i]^2 / c^2 * epsilon[0,0,i]
                if i lt nR-1 then begin
					; eT i+1/2
			    	aMat[3*i+1,3*i]	= II * nPhi * r_[i] / ( r[i]^2 * dr ) $
                                        - w[i]^2*r_[i]/(2*c^2*r[i])*epsilon_[1,0,i]
					; eZ i+1/2
			    	aMat[3*i+2,3*i]	= II * kz / dr $
                                        - w[i]^2*r_[i]/(2*c^2*r[i])*epsilon_[2,0,i]
                endif

            endelse


		;	phi component
			if i lt nR-1 then begin
           
                if bandStorage then begin
            
                    j_=3*i+1
                    i_=3*i-2
                    if i gt 0 then $
				    aMat_bandStorage[i_,nuca+j_-i_] = -r_[i-1] / ( r[i] * dr^2 ) 
                    i_=3*i
				    aMat_bandStorage[i_,nuca+j_-i_]	= -II*nPhi/(r[i]*dr) $
                                        - w[i]^2*r[i]/(2*c^2*r_[i])*epsilon[0,1,i]
                    i_=3*i+1
				    aMat_bandStorage[i_,nuca+j_-i_]	= kz^2 + r_[i]/(r[i+1]*dr^2) $
                                            + r_[i]/(r[i]*dr^2) $
                                            - w[i]^2/c^2*epsilon_[1,1,i]
                    i_=3*i+2
				    aMat_bandStorage[i_,nuca+j_-i_]	= -kz*nPhi/r_[i] - w[i]^2/(c^2)*epsilon_[2,1,i]
                    i_=3*i+3
				    aMat_bandStorage[i_,nuca+j_-i_]	= II * nPhi / (r[i+1]*dr) $
                                            - w[i]^2*r[i+1]/(2*c^2*r_[i])*epsilon[0,1,i+1]
                    i_=3*i+4
				    if i lt nR-2 then $
				    aMat_bandStorage[i_,nuca+j_-i_]	= -r_[i+1]/(r[i+1]*dr^2)

                endif else if kjInput then begin

                    if i gt 0 then $
				    aMat[3*i-2,3*i+1] = -r_[i-1] / ( r[i] * dr^2 ) 
				    aMat[3*i,3*i+1]	= -II*nPhi/(r[i]*dr) 
				    aMat[3*i+1,3*i+1]	= kz^2 + r_[i]/(r[i+1]*dr^2) $
                                            + r_[i]/(r[i]*dr^2) $
                                            - w[i]^2/c^2
				    aMat[3*i+2,3*i+1]	= -kz*nPhi/r_[i] 
				    aMat[3*i+3,3*i+1]	= II * nPhi / (r[i+1]*dr) 
				    if i lt nR-2 then $
				    aMat[3*i+4,3*i+1]	= -r_[i+1]/(r[i+1]*dr^2)
 
				endif else begin

                    if i gt 0 then $
				    aMat[3*i-2,3*i+1] = -r_[i-1] / ( r[i] * dr^2 ) 
				    aMat[3*i,3*i+1]	= -II*nPhi/(r[i]*dr) $
                                        - w[i]^2*r[i]/(2*c^2*r_[i])*epsilon[0,1,i]
				    aMat[3*i+1,3*i+1]	= kz^2 + r_[i]/(r[i+1]*dr^2) $
                                            + r_[i]/(r[i]*dr^2) $
                                            - w[i]^2/c^2*epsilon_[1,1,i]
				    aMat[3*i+2,3*i+1]	= -kz*nPhi/r_[i] - w[i]^2/(c^2)*epsilon_[2,1,i]
				    aMat[3*i+3,3*i+1]	= II * nPhi / (r[i+1]*dr) $
                                            - w[i]^2*r[i+1]/(2*c^2*r_[i])*epsilon[0,1,i+1]
				    if i lt nR-2 then $
				    aMat[3*i+4,3*i+1]	= -r_[i+1]/(r[i+1]*dr^2)
                
                endelse

		;	z component	


                if bandStorage then begin
                    j_=3*i+2
                    i_=3*i-1
                    if i gt 0 then $
				    aMat_bandStorage[i_,nuca+j_-i_] = -r[i]/(r_[i]*dr^2)
                    i_=3*i
				    aMat_bandStorage[i_,nuca+j_-i_] = -II*kz*r[i]/(r_[i]*dr)$
                                        -w[i]^2*r[i]/(2*c^2*r_[i])*epsilon[0,2,i]
                    i_=3*i+1
				    aMat_bandStorage[i_,nuca+j_-i_] = -nPhi*kz/r_[i] - w[i]^2/c^2*epsilon_[1,2,i]
                    i_=3*i+2
				    aMat_bandStorage[i_,nuca+j_-i_] = (r[i+1]+r[i])/(r_[i]*dr^2) + nPhi^2/r_[i]^2 $
                                            - w[i]^2/c^2*epsilon_[2,2,i]
                    i_=3*i+3
				    aMat_bandStorage[i_,nuca+j_-i_] = II * kz * r[i+1] / (r_[i]*dr) $
                                            - w[i]^2*r[i+1]/(2*c^2*r_[i])*epsilon[0,2,i+1]
                    i_=3*i+5
				    if i lt nR-2 then $
				    aMat_bandStorage[i_,nuca+j_-i_] = -r[i+1] / ( r_[i] * dr^2 ) 

                endif else if kjInput then begin

					t1	= -r[i]/(r_[i]*dr^2)
					t2	= -II*kz*r[i]/(r_[i]*dr) ; Propagate this "dr" correction
					t3	= -nPhi*kz/r_[i] ; Also propagate this "r_[i]" correction
					t4	= (r[i+1]+r[i])/(r_[i]*dr^2) + nPhi^2/r_[i]^2 - w[i]^2/c^2 ; Also propagate this "r_[i]^2" correction
					t5	= II * kz * r[i+1] / (r_[i]*dr) ; And this "dr" correction 
					t6	= -r[i+1] / ( r_[i] * dr^2 ) 

	                if i gt 0 then $
				    aMat[3*i-1,3*i+2] = t1
				    aMat[3*i,3*i+2]   = t2 
				    aMat[3*i+1,3*i+2] = t3
				    aMat[3*i+2,3*i+2] = t4
				    aMat[3*i+3,3*i+2] = t5
				    if i lt nR-2 then $
				    aMat[3*i+5,3*i+2] = t6 

				endif else begin
               
					t1	= -r[i]/(r_[i]*dr^2)
					t2	= -II*kz*r[i]/(r_[i]*dr)$
                                        -w[i]^2*r[i]/(2*c^2*r_[i])*epsilon[0,2,i]
					t3	= -nPhi*kz/r_[i] - w[i]^2/c^2*epsilon_[1,2,i]
					t4	= (r[i+1]+r[i])/(r_[i]*dr^2) + nPhi^2/r_[i]^2 $
                                            - w[i]^2/c^2*epsilon_[2,2,i]
					t5	= II * kz * r[i+1] / (r_[i]*dr) $
                                            - w[i]^2*r[i+1]/(2*c^2*r_[i])*epsilon[0,2,i+1]
					t6	= -r[i+1] / ( r_[i] * dr^2 ) 

	                if i gt 0 then $
				    aMat[3*i-1,3*i+2] = t1;-r[i]/(r_[i]*dr^2)
				    aMat[3*i,3*i+2] = t2;-II*kz*r[i]/(r_[i]*dr)$
                                        ;-w[i]^2*r[i]/(2*c^2*r_[i])*epsilon[0,2,i]
				    aMat[3*i+1,3*i+2] = t3;-nPhi*kz/r_[i] - w[i]^2/c^2*epsilon_[1,2,i]
				    aMat[3*i+2,3*i+2] = t4;(r[i+1]+r[i])/(r_[i]*dr^2) + nPhi^2/r_[i]^2 $
                                           ; - w[i]^2/c^2*epsilon_[2,2,i]
				    aMat[3*i+3,3*i+2] = t5;II * kz * r[i+1] / (r_[i]*dr) $
                                          ;  - w[i]^2*r[i+1]/(2*c^2*r_[i])*epsilon[0,2,i+1]
				    if i lt nR-2 then $
				    aMat[3*i+5,3*i+2] = t6;-r[i+1] / ( r_[i] * dr^2 ) 

                endelse
	
			endif

	endfor
	
	if bandStorage then aMat = aMat_bandStorage

end
