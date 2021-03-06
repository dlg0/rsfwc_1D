pro stixVariables, w, w_, wc, wc_, specData, $
	stixVars = stixVars, $
    noHalf = noHalf

    sumDim  = size ( specData.wp, /n_dim )
	nSpec = n_elements ( specData )

    @dlg_constants

	stixR	= 1d0 - specData.wp^2 / ( w * ( wc + specData.wc ) )
	stixL	= 1d0 - specData.wp^2 / ( w * ( wc - specData.wc ) )

    sig1 = _II * _e0 * wc * specData.wp^2 / (wc^2 - specData.wc^2 ) ; Stix S 
    sig2 = _e0 * specData.wc * specData.wp^2 / (specData.wc^2 - wc^2 ) ; Stix D
    sig3 = _II * _e0 * specData.wp^2 / wc ; Stix P
    
	stixS	= 0.5d0 * ( stixR + stixL )
	stixD	= 0.5d0 * ( stixR - stixL )

	stixP	= 1d0 - specData.wp^2 / (w*wc)

    if not keyword_set ( noHalf ) then begin 

		stixR_	= 1d0 - specData.wp_^2 / ( w_ * ( wc_ + specData.wc_ ) )
		stixL_	= 1d0 - specData.wp_^2 / ( w_ * ( wc_ - specData.wc_ ) )

        sig1_ = _II * _e0 * wc_ * specData.wp_^2 / (wc_^2 - specData.wc_^2 ) ; Stix S 
        sig2_ = _e0 * specData.wc_ * specData.wp_^2 / (specData.wc_^2 - wc_^2 ) ; Stix D
        sig3_ = _II * _e0 * specData.wp_^2 / wc_ ; Stix P
	
		stixS_	= 0.5d0 * ( stixR_ + stixL_ )
		stixD_	= 0.5d0 * ( stixR_ - stixL_ )

		stixP_	= 1d0 - specData.wp_^2 / (w_*wc_)

    	stixVars = { stixR : stixR, $
						stixL : stixL, $
						stixS : stixS, $
						stixD : stixD, $
						stixP : stixP, $
						stixR_ : stixR_, $
						stixL_ : stixL_, $
						stixS_ : stixS_, $
						stixD_ : stixD_, $
						stixP_ : stixP_, $
                        sig1 : sig1, $
                        sig2 : sig2, $
                        sig3 : sig3, $
                        sig1_ : sig1_, $
                        sig2_ : sig2_, $
                        sig3_ : sig3_  }

    endif else begin

        stixVars = { stixR : stixR, $
				stixL : stixL, $
				stixS : stixS, $
				stixD : stixD, $
				stixP : stixP, $
                sig1 : sig1, $
                sig2 : sig2, $
                sig3 : sig3 }

    endelse

end
