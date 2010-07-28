pro stixVariables, wReal, specData, $
	stixVars = stixVars, $
    noHalf = noHalf

    sumDim  = size ( specData.wp, /n_dim )

	stixR	= 1d0 - total ( specData.wp^2 / ( wReal * ( wReal + specData.wc ) ),sumDim)
	stixL	= 1d0 - total ( specData.wp^2 / ( wReal * ( wReal - specData.wc ) ),sumDim)
	stixS	= 0.5d0 * ( stixR + stixL )
	stixD	= 0.5d0 * ( stixR - stixL )
	stixP	= 1d0 - total ( specData.wp^2 / wReal^2, sumDim )


    if not keyword_set ( noHalf ) then begin 

	stixR_	= 1d0 - total ( specData.wp_^2 / ( wReal * ( wReal + specData.wc_ ) ),sumDim)
	stixL_	= 1d0 - total ( specData.wp_^2 / ( wReal * ( wReal - specData.wc_ ) ),sumDim)
	stixS_	= 0.5d0 * ( stixR_ + stixL_ )
	stixD_	= 0.5d0 * ( stixR_ - stixL_ )
	stixP_	= 1d0 - total ( specData.wp_^2 / wReal^2, sumDim )

    stixVars = { stixR : stixR, $
					stixL : stixL, $
					stixS : stixS, $
					stixD : stixD, $
					stixP : stixP, $
					stixR_ : stixR_, $
					stixL_ : stixL_, $
					stixS_ : stixS_, $
					stixD_ : stixD_, $
					stixP_ : stixP_ }

    endif else begin

            stixVars = { stixR : stixR, $
					stixL : stixL, $
					stixS : stixS, $
					stixD : stixD, $
					stixP : stixP }

    endelse

end
