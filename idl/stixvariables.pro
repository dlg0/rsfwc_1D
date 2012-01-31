pro stixVariables, wReal, specData, $
	stixVars = stixVars, $
    noHalf = noHalf, $
	electronsOnly = electronsOnly

    sumDim  = size ( specData.wp, /n_dim )
	nSpec = n_elements ( specData )

	if keyword_set(electronsOnly) then begin
		stixR	= 1d0 - specData[nSpec-1].wp^2 / ( wReal * ( wReal + specData[nSpec-1].wc ) )
		stixL	= 1d0 - specData[nSpec-1].wp^2 / ( wReal * ( wReal - specData[nSpec-1].wc ) )
	endif else begin
		stixR	= 1d0 - total ( specData.wp^2 / ( wReal * ( wReal + specData.wc ) ),sumDim)
		stixL	= 1d0 - total ( specData.wp^2 / ( wReal * ( wReal - specData.wc ) ),sumDim)
	endelse
	stixS	= 0.5d0 * ( stixR + stixL )
	stixD	= 0.5d0 * ( stixR - stixL )

	if keyword_set(electronsOnly) then begin
		stixP	= 1d0 - specData[nSpec-1].wp^2 / wReal^2
	endif else begin
		stixP	= 1d0 - total ( specData.wp^2 / wReal^2, sumDim )
	endelse

    if not keyword_set ( noHalf ) then begin 

		if keyword_set(electronsOnly) then begin
			stixR_	= 1d0 - specData[nSpec-1].wp_^2 / ( wReal * ( wReal + specData[nSpec-1].wc_ ) )
			stixL_	= 1d0 - specData[nspec-1].wp_^2 / ( wReal * ( wReal - specData[nSpec-1].wc_ ) )
		endif else begin
			stixR_	= 1d0 - total ( specData.wp_^2 / ( wReal * ( wReal + specData.wc_ ) ),sumDim)
			stixL_	= 1d0 - total ( specData.wp_^2 / ( wReal * ( wReal - specData.wc_ ) ),sumDim)
		endelse

		stixS_	= 0.5d0 * ( stixR_ + stixL_ )
		stixD_	= 0.5d0 * ( stixR_ - stixL_ )

		if keyword_set(electronsOnly) then begin
			stixP_	= 1d0 - specData[nspec-1].wp_^2 / wReal^2
		endif else begin
			stixP_	= 1d0 - total ( specData.wp_^2 / wReal^2, sumDim )
		endelse

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
