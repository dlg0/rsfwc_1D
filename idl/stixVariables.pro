pro stixVariables, wReal, runData, specData, $
	stixVars = stixVars

	stixR	= 1d0 - total ( specData.wp^2 / ( wReal * ( wReal + specData.wc ) ),2)
	stixL	= 1d0 - total ( specData.wp^2 / ( wReal * ( wReal - specData.wc ) ),2)
	stixS	= 0.5d0 * ( stixR + stixL )
	stixD	= 0.5d0 * ( stixR - stixL )
	stixP	= 1d0 - total ( specData.wp^2 / wReal^2, 2 )

	stixR_	= 1d0 - total ( specData.wp_^2 / ( wReal * ( wReal + specData.wc_ ) ),2)
	stixL_	= 1d0 - total ( specData.wp_^2 / ( wReal * ( wReal - specData.wc_ ) ),2)
	stixS_	= 0.5d0 * ( stixR_ + stixL_ )
	stixD_	= 0.5d0 * ( stixR_ - stixL_ )
	stixP_	= 1d0 - total ( specData.wp_^2 / wReal^2, 2 )

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

end
