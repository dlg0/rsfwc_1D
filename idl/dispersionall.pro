pro dispersionAll, wReal, epsilon, stixVars, runData, specData, $
	kR = kR__, kPhi = kPhi, kz = kz__

	common switches
	@dlg_constants
	common plotSwitches
	common dlg_colors
	common writeSwitches

    if dispersion_freeSpace then begin

        kPar	= runData.nPhi / runData.r
	    nPar	= kPar * _c / wReal
        k   = wReal / _c
        kPer1   = sqrt ( k^2 - kPar^2 )
        kPer2   = kPer1 * 0
        iPlot, r, kPer1, $
            title = 'free space dispersion', $
            /view_next

    endif 
    

	if dispersion_generalised then begin

;	Generalised dispersion analysis

	kR__	= dcomplexArr ( runData.nR, 4 )
	kPhi	= runData.nPhi / runData.r

	;print, 'Running generalised dispersion calculation ...'
	for i=0L,runData.nR-1L do begin
		
		k4	= wReal^2/_c^2 * epsilon[0,0,i]
		k3	= wReal^2/_c^2 $
				* ( kPhi[i] * ( epsilon[0,1,i] + epsilon[1,0,i] ) $
						+ runData.kz * ( epsilon[0,2,i] + epsilon[2,0,i] ) )
		k2	= wReal^2/_c^4 * ( _c^2 * kPhi[i]^2 * ( epsilon[0,0,i] + epsilon[1,1,i] ) $
						+ _c^2 * runData.kz * kPhi[i] * ( epsilon[1,2,i] + epsilon[2,1,i] ) $
						+ _c^2 * runData.kz^2 * ( epsilon[0,0,i] + epsilon[2,2,i] ) $
						+ wReal^2 * ( epsilon[0,1,i] * epsilon[1,0,i] $
										+ epsilon[0,2,i] * epsilon[2,0,i] $
										- epsilon[0,0,i] * ( epsilon[1,1,i] + epsilon[2,2,i] ) ) )
		k1	= wReal^2/_c^4 * ( _c^2 * runData.kz^2 * kPhi[i] * ( epsilon[0,1,i] + epsilon[1,0,i] ) $
						+ _c^2 * runData.kz^3 * ( epsilon[0,2,i] + epsilon[2,0,i] ) $
						+ runData.kz * ( _c^2 * kPhi[i]^2 * ( epsilon[0,2,i] + epsilon[2,0,i] ) $
							+ wReal^2 * ( epsilon[0,1,i] * epsilon[1,2,i] $
										- epsilon[1,1,i] * ( epsilon[0,2,i] + epsilon[2,0,i] ) $
										+ epsilon[1,0,i] * epsilon[2,1,i] ) ) $
						+ kPhi[i] * ( wReal^2 * ( epsilon[1,2,i] * epsilon[2,0,i] $
												+ epsilon[0,2,i] * epsilon[2,1,i] ) $
									+ ( epsilon[0,1,i] + epsilon[1,0,i] ) $
										* ( _c^2 * kPhi[i]^2 - wReal^2 * epsilon[2,2,i] ) ) )
		k0	= wReal^2/_c^6 * ( _c^4 * kPhi[i]^4 * epsilon[1,1,i] $
							+ _c^4 * runData.kz * kPhi[i]^3 * ( epsilon[1,2,i] + epsilon[2,1,i] ) $
							+ _c^2 * runData.kz * kPhi[i] $
                                * ( _c^2 * runData.kz^2 * ( epsilon[1,2,i] + epsilon[2,1,i] ) $
									+ wReal^2 * ( epsilon[0,2,i] * epsilon[1,0,i] $
													+ epsilon[0,1,i] * epsilon[2,0,i] $
													- epsilon[0,0,i] * $
														( epsilon[1,2,i] + epsilon[2,1,i] ) ) ) $
							+ _c^4 * runData.kz^4 * epsilon[2,2,i] $
							+ _c^2 * wReal^2 * runData.kz^2 $
								* ( epsilon[0,2,i] * epsilon[2,0,i] $
									+ epsilon[1,2,i] * epsilon[2,1,i] $
									- ( epsilon[0,0,i] + epsilon[1,1,i] ) * epsilon[2,2,i] ) $
							+wReal^4 * ( epsilon[0,2,i] * ( -epsilon[1,1,i] * epsilon[2,0,i] $
															+ epsilon[1,0,i] * epsilon[2,1,i] ) $
								+ epsilon[0,1,i] * ( epsilon[1,2,i] * epsilon[2,0,i] $
													- epsilon[1,0,i] * epsilon[2,2,i] ) $
								+ epsilon[0,0,i] * ( -epsilon[1,2,i] * epsilon[2,1,i] $
													+ epsilon[1,1,i] * epsilon[2,2,i] ) ) $ 
							+ _c^2 * kPhi[i]^2 * ( _c^2 * runData.kz^2 * $
								( epsilon[1,1,i] + epsilon[2,2,i] ) $
								+ wReal^2 * ( epsilon[0,1,i] * epsilon[1,0,i] $
											+ epsilon[1,2,i] * epsilon[2,1,i] $
											- epsilon[1,1,i] * $
												( epsilon[0,0,i] + epsilon[2,2,i] ) ) ) ) 

		c_	= [ k0, k1, k2, k3, k4 ]

		if noIMSL then begin
			kR__[i,*]	= fz_roots ( c_, /double ) 
		endif else begin
			kR__[i,*]	= imsl_zeroPoly ( c_, /double, /jenkins ) 
		endelse

	endfor

	kz__	= fltArr ( runData.nR ) + runData.kz

	if plotDispersionGeneral then begin

		kRPlot	= real_part ( kR__ )

		p=plot(runData.r,kRPlot[*,0],symbol="D",lineStyle='',sym_size=0.5)
		p=plot(runData.r,kRPlot[*,1],symbol="D",lineStyle='',/over,sym_size=0.5)
		p=plot(runData.r,kRPlot[*,2],symbol="D",lineStyle='',/over,sym_size=0.5)
		p=plot(runData.r,kRPlot[*,3],symbol="D",lineStyle='',/over,sym_size=0.5)

		kRPlot	= imaginary ( kR__ )

		p=plot(runData.r,kRPlot[*,0],symbol="D",lineStyle='',/over,color='b',sym_size=0.5)
		p=plot(runData.r,kRPlot[*,1],symbol="D",lineStyle='',/over,color='b',sym_size=0.5)
		p=plot(runData.r,kRPlot[*,2],symbol="D",lineStyle='',/over,color='b',sym_size=0.5)
		p=plot(runData.r,kRPlot[*,3],symbol="D",lineStyle='',/over,color='b',sym_size=0.5)

		if plotKdotB then begin	
			kMag1	= sqrt ( real_part ( kR__[*,0])^2 + kPhi^2 + kz__^2 )
			kDotB1	= real_part(kR__[*,0]) * runData.bField[*,0] $
						+ kPhi * runData.bField[*,1] $
						+ kz__ * runData.bField[*,2]
			theta1	= aCos ( kDotB1 / ( kMag1 * runData.bMag ) )

			kMag1_	= sqrt ( imaginary ( kR__[*,0])^2 + kPhi^2 + kz__^2 )
			kDotB1_	= imaginary(kR__[*,0]) * runData.bField[*,0] $
						+ kPhi * runData.bField[*,1] $
						+ kz__ * runData.bField[*,2]
			theta1_	= aCos ( kDotB1_ / ( kMag1_ * runData.bMag ) )


			kMag2	= sqrt ( real_part ( kR__[*,1])^2 + kPhi^2 + kz__^2 )
			kDotB2	= real_part(kR__[*,1]) * runData.bField[*,0] $
						+ kPhi * runData.bField[*,1] $
						+ kz__ * runData.bField[*,2]
			theta2	= aCos ( kDotB2 / ( kMag2 * runData.bMag ) )

			kMag2_	= sqrt ( imaginary ( kR__[*,1])^2 + kPhi^2 + kz__^2 )
			kDotB2_	= imaginary(kR__[*,1]) * runData.bField[*,0] $
						+ kPhi * runData.bField[*,1] $
						+ kz__ * runData.bField[*,2]
			theta2_	= aCos ( kDotB2_ / ( kMag2_ * runData.bMag ) )


			kMag3	= sqrt ( real_part ( kR__[*,2])^2 + kPhi^2 + kz__^2 )
			kDotB3	= real_part(kR__[*,2]) * runData.bField[*,0] $
						+ kPhi * runData.bField[*,1] $
						+ kz__ * runData.bField[*,2]
			theta3	= aCos ( kDotB3 / ( kMag3 * runData.bMag ) )

			kMag3_	= sqrt ( imaginary ( kR__[*,2])^2 + kPhi^2 + kz__^2 )
			kDotB3_	= imaginary(kR__[*,2]) * runData.bField[*,0] $
						+ kPhi * runData.bField[*,1] $
						+ kz__ * runData.bField[*,2]
			theta3_	= aCos ( kDotB3_ / ( kMag3_ * runData.bMag ) )


			kMag4	= sqrt ( real_part ( kR__[*,3])^2 + kPhi^2 + kz__^2 )
			kDotB4	= real_part(kR__[*,3]) * runData.bField[*,0] $
						+ kPhi * runData.bField[*,1] $
						+ kz__ * runData.bField[*,2]
			theta4	= aCos ( kDotB4 / ( kMag4 * runData.bMag ) )

			kMag4_	= sqrt ( imaginary ( kR__[*,3])^2 + kPhi^2 + kz__^2 )
			kDotB4_	= imaginary(kR__[*,3]) * runData.bField[*,0] $
						+ kPhi * runData.bField[*,1] $
						+ kz__ * runData.bField[*,2]
			theta4_	= aCos ( kDotB4_ / ( kMag4_ * runData.bMag ) )

			p=plot(runData.r, theta1 * !radeg, color = blue, $
					yRange = [0,180], title='k dot B')
			p=plot(runData.r, theta1_ * !radeg, /over, color = blue)

			p=plot(runData.r, theta2 * !radeg , /over, color = green)
			p=plot(runData.r, theta2_ * !radeg, /over, color = green)
	
			p=plot(runData.r, theta3 * !radeg , /over, color = red)
			p=plot(runData.r, theta3_ * !radeg, /over, color = red)

			p=plot(runData.r, theta4 * !radeg , /over, color = purple)
			p=plot(runData.r, theta4_ * !radeg, /over, color = purple)
		endif

	endif	

	endif

end
