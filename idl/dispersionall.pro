pro dispersionAll, wReal, epsilon, stixVars, runData, specData, $
	kR = kR__, kPhi = kPhi, kz = kz__

	common switches
	@dlg_constants
	common plotSwitches
	common dlg_colors
	common writeSwitches

	if dispersion_generalised then begin

;	Generalised dispersion analysis

	kR__	= dcomplexArr ( runData.nR, 4 )
	kPhi	= runData.nPhi / runData.r

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

	if plotDispersionGeneral then begin

		kRPlot	= real_part ( kR__ )

		p=plot(runData.r,kRPlot[*,0],symbol="D",lineStyle='',sym_size=0.5)
		p=plot(runData.r,kRPlot[*,1],symbol="D",lineStyle='',/over,sym_size=0.5)
		p=plot(runData.r,kRPlot[*,2],symbol="D",lineStyle='',/over,sym_size=0.5)
		p=plot(runData.r,kRPlot[*,3],symbol="D",lineStyle='',/over,sym_size=0.5)

		kRPlot	= imaginary ( kR__ )

		p=plot(runData.r,kRPlot[*,0],symbol="D",lineStyle='',/over,color='r',sym_size=0.5)
		p=plot(runData.r,kRPlot[*,1],symbol="D",lineStyle='',/over,color='r',sym_size=0.5)
		p=plot(runData.r,kRPlot[*,2],symbol="D",lineStyle='',/over,color='r',sym_size=0.5)
		p=plot(runData.r,kRPlot[*,3],symbol="D",lineStyle='',/over,color='r',sym_size=0.5)

	endif	

	endif

end
