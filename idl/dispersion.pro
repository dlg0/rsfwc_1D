pro dispersion, wReal, epsilon, stixVars, runData, specData

	common switches
	common constants
	common plotSwitches

	kPar	= runData.nPhi / runData.r
	nPar	= kPar * c / wReal

    if dispersion_freeSpace then begin

        k   = wReal / c
        kPer1   = sqrt ( k^2 - kPar^2 )
        kPer2   = kPer1 * 0
        iPlot, r, kPer1, $
            title = 'free space dispersion', $
            /view_next

    endif else begin

		kPerp0	= complexArr ( runData.nR, 4 )
		for i=0L,runData.nR-1L do begin
			c0	= -stixVars.stixD[i]^2 * stixVars.stixP[i] + nPar[i]^4 * stixVars.stixP[i] $
					- 2.0 * nPar[i]^2 * stixVars.stixP[i] * stixVars.stixS[i] $
					+ stixVars.stixP[i] * stixVars.stixS[i]^2
			c1	= 0
			c2	= stixVars.stixD[i]^2 + nPar[i]^2 * stixVars.stixP[i] + nPar[i]^2 * stixVars.stixS[i] $
					- stixVars.stixP[i] * stixVars.stixS[i] - stixVars.stixS[i]^2
			c3	= 0
			c4	= stixVars.stixS[i]
			c_	= [ c0, c1, c2, c3, c4 ]
			kPerp0[i,*]	= imsl_zeroPoly ( c_, /double ) * wReal / c

		endfor
	

		;	poloidal field dispersion analysis

		kPhi	= runData.nPhi / runData.r

		nPhi_	= kPhi * c / wReal	
		bPol	= sqrt ( runData.bField[*,0]^2 + runData.bField[*,2]^2 ) / runData.bMag
		bTor	= runData.bField[*,1] / runData.bMag

		c0		= -stixVars.stixD^2 * stixVars.stixP $
					+ nPhi_^4 * stixVars.stixP / bTor^4 $
					- 2 * nPhi_^2 * stixVars.stixP * stixVars.stixS / bTor^2 $
					+ stixVars.stixP * stixVars.stixS^2
		c1		= 4 * bPol * nPhi_^3 * stixVars.stixP / bTor^4 $
					- 4 * bPol * nPhi_ * stixVars.stixP * stixVars.stixS / bTor^2
		c2		= stixVars.stixD^2 $
					+ 6 * bPol^2 * nPhi_^2 * stixVars.stixP / bTor^4 $
					+ nPhi_^2 * stixVars.stixP / bTor^2 $
					+ nPhi_^2 * stixVars.stixS / bTor^2 $
					- stixVars.stixP * stixVars.stixS $
					- 2 * bPol^2 * stixVars.stixP * stixVars.stixS / bTor^2 $
					- stixVars.stixS^2
		c3		= 4 * bPol^3 * nPhi_ * stixVars.stixP / bTor^4 $
					+ 2 * bPol * nPhi_ * stixVars.stixP / bTor^2 $
					+ 2 * bPol * nPhi_ * stixVars.stixS / bTor^2
		c4		= bPol^4 * stixVars.stixP / bTor^4 $
					+ bPol^2 * stixVars.stixP / bTor^2 $
					+ stixVars.stixS $
					+ bPol^2 * stixVars.stixS / bTor^2

		kPerp__	= complexArr ( runData.nR, 4 )	
		for i = 0L, runData.nR - 1L do begin	
	
			c_	= [ c0[i], c1[i], c2[i], c3[i], c4[i] ]
			kPerp__[i,*]	= imsl_zeroPoly ( c_, /double ) * wReal / c

		endfor

		sqrtArg	= kPerp__^2 - runData.kz^2
		kR	= sqrt ( sqrtArg )
		kR2	= -sqrt ( -sqrtArg)


    endelse

;	Generalised dispersion analysis

	kR__	= complexArr ( runData.nR, 4 )
	print, 'Running generalised dispersion calculation ...'
	for i=0L,runData.nR-1L do begin
		
		k4	= wReal^2/c^2 * epsilon[0,0,i]
		k3	= wReal^2/c^2 $
				* ( kPhi[i] * ( epsilon[0,1,i] + epsilon[1,0,i] ) $
						+ runData.kz * ( epsilon[0,2,i] + epsilon[2,0,i] ) )
		k2	= wReal^2/c^4 * ( c^2 * kPhi[i]^2 * ( epsilon[0,0,i] + epsilon[1,1,i] ) $
						+ c^2 * runData.kz * kPhi[i] * ( epsilon[1,2,i] + epsilon[2,1,i] ) $
						+ c^2 * runData.kz^2 * ( epsilon[0,0,i] + epsilon[2,2,i] ) $
						+ wReal^2 * ( epsilon[0,1,i] * epsilon[1,0,i] $
										+ epsilon[0,2,i] * epsilon[2,0,i] $
										- epsilon[0,0,i] * ( epsilon[1,1,i] + epsilon[2,2,i] ) ) )
		k1	= wReal^2/c^4 * ( c^2 * runData.kz^2 * kPhi[i] * ( epsilon[0,1,i] + epsilon[1,0,i] ) $
						+ c^2 * runData.kz^3 * ( epsilon[0,2,i] + epsilon[2,0,i] ) $
						+ runData.kz * ( c^2 * kPhi[i]^2 * ( epsilon[0,2,i] + epsilon[2,0,i] ) $
							+ wReal^2 * ( epsilon[0,1,i] * epsilon[1,2,i] $
										- epsilon[1,1,i] * ( epsilon[0,2,i] + epsilon[2,0,i] ) $
										+ epsilon[1,0,i] * epsilon[2,1,i] ) ) $
						+ kPhi[i] * ( wReal^2 * ( epsilon[1,2,i] * epsilon[2,0,i] $
												+ epsilon[0,2,i] * epsilon[2,1,i] ) $
									+ ( epsilon[0,1,i] + epsilon[1,0,i] ) $
										* ( c^2 * kPhi[i]^2 - wReal^2 * epsilon[2,2,i] ) ) )
		k0	= wReal^2/c^6 * ( c^4 * kPhi[i]^4 * epsilon[1,1,i] $
							+ c^4 * runData.kz * kPhi[i]^3 * ( epsilon[1,2,i] + epsilon[2,1,i] ) $
							+ c^2 * runData.kz * kPhi[i] * ( c^2 * runData.kz^2 * ( epsilon[1,2,i] + epsilon[2,1,i] ) $
									+ wReal^2 * ( epsilon[0,2,i] * epsilon[1,0,i] $
													+ epsilon[0,1,i] * epsilon[2,0,i] $
													- epsilon[0,0,i] * $
														( epsilon[1,2,i] + epsilon[2,1,i] ) ) ) $
							+ c^4 * runData.kz^4 * epsilon[2,2,i] $
							+ c^2 * wReal^2 * runData.kz^2 $
								* ( epsilon[0,2,i] * epsilon[2,0,i] $
									+ epsilon[1,2,i] * epsilon[2,1,i] $
									- ( epsilon[0,0,i] + epsilon[1,1,i] ) * epsilon[2,2,i] ) $
							+wReal^4 * ( epsilon[0,2,i] * ( -epsilon[1,1,i] * epsilon[2,0,i] $
															+ epsilon[1,0,i] * epsilon[2,1,i] ) $
								+ epsilon[0,1,i] * ( epsilon[1,2,i] * epsilon[2,0,i] $
													- epsilon[1,0,i] * epsilon[2,2,i] ) $
								+ epsilon[0,0,i] * ( -epsilon[1,2,i] * epsilon[2,1,i] $
													+ epsilon[1,1,i] * epsilon[2,2,i] ) ) $ 
							+ c^2 * kPhi[i]^2 * ( c^2 * runData.kz^2 * $
								( epsilon[1,1,i] + epsilon[2,2,i] ) $
								+ wReal^2 * ( epsilon[0,1,i] * epsilon[1,0,i] $
											+ epsilon[1,2,i] * epsilon[2,1,i] $
											- epsilon[1,1,i] * $
												( epsilon[0,0,i] + epsilon[2,2,i] ) ) ) ) 

		c_	= [ k0, k1, k2, k3, k4 ]
		kR__[i,*]	= imsl_zeroPoly ( c_, /double ) 


	endfor

	if plotDispersionGeneral then begin

		iPlot, runData.r, kR__[*,0], sym_index = 4, lineStyle = 6, yRange = [-500,500]
		iPlot, runData.r, kR__[*,1], sym_index = 4, lineStyle = 6, /over
		iPlot, runData.r, kR__[*,2], sym_index = 4, lineStyle = 6, /over
		iPlot, runData.r, kR__[*,3], sym_index = 4, lineStyle = 6, /over
		iPlot, runData.r, imaginary(kR__[*,0]), sym_index = 3, lineStyle = 6, /over
		iPlot, runData.r, imaginary(kR__[*,1]), sym_index = 3, lineStyle = 6, /over
		iPlot, runData.r, imaginary(kR__[*,2]), sym_index = 3, lineStyle = 6, /over
		iPlot, runData.r, imaginary(kR__[*,3]), sym_index = 3, lineStyle = 6, /over
	
			range	= max ( abs( abs ( kR__ ) ) )
			for i = 0, runData.nSpec - 1 do begin
				for harm = 1, 2 do begin
	
				iiRes	= where ( abs ( specData[i].wc*harm - wReal ) $
							eq min ( abs ( specData[i].wc*harm - wReal ) ) )
	
				if iiRes ne 0 and iiRes ne runData.nR-1 then begin
	
					print, 'Cyclotron resonance found'
					iPlot, [ runData.r[iiRes], runData.r[iiRes] ], [ -range, range ],$
						thick = 10, trans = 50, /over
	
				endif			
				endfor
			endfor

	endif	

	if  plotDispersionJaeger then begin

		iPlot, runData.r, kPerp__[*,0], sym_index = 4, lineStyle = 6, yRange = [-500, 500]
		iPlot, runData.r, kPerp__[*,1], sym_index = 4, lineStyle = 6, /over
		iPlot, runData.r, kPerp__[*,2], sym_index = 4, lineStyle = 6, /over
		iPlot, runData.r, kPerp__[*,3], sym_index = 4, lineStyle = 6, /over
		iPlot, runData.r, imaginary(kPerp__[*,0]), sym_index = 3, lineStyle = 6, /over
		iPlot, runData.r, imaginary(kPerp__[*,1]), sym_index = 3, lineStyle = 6, /over
		iPlot, runData.r, imaginary(kPerp__[*,2]), sym_index = 3, lineStyle = 6, /over
		iPlot, runData.r, imaginary(kPerp__[*,3]), sym_index = 3, lineStyle = 6, /over

	endif

	if plotDispersionNoPol then begin

		iPlot, runData.r, kPerp0[*,0], sym_index = 4, lineStyle = 6
		iPlot, runData.r, kPerp0[*,1], sym_index = 4, lineStyle = 6, /over
		iPlot, runData.r, kPerp0[*,2], sym_index = 4, lineStyle = 6, /over
		iPlot, runData.r, kPerp0[*,3], sym_index = 4, lineStyle = 6, /over
		iPlot, runData.r, imaginary(kPerp0[*,0]), sym_index = 3, lineStyle = 6, /over
		iPlot, runData.r, imaginary(kPerp0[*,1]), sym_index = 3, lineStyle = 6, /over
		iPlot, runData.r, imaginary(kPerp0[*,2]), sym_index = 3, lineStyle = 6, /over
		iPlot, runData.r, imaginary(kPerp0[*,3]), sym_index = 3, lineStyle = 6, /over

	endif

end
