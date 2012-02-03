pro dispersionAll, wReal, epsilon, stixVars, runData, specData, $
	kR = kR__, kPhi = kPhi, kz = kz__

	common switches
	common constants
	common plotSwitches
	common dlg_colors
	common writeSwitches

    if dispersion_freeSpace then begin

        kPar	= runData.nPhi / runData.r
	    nPar	= kPar * c / wReal
        k   = wReal / c
        kPer1   = sqrt ( k^2 - kPar^2 )
        kPer2   = kPer1 * 0
        iPlot, r, kPer1, $
            title = 'free space dispersion', $
            /view_next

    endif 
    
    if dispersion_noPoloidal then begin

        kR__	= complexArr ( runData.nR, 4 )
        kPar	= runData.nPhi / runData.r
	    nPar	= kPar * c / wReal

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
			kR__[i,*]	= imsl_zeroPoly ( c_, /double ) * wReal / c

		endfor

    endif    

    if dispersion_jaegerPRL then begin

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

			if noIMSL then begin
				kPerp__[i,*]	= fz_roots ( c_, /double ) * wReal / c
			endif else begin
				kPerp__[i,*]	= imsl_zeroPoly ( c_, /double ) * wReal / c
			endelse

		endfor

		sqrtArg	= kPerp__^2 - runData.kz^2
		kR	= sqrt ( sqrtArg )
		kR2	= -sqrt ( -sqrtArg)


    endif


	if dispersion_generalised then begin

;	Generalised dispersion analysis

	kR__	= dcomplexArr ( runData.nR, 4 )
	kPhi	= runData.nPhi / runData.r

	;print, 'Running generalised dispersion calculation ...'
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
							+ c^2 * runData.kz * kPhi[i] $
                                * ( c^2 * runData.kz^2 * ( epsilon[1,2,i] + epsilon[2,1,i] ) $
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

        ;print, k0, k1, k2, k3, k4
		c_	= [ k0, k1, k2, k3, k4 ]

		if noIMSL then begin
			kR__[i,*]	= fz_roots ( c_, /double ) 
		endif else begin
			kR__[i,*]	= imsl_zeroPoly ( c_, /double, /jenkins ) 
		endelse

	endfor

	kz__	= fltArr ( runData.nR ) + runData.kz

	if plotDispersionGeneral then begin

        nonLinearSF = 1.0
		kRPlot	= real_part ( kR__ )
		iiNeg	= where ( kRPlot lt 0, iiNegCnt )
		if iiNegCnt gt 0 then kRPlot[iiNeg] = -1.0 * kRPlot[iiNeg] 
		kRPlot	= ( kRPlot )^(1d0/nonLinearSF)
		if iiNegCnt gt 0 then kRPlot[iiNeg] = -1.0 * kRPlot[iiNeg]

		;iPlot, runData.r, kRPlot[*,0], sym_index = 4, lineStyle = 6, $
		;	view_grid = [2,1], color = blue, $
		;	/strech_to_fit, /zoom_on_resize
		;iPlot, runData.r, kRPlot[*,1], sym_index = 4, lineStyle = 6, /over, color = green
		;iPlot, runData.r, kRPlot[*,2], sym_index = 4, lineStyle = 6, /over, color = red
		;iPlot, runData.r, kRPlot[*,3], sym_index = 4, lineStyle = 6, /over, color = purple

		p1=plot(runData.r, kRPlot[*,0], thick = 2.0, color = blue, title='dispersion')
		p2=plot(runData.r, kRPlot[*,1], thick = 2.0, /over, color = green)
		p3=plot(runData.r, kRPlot[*,2], thick = 2.0, /over, color = red)
		p4=plot(runData.r, kRPlot[*,3], thick = 2.0, /over, color = purple)

		kRPlot_	= kRPlot
		kRPlot	= imaginary ( kR__ )
		iiNeg	= where ( kRPlot lt 0, iiNegCnt )
		if iiNegCnt gt 0 then kRPlot[iiNeg]	= -1.0 * kRPlot[iiNeg] 
		kRPlot	= ( kRPlot )^(1d0/nonLinearSF)
		if iiNegCnt gt 0 then kRPlot[iiNeg]	= -1.0 * kRPlot[iiNeg]

		p5=plot(runData.r, kRPlot[*,0], lineStyle = 5, /over, color = blue)
		p6=plot(runData.r, kRPlot[*,1], lineStyle = 5, /over, color = green)
		p7=plot(runData.r, kRPlot[*,2], lineStyle = 5, /over, color = red)
		p8=plot(runData.r, kRPlot[*,3], lineStyle = 5, /over, color = purple)
	
			range	= max ( abs( [kRPlot,kRPlot_] ) )
			for i = 0, runData.nSpec do begin
				for harm = 1, 15 do begin
	
				iiRes	= where ( abs ( abs(specData[i].wc*harm) - wReal ) $
							eq min ( abs ( abs(specData[i].wc*harm) - wReal ) ) )

					for jj=0,n_elements(iiRes)-1 do begin	
						if iiRes[jj] ne 0 and iiRes[jj] ne runData.nR-1 $
						and iiRes[jj] ne runData.nR/2 $
						and iiRes[jj] ne runData.nR/2-1 $
						then begin
	
							print, harm, 'th Cyclotron resonance found at ', runData.r[iiRes[jj]]
							p9=plot( [ runData.r[iiRes[jj]], runData.r[iiRes[jj]] ], [ -range, range ],$
								thick = 6, trans = 80, /over )
	
						endif			
					endfor
				endfor
			endfor

	
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

	if  plotDispersionJaeger then begin

		iPlot, runData.r, kPerp__[*,0], sym_index = 4, lineStyle = 6;, yRange = [-500, 500]
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


	if writeDispersionTxt then begin

		;	Write text file for comparison with GCC

		outData	= replicate ( { r : 0.0 , $
					ne_ : 0.0, $
					n1 : 0.0, $
					n2 : 0.0, $
					n3 : 0.0, $
					bField : [ 0.0, 0.0, 0.0 ] , $
					kR : complexArr ( 4 ) , $
					wpe : 0.0, $
					wce : 0.0, $
					wp1 : 0.0, $
					wc1 : 0.0, $
					wp2 : 0.0, $
					wc2 : 0.0 }, n_elements ( runData.r ) )

		outData.r = runData.r
		outData.ne_	= specData[2].n
		outData.n1	= specData[0].n
		outData.n2	= specData[1].n
		;outData.n3	= specData[2].n
		outData.bField = transpose(runData.bField)
		outData.kR	= transpose ( kR__ )
		outData.wpe	= specData[2].wp	
		outData.wce	= specData[2].wc	
		outData.wp1	= specData[0].wp	
		outData.wc1	= specData[0].wc	
		outData.wp2	= specData[1].wp	
		outData.wc2	= specData[1].wc	


		openw, lun, 'dlg_dispersion.txt', /get_lun
		for i=0,n_elements(runData.r)-1 do $
			printf, lun, outData[i], $
				format = '(f8.5,2x,4(e12.4,2x),3(f8.5,2x),8(e12.4,2x),6(e12.4,2x))'
		close, lun
					
	endif

end