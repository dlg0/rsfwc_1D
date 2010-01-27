pro grad_freq_scan

	@load_colors

	nn	= 50.0
	nkz	= 10.0
	ng	= 10.0

	;mag	= fltArr ( nn, nkz )
	mag	= fltArr ( ng, nkz )


	kzStart	= 63.0
	kzEnd	= 66.0
	kzStep	= ( kzEnd - kzStart ) / nkz

	nStart	= 1.0 * 1d17
	nEnd	= 20.0 * 1d17
	nStep	= ( nEnd - nStart ) / nn 

	gStart	= 0.1
	gEnd	= 0.3
	gStep	= ( gEnd - gStart ) / ng 

	kz	= fIndGen ( nkz ) * kzStep + kzStart
	n	= fIndGen ( nn ) * nStep + nStart
	g	= fIndGen ( ng ) * gStep + gStart

	cnt	= 0
	;for i=0.0,nn-1 do begin
	;for j=0.0,ng-1 do begin
	for k=0.0,nkz-1 do begin

		rsfwc_1d, eR = eR, rFull = r, rHalf = r_, $
			ez = ez, ePhi = ePhi, $
			;nMax = i * nStep + nStart, $
			kz = k * kzStep + kzStart,$
			;gradSize = j * gStep + gStart, $
			iiAnt = iiAnt

		;iiRange	= where ( abs ( r_ - 1.17 ) eq min ( abs ( r_ - 1.17 ) ) ) 
		;;mag[i,k]	= max((abs(ePhi))[iiRange]) ;/ abs(ePhi[iiAnt])
		;mag[j,k]	= max((abs(ePhi))[iiRange[0]]) / abs(ePhi[iiAnt])

		if cnt eq 0 then begin
			iPlot, r_, abs(ePhi) / abs ( ePhi[iiAnt[0]] )
		endif else begin
			iPlot, r_, abs(ePhi) / abs ( ePhi[iiAnt[0]] ), $
				trans = 70, /over
		endelse
		cnt++
	;endfor
	;endfor
	endfor

stop

	;levels	= fIndGen ( 20 ) * 0.005;10.0e-5
	nLevs	= 10
	levels	= 2.0^fIndGen ( nLevs ) * 1e-3
	colors	= bytScl ( fIndGen ( nLevs ), top = 253 ) + 1
	device, decomposed = 0
	loadct, 1
	;contour, mag, n, kz, levels = levels, c_colors = colors, /fill
	contour, mag, g, kz, levels = levels, c_colors = colors, /fill

stop
end
