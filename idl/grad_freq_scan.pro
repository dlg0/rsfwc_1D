pro grad_freq_scan

	@load_colors

	for i=5.0,7.0,0.1 do begin

		print, i
		rsfwc_1d, eR = eR, rFull = r, rHalf = r_, $
			ez = ez, ePhi = ePhi, nMax = i*1e18

		if i eq 0 then begin
			iPlot, r_, abs(ePhi)
		endif else begin
			iPlot, r_, abs(ePhi), $
				trans = 70, /over

		endelse

	endfor

stop
end
