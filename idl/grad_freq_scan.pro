pro grad_freq_scan

	nP = 1024 
	wArray = 0
	eRArray = fltArr ( nP )
	i	= 0		
	for w=20e6*2*!pi,40e6*2*!pi,0.5e6*2*!pi do begin

		print, w
		rsfwc_1d, /band, nR = nP, plot = 0, nPhi = -22, w = w, kz = 20, $
			eR = eRTmp, rFull = r, rHalf = r_, $
			ez = ezTmp, ePhi = ePhiTmp, $
			kR1 = kR1, kR2 = kR2, damping = 0.00, $
			/useEqds, /useProf

		wArray	= [ wArray, w ]
		eRArray	= [ [eRArray], [eRTmp ] ]
	
		if i eq 0 then begin
			iPlot, r[1:*], ezTmp, id = id, $
			view_grid = [2,1]	
			;iPlot, r, real_part(kR1), /view_next, /yLog, yRange = [1,1e3]
			;iPlot, r, real_part(kR2), /over
			;iPlot, r, real_part(kR1), /over, lineType = 2
			;iPlot, r, real_part(kR2), /over, lineType = 2
		endif else begin
			iPlot, r[1:*], ezTmp, overPlot = id, $
				trans = 70, view_number = 1
			;iPlot, r, real_part(kR1), view_number = 2, /yLog, yRange = [1,1e3]
			;iPlot, r, real_part(kR2), /over
			;iPlot, r, imaginary(kR1), /over, lineType = 2
			;iPlot, r, imaginary(kR2), /over, lineType = 2

		endelse

		i++	

	endfor

	wArrary	= wArray[1:*]
	eRArray	= eRArray[*,1:*]

stop
end
