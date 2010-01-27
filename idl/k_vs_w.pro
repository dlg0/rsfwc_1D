pro k_vs_w

	@load_colors

	nf	= 400
	fStart	= 1d6
	fEnd	= 800d6
	df	= ( fEnd - fStart ) / (nf - 1)

	kRAll	= complexArr ( nf, 4 )

	rSlice	= 1.11

	w	= ( fIndGen ( nf ) * df + fStart ) * 2 * !pi


	freq	= 10d0^(fIndGen(nf)/(nf-1)*6+5)
	for i=0,nf-1 do begin
	
		rsfwc_1d, $
			/dispersionOnly, $
			freq = freq[i], $
			kR = kR, rOut = r, $
			specData = specData

		iiSlice	= where ( abs ( r - rSlice ) eq min ( abs ( r - rSlice ) ) )

		kRAll[i,*]	= kR[iiSlice[0],*]
		wpe	= specData[1].wp[iiSlice[0]]
		wce	= specData[1].wc[iiSlice[0]]
		wpi	= specData[0].wp[iiSlice[0]]
		wci	= specData[0].wc[iiSlice[0]]
		wLH	= specData[0].wLH[iiSlice[0]]
		wUH	= specData[0].wUH[iiSlice[0]]

		close, /all

	endfor

	iPlot, freq*2*!Pi/wci, kRAll[*,0]^2, sym_index = 4, lineStyle = 6, color = blue, $
		/xLog, /zoom_on_resize, /stretch_to_fit, symSize = 0.5
	iPlot, freq*2*!Pi/wci, kRAll[*,1]^2, sym_index = 4, lineStyle = 6, /over, color = green, symSize = 0.5
	iPlot, freq*2*!Pi/wci, kRAll[*,2]^2, sym_index = 4, lineStyle = 6, /over, color = red, symSize = 0.5
	iPlot, freq*2*!Pi/wci, kRAll[*,3]^2, sym_index = 4, lineStyle = 6, /over, color = purple, symSize = 0.5

	iPlot, freq*2*!Pi/wci, imaginary(kRAll[*,0])^2, sym_index = 3, lineStyle = 6, color = blue, /over
	iPlot, freq*2*!Pi/wci, imaginary(kRAll[*,1])^2, sym_index = 3, lineStyle = 6, /over, color = green
	iPlot, freq*2*!Pi/wci, imaginary(kRAll[*,2])^2, sym_index = 3, lineStyle = 6, /over, color = red
	iPlot, freq*2*!Pi/wci, imaginary(kRAll[*,3])^2, sym_index = 3, lineStyle = 6, /over, color = purple

	yPlot	= [ -max(abs(kRAll)),max(abs(kRAll)) ]
	iPlot, [wLH,wLH]/wci, yPlot, /over, thick = 4, trans = 70, name = 'wLH', /insert_legend 
	iPlot, [wUH,wUH]/wci, yPlot, /over, thick = 4, trans = 70, name = 'wUH', /insert_legend, lineStyle = 2
	iPlot, [wpe,wpe]/wci, yPlot, /over, thick = 4, trans = 70, name = 'wpe', /insert_legend, color = blue
	iPlot, abs([wce,wce])/wci, yPlot, /over, thick = 4, trans = 70, name = 'wce', /insert_legend, color = green
	iPlot, [wci,wci]/wci, yPlot, /over, thick = 4, trans = 70, name = 'wci', /insert_legend, color = red
	iPlot, [wpi,wpi]/wci, yPlot, /over, thick = 4, trans = 70, name = 'wpi', /insert_legend, color = purple

	wNSTX	= 30.0d6*2*!pi
	iPlot, [wNSTX,wNSTX]/wci, yPlot, /over, thick = 8, trans = 70, name = 'wNSTX', /insert_legend 


stop
end
