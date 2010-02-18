pro k_vs_w, $
	krSq_plot = krSq_plot, $
	freq = freq, $
	wci = wci, $
	kRAll = kRAll, $
	plot_kR = plot_kR, $
	plot_over = plot_over, $
	poloidalScale = poloidalScale, $
	nPhi = nPhi, $
	kz = kz, $
	color_no = color_no

	@load_colors

	plot_vG	= 1
	plot_vGvP	= 1
	plot_thRC	= 1
	plot_thKB	= 1

	if not keyword_set ( plot_kR ) then plot_kR = 1
	if not keyword_set ( plot_over ) then plot_over = 0
	if not keyword_set ( kz ) then kz = 0.0
	if not keyword_set ( nPhi ) then nPhi = 0.0
	if not keyword_set ( poloidalScale ) then poloidalScale = 1.0
	if not keyword_set ( color_no ) then color_no = 0

	if color_no eq 0 then color = black
	if color_no eq 1 then color = red
	if color_no eq 2 then color = green
	if color_no eq 3 then color = purple
	if color_no eq 4 then color = blue

	nf	= 500
	fStart	= 1d6
	fEnd	= 800d6
	df	= ( fEnd - fStart ) / (nf - 1)

	kRAll	= complexArr ( nf, 4 )
	S	= complexArr ( nf )
	P	= complexArr ( nf )

	rSlice	= 1.175;1.11

	;w	= ( fIndGen ( nf ) * df + fStart ) * 2 * !pi

	freq	= 10d0^(fIndGen(nf)/(nf-1)*5+6)

	stixS_prev = 1
	stixR_prev = 1
	stixL_prev = 1
	stixP_prev = 1

	for i=0,nf-1 do begin
	
		rsfwc_1d, $
			/dispersionOnly, $
			freq = freq[i], $
			kR = kR, rOut = r, $
			specData = specData, $
			stixVars = stixVars, $
			runData = runData, $
			kPhi = kPhi, $
			in_kz = kz, $
			poloidalScale = poloidalScale, $
			nPhi = nPhi


		iiSlice	= where ( abs ( r - rSlice ) eq min ( abs ( r - rSlice ) ) )

		if i eq 0 then begin
		
			bR		=	(runData.bField[*,0])[iiSlice[0]]
			bPhi	=	(runData.bField[*,1])[iiSlice[0]]
			bz		=	(runData.bField[*,2])[iiSlice[0]]
			bMag	= 	(runData.bMag)[iiSlice[0]]

		endif

		kPhi	= kPhi[iiSlice[0]]
		kRAll[i,*]	= kR[iiSlice[0],*]
		wpe	= specData[1].wp[iiSlice[0]]
		wce	= specData[1].wc[iiSlice[0]]
		wpi	= specData[0].wp[iiSlice[0]]
		wci	= specData[0].wc[iiSlice[0]]
		wLH	= specData[0].wLH[iiSlice[0]]
		wUH	= specData[0].wUH[iiSlice[0]]

		stixS	= stixVars.stixS[iiSlice[0]]
		stixP	= stixVars.stixP[iiSlice[0]]
		stixR	= stixVars.stixR[iiSlice[0]]
		stixL	= stixVars.stixL[iiSlice[0]]

		S[i]	= stixS
		P[i]	= stixP

		if i ne 0 then begin
			
			if stixS / stixS_prev lt 0 then begin
				print, 'Found S=0 resonance'
				if size(S0,/type) eq 0 then begin
					S0	= freq[i] * 2 * !pi 
				endif else begin
					S0	= [ S0, freq[i] * 2 * !pi ]
				endelse
			endif 

			stixS_prev	= stixS

			if stixR / stixR_prev lt 0 then begin
				print, 'Found R=0 cutoff'
				if size(R0,/type) eq 0 then begin
					R0	= freq[i] * 2 * !pi 
				endif else begin
					R0	= [ R0, freq[i] * 2 * !pi ]
				endelse
			endif 

			stixR_prev	= stixR

			if stixL / stixL_prev lt 0 then begin
				print, 'Found L=0 cutoff'
				if size(L0,/type) eq 0 then begin
					L0	= freq[i] * 2 * !pi 
				endif else begin
					L0	= [ L0, freq[i] * 2 * !pi ]
				endelse
			endif 

			stixL_prev	= stixL

			if stixP / stixP_prev lt 0 then begin
				print, 'Found P=0 cutoff'
				if size(P0,/type) eq 0 then begin
					P0	= freq[i] * 2 * !pi 
				endif else begin
					P0	= [ P0, freq[i] * 2 * !pi ]
				endelse
			endif 

			stixP_prev	= stixP

		endif

		close, /all

	endfor

	kPhi	= fltArr ( size ( kRAll, /dim ) ) + kPhi
	kz		= fltArr ( size ( kRAll, /dim ) ) + kz

	kdotB	= fltarr ( size ( krall, /dim ) )
	kPerp	= fltArr ( size ( kRAll, /dim ) )
	kMag	= sqrt ( real_part(kRAll)^2+kPhi^2+kz^2 )
	kDotB	= ( kRAll * bR + kPhi * bPhi + kz * bz ) / bMag
	kPerp	= sqrt ( kMag - kDotB )

	kDotB_	= kRAll * bR + kPhi * bPhi + kz * bz 
	theta_k_B	= aCos ( kDotB_ / ( kMag * bMag ) )

	dw	= (freq[1]-freq[0])*2*!pi
	vgPerp	= dw / (kPerp[1:*,*] -kPerp[0:n_elements(kPerp[*,0])-2,*])
	vpPerp	= rebin(freq,n_elements(freq),4)*2*!pi/kPerp

	thetaRC	= aTan ( sqrt ( complex ( -P, P*0 ) ), sqrt ( complex ( S, S*0 ) ) )
stop
	;	Resonance cone angle

	if plot_thRC then begin	
		iPlot, freq*2*!pi/wci, thetaRC * !radeg, /xLog, $
			window_title = 'Resonance cone angle'
		iPlot, freq*2*!pi/wci, imaginary(thetaRC) * !radeg, /over, lineStyle = 2
	endif

	;	Angle between k and B

	if plot_thKB then begin
		iPlot, freq*2*!pi/wci, theta_k_B[*,0] * !radeg, /xLog, sym_index = 4, lineStyle = 6, $
			window_title = 'k dot B angle'
		iPlot, freq*2*!pi/wci, theta_k_B[*,1] * !radeg, /over, sym_index = 4, lineStyle = 6
		iPlot, freq*2*!pi/wci, theta_k_B[*,2] * !radeg, /over, sym_index = 4, lineStyle = 6
		iPlot, freq*2*!pi/wci, theta_k_B[*,3] * !radeg, /over, sym_index = 4, lineStyle = 6
	endif

	;	vG / vP	

	if plot_vGvP then begin
		iPlot, freq[1:*]*2*!Pi/wci, vgPerp[*,0]/vpPerp[1:*,0], $
			sym_index = 4, lineStyle = 6, color = blue, $
			/xLog, /zoom_on_resize, /stretch_to_fit, symSize = 0.5, $
			window_title = 'vG / vP'
		iPlot, freq[1:*]*2*!Pi/wci, vgPerp[*,1]/vpPerp[1:*,1], $
			sym_index = 4, lineStyle = 6, /over, color = green, symSize = 0.5
		iPlot, freq[1:*]*2*!Pi/wci, vgPerp[*,2]/vpPerp[1:*,2], $
			sym_index = 4, lineStyle = 6, /over, color = red, symSize = 0.5
		iPlot, freq[1:*]*2*!Pi/wci, vgPerp[*,3]/vpPerp[1:*,3], $
			sym_index = 4, lineStyle = 6, /over, color = purple, symSize = 0.5
	endif

	;	vG - group velocity

	if plot_vG then begin
		iPlot, freq[1:*]*2*!Pi/wci, vgPerp[*,0], sym_index = 4, lineStyle = 6, color = blue, $
			/xLog, /zoom_on_resize, /stretch_to_fit, symSize = 0.5, window_title = 'Perp group velocity'
		iPlot, freq[1:*]*2*!Pi/wci, vgPerp[*,1], $
			sym_index = 4, lineStyle = 6, /over, color = green, symSize = 0.5
		iPlot, freq[1:*]*2*!Pi/wci, vgPerp[*,2], $
			sym_index = 4, lineStyle = 6, /over, color = red, symSize = 0.5
		iPlot, freq[1:*]*2*!Pi/wci, vgPerp[*,3], $
			sym_index = 4, lineStyle = 6, /over, color = purple, symSize = 0.5
	endif


	if plot_kR then begin

		;	kR vs w

		kRPlot	= real_part ( kRAll^2 )
		iiNeg	= where ( kRPlot lt 0 )
		kRPlot[iiNeg]	= -1.0 * kRPlot[iiNeg] 
		kRPlot	= ( kRPlot )^(1d0/5d0)
		kRPlot[iiNeg]	= -1.0 * kRPlot[iiNeg]

		if not plot_over then begin
		iPlot, freq*2*!Pi/wci, kRPlot[*,0], sym_index = 4, lineStyle = 6, color = color, $
			/xLog, /zoom_on_resize, /stretch_to_fit, sym_size = 0.4, id = kRSq_plot
		endif else begin
		iPlot, freq*2*!Pi/wci, kRPlot[*,0], $
			sym_index = 4, lineStyle = 6, over = kRSq_plot, color = color, sym_size = 0.4
		endelse
		iPlot, freq*2*!Pi/wci, kRPlot[*,1], $
			sym_index = 4, lineStyle = 6, over = kRSq_plot, color = color, sym_size = 0.4
		iPlot, freq*2*!Pi/wci, kRPlot[*,2], $
			sym_index = 4, lineStyle = 6, over = kRSq_plot, color = color, sym_size = 0.4
		iPlot, freq*2*!Pi/wci, kRPlot[*,3], $
			sym_index = 4, lineStyle = 6, over = kRSq_plot, color = color, sym_size = 0.4

		;iPlot, freq*2*!Pi/wci, imaginary(kRAll[*,0]), $
		;	sym_index = 3, lineStyle = 6, color = blue, /over
		;iPlot, freq*2*!Pi/wci, imaginary(kRAll[*,1]), $
		;	sym_index = 3, lineStyle = 6, /over, color = green
		;iPlot, freq*2*!Pi/wci, imaginary(kRAll[*,2]), $
		;	sym_index = 3, lineStyle = 6, /over, color = red
		;iPlot, freq*2*!Pi/wci, imaginary(kRAll[*,3]), $
		;	sym_index = 3, lineStyle = 6, /over, color = purple

		;	overplot frequencies of interest

		if not plot_over then begin
	
			yPlot	= [ -max(abs(kRPlot)^2),max(abs(kRPlot)^2) ]
			iPlot, [wLH,wLH]/wci, yPlot, $
				/over, thick = 4, trans = 20, name = 'wLH', /insert_legend 
			iPlot, [wUH,wUH]/wci, yPlot, $
				/over, thick = 4, trans = 20, name = 'wUH', /insert_legend, lineStyle = 2
			iPlot, [wpe,wpe]/wci, yPlot, $
				/over, thick = 4, trans = 20, name = 'wpe', /insert_legend, color = blue
			iPlot, abs([wce,wce])/wci, yPlot, $
				/over, thick = 4, trans = 20, name = 'wce', /insert_legend, color = green
			iPlot, [wci,wci]/wci, yPlot, $
				/over, thick = 4, trans = 20, name = 'wci', /insert_legend, color = red
			iPlot, [wpi,wpi]/wci, yPlot, $
				/over, thick = 4, trans = 20, name = 'wpi', /insert_legend, color = purple
	
			wNSTX	= 30.0d6*2*!pi
			iPlot, [wNSTX,wNSTX]/wci, yPlot, $
				/over, thick = 8, trans = 70, name = 'wNSTX', /insert_legend 
	
			;	overplot the S,P,R and L = 0 surfaces
	
			for i=0,n_elements(S0)-1 do begin
				iPlot, [S0[i],S0[i]]/wci, yPlot, /over, thick = 2, trans = 40, $
					name = 'S=0', /insert_legend, color = blue, lineStyle = 2
			endfor
			for i=0,n_elements(R0)-1 do begin
				iPlot, [R0[i],R0[i]]/wci, yPlot, /over, thick = 2, trans = 40, $
					name = 'R=0', /insert_legend, color = red, lineStyle = 2
			endfor
			for i=0,n_elements(L0)-1 do begin
				iPlot, [L0[i],L0[i]]/wci, yPlot, /over, thick = 2, trans = 40, $
					name = 'L=0', /insert_legend, color = green, lineStyle = 2
			endfor
			for i=0,n_elements(L0)-1 do begin
				iPlot, [P0[i],P0[i]]/wci, yPlot, /over, thick = 2, trans = 40, $
					name = 'P=0', /insert_legend, color = purple, lineStyle = 2
			endfor

			iPlot, [min(freq)*2*!pi/wci,max(freq)*2*!pi/wci], $
				[0,0], /over

		endif
	endif

end
