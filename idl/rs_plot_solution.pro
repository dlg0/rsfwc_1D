pro rs_plot_solution, antLoc, dR, nR, $
	e1, e2, e3, $
	kR = kR, r_kR = r_kR, $
	r1 = r1, r2 = r2, r3 = r3, $
    kjIn=kj

	common dlg_colors
	common plotSwitches
    common switches
   
    ;	Visualise solution

		eRange	= max(sqrt ( abs(e1)^2 + abs(e2)^2 + abs(e3)^2 ))
   		p_r = plot ( r1, e1, layout=[1,3,1],$
				title='Er',ytitle='Er [V/m]',name='Re',window_title='rsfwc_1d',$
                yRange = [-eRange,eRange])
		p_i = plot ( r1, imaginary(e1), color='red',/over,name='Im')

        if ar2EField then begin
            p = plot(r1,kj.eR*kj.replace,/over,color='blue',linestyle='--')
            p = plot(r1,imaginary(kj.eR)*kj.replace,/over,color='red',linestyle='--')
            ;p = plot(r1,e1-(kj.er*kj.replace),/over,color='g',thick=2)
            ;p = plot(r1,imaginary(e1-(kj.er*kj.replace)),/over,color='g',thick=1)
        endif

		l = legend(target=[p_r,p_i],position=[0.98,0.95],/norm,font_size=10,horizontal_alignment='RIGHT')

		p_r = plot ( r2, e2, layout=[1,3,2],/current,$
				title='Et',ytitle='Et [V/m]',name='Re')
		p_i = plot ( r2, imaginary(e2), color='red',/over,name='Im')
        if ar2EField then begin
            p = plot(r1,kj.et*kj.replace,/over,color='blue',linestyle='--')
            p = plot(r1,imaginary(kj.et)*kj.replace,/over,color='red',linestyle='--')
            ;p = plot(r1,e2-(kj.et*kj.replace),/over,color='g',thick=2)
            ;p = plot(r1,imaginary(e2-(kj.et*kj.replace)),/over,color='g',thick=1)
        endif


		l = legend(target=[p_r,p_i],position=[0.98,0.62],/norm,font_size=10,horizontal_alignment='RIGHT')

		p_r = plot ( r3, e3, layout=[1,3,3],/current,$
				title='Ez',ytitle='Ez [V/m]',name='Re')
		p_i = plot ( r3, imaginary(e3), color='red',/over,name='Im')
        if ar2EField then begin
            p = plot(r1,kj.ez*kj.replace,/over,color='blue',linestyle='--')
            p = plot(r1,imaginary(kj.ez)*kj.replace,/over,color='red',linestyle='--')
            ;p = plot(r1,e3-(kj.ez*kj.replace),/over,color='g',thick=2)
            ;p = plot(r1,imaginary(e3-(kj.ez*kj.replace)),/over,color='g',thick=1)
        endif


		l = legend(target=[p_r,p_i],position=[0.98,0.28],/norm,font_size=10,horizontal_alignment='RIGHT')

        ;loadct, 12, /sil
    
        ;iPlot, r1, e1, $
        ;    title = 'e1', $
		;	id = solution, $
		;	view_grid = [1,3], $
		;	/zoom_on_resize, $
		;	/stretch_to_fit
        ;iPlot, r1, imaginary ( e1 ), $
		;	view_number = 1, $
        ;    over = solution, $
        ;    color = transpose ( ct12[8*16-1,*] )
        ;;iPlot, r1, abs ( e1 ), $
		;;	view_number = 1, $
        ;;    over = solution, $
        ;;    color = red
	
		;iPlot, [antLoc,antLoc], [-max(abs(e1)),max(abs(e1))], $
		;	transparency = 80, $
		;	thick = 10, $
		;	over = solution, $
		;	yRange = [ -max ( eMag ), max ( eMag ) ]
	
    	;iPlot, r2, e2, $
        ;    psym = -4, $
        ;    view_number = 2, $
        ;    title = 'e2', $
		;	over = solution, $
		;	/stretch_to_fit

    	;iPlot, r2, imaginary(e2), $
        ;    psym = -4, $
        ;    over = solution, $
		;	view_number = 2, $
        ;    color = transpose ( ct12[8*16-1,*] )
    	;;iPlot, r2, abs(e2), $
        ;;    psym = -4, $
        ;;    over = solution, $
		;;	view_number = 2, $
        ;;    color = red 
		;iPlot, [antLoc,antLoc], [-max(abs(e2)),max(abs(e2))], $
		;	transparency = 80, $
		;	thick = 10, $
		;	over = solution, $
		;	yRange = [ -max ( eMag ), max ( eMag ) ] 

    	;iPlot, r3, e3, $
        ;    psym = -4, $
        ;    view_number = 3, $
        ;    title = 'e3', $
		;	over = solution, $
		;	/stretch_to_fit

    	;iPlot, r3, imaginary(e3), $
        ;    psym = -4, $
        ;    over = solution, $
        ;    color = transpose ( ct12[8*16-1,*] )
    	;;iPlot, r3, abs(e3), $
        ;;    psym = -4, $
        ;;    over = solution, $
        ;;    color = red 
		;iPlot, [antLoc,antLoc], [-max(abs(e3)),max(abs(e3))], $
		;	transparency = 80, $
		;	thick = 10, $
		;	over = solution, $
		;	yRange = [ -max ( eMag ), max ( eMag ) ] 

		print, 'FFTing solution ...'
		nFFT	= 32 
		fftData_e1	= complexArr ( nFFT )
		fftData_e2	= complexArr ( nFFT )
		fftData_e3	= complexArr ( nFFT )
		rFFTData	= fltArr ( nFFT )
		fftStep	= nR / nFFT 
		fftWidth	= nR / 10 
		fftWidth	+= (fftWidth mod 2)
		osf	= 10

		fftData2d	= complexArr ( nFFT, fftWidth*osf )
		hannWin	= hanning ( fftWidth )
		kRFFT	= findgen ( fftWidth * osf ) / ( fftWidth * osf * dr ) * ( 2 * !pi )
		eTmp	= complexArr ( fftWidth * osf )

		for i=0L,nFFT-1L do begin
		
			iiLow	= i*fftStep-fftWidth/2
			iiHig	= i*fftStep+fftWidth/2-1

			rFFTData[i]	= r1[i*fftStep]

			if iiLow gt 0 and iiHig lt nR-2 then begin

				e3Tmp	= ( e3[iiLow:iiHig] + e3[iiLow+1:iiHig+1] )	/ 2
				e2Tmp	= ( e2[iiLow:iiHig] + e2[iiLow+1:iiHig+1] )	/ 2
				e1Tmp	= e1[iiLow:iiHig] 
				eTmp[fftWidth*osf/2-fftWidth/2:fftWidth*osf/2+fftWidth/2-1]	= $
					e1Tmp * hannWin
				fftTmp	= abs ( fft ( eTmp ) )^2
				fftData2D[i,*]	= fftTmp

			endif

		endfor
		print, 'DONE'
	
		nLevs	= 21	
		levels	= 10.0^fIndGen ( nLevs )/1e20
		colors	= 255 - ( bytScl ( fIndGen(nLevs), top = 253 ) + 1 )

		fftData2D_pm	= fftData2D
		nqII	= n_elements(fftData2D[0,*])/2
		fftData2D_pm[*,nqII:*]	= fftData2D[*,0:nqII-1]
		fftData2D_pm[*,0:nqII-1]	= fftData2D[*,nqII:*]
		kRFFT_pm	= kRFFT - kRFFT[nqII]	

		kMax	= 2000 
		iiPlotk	= where ( abs ( kRFFT_pm ) le kMax )
		;iContour, fftData2D_pm[*,iiPlotk], $
		;	rFFTData, kRFFT_pm[iiPlotk], $
		;	yRange = [-kMax,kMax], $
		;	xRange = [min(r1),max(r1)], $
		;	rgb_indices = colors, $
		;	rgb_table = 1, $
		;	/fill, $
		;	c_value = levels, $
		;	/zoom_on_resize, $
		;	id = spectrum, $
		;	/stretch_to_fit

		if plotFFTSolution then $
		c=contour(fftData2D_pm,rfftdata, krfft_pm, yrange=[-kmax,kmax], rgb_table=1,n_levels=20)

		;if keyword_set ( kR ) then begin

		;	p=plot(r_kR, real_part ( kR[*,0] ), /over)
		;	p=plot(r_kR, real_part ( kR[*,1] ), /over)
		;	p=plot(r_kR, real_part ( kR[*,2] ), /over)
		;	p=plot(r_kR, real_part ( kR[*,3] ), /over)
		;	
		;endif

		if keyword_set ( divD ) then begin

    		;   plot the nume1ical and real charge density
 
        	iPlot, r_, rho_, $
        	    color = transpose ( ct12[8*16-1,*] ), $
        	    /over, $
        	    title = 'rho', $
        	    /view_next
        	iPlot, r_, divD_, $
        	    /over
        	iPlot, r_, imaginary ( divD_ ), $
        	    /over, $
        	    lineStyle = 1

        	iPlot, r_, divD_, $
        	    title = 'div D', $
        	    /over, $
        	    /view_next
        	iPlot, r_, imaginary ( divD_ ), $
        	    /over, $
        	    lineStyle = 1

		endif
end
