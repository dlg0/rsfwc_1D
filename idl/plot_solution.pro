pro plot_solution, r, r_, antLoc, dR, nR, $
	eR, ePhi, ez

	common dlg_colors
   
    ;	Visualise solution
    
        loadct, 12, /sil
		eMag	= sqrt ( abs(er)^2 + abs(ePhi)^2 + abs(ez)^2 )
    
        iPlot, r, eR, $
            view_number = 4, $
            title = 'eR', $
			over = panel
        iPlot, r, imaginary ( eR ), $
			view_number = 4, $
            over = panel, $
            color = transpose ( ct12[8*16-1,*] )
		iPlot, [antLoc,antLoc], [-max(abs(eR)),max(abs(eR))], $
			transparency = 80, $
			thick = 10, $
			over = panel, $
			yRange = [ -mean ( eMag ), mean ( eMag ) ] * 2
	
    	iPlot, [r_[0]-dr,r_,r_[nR-2]+dr], [0,ePhi,0], $
            psym = -4, $
            view_number = 5, $
            title = 'ePhi', $
			over = panel
    	iPlot, r_, imaginary(ePhi), $
            psym = -4, $
            over = panel, $
			view_number = 5, $
            color = transpose ( ct12[8*16-1,*] )
		iPlot, [antLoc,antLoc], [-max(abs(ePhi)),max(abs(ePhi))], $
			transparency = 80, $
			thick = 10, $
			over = panel, $
			yRange = [ -mean ( eMag ), mean ( eMag ) ] * 2

    	iPlot, [r_[0]-dr,r_,r_[nR-2]+dr], [0,ez,0], $
            psym = -4, $
            view_number = 6, $
            title = 'ez', $
			over = panel
    	iPlot, [r_[0]-dr,r_,r_[nR-2]+dr], [0,imaginary(ez),0], $
            psym = -4, $
            over = panel, $
            color = transpose ( ct12[8*16-1,*] )
		iPlot, [antLoc,antLoc], [-max(abs(ez)),max(abs(ez))], $
			transparency = 80, $
			thick = 10, $
			over = panel, $
			yRange = [ -mean ( eMag ), mean ( eMag ) ] * 2

		print, 'FFTing solution ...'
		nFFT	= 128 
		fftData_eR	= complexArr ( nFFT )
		fftData_ePhi	= complexArr ( nFFT )
		fftData_ez	= complexArr ( nFFT )
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

			rFFTData[i]	= r[i*fftStep]

			if iiLow gt 0 and iiHig lt nR-2 then begin

				ezTmp	= ( ez[iiLow:iiHig] + ez[iiLow+1:iiHig+1] )	/ 2
				ePhiTmp	= ( ePhi[iiLow:iiHig] + ePhi[iiLow+1:iiHig+1] )	/ 2
				eRTmp	= eR[iiLow:iiHig] 
				eTmp[fftWidth*osf/2-fftWidth/2:fftWidth*osf/2+fftWidth/2-1]	= $
					eRTmp * hannWin
				fftTmp	= abs ( fft ( eTmp ) )^2
				fftData2D[i,*]	= fftTmp

			endif

		endfor
		print, 'DONE'
	
		nLevs	= 21	
		levels	= 10.0^fIndGen ( nLevs )/1e10
		colors	= 255 - ( bytScl ( fIndGen(nLevs), top = 253 ) + 1 )
		kMax	= 1000 
		iiPlotk	= where ( kRFFT gt 0 and kRFFT le kMax )
		iContour, fftData2D[*,iiPlotk]/max(fftData2D), rFFTData, kRFFT[iiPlotk], $
			yRange = [1,kMax], $
			view_number = 3, $
			xRange = [min(r),max(r)], $
			rgb_indices = colors, $
			rgb_table = 1, $
			/fill, $
			c_value = levels, $
			/zoom_on_resize, $
			over = panel


		iPlot, [antLoc,antLoc], [1,1e3], $
			transparency = 80, $
			thick = 10, $
			over = panel
	

		if keyword_set ( divD ) then begin

    		;   plot the numerical and real charge density
 
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
