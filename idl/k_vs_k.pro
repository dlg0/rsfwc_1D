pro k_vs_k

	@load_colors


	rSlice	= 1.175
	freq	= 30d6 

    n_nPhi  = 99 
    s_nPhi  = -10.0
    e_nPhi  = 10.0
    d_nPhi  = ( e_nPhi - s_nPhi ) / n_nPhi 
    nPhi    = fIndGen ( n_nPhi ) * d_nPhi + s_nPhi
	;nPhi	= 10d0^(fIndGen(n_nPhi)/(n_nPhi-1)*10-9)

    n_kz  = 9 
    s_kz  =  -100.0
    e_kz  =  100.0
    d_kz  = ( e_kz - s_kz ) / n_kz
    kz    = fIndGen ( n_kz ) * d_kz + s_kz

    kRAll   = complexArr ( n_nPhi, n_kz, 4 )
    kPar_re   = complexArr ( n_nPhi, n_kz, 4 )
    kPer_re   = complexArr ( n_nPhi, n_kz, 4 )
    kPar_im   = complexArr ( n_nPhi, n_kz, 4 )
    kPer_im   = complexArr ( n_nPhi, n_kz, 4 )
    kPar   = complexArr ( n_nPhi, n_kz, 4 )
    kPer   = complexArr ( n_nPhi, n_kz, 4 )

    kPhiAll = fltArr ( n_nPhi, n_kz, 4 )
    kzAll = fltArr ( n_nPhi, n_kz, 4 )

	for i=0,n_nPhi-1 do begin
        print, i, n_nPhi
        for j=0,n_kz-1 do begin
	
		    rsfwc_1d, $
		    	/dispersionOnly, $
		    	freq = freq, $
		    	kR = kR, rOut = r, $
		    	specData = specData, $
		    	stixVars = stixVars, $
		    	runData = runData, $
		    	kPhi = kPhi, $
		    	in_kz = kz[j], $
		    	poloidalScale = 1d-6, $
		    	nPhi = nPhi[i]

		    iiSlice	= where ( abs ( r - rSlice ) eq min ( abs ( r - rSlice ) ) )

		    if i eq 0 then begin
		    
		    	bR		=	(runData.bField[*,0])[iiSlice[0]]
		    	bPhi	=	(runData.bField[*,1])[iiSlice[0]]
		    	bz		=	(runData.bField[*,2])[iiSlice[0]]
		    	bMag	= 	(runData.bMag)[iiSlice[0]]

		    endif

		    kPhiAll[i,j,*]	= nPhi[i]/rSlice;kPhi[iiSlice[0]]
		    kRAll[i,j,*]	= kR[iiSlice[0],*]
		    kzAll[i,j,*]	= kz[j]

            kMag	= sqrt ( kRAll[i,j,*]^2+kPhiAll[i,j,*]^2+kz[j]^2 )
	        kPar[i,j,*]	= ( kRAll[i,j,*] * bR + kPhiAll[i,j,*] * bPhi + kz[j] * bz ) / bMag
	        kPer[i,j,*]	= sqrt ( kMag - kPar[i,j,*] )
	        theta_k_B	= aCos ( kPar[i,j,*] / ( kMag * bMag ) )

            close, /all

        endfor
	endfor

    !p.multi = [0,2,2]
    plot, 1/kpar[*,*,0], 1/kper[*,*,0], psym =4;, xRange = [-1,1], yRange = [-0.5,0.5]

    plot, 1/kpar[*,*,1], 1/kper[*,*,1], psym =4;, xRange = [-1,1], yRange = [-0.5,0.5]

    plot, 1/kpar[*,*,2], 1/kper[*,*,2], psym =4;, xRange = [-1,1], yRange = [-0.5,0.5]

    plot, 1/kpar[*,*,3], 1/kper[*,*,3], psym =4;, xRange = [-1,1], yRange = [-0.5,0.5]

stop
end
