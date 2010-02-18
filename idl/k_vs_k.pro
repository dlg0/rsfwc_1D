pro k_vs_k

	@load_colors
    @constants
    @switches

    dispersion_generalised = 1
    dispersion_noPoloidal = 0

    ionSpecZ    = [ 2 ]
    ionSpecAmu  = [ 4 ]

    nMax     = [ 5.7e18 ]

    bR  = 0.0
    bPhi    = 0.53
    bz  = 0.0
    bMag    = sqrt ( bR^2 + bPhi^2 + bz^2 )

    create_specData, ionSpecZ, ionSpecAmu, nMax, bMag, $
        specData = specData 

	r   	= 1.175
	freq	= 30d6 
    wReal   = 2.0 * !pi * freq

    stixVariables, wReal, specData, $
	    stixVars = stixVars, /noHalf

	runData = { bField : [ [bR], [bPhi], [bz] ], $
				bMag : bMag, $
				nR : 1, $
				r : r, $
                kz : 0.0, $
				nPhi : 0 }

    dielectric, runData, stixVars, $
	    epsilonFull = epsilon

    n_nPhi  = 401 
    s_nPhi  = -100.0
    e_nPhi  = 100.0
    d_nPhi  = ( e_nPhi - s_nPhi ) / n_nPhi 
    nPhi    = fIndGen ( n_nPhi ) * d_nPhi + s_nPhi
	;nPhi	= 10d0^(fIndGen(n_nPhi)/(n_nPhi-1)*10-9)

    n_kz  = 1 
    s_kz  =  00.0
    e_kz  =  00.0
    d_kz  = ( e_kz - s_kz ) / n_kz
    kz    = fIndGen ( n_kz ) * d_kz + s_kz

    kRAll   = complexArr ( n_nPhi, n_kz, 4 )
    kPar_re   = complexArr ( n_nPhi, n_kz, 4 )
    kPer_re   = complexArr ( n_nPhi, n_kz, 4 )
    kPar_im   = complexArr ( n_nPhi, n_kz, 4 )
    kPer_im   = complexArr ( n_nPhi, n_kz, 4 )
    kMag_re   = complexArr ( n_nPhi, n_kz, 4 )
    theta_re   = complexArr ( n_nPhi, n_kz, 4 )
    kMag_im   = complexArr ( n_nPhi, n_kz, 4 )
    theta_im   = complexArr ( n_nPhi, n_kz, 4 )

    kPar   = complexArr ( n_nPhi, n_kz, 4 )
    kPer   = complexArr ( n_nPhi, n_kz, 4 )
    kMag   = complexArr ( n_nPhi, n_kz, 4 )
    theta   = complexArr ( n_nPhi, n_kz, 4 )
    theta2   = complexArr ( n_nPhi, n_kz, 4 )

    kPhiAll = fltArr ( n_nPhi, n_kz, 4 )
    kzAll = fltArr ( n_nPhi, n_kz, 4 )

	for i=0,n_nPhi-1 do begin
        ;print, i, n_nPhi
        for j=0,n_kz-1 do begin

            runData.kz  = kz[j]    
            runData.nPhi    = nPhi[i]
            dispersion, wReal, epsilon, stixVars, runData, specData, $
	            kR = kR__, kPhi = nPhi[i] / r 

		    kPhiAll[i,j,*]	= nPhi[i]/r
		    kRAll[i,j,*]	= kR__
		    kzAll[i,j,*]	= kz[j]

            kMag_re[i,j,*]	= sqrt ( real_part(kRAll[i,j,*])^2+kPhiAll[i,j,*]^2+kz[j]^2 )
	        kDotB	= ( real_part(kRAll[i,j,*]) * bR + kPhiAll[i,j,*] * bPhi + kz[j] * bz )
	        theta_re[i,j,*]	= aCos ( kDotB / ( kMag_re[i,j,*] * bMag ) )

            kMag_im[i,j,*]	= sqrt ( imaginary(kRAll[i,j,*])^2+kPhiAll[i,j,*]^2+kz[j]^2 )
	        kDotB	= ( imaginary(kRAll[i,j,*]) * bR + kPhiAll[i,j,*] * bPhi + kz[j] * bz )
	        theta_im[i,j,*]	= aCos ( kDotB / ( kMag_im[i,j,*] * bMag ) )

            kMag[i,j,*]	= sqrt ( kRAll[i,j,*]^2+kPhiAll[i,j,*]^2+kz[j]^2 )
            for mm=0,3 do begin
	            kDotB	= transpose ( [ kRAll[i,j,mm], kPhiAll[i,j,mm], kz[j] ] ) # [ bR, bPhi, bz ]
	            theta[i,j,mm]	= aCos ( kDotB / ( (kMag[i,j,mm]) * bMag ) )
	            theta2[i,j,mm]	= aCos ( kDotB / ( imaginary(kMag[i,j,mm]) * bMag ) )
            endfor

            close, /all

        endfor
	endfor

    xRange = [-1.5,1.5];[-0.04,0.04];
    yRange = [-1.5,1.5];[-0.04,0.04];
    xRange = [-0.04,0.04];
    yRange = [-0.04,0.04];


    window, 0, xSize = 600, ySize = 600
    plot, 1/real_part(kMag[*,*,0]), theta[*,*,0], psym =4, /polar, /iso, xRange = xRange, yRange = yRange
    oplot, 1/real_part(kMag[*,*,1]), theta[*,*,1], psym =4, /polar
    oplot, 1/real_part(kMag[*,*,2]), theta[*,*,2], psym =4, /polar
    oplot, 1/real_part(kMag[*,*,3]), theta[*,*,3], psym =4, /polar
    window, 1, xSize = 600, ySize = 600
    plot, 1/imaginary(kMag[*,*,0]), (theta2[*,*,0]), psym =4, /polar, /iso, xRange = xRange, yRange = yRange
    oplot, 1/imaginary(kMag[*,*,1]),(theta2[*,*,1]), psym =4, /polar
    oplot, 1/imaginary(kMag[*,*,2]), (theta2[*,*,2]), psym =4, /polar
    oplot, 1/imaginary(kMag[*,*,3]), (theta2[*,*,3]), psym =4, /polar
 
    window, 2, xSize = 600, ySize = 600
    plot, 1/kMag_im[*,*,0], theta_im[*,*,0], psym =4, /polar, /iso, xRange = xRange, yRange = yRange
    oplot, 1/kMag_im[*,*,1], theta_im[*,*,1], psym =4, /polar
    oplot, 1/kMag_im[*,*,2], theta_im[*,*,2], psym =4, /polar
    oplot, 1/kMag_im[*,*,3], theta_im[*,*,3], psym =4, /polar

    oplot, 1/kMag_re[*,*,0], theta_re[*,*,0], psym =4, /polar;, xRange = xRange, yRange = yRange
    oplot, 1/kMag_re[*,*,1], theta_re[*,*,1], psym =4, /polar
    oplot, 1/kMag_re[*,*,2], theta_re[*,*,2], psym =4, /polar
    oplot, 1/kMag_re[*,*,3], theta_re[*,*,3], psym =4, /polar




stop
end
