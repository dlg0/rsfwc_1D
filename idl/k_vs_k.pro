pro k_vs_k

	@load_colors
    @constants
    @switches

    dispersion_generalised = 1
    dispersion_noPoloidal = 0

    ionSpecZ    = [ 1.0 ]
    ionSpecAmu  = [ 2.0 ]

    r   	= 1.00
    z       = 0.00
	freq	= 30d6 
    wReal   = 2.0 * !pi * freq

    nstx_eqdsk	= '../eqdsk/g120740.00275.EFIT02.mds.uncorrected.qscale_1.00000.dlgMod_1.67'
	eqdsk	= readgeqdsk ( nstx_eqdsk )
	
	bR  = interpolate ( eqdsk.bR, $
			( r - eqdsk.rleft ) / eqdsk.rdim * (eqdsk.nW-1.0), $
        ( z - min ( eqdsk.z ) ) / eqdsk.zdim * (eqdsk.nH-1.0), $
		cubic = -0.5 )
    bPhi  = interpolate ( eqdsk.bPhi, $
			( r - eqdsk.rleft ) / eqdsk.rdim * (eqdsk.nW-1.0), $
        ( z - min ( eqdsk.z ) ) / eqdsk.zdim * (eqdsk.nH-1.0), $
		cubic = -0.5 )
    bz  = interpolate ( eqdsk.bz, $
			( r - eqdsk.rleft ) / eqdsk.rdim * (eqdsk.nW-1.0), $
        ( z - min ( eqdsk.z ) ) / eqdsk.zdim * (eqdsk.nH-1.0), $
		cubic = -0.5 )

    r0  = 1.0
    b0  = 0.53
    bPhi    = b0 / r * r0
    bR  = 0.5 * bPhi
    bz  = 0.5 * bPhi

    bMag    = sqrt ( bR^2 + bPhi^2 + bz^2 )

    nstx_profile	= '../profiles/dlg_profiles_x1.00.nc' 
	cdfId = ncdf_open ( nstx_profile, /noWrite ) 
    ncdf_varGet, cdfId, 'ne', xMap_ne
    ncdf_varGet, cdfId, 'r2d', xMap_R2D
    ncdf_varGet, cdfId, 'z2d', xMap_z2D
    nCdf_close, cdfId
         
    nX  = n_elements ( xMap_R2D[*,0] )
    nY  = n_elements ( xMap_R2D[0,*] )

	ne_  = interpolate ( xMap_ne, $
			( r - min(xMap_R2d) ) / (max(xMap_R2d)-min(xMap_R2d)) * (nX-1.0), $
        ( z - min ( xMap_z2d ) ) / (max(xMap_z2d)-min(xMap_z2d)) * (nY-1.0), $
		cubic = -0.5 )

    nMax    = [ ne_ / ionSpecZ[0] ]
    nMax    = 2.0d18

    create_specData, ionSpecZ, ionSpecAmu, nMax, bMag, $
        specData = specData 

    stixVariables, wReal, specData, $
	    stixVars = stixVars, /noHalf

	runData = { bField : [ [bR], [bPhi], [bz] ], $
				bMag : bMag, $
				nR : 1, $
				r : r, $
                kz : 0.0, $
				nPhi : 0.0 }

    dielectric, runData, stixVars, $
	    epsilonFull = epsilon, /noHalfGrid

    n_nPhi  = 2.0 
    s_nPhi  = 0.0
    e_nPhi  = 0.0
    d_nPhi  = ( e_nPhi - s_nPhi ) / n_nPhi 
    nPhi    = fIndGen ( n_nPhi ) * d_nPhi + s_nPhi
	;nPhi	= 10d0^(fIndGen(n_nPhi)/(n_nPhi-1)*10-9)

    n_kz  = 2 
    s_kz  =  0.0
    e_kz  =  0.0
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
    thetaC   = complexArr ( n_nPhi, n_kz, 4 )

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

            kMag_im[i,j,*]	= sqrt( imaginary(kRAll[i,j,*])^2+kPhiAll[i,j,*]^2+kz[j]^2 )
	        kDotB	= ( imaginary(kRAll[i,j,*]) * bR + kPhiAll[i,j,*] * bPhi + kz[j] * bz )
	        theta_im[i,j,*]	= aCos ( kDotB / ( kMag_im[i,j,*] * bMag ) )

            for mm=0,3 do begin

                kMag[i,j,mm]    = sqrt ( $
                        conj ( kRAll[i,j,mm] ) * kRAll[i,j,mm] $
                        + kPhiAll[i,j,mm]^2 $
                        + kz[j]^2 )

                innerProdC_k_B  = conj ( kRAll[i,j,mm] ) * bR $
                        + kPhiAll[i,j,mm] * bPhi $
                        + kz[j] * bz

                thetaC[i,j,mm]	= aCos ( innerProdC_k_B / ( kMag[i,j,mm] * bMag ) )

                AA  = fltArr ( 3*2 )
                BB  = fltArr ( 3*2 )

                AA[0]   = real_part ( kRAll[i,j,mm] )
                AA[1]   = imaginary ( kRAll[i,j,mm] )
                AA[2]   = kPhiAll[i,j,mm]
                AA[3]   = 0.0
                AA[4]   = kz[j]
                AA[5]   = 0.0

                BB[0]   = bR
                BB[1]   = 0.0
                BB[2]   = bPhi
                BB[3]   = 0.0
                BB[4]   = bz
                BB[5]   = 0.0

                innerProdR_k_B  = transpose ( AA ) # BB
                innerProdR_k_B  = total ( AA * BB )

                AAMag   = sqrt ( total ( AA^2 ) )
                BBMag   = sqrt ( total ( BB^2 ) )
                BBMag   = bMag

                theta[i,j,mm]   = aCos ( innerProdR_k_B / ( AAMag * BBMag ) ) 
 
            endfor

            close, /all

        endfor
	endfor
   
    restore, 'u1u2.sav'
    magRange  = max ( abs ( [u1[iiNSTX,*], u2[iiNSTX,*] ] ) )
    xRange  = [-magRange, magRange]*1.2
    yRange  = xRange
    device, decomposed = 0
    loadct, 12
    !p.background = 255
    window, 0, xSize = 800, ySize = 800
    plot, real_part(u1[iiNSTX,*]), thetaArr*!dtor, $
            /polar, /iso, xRange = xRange, yRange = yRange, $
            xStyl = 1, yStyle=1, color = 0, $
            xTickLen = 1, yTickLen = 1, xTickInterval = 0.02, yTickInterval=0.02
    oPlot, imaginary(u1[iiNSTX,*]), thetaArr*!dtor, /polar, lineStyle = 2, color = 0
    oPlot, real_part(u2[iiNSTX,*]), thetaArr*!dtor, /polar, color = 0
    oPlot, imaginary(u2[iiNSTX,*]), thetaArr*!dtor, /polar, lineStyle = 2, color = 0
   
    oplot, real_part(1/kMag), theta, /polar, psym = 4, symSize = 1, color = 12*16-1, thick = 3
    oplot, imaginary(1/kMag), imaginary(theta), /polar, psym = 6, symSize = 1, color = 12*16-1, thick = 3
    oPlot, 1/kMag_re, real_part(theta_re), /polar, psym = 5, color = 8*16-1, thick = 2
    oPlot, 1/kMag_im, real_part(theta_im), /polar, psym = 6, color = 3*16-1, thick = 2

    ;   Calculate the nPhi/kz enforcement surface

    kRScan  = (fIndGen(10001)/10000-0.5)*2000
    kMagScan    = sqrt(kRScan^2+(s_nPhi/r)^2+s_kz^2)
    kDotB   = kRScan * bR + (s_nPhi/r) * bPhi + s_kz * bz
    thetaScan   = aCos ( kDotB / ( kMagScan * bMag ) )
    kPerScan    = sqrt(kMagScan^2-(kDotB/bMag)^2)

    oPlot, 1/kMagScan, thetaScan, /polar, color = 0, lineStyle = 1
    oPlot, 1/kMagScan, -thetaScan, /polar, color = 0, lineStyle = 1

    ;window, 1, xSize = 600, ySize = 600
    ;plot, 1/imaginary(kMag[*,*,0]), (theta2[*,*,0]), psym =4, /polar, /iso, xRange = xRange, yRange = yRange
    ;oplot, 1/imaginary(kMag[*,*,1]),(theta2[*,*,1]), psym =4, /polar
    ;oplot, 1/imaginary(kMag[*,*,2]), (theta2[*,*,2]), psym =4, /polar
    ;oplot, 1/imaginary(kMag[*,*,3]), (theta2[*,*,3]), psym =4, /polar
 
    ;;window, 2, xSize = 600, ySize = 600

    ;iPlot, real_part(u1[iiNSTX,*]), thetaArr*!dtor, /polar, color = red, /iso
    ;iPlot, imaginary(u1[iiNSTX,*]), thetaArr*!dtor, /polar, /over, color = red, lineStyle = 2
    ;iPlot, real_part(u2[iiNSTX,*]), thetaArr*!dtor, /polar, /over, color = blue
    ;iPlot, imaginary(u2[iiNSTX,*]), thetaArr*!dtor, /polar, /over, color = blue, lineStyle = 2

    ;iplot, 1/kMag_re[*,*,0], theta_re[*,*,0],sym_index =4, /polar, lineStyle = 6, /over, color = blue
    ;iplot, 1/kMag_re[*,*,1], theta_re[*,*,1],sym_index =4, /polar, /over, lineStyle = 6, color = red
    ;iplot, 1/kMag_re[*,*,2], theta_re[*,*,2],sym_index =4, /polar, /over, lineStyle = 6, color = green
    ;iplot, 1/kMag_re[*,*,3], theta_re[*,*,3],sym_index =4, /polar, /over, lineStyle = 6, color = purple

    ;iplot, 1/kMag_im[*,*,0], theta_im[*,*,0],sym_index =5, /polar, /over, lineStyle = 6, color = blue
    ;iplot, 1/kMag_im[*,*,1], theta_im[*,*,1],sym_index =5, /polar, /over, lineStyle = 6, color = red
    ;iplot, 1/kMag_im[*,*,2], theta_im[*,*,2],sym_index =5, /polar, /over, lineStyle = 6, color = green
    ;iplot, 1/kMag_im[*,*,3], theta_im[*,*,3],sym_index =5, /polar, /over, lineStyle = 6, color = purple

stop
end
