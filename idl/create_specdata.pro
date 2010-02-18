pro create_specdata, $
    ionSpecZ, $
    ionSpecAmu, $
    nMax, $
    bMag, $
    specData = specData

    @constants

    nIonSpec	= n_elements ( ionSpecZ )

	specData	= replicate ( $
			{ 	q : 0d0, $
				m : 0d0, $
				wp : 0.0, $
				wc : 0.0, $
				n : 0.0, $
				wLH : 0.0, $
				wUH : 0.0 }, nIonSpec + 1 )

	for i=0,nIonSpec do begin

		; ions
		if i lt nIonSpec then begin

			specData[i].q 		= ionSpecZ[i] * e
			specData[i].m 		= ionSpecAmu[i] * mi
			specData[i].n 		= nMax[i] 

		; electrons for charge neutrality
		endif else begin

			specData[i].q 		= -e
			specData[i].m 		= me
			if nIonSpec gt 1 then begin
				specData[i].n 		= total ( specData[0:nIonSpec-1].n $
										* ionSpecZ ) 
			endif else begin
				specData[i].n 		= specData[0:nIonSpec-1].n $
										* ionSpecZ 
			endelse

		endelse

	endFor


    for i = 0, nIonSpec do begin

		specData[i].wp	= sqrt ( specData[i].n * specData[i].q^2 $
                            / ( specData[i].m*e0 ))
		specData[i].wc	= specData[i].q * bMag / specData[i].m
	
	endfor

	for i = 0, nIonSpec - 1 do begin

		specData[i].wLH	= sqrt ( specData[i].wc^2 $
			+ specData[i].wp^2 / ( 1.0 + specData[nIonSpec].wp^2 / specData[nIonSpec].wc^2 ) )
		specData[i].wUH	= sqrt ( specData[nIonSpec].wp^2 + specData[nIonSpec].wc^2 ) 
	
	endfor


end
