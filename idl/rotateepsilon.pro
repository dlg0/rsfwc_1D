function rotateEpsilon, epsilonIn, bUnit_car

;	;	y towards z
;	rotTh1	=	-aTan ( bUnit_car[1], bUnit_car[2] )
;	rot_x	= [ [ 1, 0, 0], $
;				[ 0, cos ( rotTh1 ), -sin ( rotTh1 ) ], $
;				[ 0, sin ( rotTh1 ), cos ( rotTh1 ) ] ]
;	inv_rot_x	= transpose ( rot_x )
;
;	;	z towards x
;	rotTh2	=	0.0;-aTan ( bUnit_car[0], bUnit_car[2] )
;	rot_y	= [ [ cos ( rotTh2 ), 0, -sin ( rotTh2 ) ], $
;				[ 0, 1, 0 ], $
;				[ sin ( rotTh2 ), 0, cos ( rotTh2 ) ] ]
;	inv_rot_y	= transpose ( rot_y )
;
;	;	x towards y
;	rotTh3	=	-aTan ( bUnit_car[0], bUnit_car[1] ) 
;	rot_z	= [ [ cos ( rotTh3 ), -sin ( rotTh3 ), 0 ], $
;				[ sin ( rotTh3 ), cos ( rotTh3 ), 0 ], $
;				[ 0, 0, 1 ] ]
;	inv_rot_z	= transpose ( rot_z )
;
;	print, rotth1*!radeg, rotth2*!radeg, rotth3*!radeg
;
;	;	test rotations
;	;	There seems to be the requirement to rotate z towards x
;	;	by -!pi/2 to get everything to work correctly. I think
;	;	this is due to the ambiguity in the orientation of the x/y
;	;	axis in the stix frame. So I'm still not 100% on this, but
;	;	without the extra rotation the dispersion calculation for
;	;	any case with a z-component of the poloidal field does not 
;	;	agree (nor show any +/- asymmetry in kR) with the dispersion 
;	;	calculation for kPerp based on the PRL of Jaeger, V90 N19, May 2003. 
;	;	Update: It seems that perhaps the dispersion calculation of 
;	;	Jaeger is not particularly useful for seperating the components
;	;	of the poloidal field and that the additional rotation that
;	;	I was unsure of may not be needed at all. The wave solution
;	;	seems to show not much change at all after including a vertical 
;	;	poloidal field, perhaps because the fast wave is perpendicular?
;	;	Anyway, the results seem to agree with AORSA1D within a minus 
;	;	sign.
;

;   Try a different rotation approach
;   ---------------------------------

    ;   get vector perp to both z axis and b

    zaxis   = [ [0.0], [0.0], [1.0] ]
    perp    = [ [zaxis[1]*bUnit_car[2]-zaxis[2]*bUnit_car[1]], $
                [-(zaxis[0]*bUnit_car[2]-zaxis[2]*bUnit_car[0])], $
                [zaxis[0]*bUnit_car[1]-zaxis[1]*bUnit_car[0]] ] 

    ;   check angle between perp and b

    perpDotB    = perp[0]*bUnit_car[0]+perp[1]*bUnit_car[1]+perp[2]*bUnit_car[2]
    perpMag = sqrt ( total ( perp^2 ) )
    bUnitMag    = sqrt ( total ( bUnit_car^2 ) )
    checkTh = aCos ( perpDotB / ( perpMag * bUnitMag ) ) * !radeg

    perpUnit    = perp / perpMag

    ;   get angle between z axis and b

    theta   = aCos ( bUnit_car[2] )

    ;   calculate the quaternions

    q0  = cos ( theta / 2.0 )
    q1  = sin ( theta / 2.0 ) * perpUnit[0]
    q2  = sin ( theta / 2.0 ) * perpUnit[1] 
    q3  = sin ( theta / 2.0 ) * perpUnit[2]

    ;   construct the rotation matrix

    rotQ    = [ [ q0^2+q1^2-q2^2-q3^2, 2*(q1*q2-q0*q3), 2*(q1*q3+q0*q2) ], $
                [ 2*(q2*q1+q0*q3), q0^2-q1^2+q2^2-q3^2, 2*(q2*q3-q0*q1) ], $
                [ 2*(q3*q1-q0*q2), 2*(q3*q2+q0*q1), q0^2-q1^2-q2^2+q3^2 ] ]

    inv_rotQ    = transpose ( rotQ )

    testRot = rotQ ## zaxis

    ;@load_colors
    ;iPlot, [0,bUnit_car[0]],[0,bUnit_car[1]],[0,bUnit_car[2]], $
    ;       thick = 4, color = blue, /iso
    ;iPlot, [0,perp[0]],[0,perp[1]],[0,perp[2]], thick = 4, /over, color = green
    ;iPlot, [0,testRot[0]],[0,testRot[1]],[0,testRot[2]], thick = 4, /over




;	testVec	= [ [ 0d0 ], [ 0d0 ], [ 1d0 ] ]
;	print, 'Test0', testVec[*]
;    iPlot, [0,bUnit_car[0]],[0,bUnit_car[1]],[0,bUnit_car[2]], $
;           thick = 4, color = blue, /iso
;    iPlot, [0,testVec[0]],[0,testVec[1]],[0,testVec[2]], thick = 4, /over
;
;	testVec	= rot_x ## testVec
;	print, 'Test1', testVec[*]
;    iPlot, [0,testVec[0]],[0,testVec[1]],[0,testVec[2]], /over, thick = 2, color = red
;
;	testVec	= rot_y ## testVec
;	print, 'Test2', testVec[*]
;    iPlot, [0,testVec[0]],[0,testVec[1]],[0,testVec[2]], /over, thick = 2, color = green
;
;	testVec	= rot_z ## testVec
;    iPlot, [0,testVec[0]],[0,testVec[1]],[0,testVec[2]], /over, thick = 2, color = purple
;
;	print, 'Test3', testVec[*]
;	print, 'BUnit', (bUnit_car[*])[*]
;
;	epsilonOut	= rot_x ## epsilonIn ## inv_rot_x 
;	epsilonOut	= rot_y ## epsilonOut ## inv_rot_y
;	epsilonOut	= rot_z ## epsilonOut ## inv_rot_z

    epsilonOut  = rotQ ## epsilonIn ## inv_rotQ

	return, epsilonOut
end


