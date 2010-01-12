function rotateEpsilon, epsilonIn, bUnit_car

	;	y towards z
	rotTh1	=	-(!dpi/2d0-aTan ( bUnit_car[2], bUnit_car[1] ))
	rot_x	= [ [ 1, 0, 0], $
				[ 0, cos ( rotTh1 ), -sin ( rotTh1 ) ], $
				[ 0, sin ( rotTh1 ), cos ( rotTh1 ) ] ]
	inv_rot_x	= transpose ( rot_x )

	;	z towards x
	rotTh2	=	0d0;-!dpi/2d0;aTan ( bUnit_car[2], bUnit_car[0] )
	rot_y	= [ [ cos ( rotTh2 ), 0, -sin ( rotTh2 ) ], $
				[ 0, 1, 0 ], $
				[ sin ( rotTh2 ), 0, cos ( rotTh2 ) ] ]
	inv_rot_y	= transpose ( rot_y )

	;	x towards y
	rotTh3	=	-aTan ( bUnit_car[0], bUnit_car[1] ) 
	rot_z	= [ [ cos ( rotTh3 ), -sin ( rotTh3 ), 0 ], $
				[ sin ( rotTh3 ), cos ( rotTh3 ), 0 ], $
				[ 0, 0, 1 ] ]
	inv_rot_z	= transpose ( rot_z )

	;print, rotth1*!radeg, rotth2*!radeg, rotth3*!radeg

	;	test rotations
	;	There seems to be the requirement to rotate z towards x
	;	by -!pi/2 to get everything to work correctly. I think
	;	this is due to the ambiguity in the orientation of the x/y
	;	axis in the stix frame. So I'm still not 100% on this, but
	;	without the extra rotation the dispersion calculation for
	;	any case with a z-component of the poloidal field does not 
	;	agree (nor show any +/- asymmetry in kR) with the dispersion 
	;	calculation for kPerp based on the PRL of Jaeger, V90 N19, May 2003. 
	;	Update: It seems that perhaps the dispersion calculation of 
	;	Jaeger is not particularly useful for seperating the components
	;	of the poloidal field and that the additional rotation that
	;	I was unsure of may not be needed at all. The wave solution
	;	seems to show not much change at all after including a vertical 
	;	poloidal field, perhaps because the fast wave is perpendicular?
	;	Anyway, the results seem to agree with AORSA1D within a minus 
	;	sign.

	;testVec	= [ [ 0d0 ], [ 0d0 ], [ 1d0 ] ]
	;print, 'Test0', testVec[*]

	;testVec	= rot_x ## testVec
	;print, 'Test1', testVec[*]

	;testVec	= rot_y ## testVec
	;print, 'Test2', testVec[*]

	;testVec	= rot_z ## testVec

	;print, 'Test3', testVec[*]
	;print, 'BUnit', (bUnit_car[*])[*]

	epsilonOut	= rot_x ## epsilonIn ## inv_rot_x 
	epsilonOut	= rot_y ## epsilonOut ## inv_rot_y
	epsilonOut	= rot_z ## epsilonOut ## inv_rot_z

	return, epsilonOut
end


