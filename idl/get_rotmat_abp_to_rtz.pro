function cross, v1,v2

    r1 = v1[1]*v2[2]-v1[2]*v2[1]
    r2 = -(v1[0]*v2[2]-v1[2]*v2[0])
    r3 = v1[0]*v2[1]-v1[1]*v2[0]

    return, transpose([r1,r2,r3])

end

function mag, v

    return, sqrt(v[0]^2+v[1]^2+v[2]^2)

end 

function dot, v1,v2

    return, v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]

end 

function get_RotMat_abp_to_rtz, bu_rtz

;   Try a different rotation approach - COPY OF THE AORSA CODE
;   ---------------------------------

    ;   get vector perp to both z axis and b

    ru_rtz = [[1],[0],[0]] ; these are 1 column x 3 row vectors
    tu_rtz = [[0],[1],[0]]
    zu_rtz = [[0],[0],[1]]

    pu_rtz = transpose(bu_rtz[*])

    a_rtz = cross(zu_rtz,pu_rtz)
    au_rtz = a_rtz/mag(a_rtz)

    b_rtz = cross(pu_rtz,au_rtz)
    bu_rtz = b_rtz/mag(b_rtz)

if keyword_set(debug) then begin

        au_rtz2 = au_rtz
        bu_rtz2 = bu_rtz
        pu_rtz2 = pu_rtz

        print, 'au_rtz:', au_rtz
        print, 'bu_rtz:', bu_rtz
        print, 'pu_rtz:', pu_rtz

        epsilon0 = epsilonIn

endif

    ; Rotation 1

    theta = acos( dot(ru_rtz,au_rtz) )

if keyword_set(debug) then begin
        print, 'theta: ', theta*!radeg
endif

    q0  = cos ( theta / 2.0 )
    q1  = sin ( theta / 2.0 ) * (-zu_rtz[0]) 
    q2  = sin ( theta / 2.0 ) * (-zu_rtz[1])  
    q3  = sin ( theta / 2.0 ) * (-zu_rtz[2]) 

    ;   construct the rotation matrix

    rot1    = [ [ q0^2+q1^2-q2^2-q3^2, 2*(q1*q2-q0*q3), 2*(q1*q3+q0*q2) ], $
                [ 2*(q2*q1+q0*q3), q0^2-q1^2+q2^2-q3^2, 2*(q2*q3-q0*q1) ], $
                [ 2*(q3*q1-q0*q2), 2*(q3*q2+q0*q1), q0^2-q1^2-q2^2+q3^2 ] ]

if keyword_set(debug) then begin
    print, rot1
endif

    au_rtz = rot1 ## au_rtz
    bu_rtz = rot1 ## bu_rtz
    pu_rtz = rot1 ## pu_rtz

if keyword_set(debug) then begin
    print, 'au_rtz 1: ', au_rtz[*]
    print, 'bu_rtz 1: ', bu_rtz[*]
    print, 'pu_rtz 1: ', pu_rtz[*]

    epsilon1 = rot1 ## epsilon0 ## transpose(rot1)
    print, 'epsilon1: '
    print, epsilon1
endif

    ; Rotation 2

    theta = acos( dot(zu_rtz,pu_rtz) )

if keyword_set(debug) then begin
    print, 'theta: ', theta * !radeg
endif
    q0  = cos ( theta / 2.0 )
    q1  = sin ( theta / 2.0 ) * (-ru_rtz[0]) 
    q2  = sin ( theta / 2.0 ) * (-ru_rtz[1])  
    q3  = sin ( theta / 2.0 ) * (-ru_rtz[2]) 

    ;   construct the rotation matrix

    rot2    = [ [ q0^2+q1^2-q2^2-q3^2, 2*(q1*q2-q0*q3), 2*(q1*q3+q0*q2) ], $
                [ 2*(q2*q1+q0*q3), q0^2-q1^2+q2^2-q3^2, 2*(q2*q3-q0*q1) ], $
                [ 2*(q3*q1-q0*q2), 2*(q3*q2+q0*q1), q0^2-q1^2-q2^2+q3^2 ] ]

if keyword_set(debug) then begin
    print, rot2
endif

    au_rtz = rot2 ## au_rtz
    bu_rtz = rot2 ## bu_rtz
    pu_rtz = rot2 ## pu_rtz

if keyword_set(debug) then begin
    print, 'au_rtz 2: ', au_rtz[*]
    print, 'bu_rtz 2: ', bu_rtz[*]
    print, 'pu_rtz 2: ', pu_rtz[*]
    epsilon2 = rot2 ## epsilon1 ## transpose(rot2)
    print, 'epsilon2: '
    print, epsilon2
endif

    rot_ = rot2 ## rot1
    ;rot_ = transpose(rot_)

if keyword_set(debug) then begin
    print, rot_

    print, 'ru_abp: ', (rot_ ## ru_rtz)[*]
    print, 'tu_abp: ', (rot_ ## tu_rtz)[*]
    print, 'zu_abp: ', (rot_ ## zu_rtz)[*]

    print, 'rot_ ## au_rtz2: ', (rot_ ## au_rtz2)[*]
    print, 'rot_ ## bu_rtz2: ', (rot_ ## bu_rtz2)[*]
    print, 'rot_ ## pu_rtz2: ', (rot_ ## pu_rtz2)[*]
endif

    rot_abp_to_rtz = transpose(rot_)

	return, rot_abp_to_rtz
end


