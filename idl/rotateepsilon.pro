function rotateEpsilon, epsilonIn, bu_rtz, w, debug = debug, i = i

    rot_abp_to_rtz = get_RotMat_abp_to_rtz(bu_rtz)

    if keyword_set(i) then begin
            if i eq 256 then print, rot_abp_to_rtz
    endif

    epsilonOut  = rot_abp_to_rtz ## epsilonIn ## transpose(rot_abp_to_rtz)

if keyword_set(debug) then begin

    @constants
    print, 'epsilonCheck against rot1 then rot2 (epsilon2 above):'
    epsilonCheck  = rot_ ## epsilonIn ## transpose(rot_)
    print, epsilonCheck
    print, ''
    print, 'epsilonOut:'
    print, epsilonOut
    thisSigma = (epsilonOut-identity(3))*w*e0/II
    print, 'sigmaOut [rotated to RTZ]:'
    print, thisSigma
stop 
endif

	return, epsilonOut
end


