function rotateEpsilon, epsilonIn, bu_rtz

    rot_abp_to_rtz = get_RotMat_abp_to_rtz(bu_rtz)

    epsilonOut  = rot_abp_to_rtz ## epsilonIn ## transpose(rot_abp_to_rtz)

	return, epsilonOut

end


