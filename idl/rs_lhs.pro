function rs_lhs

    ; Load the E field from file

    rs = rs_read_solution('./') 

    ; Load the A matrix

    restore, 'rs-amat.sav'

    ; Put the Et and Ez components on the half grid

    E_r = rs('E_r')
    E_t = rs('E_t')
    E_z = rs('E_z')

    E_t_re_ = interpol(real_part(E_t), r, r_, /spline)
    E_z_re_ = interpol(real_part(E_z), r, r_, /spline)
    E_t_im_ = interpol(imaginary(E_t), r, r_, /spline)
    E_z_im_ = interpol(imaginary(E_z), r, r_, /spline)
 
    E_t_ = complex(E_t_re_,E_t_im_)
    E_z_ = complex(E_z_re_,E_z_im_)

    ; Reconstruct the E column vector
    ; This is complicated by the half grid and BCs used in RS.

    nR = n_elements(rs('r'))

    E = complexArr(nR*3+2)

	ii_eR	= 2 + lIndGen(nR)*3
	ii_eT	= ii_eR[0:nR-2]+1
	ii_eZ	= ii_eT+1

    E[ii_eR] = E_r
    E[ii_eT] = E_t_
    E[ii_eZ] = E_z_

    E[0] = eField[0]
    E[1] = eField[1]
    E[-2] = eField[-2]
    E[-1] = eField[-1]

    ; Apply A

    ;LHS = Amat ## Efield
    LHS = Amat ## E

    return, reform(LHS)

end
