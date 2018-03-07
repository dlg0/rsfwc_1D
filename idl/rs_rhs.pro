function rs_rhs

    @dlg_constants

    ; Load the E field from file

    rs = rs_read_solution('./') 
    f = rs('freq')

    ; Cold jP_c

    jP_c_r = rs('jP_r')
    jP_c_t = rs('jP_t')
    jP_c_z = rs('jP_z')

    ; Load the A matrix

    restore, 'rs-amat.sav'

    RHS0 = RHS

    ; Hot jP_h

    kj_stix_current, hot=1, useRS=1, jr=jP_h_r, jt=jP_h_t, jz=jP_h_z

    ; Sum over species

    jP_h_r = total(jP_h_r,2)
    jP_h_t = total(jP_h_t,2)
    jP_h_z = total(jP_h_z,2)

    delta_jP_r = jP_h_r - jP_c_r
    delta_jP_t = jP_h_t - jP_c_t
    delta_jP_z = jP_h_z - jP_c_z

    ; Put the Jt and Jz components on the half grid

    delta_jP_t_re_ = interpol(real_part(delta_jP_t), r, r_, /spline)
    delta_jP_z_re_ = interpol(real_part(delta_jP_z), r, r_, /spline)
    delta_jP_t_im_ = interpol(imaginary(delta_jP_t), r, r_, /spline)
    delta_jP_z_im_ = interpol(imaginary(delta_jP_z), r, r_, /spline)
 
    delta_jP_t_ = dcomplex(delta_jP_t_re_,delta_jP_t_im_)
    delta_jP_z_ = dcomplex(delta_jP_z_re_,delta_jP_z_im_)

    ; Reconstruct the jP column vector
    ; This is complicated by the half grid and BCs used in RS.

    nR = n_elements(rs('r'))

    delta_jP = dcomplexArr(nR*3+2)

	ii_R	= 2 + lIndGen(nR)*3
	ii_T	= ii_R[0:nR-2]+1
	ii_Z	= ii_T+1

    delta_jP[ii_R] = delta_jP_r
    delta_jP[ii_T] = delta_jP_t_
    delta_jP[ii_Z] = delta_jP_z_

    delta_jP[0]  = 0 
    delta_jP[1]  = 0 
    delta_jP[-2] = 0 
    delta_jP[-1] = 0 

    ; Add into the RHS

    w = 2*!pi*f

    ;RHS = RHS0 
    RHS = RHS0 +_ii*w*_u0*delta_Jp

    return, RHS

end
