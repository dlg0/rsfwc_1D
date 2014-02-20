pro rsfwc_read_iterations, nIterations, MPE_FileName = MPE_FileName

    RightName = '~/scratch/rsfwc_1d/ar2_vorpal/right_simple_'
    CompareDir = '~/scratch/rsfwc_1d/ar2_vorpal/simple_full'

    for n=0,nIterations-1 do begin

        right_ThisDir =  RightName+StrCompress(string(n),/rem)
        print, 'Reading ... ', right_ThisDir
        right_s = rsfwc_read_solution (right_ThisDir)

        if jRAll eq !null then begin
            nR = n_elements(right_s.jP_r)
            jRAll = ComplexArr(nR,nIterations)
            jTAll = ComplexArr(nR,nIterations)
            jZAll = ComplexArr(nR,nIterations)
        endif

        jRAll[*,n] = right_s.jP_r
        jTAll[*,n] = right_s.jP_t
        jZAll[*,n] = right_s.jP_z

    endfor

    compare_s = rsfwc_read_solution (CompareDir)

    jR_mpe = kj_mpe(jRAll)
    jT_mpe = kj_mpe(jTAll)
    jZ_mpe = kj_mpe(jZAll)

    range = max(abs([jR_mpe,jT_mpe,jZ_mpe]))

    p=plot(right_s.r, jRAll[*,0],layout=[1,3,1], window_title=MPE_FileName)
    for n=1,nIterations-1 do p=plot(right_s.r,jRAll[*,n],/over)
    p=plot(compare_s.r, compare_s.jP_r, color='g',thick=3,/over)
    p=plot(right_s.r, jR_mpe,/over,color='b',thick=3,yRange=[-range,range])

    p=plot(right_s.r, jTAll[*,0],layout=[1,3,2],/current)
    for n=1,nIterations-1 do p=plot(right_s.r,jTAll[*,n],/over)
    p=plot(compare_s.r, compare_s.jP_t, color='g',thick=2,/over)
    p=plot(right_s.r, jT_mpe,/over,color='b',thick=2,yRange=[-range,range])
    p=plot(compare_s.r, imaginary(compare_s.jP_t), color='g',thick=2,/over,linestyle='--')
    p=plot(right_s.r, imaginary(jT_mpe),/over,color='b',thick=2,linestyle='--')

    p=plot(right_s.r, jZAll[*,0],layout=[1,3,3],/current)
    for n=1,nIterations-1 do p=plot(right_s.r,jZAll[*,n],/over)
    p=plot(compare_s.r, compare_s.jP_z, color='g',thick=2,/over)
    p=plot(right_s.r, jZ_mpe,/over,color='b',thick=2,yRange=[-range,range])
    p=plot(compare_s.r, imaginary(compare_s.jP_z), color='g',thick=2,/over,linestyle='--')
    p=plot(right_s.r, imaginary(jZ_mpe),/over,color='b',thick=2,linestyle='--')

    fName = 'restart_mpe.nc'
    if keyword_set(MPE_FileName) then fName = MPE_FileName
	nc_id = nCdf_create (fName, /clobber )

	nCdf_control, nc_id, /fill
	
	nr_id = nCdf_dimDef ( nc_id, 'nR', n_elements(right_s.R) )
	nrH_id = nCdf_dimDef ( nc_id, 'nR_', n_elements(right_s.R_) )
	scalar_id = nCdf_dimDef ( nc_id, 'scalar', 1 )

	r_id = nCdf_varDef ( nc_id, 'r', nr_id, /float )
	rH_id = nCdf_varDef ( nc_id, 'r_', nrH_id, /float )
	z_id = nCdf_varDef ( nc_id, 'z', nr_id, /float )
	zH_id = nCdf_varDef ( nc_id, 'z_', nrH_id, /float )

	jP_r_re_id = nCdf_varDef ( nc_id, 'jP_r_re', nr_id, /float )
	jP_r_im_id = nCdf_varDef ( nc_id, 'jP_r_im', nr_id, /float )
	jP_p_re_id = nCdf_varDef ( nc_id, 'jP_p_re', nr_id, /float )
	jP_p_im_id = nCdf_varDef ( nc_id, 'jP_p_im', nr_id, /float )
	jP_z_re_id = nCdf_varDef ( nc_id, 'jP_z_re', nr_id, /float )
	jP_z_im_id = nCdf_varDef ( nc_id, 'jP_z_im', nr_id, /float )

	nCdf_control, nc_id, /enDef

	nCdf_varPut, nc_id, r_id, right_s.r
	nCdf_varPut, nc_id, rH_id, right_s.r_

	nCdf_varPut, nc_id, z_id, right_s.z
	nCdf_varPut, nc_id, zH_id, right_s.z_

	nCdf_varPut, nc_id, jP_r_re_id, real_part(jr_mpe)
	nCdf_varPut, nc_id, jP_r_im_id, imaginary(jr_mpe) 
	nCdf_varPut, nc_id, jP_p_re_id, real_part(jt_mpe) 
	nCdf_varPut, nc_id, jP_p_im_id, imaginary(jt_mpe) 
	nCdf_varPut, nc_id, jP_z_re_id, real_part(jz_mpe) 
	nCdf_varPut, nc_id, jP_z_im_id, imaginary(jz_mpe) 

	nCdf_close, nc_id

end
