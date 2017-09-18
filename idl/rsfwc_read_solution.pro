function rsfwc_read_solution, runFolderName

    SolutionFile = expand_path(runFolderName)+'/output/rs_solution.nc'

    cdfId = nCdf_open(SolutionFile)

        ;nCdf_varGet, cdfId, 'freq', freq

        nCdf_varGet, cdfId, 'r', r
        nCdf_varGet, cdfId, 'r_', r_
        nCdf_varGet, cdfId, 'z', z
        nCdf_varGet, cdfId, 'z_', z_

        ;nCdf_varGet, cdfId, 'B0_r', b_r 
        ;nCdf_varGet, cdfId, 'B0_p', b_t 
        ;nCdf_varGet, cdfId, 'B0_z', b_z 

        nCdf_varGet, cdfId, 'e_r_re', e_r_re 
        nCdf_varGet, cdfId, 'e_p_re', e_t_re 
        nCdf_varGet, cdfId, 'e_z_re', e_z_re 
        nCdf_varGet, cdfId, 'e_r_im', e_r_im 
        nCdf_varGet, cdfId, 'e_p_im', e_t_im 
        nCdf_varGet, cdfId, 'e_z_im', e_z_im 

        e_r = complex(e_r_re,e_r_im)
        e_t = complex(e_t_re,e_t_im)
        e_z = complex(e_z_re,e_z_im)

        nCdf_varGet, cdfId, 'jP_r_re', jP_r_re 
        nCdf_varGet, cdfId, 'jP_p_re', jP_t_re 
        nCdf_varGet, cdfId, 'jP_z_re', jP_z_re 
        nCdf_varGet, cdfId, 'jP_r_im', jP_r_im 
        nCdf_varGet, cdfId, 'jP_p_im', jP_t_im 
        nCdf_varGet, cdfId, 'jP_z_im', jP_z_im 

        nCdf_varGet, cdfId, 'jP_r_re_spec', jP_r_re_spec 
        nCdf_varGet, cdfId, 'jP_p_re_spec', jP_t_re_spec 
        nCdf_varGet, cdfId, 'jP_z_re_spec', jP_z_re_spec 
        nCdf_varGet, cdfId, 'jP_r_im_spec', jP_r_im_spec 
        nCdf_varGet, cdfId, 'jP_p_im_spec', jP_t_im_spec 
        nCdf_varGet, cdfId, 'jP_z_im_spec', jP_z_im_spec 

        jP_r = complex(jP_r_re,jP_r_im)
        jP_t = complex(jP_t_re,jP_t_im)
        jP_z = complex(jP_z_re,jP_z_im)

        jP_r_spec = complex(jP_r_re_spec,jP_r_im_spec)
        jP_t_spec = complex(jP_t_re_spec,jP_t_im_spec)
        jP_z_spec = complex(jP_z_re_spec,jP_z_im_spec)

        ;nCdf_varGet, cdfId, 'jA_r_re', jA_r_re 
        ;nCdf_varGet, cdfId, 'jA_p_re', jA_t_re 
        ;nCdf_varGet, cdfId, 'jA_z_re', jA_z_re 
        ;nCdf_varGet, cdfId, 'jA_r_im', jA_r_im 
        ;nCdf_varGet, cdfId, 'jA_p_im', jA_t_im 
        ;nCdf_varGet, cdfId, 'jA_z_im', jA_z_im 

        ;jA_r = complex(jA_r_re,jA_r_im)
        ;jA_t = complex(jA_t_re,jA_t_im)
        ;jA_z = complex(jA_z_re,jA_z_im)

        nCdf_close, cdfId

        solution = { $
                r: r, $
                r_: r_, $
                z: z, $
                z_: z_, $
                x: r, $
                y: z, $
                E_r: E_r, $
                E_t: E_t, $
                E_z: E_z, $
                jP_r: jP_r_spec, $
                jP_t: jP_t_spec, $
                jP_z: jP_z_spec } ;, $
                ;jP_r_spec: jP_r_spec, $
                ;jP_t_spec: jP_t_spec, $
                ;jP_z_spec: jP_z_spec  }
                
        return, solution
end
