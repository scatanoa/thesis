module hg_model_fldbin_ubc_qt
    !13dec2019: subroutines are without modularize. Only its future MODULE is commented above in capital letters. Final subroutine is the main: model_sc.
!	implicit none
	use ddt_class_fld
	use balance_fldbin_ubc_qt

!!	in newton_phd.f90:
!    use newton_raphson
!    use functions
	
!	use bin_io
	
!	integer i,j,k	!index for loops in 4D arrays: time, space, grainSize and widthPortion
	integer i1, j1, k1	!i and j are used by ddt_class.
	

contains

	subroutine write_FloatArrIJ(fileunit,path,vect,record,N_x)  	!read function is in other f90
		!Variables de entrada
		integer, intent(in) :: fileunit,record, N_x
		character*255, intent(in) :: path
		real, intent(in) :: vect(N_x)
		!Variables internas
		character*10 stt
		!Escritura     
		stt='old'
		if (record.eq.1) stt='replace'
		open(fileunit,file=path,form='unformatted',status=stt,access='direct',RECL=4*N_x)
			write(fileunit,rec=record) vect
		close(fileunit)
	end subroutine
    
    subroutine update_geom(nR, jp, ju, jd, bl, fs, SoMin, Sor, dx, zx, So)
        integer, intent(in) :: nR, jp, ju, jd, fs
            !spatially smoothed (central difference?) as Parker ebook, ch17, pg11
        real, intent(in) :: bl, SoMin, Sor, dx, zx(nR)
        real, intent(out) :: So
        integer kdx        
        real z_up, z_down
                                       
        if (fs.eq.1) then!                    
            if (jd.gt.0 .and. ju.gt.0) then  
            !intermediate nodes, with neighbors both up and downstream    
                z_up = zx(ju)
                z_down = zx(jd)
                kdx = 2
		    elseif (jd.eq.0) then  !final node= basin outlet		        
	            z_up = zx(jp)
		        z_down = bl  !weir?
		        kdx = 1
		    else   !initial node= headwater			        
	            z_up = zx(jp)
                z_down = zx(jd)
	            kdx = 1
		    endif
        else
            kdx = 1
            z_up = zx(jp)
            if (jp.lt.nR) then              
                z_down = zx(jd)
            else
                z_down = bl
            endif
        endif 
		    So = Sor + (z_up - z_down)/(kdx*dx)		    
		    So = max(SoMin, So)
    end subroutine
    
    
    integer function iplot(flag_day_else_tStep, i1, day)
        integer, intent(in) :: flag_day_else_tStep, i1, day   
	    if (flag_day_else_tStep.eq.0) then
		    iplot = i1
		else
		    iplot = day
	    endif
    end function
          
    function f_acu(C_else_F, nD, f) result(facu)
        logical :: C_else_F
        integer, intent(in) :: nD
        real, intent(in) :: f(nD)
        integer :: ii
        real :: facu(nD)
        do ii=1, nD
            facu(ii) = sum(f(1:ii))
        end do
        if (abs(facu(nD)-1) .gt. .05) then
            if (C_else_F) then
                print*, 'facu BUG in channel'
            else
                print*, 'facu BUG in floodplain'            
            endif                
            print*, 'facu(nD) = ', facu(nD)
            print*, 'take this t,x to edit print_range for DEBUG: i,j:', i1,j1
            stop 'check f_acu, whose greater is not 1'
        endif           
    end function
    
    function f_coefs_wqsy(p) result(c)
        !see SCphd notes p148dr.
        real, intent(in) :: p(5)
        real c(8)
        real beta, u, u2, r2, ur, u2r2, ka, fi
        beta = p(1);  fi = p(2); ka = p(5)
        r = p(3);  u = p(4)
        u2 = u**2;  r2 = r**2;  ur = u*r; u2r2 = u2*r2
        c(1) = beta;  c(2) = fi; c(3) = ka        
        c(4) = 1 - u2r2
        c(5) = (1+ur)**2
        c(6) = 4*u2r2
        c(7) = -2*ur*(1+ur)        
        c(8) = u2
    end function        
        
    subroutine get_neighbors(Rd, nR, jup, jdown)
        !Rd: reach downstream of each one, from topology array X.
        integer, intent(in) :: nR
        real, intent(in) :: Rd(nR)
        integer, intent(out) :: jdown(nR), jup(nR, 2)
        integer Rdi(nR), jdd

        jup = 0;  jdown = 0
        Rdi = int(Rd)
		do j1=1, nR
		    jdd = Rdi(j1)
		    if (jdd .gt. 0) then  !if it is not basin outlet
		        if (jup(jdd, 1) .eq. 0) then   !NO upstream reach had been assigned to reach jdd yet.
    		        jup(jdd, 1) = j1
                else     !upstream reach had been assigned to reach jdd, so fill the 2nd column, as jdd results from confluence.
                    jup(jdd, 2) = j1
                endif                    
    		else
		        if (jup(jdd, 1) .eq. 0) then
    		        jup(jdd, 1) = nR-1
                else
                    jup(jdd, 2) = nR-1
                endif                    
    		endif    		
    		jdown(j1) = Rdi(j1)
		end do
    end subroutine

	subroutine model_sc(iref_ini, iref_fin, jref_ini, jref_fin, &   !for printing.
	output_folder, pos_ini_var_py, &    !for write.
	which_flume, &
    flume_else_basin, flagSupplyWash, flagSupplySlide, flagSo_smooth, & !flags
    flag_day_else_tStep, flag_printDetails, flg_backwat,  & 	!flags
	ddt_args_py, dt_dQ_crs_arr_py, dt_py, tb_hill, & !temporal resolution.
    path_list_py, &   !for write; here after declaring size via sizers.
	h_s_ini, w_b_ini, fl_fg_py, fl_fh_py, fl_fqs_py, fl_fw_py, fl_fwqsy_py, fl_gtm_py, calib, calibv, &
	X_py, fb_py, Dfbr_py, slid_sh_py, wash_sh_py, Qhill_py, & !physical.
	Qh, Q, wash_sh, slid_sh, v, h, vf, w, c, D50c, D84c, D50c_ss, D84c_ss, So, & 	   !physical, intent(out) to python.
	h_s, h_sf, w_c, w_b, w_v, ts_c, tc_c, tsm_c, Fr, aspect, subm, Mgtm, gam, pi, & 	   !physical, intent(out) to python.
	nD, nt, nR, n_topol_vars, num_paths, nvars, num_char, n_coarse_dt)
	    !vars "n.." (sizers) are optional in f2py .__doc__. so they come at last.
	
	    !VARIABLE DEFINITION:
	    
!		!Flags/parameters related to model type is 1st, do not move following line, so it can be visible:
!		flagSupply=0; flagSo_smooth=1; 
!		!output printing options:
!		flag_printDetails=1; flag_day_else_tStep=0;
!		flag_trench=1;
		
		!input variables:
        integer, intent(in) :: nD, nt, nR, n_topol_vars, num_paths, nvars, num_char, n_coarse_dt
		integer, intent(in) :: flume_else_basin, flagSupplyWash, flagSupplySlide, flagSo_smooth, flag_day_else_tStep
		integer, intent(in) :: flag_printDetails, flg_backwat, which_flume
        integer, intent(in) :: ddt_args_py(0:6-1)
		integer, intent(in) :: iref_ini, iref_fin, jref_ini, jref_fin, pos_ini_var_py(0:nvars-1)
		real, intent(in) :: fl_fg_py(0:6-1), fl_fh_py(0:9-1), fl_fqs_py(0:3-1)
		real, intent(in) :: fl_fw_py(0:6-1), fl_fwqsy_py(0:6-1), fl_gtm_py(0:7-1)
		real, intent(in) :: X_py(0:nR-1, 0:n_topol_vars-1), fb_py(0:nD-1), Dfbr_py(0:nD-1)
		real, intent(in) :: tb_hill, h_s_ini, w_b_ini, calib, calibv, dt_dQ_crs_arr_py(0:n_coarse_dt-1), dt_py(0:nt-1)
		real, intent(in), dimension(0:n_coarse_dt-1, 0:nR-1) :: Qhill_py, wash_sh_py, slid_sh_py
		character, intent(in), dimension(0:num_char-1, 0:num_paths-1) :: path_list_py
		
		!17apr: this block is intent(out) but this clause was erased after p6 in:
    		!http://folk.uio.no/inf3330/scripting/doc/python/fc/f2py-tutorial/tutorial.pdf
		real, intent(out), dimension(nt, nR) :: Qh, Q, wash_sh, slid_sh, aspect, subm, Mgtm, gam
		real, intent(out), dimension(nt, nR) :: v, h, vf, w, c, So, h_s, h_sf, w_c, w_b, ts_c, tc_c, tsm_c, Fr
		real, intent(out), dimension(nt, nR) :: D50c, D84c, D50c_ss, D84c_ss
		real, intent(out), dimension(nt, nR, nD) :: pi
		real, intent(out)  :: w_v(nR)
		
		!parameters:
		real :: base_level, SoMin, D84r, W_valley, dry_repos_ang
		real :: g = 9.8
		real :: a_s, b_s, a_w, b_w, a_1, a_2, b_1, b_2, b_3, a_qs, b1_qs, b2_qs
		real :: ta, ta4, c_w_Qf, ang_rad, ca, sa, kLa, ta_inv_sqd
		
!		dry_repos_ang = 30*3.1416/180
!        poros=.3	!cui07
!        D84r=.03	!Mar23, 2018: bedrock roughness. May27, 2018: 
!        tb_hill = 6 !!default: 6 !hours. Must be related to hillslope area. 6 hours is assumed for 4km2	

		!input-receiving variables:
		character, dimension(num_char, num_paths) :: path_list
		character*255  :: output_folder
		integer :: pos_ini_var(nvars), ddt_args(6)
		real :: X(nR, n_topol_vars), dx(nR), fl_fg(6), fl_fh(9), fl_fqs(3), fl_fw(6), fl_fwqsy(6), fl_gtm(7)
		real, dimension(n_coarse_dt, nR) :: Qhill_in, wash_sh_in, slid_sh_in    !to pass simple fortran indexes >=1 in ddt subrout.
		real :: fb(nD), Dfbr(nD), dt(nt), dt_dQ_crs_arr(n_coarse_dt)
		    !26mar20. Rigth side of slice of GSD whose probability is fb.
		
		!main result variables:
		real, dimension(nt, nR) :: Qf, D84f
            
        !module setup:
        integer wv, range2zoomDdt(3), day_of_calc(nt)    !for ddt.
        real Qhill_mean(nR)  !for ddt.
        real Sor(nR), psolid !for geom subroutine.
            !for boundaries subroutine.
            !for init_conds subroutine.
        logical :: fww(n_coarse_dt, nR), print_cnd, C_else_F
        
        integer :: contAdverseSlope = 0  !see comment below.
		
		!module travel_hgy:		
!		real roug, yp, discr, sum_y, sum_yz, sum_yz2, sum_w, sum_wz, aa, bb, cc, zw1, zw2,d1,d2,v1,v2,Fr1,Fr2,yy,zz, hm
            !for iterD84hco subroutine. Dependerá de nuevo método.
		integer :: prt_cd, minifld, repeater, jdown(nR),jup(nR, 2)     !if reach receives confluence so it has 2 upstream reaches
		real :: errQadm(nR), Qm_appr(nR), q_prev(nR), Rw_prev(nR), Qup, Rw_ppr, cumWatErr_m3, cumWat_m3, err_m3
		real :: HS_ord(nR), Aw, sg8, Sw_prev, h_prev(nR), Iw_prev(nR), Is_prev(nR, nD), Iw(nR), Is(nR, nD)
		    !for transit subroutine(s), hgy and sed.       
        real :: coef_fu(12), coef_fu2(16), coef_fu3(14), coef_fu4(17), morp(9), hc(7)
        real :: Is_pr_allD(nt, nR), Qso_c_k(nD), Qso_f_k(nD), Qsi_k(nD), Qsi_c_k(nD), Qsi_f_k(nD)
            !for newton nonlinear solver.
        real :: Qf_try0, a12_2, a1_2, a2_2, q_dl, v_prev, hbf_n, ta_n, deepen_indx, errwc
        logical :: isFlood, lowQ_CH_narrower
        !module travel_sed:
        integer trench, np, keyDpct(3), indDp(3), dhsf_type, colm0, xtcum   !trench is for widenByIncision subroutine.
        real indDp_r(3), wc_ini(nR), Adr(nR), qfmx(nR), wb_thw(nR), hbf_thw(nR)
        real Rc, Rf, Qs, Qs_fl, Qs_c, Qspp, w_sf, net_Qs_flow, pbl(nD), pblf(nD), fac_acu(nD), faf_acu(nD)
        real fsc_acu(nD), qsbx, pf, pc, cmx, cveg, pWreveg, wc_ambit
		real :: Qsb, port(21), D90c, D90f, D50f, D90c_ss, dhs_dt, dhsf2_dt, dwb_dt, dwf_dt, La_ini 
		    !28mar20: make size(port) > nD to avoid bugs in gtm (fractional sediment transport).
		real, dimension(nt, nR) :: Lac, Laf, wproy, colm, Qbf, Qsbf, h_bf, dvol
		real, dimension(nt, nR, nD) :: fac, faf, fsc, fsf
		real :: tcum, xnp, ynp, coefs_wqsy(8), ss_trg(3)
        real :: supply_Di(nt, nR, nD), supply_tot(nt, nR)
 
        integer :: i_plot !print control

        !module writeBin
        integer fileunit, nvar
        character*255 :: path_Q, path_v, path_h, path_W, path_z, path_c, path_D84, path_D50
        character*255 :: path_zw, xx, xy(nvars)
        
        print*, '**********************************from py, in hg_model_class_fld.f90*********************************************'
        print*, 'calib, calibv=', calib, calibv
        print*, 'output_path: ', output_folder
        print*, 'nt,nR',nt,nR
        print*, 'ii,if,ji,jf', iref_ini, iref_fin, jref_ini, jref_fin
        if (flume_else_basin.eq.1 .and. flag_day_else_tStep.eq.1) &
            stop 'if flume then set NULL flag_day_else_tStep, for easy plot'
               
        prt_cd = 0; repeater = 0; minifld = 0
        Qh = 0.0; wash_sh = 0.0; slid_sh = 0.0; Q = 0.0; Qf = 0.0
        v = 0.0; h = 0.0; W = 0.0; c = 0.0; D84f = 0.0; wproy=0.0; colm=0; dvol = 0; v_prev=0
        D50c= 0.0; D84c= 0.0; D50c_ss= 0.0; D84c_ss = 0.0
        h_s = 0.0; h_sf = 0.0; w_b = 0.0; w_c = 0.0; ts_c = 0; tc_c=0; tsm_c = 0
        Lac = 0.0; Laf = 0.0; fac = 0.0; faf = 0.0; fsc = 0.0; fsf = 0.0
        Is_pr_allD = 0
        So = 0.0
        Iw_prev = 0.0; Is_prev = 0.0
        pi = 0; supply_Di = 0; colm0 = 0; xtcum = 0; cumWatErr_m3 = 0; cumWat_m3 = 0; err_m3 = 0
!        h_prev = 0.0
        		
		ddt_args = ddt_args_py;  X = X_py;  dt = dt_py; dt_dQ_crs_arr = dt_dQ_crs_arr_py  ![dt] = sec.
		fb = fb_py; Dfbr = Dfbr_py; 
		fl_fg = fl_fg_py; fl_fh = fl_fh_py; fl_fqs = fl_fqs_py; fl_fw = fl_fw_py; fl_fwqsy = fl_fwqsy_py
		fl_gtm = fl_gtm_py		
	
		Qhill_in = Qhill_py;  wash_sh_in = wash_sh_py;  slid_sh_in = slid_sh_py
        path_list = path_list_py(0:num_char-1, 0:num_paths-1)
		pos_ini_var=pos_ini_var_py		

        print*, 'max(dt),min(dt),len(dt)', maxval(dt),minval(dt),size(dt)		
		
        wc_ini = X(:,10); w_v = X(:,9); Sor= X(:,6);  dx= X(:,7);  base_level = h_s_ini
        h_s(1, 1:nR-1) = h_s_ini;  
        
        dry_repos_ang = fl_fg(1); thw_seed = fl_fg(3); SoMin = fl_fg(4); 
        kLa = fl_fg(5)
        a_s = fl_fh(1) ; b_s = fl_fh(2); a_w = fl_fh(3); b_w = fl_fh(4); 
        a_1 = calibv*fl_fh(5); a_2 = fl_fh(6); b_1 = fl_fh(7); b_2 = fl_fh(8); b_3= fl_fh(9)
        a_qs = fl_fqs(1); b1_qs = fl_fqs(2); b2_qs = fl_fqs(3)
        coefs_wqsy = f_coefs_wqsy(fl_fwqsy); psolid = 1 - .3
        
        np = size(port)
        xnp = np-1
        ynp = 1./(np-1)
!        print*, '28mar, np=', np
        if (.not.( size(port) .eq. np)) stop '28mar, in hg_model_class.f90; repair, so size(port) = np'
!        print*, '28mar, after cnd port'
        keyDpct= (/50, 84, 90/)
!        print*, '28mar, keyDpct=', keyDpct
        indDp_r = 1 + keyDpct/100.*(np-1)
!        print*, '28mar, indDp_r=', indDp_r
        indDp = nint(indDp_r)
!        print*, '28mar, indDp=', indDp
        ang_rad = dry_repos_ang*3.14/180. 
        ta_orig = tan(ang_rad); ca = cos(ang_rad); sa = sin(ang_rad); ang_eff = coefs_wqsy(3)*ang_rad
        ta = ta_orig
        ta_inv_sqd = 1/ta**2
        print*, '1may20; degrees dry_repos_ang, ka: ', dry_repos_ang, coefs_wqsy(3)
        tas = tan(ang_eff); cas = cos(ang_eff); sas = sin(ang_eff)
        ss_trg = (/sas, cas, tas/)  !subscript s: side = near trapeze corner, slope ka*ang_rad exist to get qsy to widen channel.
        ta4 = ta/4; ta2 = 2/ta; ta22 = ta2**2
        a1_2 = 6.5*6.5; a2_2 = 2.5*2.5; a12_2 = a1_2*a2_2 !for D=D84. see p207ur for ferguson ranges.
        diag_un = sqrt(1+ta**2)
        ex_w_Qf = 2.5*b_w
        cveg = fl_fwqsy(6)/3.1536e7 !26may20; [1/s]
        port = get_port(np)
                
            !22mar: even in basin, fluvial reaches start evolution without floodplain, i.e. in trapeze trenched config. 
            !This will develop height via colmatation and thalwegization.
        Adr = X(:, 5)
        if (flume_else_basin .eq. 1) then
            w_b(1, 1:nR-1) = w_b_ini
        else
            w_b(1, 1:nR-1) = .1 !17may20; wbmin=.1, see get_hs() & deprec 1.5*Adr(1:nR-1)**.42  
                                                !wb~5/6*wc: note wc defined in topol has coeficient 1.8.
        endif
        h_sf(1, 1:nR-1) = h_s(1, 1:nR-1) + .5*(wc_ini - w_b(1, 1:nR-1))*ta   !22mar, see notes p112ur.
!        print*, '7may20; wc_ini: ', wc_ini
!        print*, '7may20; w_b(1, 1:nR-1): ', w_b(1, 1:nR-1)
!        print*, '7may20; h_sf(1, 1:nR-1): ', h_sf(1, 1:nR-1)
!        print*, '7may20; h_s(1, 1:nR-1): ', h_s(1, 1:nR-1)
        w_c(1, 1:nR-1) = wc_ini  !wc_ini is w_v for trench exp and database w_c for flood exp.
        do j1=1,nR    !7jun20; p246ul; make equal to avoid wcini>wv leading to false widening need in colm=-2 at get_hs()
            if ( abs( w_c(1, j1)-w_v(j1) ) .lt. 1e-3 )  w_c(1, j1) = w_v(j1)
        enddo

        print*, 'wcini:', w_c(1, :)
        print*, 'wv:', w_v
        
        !static outlet geometry (critical section):
        w_c(:, nR) = w_v(nR)    !last reach is only rigid critical rectangular conveyor that only controls hydraulics.
                                !In jR, any sediment dynamics on nulled as wf=wv-wc=0 given wc=wv.                                
        w_b(:, nR) = w_v(nR)
        h_s(:, nR) = h_s_ini 
        h_sf(:, nR) = h_s(:, nR)
        
!        if (flume_else_basin .eq. 1) then !21may20, p226r.
        wb_thw = .1 !to unalter trench calid
!        else            
!            wb_thw = thw_seed*w_b(1, :)  !max(.1,thw_seed*w_b(1, :))
!        endif
        hbf_thw = thw_seed*(h_sf(1, :) - h_s(1, :)) 
            !given same thw_seed factor for w and h, same aspect ratio as initial crosssection is assumed.
        
        fac_acu = f_acu(.true., nD, fb)
        call get_GSDpctiles(print_cnd, nD, np, port, fac_acu, Dfbr, D50c(i1, j1), D84c(i1, j1), D90c, indDp)
        La_ini = kLa*D90c   !initially both channel and floodplain active layer is well mixed with alluvium.
        Lac(1, :) = La_ini; Laf(1, :) = La_ini
        do k1=1,nD
            fac(1, :, k1) = fb(k1)
            fsc(1, :, k1) = fb(k1)
            faf(1, :, k1) = fb(k1)            
            fsf(1, :, k1) = fb(k1)            
        end do
!        print*, '14may20; at j5 (rnd); ini fac: ', fac(1,5,:)
!        print*, '14may20; at j5 (rnd); ini fsc: ', fsc(1,5,:)
!        print*, '14may20; at j5 (rnd); ini faf: ', faf(1,5,:)
!        print*, '14may20; at j5 (rnd); ini fsf: ', fsf(1,5,:)
		range2zoomDdt = [flag_day_else_tStep, iref_ini, iref_fin]
        !apply ddt to Qhill:
        Qhill_mean = temporal_averager_2d(nR, n_coarse_dt, Qhill_in)        
        fww = FloodWhenWhere(n_coarse_dt, nR, Qhill_mean, Qhill_in)
        wv = 1  !wv: which variable (1: Q, 2: wash, 3: slide)
        call applyDdt(flag_printDetails, wv, .true., flume_else_basin, nt, nR, ddt_args, n_coarse_dt, & 
            tb_hill, dt_dQ_crs_arr, dt, range2zoomDdt, fww, Qhill_in, Qh, day_of_calc)
        wv = 2
        call applyDdt(flag_printDetails, wv, .true., flume_else_basin, nt, nR, ddt_args, n_coarse_dt, &
            tb_hill, dt_dQ_crs_arr, dt, range2zoomDdt, fww, wash_sh_in, wash_sh, day_of_calc)  !get slid_sh
        wv = 3                
        call applyDdt(flag_printDetails, wv, .false., flume_else_basin, nt, nR, ddt_args, n_coarse_dt, &
            tb_hill, dt_dQ_crs_arr, dt, range2zoomDdt, fww, slid_sh_in, slid_sh, day_of_calc)   !get slid_sh
        supply_tot = flagSupplyWash*wash_sh + flagSupplySlide*slid_sh
        
!        print*, '29ago20; check pp13 suply. 3rd time supply_tot all j: ', supply_tot(3,:)
!        stop
        
        do k1=1,nD
            supply_Di(:,:,k1) = fb(k1) * supply_tot
        enddo
!        print*, '10may20; sum(fb):', sum(fb)
        call get_neighbors(X(:, 4), nR, jup, jdown)
        
        !Water balance, for constant or transient Q according to process happening in flume or basin:        
        if (flume_else_basin.eq.1) Q = Qh
            !6may20; else Iw_prev & water balance solved with friction in balance.f90 numSolver for f3&f4. See notes p194.
            !I imagined 2 packs i,j; for flume and basin; but too long ugly duplicated code will be used only if slow PCperform.
        tcum = 0        
        hprev = .05 !asummed. used only for flume
        Qm_appr = .05*Adr        
        q_prev =  Qm_appr/ w_c(1, :) !mean water yield: 50l/s/km2; 5th column of X is A[km2]. Assume wet width is ~ bankfull.
        errQadm = .02*Qm_appr
!        qfmx = 1*Adr/ w_c(1, :) !7may20; 1m3/s/km2 flowing ALL on floodplain. Assume init w_c is near expected dyn average.
        Rw_prev = 0 !used only for basin
        HS_ord = X(:,3)

        lowQ_CH_narrower=.false.
        
        do i1= 1, nt-1
            if ((mod(i1, 100) .eq. 0) .or. (i1 .eq. nt-1)) print*, 't step #', i1, 'of ', nt
            i_plot = iplot(flag_day_else_tStep, i1, day_of_calc(i1))
            Iw = 0.0; Is = 0.0; tcum = tcum + dt(i1)
            pWreveg = cveg*dt(i1)
!            do j1= 1, nR !17apr20; at j1=nR=static outlet, only Fr=1 flow (in flood solver) and sediment bypass is calculated.

            j1=0 !10may20; to repeat j1 if prt_cd to study flood bug in 'random' locations.
            do while (j1.le.nR-1)
                j1 = j1+1
                dx_dt = dx(j1)/dt(i1)
                print_cnd = (prt_cd .eq. 1) .or. &
                    (flag_printDetails.eq.1 .and. j1.le.jref_fin .and. j1.ge.jref_ini .and. &
                    i_plot.le.iref_fin .and. i_plot.ge.iref_ini)
                call update_geom(nR, j1, jup(j1, 1), jdown(j1), base_level, flagSo_smooth, &
                    SoMin, Sor(j1), dx(j1), h_s(i1, :), So(i1, j1)) 
                    !if reach after confluence then slope is set by any(1/2) jup.                
                        !For simplicity, energy base_level is the same as that of the channel.
                sg8 = 80*So(i1, j1)                        
                if (j1 < nR)  fac_acu = f_acu(.true., nD, fac(i1, j1, :))
                call get_GSDpctiles(print_cnd, nD, np, port, fac_acu, Dfbr, D50c(i1, j1), D84c(i1, j1), D90c, indDp)
                call get_GSDpctiles(print_cnd, nD, np, port, faf_acu, Dfbr, D50f, D84f(i1, j1), D90f, indDp)
                h_bf(i1, j1) = h_sf(i1, j1) - h_s(i1, j1)
                
                !WATER BALANCE
                isFlood = .false.
                if (j1 .eq. nR) isFlood = .true. !16apr20; to solve critical depth with rigid walls (Fr=1).                
                if (flume_else_basin .eq. 0) then
                    q_dl = 0
                    Iw_prev(j1) = Iw_prev(j1) + Qh(i1, j1)*dt(i1)
                    Acd = h_bf(i1, j1) * ( w_c(i1, j1) - h_bf(i1, j1)/ta )
                    if (.not. isFlood) then
                        if (repeater .eq. 1) then
                            Rw_prev(j1) = Rw_ppr !repeater was set to 1 as prt_cd=1 led to repeat calcs.
                            q_prev(j1) = q_ppr
                        endif
                        cp = Rw_prev(j1) + Iw_prev(j1)
                        cpp = cp/dt(i1)
                        coef_fu3 = (/ cp, dx(j1), dt(i1), D84c(i1, j1), So(i1, j1), g, w_b(i1, j1), b_1, b_3, &
                            b_2, a_1, a_2, ta2, ta22 /)
                        if (print_cnd) then
                            print*, ''; print*, ''; 
                            if (repeater .eq. 1) print*, '8jun20; REPEATED RUN FOR THIS IJ'
                            print*, 'i,dt(min),tcum(h)=', i1,dt(i1)/60,tcum/3600, &
                                ', j,HS_ord,Adr=', j1, HS_ord(j1),Adr(j1), 'Qhill=', Qh(i1, j1)
                            print*, 'coef_fu3=', coef_fu3
                            print*, 'Rw_prev(j1),Qh(i1, j1)*dt(i1),Iw_prev(j1):', Rw_prev(j1),Qh(i1, j1)*dt(i1),Iw_prev(j1)
                        endif
                        call solve_flow_basin(minifld, i1,j1,print_cnd, ta, q_prev(j1), w_c(i1, j1),h_bf(i1, j1), &
                            coef_fu3, &
                            w(i1, j1), &
                            h(i1, j1), v(i1, j1), Fr(i1, j1), aspect(i1,j1), sqc, sq3c, Q(i1, j1), Qf_try0, Aw, &
                            isFlood, prt_cd)
                        q_prev(j1) = Q(i1, j1) / w(i1, j1)
                        q_dl = q_prev(j1)/sq3c
                    endif  !needed endif as through the NOflood solution try is how the flood gets discovered.
                    if(isFlood) then
                        cpp = cp/dt(i1); ![cpp]=m3/s
                        cp3 = cpp-dx_dt*Acd; cp4 = dx_dt*w_v(j1) !p196u
                        if (print_cnd) print*, 'h_bf,ta,w_c(i1, j1),Acd, w_v(j1) = ', &
                            h_bf(i1, j1),ta,w_c(i1, j1),Acd, w_v(j1)
                        if (print_cnd) print*, 'cpp,dx_dt,cp3,cp4: ', cpp,dx_dt,cp3,cp4
                        coef_fu4 = (/ b_1, b_3, b_2, a_1, a_2, g, w_c(i1, j1), w_v(j1), &
                            D84c(i1, j1), D84f(i1, j1), ta, h_bf(i1, j1), So(i1, j1), &
                            cp3, cp4, sqc, sq3c /)
                        if (j1 .eq. nR) Qup = Q(i1, j1-1)
                        call solve_flood_basin(i1,j1,nR,print_cnd, h(i1, j1), w(i1, j1), v(i1, j1), sg8, Qup, Qf_try0, Acd, &
                            dx_dt, cpp, sa, a1_2, a2_2, a12_2, &
                            coef_fu4, errQadm(j1), & 
                            h(i1, j1), v(i1, j1), &
                            vf(i1, j1), Fr(i1, j1), aspect(i1, j1), w(i1, j1), sq3f, Qf(i1, j1), Q(i1, j1), Aw, &
                            prt_cd, minifld, repeater)
                            !this output v/Fr refers to that in channel (in w=wc).
                    endif                    
		            if (jdown(j1).gt.0 .and. prt_cd .eq. 0) then     ! .and. i.lt.nt
	                    Iw(jdown(j1)) = Iw(jdown(j1)) + Q(i1, j1)*dt(i1+1)  
	                        !30may20; future dt so it's same as courant involved in transient fut mass balance  at jdown.
		            endif
                    Rw_prev(j1) = dx(j1)*Aw
                    Iw_prev(j1) = Iw(j1)
                    if (print_cnd) then
                        print*, '  after solving transient hg&hc; Q=', Q(i1, j1)
                    endif                    
                else !flume
                    if (.not. isFlood) then
                        coef_fu = (/ Q(i1, j1), D84c(i1, j1), So(i1, j1), g, w_b(i1, j1), b_1, b_3, b_2, a_1, a_2, ta, ta4 /)
                        if (print_cnd) then
                            print*, ''; print*, ''; print*, 'i,tcum(min)=', i1,tcum/60, ', j=', j1, 'Q=', Q(i1, j1)
                            print*, 'coef_fu=', coef_fu
                        endif                    
                        call solve_flow_flume(i1,j1,print_cnd, hprev, w_c(i1, j1),h_bf(i1, j1),coef_fu, w(i1, j1), h(i1, j1), &
                            v(i1, j1), &
                            Fr(i1, j1), aspect(i1,j1), sqc, sq3c, isFlood)
                        hprev = h(i1, j1) !if (i1 .gt. 1)
                    endif  !needed endif as through the NOflood solution try is how the flood gets discovered.
                    if(isFlood) then
                        coef_fu2 = (/Q(i1, j1), D84f(i1, j1), So(i1, j1), h_bf(i1, j1), g, &
                            w_c(i1, j1), w_v(j1), b_1, b_3, b_2, a_1, a_2, ta, sqc, sq3c, D84c(i1, j1)/)

!                        !27ago20
!                        print*, 'i1, j1, which_flume', i1, j1, which_flume
!                        stop               

                        call solve_flood(which_flume, ta_inv_sqd, i1,j1,nR,print_cnd, coef_fu, coef_fu2, h(i1, j1), v(i1, j1), &
                            vf(i1, j1), &
                            Fr(i1, j1), aspect(i1, j1), w(i1, j1), sq3f, Qf(i1, j1))
                            !this output v/Fr refers to that in channel (in w=wc).
                    endif
                endif               
                    !for flume, no topological Q aggregation is needed as all nodes have Q from supply file '_flm.txt'.
                    
                !SEDIMENT BALANCE:
                if (j1<nR) then
                    if (print_cnd) print*, '10may20; flagSupplyWash, flagSupplySlide, supply_tot(i1,j1):', &
                        flagSupplyWash, flagSupplySlide, supply_tot(i1,j1)
!                    if (print_cnd) print*, &'10may20; flagSupply, supply_tot(i1,j1),sum(supply_Di(i1, j1, :)),sum(Is_prev(j1, :)):',&
!                        flagSupply, supply_tot(i1,j1),sum(supply_Di(i1, j1, :)),sum(Is_prev(j1, :))
                    Is_prev(j1, :) = Is_prev(j1, :) + supply_Di(i1, j1, :)*dt(i1)
                    
!                    !27ago20; supply_tot
!                    print*, 'i1, j1, Is_prev(j1, :): ', i1, j1, Is_prev(j1, :)
!                    stop                        
                    
                    if (print_cnd) print*, '3may20; pre get_alluv_derivs; wc,hbf: ', w_c(i1, j1),h_sf(i1, j1)-h_s(i1, j1)
                    if (i1>1) then
                        colm_prev = colm(i1-1, j1)
                    else
                        colm_prev = 0
                    endif                        
                    call get_alluv_derivs(i1,j1,print_cnd, h_bf(i1, j1), Q(i1, j1), w_c(i1, j1), w_v(j1), h_s(i1, j1), &
                        h_sf(i1, j1), w_b(i1, j1), colm_prev, ta, dx(j1), dvol_dwc, dvol_dhsf, dvol_dhs, w_sf)

!                    if (print_cnd) then
!                        print*, 'flag_day_else_tStep,flag_printDetails,i1,j1,print_cnd:', &
!                            flag_day_else_tStep,flag_printDetails,i1,j1,print_cnd
!                        stop
!                    endif
                    
                    if ( faf(i1, j1, 1) .ne. faf(i1, j1, 1) ) then !NaN is not equal to itself, says stack web forum.
                        print*, 'take this t,x to edit print_range for DEBUG: i,j:', i1,j1
                        stop 'NaN faf here'
                    endif
                    faf_acu = f_acu(.false., nD, faf(i1, j1, :))                
                    call get_Qs_capac( print_cnd, coefs_wqsy, ss_trg, calib, Q(i1, j1), Qf(i1, j1), &
                        v(i1, j1), w_sf, w_v(j1), w_c(i1, j1), w_b(i1, j1), w(i1, j1), &
                        h(i1, j1), h_bf(i1, j1), sa, ta, D50c(i1, j1), D50f, D84c(i1, j1), D84f(i1, j1), So(i1, j1), i1, j1, &
                        a_qs, b1_qs, b2_qs, flume_else_basin, Acd, a1_2, a2_2, a12_2, sg8, &
                        fl_gtm, nD, np, xnp, ynp, sq3c, sq3f, port, fac(i1, j1, :), fac_acu, faf(i1, j1, :), faf_acu, Qs_c, Qs_fl, &
                        pbl, pblf, ts_c(i1, j1), tc_c(i1, j1), tsm_c(i1, j1), subm(i1,j1), Mgtm(i1,j1), gam(i1,j1), qsbx, & 
                        fic50, fif50, fibc50, cmx, Qbf(i1, j1), Qsbf(i1, j1) )
                else !rigid output means sediment bypass.
                     Qs_c = c(i1, j1-1)*Q(i1, j1-1)
                     Qs_fl = 0                                         
                endif                        
!                if (print_cnd) print*, '12apr20: passed get_Qs_capac'
                Qs = Qs_c + Qs_fl
                if (Qs .ne. Qs) then !10jun20; to solve c= NaN bug, p255ul.
                    print*, 'take this t,x to edit print_range for Qs NaN, out of get_Qs_capac(), DEBUG: i,j:', &
                        i1,j1
                    print*,'Qs_c,Qs_fl'
                    stop
                endif
!                call get_wlog(Qs_c, Q(i1, j1), w(i1, j1), dt(i1), fl_fw, wproy)
!                if (print_cnd) print*, '12apr20: passed get_wlog'
                if (j1<nR) then !otherwise i1+1 arrays remain static and valued as set initially for all i1 (all t).
                    call get_w_qsy(print_cnd, flume_else_basin, isFlood, coefs_wqsy, ss_trg, qsbx, fibc50, psolid, h_bf(i1, j1), &
                        w_c(i1, j1), dt(i1), w(i1, j1), pWreveg, wproy(i1,j1))
                    Is_pr_allD(i1, j1) = sum(Is_prev(j1, :))
!                    print*, '10may20; Is_pr_allD:', Is_pr_allD
                    call get_wc_hsf(print_cnd, wproy(i1,j1), Q(i1, j1), Qf(i1, j1), v(i1, j1), ta, w_v(j1), h(i1, j1), &
                        h_bf(i1, j1), &
                        w_b(i1, j1), &
                        w_c(i1, j1), h_s(i1, j1), h_sf(i1, j1), w_sf, psolid, Is_pr_allD(i1, j1), Qs_fl, dt(i1), dx(j1), &
                        w_b(i1 + 1, j1), &
                        w_c(i1 + 1, j1), h_sf(i1 + 1, j1), dwc_dt, dhsf_dt, pf, dhsf_type, wc_ambit)                        
                    call get_hs(i1,j1,nR,print_cnd, h_bf(i1, j1), diag_un, ta, ca, sa, w_v(j1), dx(j1), psolid, h_s(i1, j1), &
                        h_sf(i1, j1), w_b(i1, j1), w_c(i1, j1), fic50, fif50, fibc50, &
                        Is_pr_allD(i1, j1), Q(i1, j1), Qs, dt(i1), wb_thw(j1), hbf_thw(j1), cmx, dhsf_type, wc_ambit, &
                        dvol_dwc, dvol_dhsf, dvol_dhs, dwc_dt, dhsf_dt, w_c(i1 + 1, j1), & 
                        h_sf(i1 + 1, j1), w_b(i1 + 1, j1), h_s(i1 + 1, j1), net_Qs_flow, dhs_dt, dhsf2_dt, dwb_dt, dwf_dt, &
                        colm(i1,j1), Qspp)
!                    print*, '3may20; just out of hs, inout var wcfut must be wv=1 for trenchExp; wcfut,hbf_fut,colm: ', &
!                        w_c(i1 + 1, j1), h_sf(i1 + 1, j1)-h_s(i1 + 1, j1),colm(i1,j1)

!                    print*, '9jun20; just out of hs; wcfut,h_sf(i1 + 1, j1),h_s(i1 + 1, j1),hbf_fut,colm,hbf_thw: ', &
!                        w_c(i1 + 1, j1), h_sf(i1 + 1, j1), h_s(i1 + 1, j1), h_sf(i1 + 1, j1)-h_s(i1 + 1, j1), &
!                        colm(i1,j1),hbf_thw(j1)
!                    stop
!                    print*, '9jun20; just out of hs; fut geom of wb,wc,hs,hsf: ', &
!                        w_b(i1 + 1, j1), w_c(i1 + 1, j1),  h_s(i1 + 1, j1), h_sf(i1 + 1, j1)

                    pc = 1-pf                          
                    Qs_c = Qs_c + pc*Qspp; Qs_fl = Qs_fl + pf*Qspp; Qs = Qs_c + Qs_fl
                endif
                !channel texture:
                if (j1<nR) then
                    Qso_c_k = Qs_c*pbl !pbl recking16gtm physics deviates when Qs_c grows due to (w+) Qspp from get_hs.
                            !Qspp, which was not needed for flume calibs(check in plot), explains eventual diff gtm that flume.
                            !aunq pa metodo práctico (flume&basin) no debe haber varios combos de params 
                            !(trench-pp13 interpol?).                
                    Qso_f_k = Qs_fl*pblf !if qf=0, Qs_fl is 0, but this calculation is needed to cumulate Is.
                    if ( (colm(i1,j1).lt.1) .or. (colm(i1,j1).gt.3) ) then 
                        !avoids colm1,2,3 that erases (at least)bed surface texture. !p224dr.
                        if (print_cnd) print*, 'to solve CHANNEL texture, with pbl=', pbl
                        Qsi_k = Is_prev(j1, :)/dt(i1)
                        Qsi_c_k = pc*Qsi_k
                        call update_La(print_cnd, dt(i1), kLa, D90c, Lac(i1, j1), Lac(i1 + 1, j1), dLac_dt)
                        call update_text(i1, j1, print_cnd, dt(i1), nD, psolid, h_s(i1, j1), Lac(i1, j1), w_b(i1, j1), dhs_dt, & 
                            dLac_dt, &
                            dwb_dt, dx(j1), &
                            fac(i1, j1, :), fsc(i1, j1, :), Qsi_c_k, Qso_c_k,  &
                            fac(i1+1, j1, :), fsc(i1+1, j1, :))
                        !floodplain texture, if flow on it:
    !                    if (j1 .eq. 8) then
    !                    print*, '14may20; i1:', i1
    !                    print*, '14may20; at j5 (rnd); fac: ', fac(i1,8,:)
    !                    print*, '14may20; at j5 (rnd); fsc: ', fsc(i1,8,:)
    !                    print*, '14may20; at j5 (rnd); faf: ', faf(i1,8,:)
    !                    print*, '14may20; at j5 (rnd); fsf: ', fsf(i1,8,:)                    
    !                    endif                    
                        if (pf .eq. 0) then
                            Laf(i1 + 1, j1) = Laf(i1, j1)
                            faf(i1+1, j1, :) = faf(i1, j1, :)
                            fsf(i1+1, j1, :) = faf(i1, j1, :)
                        else                
                            if (print_cnd) print*, 'to solve FLOODPLAIN texture, with Qf/Q = pf =', pf
                            if (print_cnd) print*, ' pblf=', pblf
                            call update_La(print_cnd, dt(i1), kLa, D90f, Laf(i1, j1), Laf(i1 + 1, j1), dLaf_dt)
                            Qsi_f_k = pf*Qsi_k
                            call update_text(i1, j1, print_cnd, dt(i1), nD, psolid, h_sf(i1, j1), Laf(i1, j1), &
                                w_v(j1) - w_c(i1, j1), &
                                dhsf2_dt, dLaf_dt, dwf_dt, dx(j1), &
                                faf(i1, j1, :), fsf(i1, j1, :), Qsi_f_k, Qso_f_k,  &
                                faf(i1+1, j1, :), fsf(i1+1, j1, :))                        
                        endif
                    else
!                        print*, 'fb: ', fb
!                        stop '20may20; ve q si coge colm asi sea real'
                        call update_La(print_cnd, dt(i1), kLa, D90c, Lac(i1, j1), Lac(i1 + 1, j1), dLac_dt)
                        call update_La(print_cnd, dt(i1), kLa, D90f, Laf(i1, j1), Laf(i1 + 1, j1), dLaf_dt)
                        faf(i1+1, j1, :) = fb;  fsf(i1+1, j1, :) = fb; 
                        fac(i1+1, j1, :) = fb;  fsc(i1+1, j1, :) = fb;
                        if (print_cnd) print*, 'in TEXTcalcs; hs() colm>0; sum(Qso_c_k),sum(Qso_f_k),fac(i1+1, j1, :):', &
                            sum(Qso_c_k),sum(Qso_f_k),fac(i1+1, j1, :)
                    endif                   
                    
		            if (jdown(j1).gt.0) then     ! .and. i.lt.nt
	                    Is(jdown(j1), :) = Is(jdown(j1), :) + (Qso_c_k + Qso_f_k)*dt(i1+1) !2jun20; as for ok Iw fld transit.
!	                    if ( sum(Is(jdown(j1), :)) .ne. sum(Is(jdown(j1), :)) ) then
!	                        print*, '20may20; here NaN Is jdown, so will stop, i1,j1:', i1,j1
!	                        stop
!	                    endif
		            endif
                    fsc_acu = f_acu(.true., nD, fsc(i1, j1, :))		                      
                    if ((i1>1) .and. (flume_else_basin .eq. 0)) then !as i1=1 is filling channel so brings bugs often.
!                        print*, '10jun20; just before check_wat_balance; repeater, fut geom of wb,wc,hs,hsf: ', &
!                            repeater, &
!                            w_b(i1 + 1, j1), w_c(i1 + 1, j1),  h_s(i1 + 1, j1), h_sf(i1 + 1, j1)
                        morp = (/ colm(i1, j1), w_b(i1,j1), w_c(i1,j1), h_s(i1,j1), h_sf(i1,j1), &
                            w_b(i1+1,j1), w_c(i1+1,j1), h_s(i1+1,j1), h_sf(i1+1,j1) /)
                        hc = (/ So(i1, j1), w(i1, j1), h(i1, j1), h_bf(i1, j1), v(i1, j1), vf(i1, j1), Fr(i1, j1) /)
!                        print*, '10jun20; morp:', morp
                        cumWat_m3 = cumWat_m3 + Q(i1,j1)*dt(i1)
!                        if ( (q_dl>.1) .or. lowQ_CH_narrower) then !i.e. in range of recking11 hydraulics database
                            !condition lowQ_CH_narrower to avoid #try>1.

                        !lowQ_CH_narrower deprec; see p256d; increae tol instead.
                        call check_wat_balance(i1,j1, q_dl, tcum, cumWat_m3, cumWatErr_m3, print_cnd, &
                            isFlood, Q(i1, j1), &
                            hc, morp, dt(i1), dx_dt, cpp, Aw, err_m3) !p231u
                        cumWatErr_m3 = cumWatErr_m3 + err_m3

!                            if (lowQ_CH_narrower) then
!                                !turn off indexes:
!                                print*, 'after recalc by CHnarrowing, check if v>vp; i1,j1,v,v_prev, NEW q_dl :', &
!                                    i1,j1,v(i1,j1),v_prev,q_dl
!                                prt_cd = 0; repeater = 0
!                                lowQ_CH_narrower = .false.
!                                ta = ta_orig
!                                v_prev = 1
!                            endif                                
!                        elseif (q_dl<.1 .and. v_prev .ne. 1) then !v_prev .ne. 1 to assure no try narrowing again.
!                            prt_cd = 1 !to repeat calc for this i1,j1 with narrower channel, as per p256u theory.
!                            lowQ_CH_narrower=.true. !10jun20; p256u.
!                            v_prev = v(i1,j1)
!                            deepen_indx= ( w_c(i1,j1) + w_b(i1,j1) ) / ( w_c(i1,j1) + wb_thw(j1) )
!                            hbf_n = h_bf(i1,j1) * deepen_indx
!                            ta_n = 2*hbf_n / ( w_c(i1,j1) - wb_thw(j1) )                            
!                            w_b(i1, j1) = wb_thw(j1);  h_s(i1, j1) =  h_s(i1, j1) - ( hbf_n - h_bf(i1,j1) )
!                                !updated geom at current i1,j1, as this condition will be recalc. 
!                                !Note this new geom remains stored, as after here whatever solution advances i1,j1.
!                            if (deepen_indx < 1) stop 'deepen_indx<1'
!                            print*,'10jun20; lowQ; will dwb- to increase v (p256u). i1,j1,q_dl,ta,ta_n:', &
!                                i1,j1,q_dl,ta,ta_n
!                            print*,' h_bf(i1,j1),hbf_n,w_b(i1, j1),w_c(i1, j1),h_s(i1, j1),h_sf(i1, j1):', &
!                                h_bf(i1,j1),hbf_n,w_b(i1, j1),w_c(i1, j1),h_s(i1, j1),h_sf(i1, j1)
!                            errwc = w_b(i1, j1) + 2*hbf_n/ta_n - w_c(i1, j1)
!                            print*, ' errwc, i.e. overstim new wc from wb; wb+2*hbf_n/ta_n - wc:', errwc
!                            if (errwc>.01) stop                                
!                            ta = ta_n !provisional bank slope for the recalc coming.
!                        endif
                    endif
		            if (prt_cd .eq. 1) then !10jun20; well to reprint hc sln bugs, or to narrow CH to increase v of lowQ.
		                Rw_ppr = Rw_prev(j1) !to solve Aw_NaN bug by repeating cals to print details to avoid NaN R.
		                q_ppr = q_prev(j1)
                    endif                    
                    call check_sed_balance(i1,j1,print_cnd, colm(i1, j1), w_v(j1), ta, dx(j1), net_Qs_flow*dt(i1), &
                    dt(i1), Q(i1, j1), psolid, h_sf(i1, j1), w_c(i1, j1), w_b(i1, j1), h_s(i1, j1), h_sf(i1 + 1, j1), &
                    w_c(i1 + 1, j1), w_b(i1 + 1, j1), h_s(i1 + 1, j1), &
                    dvol(i1, j1))
!                    call check_wat_balance: 7may20: will replace Qcalc 5% tol of previous versions of this file.
                endif
                Is_prev(j1, :) = Is(j1, :)
                c(i1, j1) = Qs / Q(i1, j1)
                if (c(i1, j1) .gt. 1.01*cmx) then
                    print*, 'take this t,x to edit print_range for cmxExcess DEBUG: i,j:', i1,j1
                    print*, 'c(i1, j1),cmx', c(i1, j1),cmx
                    stop
                endif
                
!                if (c(i1, j1) .ne. c(i1, j1)) then !15may20; NaN check for c, as such values damage plotter.
!                    print*, 'take this t,x to edit print_range for NaN_c DEBUG: i,j:', i1,j1
!                    print*, 'Qs,Q', Qs,Q(i1, j1)
!                    stop                    
!                endif
                
                pi(i1, j1, :) = pbl
                call get_GSDpctiles(print_cnd, nD, np, port, fsc_acu, Dfbr, D50c_ss(i1, j1), D84c_ss(i1, j1), D90c_ss, indDp)
                                
                if (repeater .eq. 1) then
                    if (lowQ_CH_narrower) then
                        print*, '10jun20; these low vals are why will repeat calcs with narrower CH: v,q_dl:', v(i1,j1),q_dl
                        stop
                    else
                        stop '13may20; above ij printed, check which was the bug'
                    endif
                end if                    
                if ( (minifld .eq. 1) .or. (prt_cd .eq. 1) ) j1 = j1-1 
                    !note that, for minifld, repeater remains off so runs will continue after this ij repetition.
                if (prt_cd .eq. 1) repeater = 1 !repeat iter to print
                
                if (colm(i1,j1) .eq. 0)  colm0 = colm0+1
                xtcum = xtcum + 1 !10jun20; spatiotemp step
                print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
                print*, '9jun20; colm0,xtcum: ',colm0,xtcum
                if ( (i1 .eq. 4000) .and. (colm0/i1<.1) ) stop 'too infreq colm0'
                print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'                
                
!                if (i1 .eq. 360) then
!                    print*, '25ago20; stop to debug hs j1'
!                    stop
!                endif
                
                
!                if ((i1.eq.10027) .and. (j1.eq.16)) stop '6jun20; check location context to see likely bug cause' !p243cl



                
!                if (j1 .eq. 5) then         
!                    print*, '----------'
!                    print*, '----------i1:', i1
!                    print*, '14may20; j1=5; dt(i1)/60,tcum/3600,Qh(i1, 5)',dt(i1)/60,tcum/3600,Qh(i1, 5)
!                    print*, '14may20; j1=5; Q(i1,5),wc(i1,5),w_v(5),v(i1,5):', &
!                        Q(i1,5),w_c(i1,5),w_v(5),v(i1,5)
!                    print*, '14may20; j1=5; c(i1,5),So((i1,5),D84c(i1,5):', c(i1,5),So(i1,5),D84c(i1,5)
!                    print*, '14may20; j1=5; h,hbf,isFlood:', h(i1,5),h_bf,isFlood
!                    print*, '----------'
!                elseif (j1 .eq. 20) then                                        
!                    print*, '14may20; j1=20; Qh(i1, 20)',Qh(i1, 20)
!                    print*, '14may20; j1=20; Q(i1,20),wc(i1,20),w_v(20),v(i1,20):', &
!                        Q(i1,20),w_c(i1,20),w_v(20),v(i1,20)
!                    print*, '14may20; j1=20; c(i1,20),So((i1,20),D84c(i1,20):', c(i1,20),So(i1,20),D84c(i1,20)
!                    print*, '14may20; j1=20; h,hbf,isFlood:', h(i1,20),h_bf,isFlood
!                endif
                
            end do            

		
			fileunit=12			

			nvar=1
			write(fmt=*, unit= xy(nvar)), path_list(1:num_char,pos_ini_var(nvar))
            xx=trim(output_folder)//'/'//adjustl(xy(nvar))
			call write_FloatArrIJ(fileunit,xx,Qh(i1,:),i1,nR)
			
!			nvar=2
!			write(fmt=*, unit= xy(nvar)), path_list(1:num_char,pos_ini_var(nvar))
!            xx=trim(output_folder)//'/'//adjustl(xy(nvar))
!			call write_FloatArrIJ(fileunit,xx,Q(i1,:),i1,nR)
			
			nvar=2
			write(fmt=*, unit= xy(nvar)), path_list(1:num_char,pos_ini_var(nvar))
            xx=trim(output_folder)//'/'//adjustl(xy(nvar))
			call write_FloatArrIJ(fileunit,xx,wash_sh(i1,:),i1,nR)
			
			nvar=3
			write(fmt=*, unit= xy(nvar)), path_list(1:num_char,pos_ini_var(nvar))
            xx=trim(output_folder)//'/'//adjustl(xy(nvar))
			call write_FloatArrIJ(fileunit,xx,slid_sh(i1,:),i1,nR)			
			
!			
!			nvar=7
!			write(fmt=*, unit= xy(nvar)), path_list(1:num_char,pos_ini_var(nvar))
!            xx=trim(output_folder)//'/'//adjustl(xy(nvar))
!			call write_FloatArrIJ(fileunit,xx,Mgtm(i1,:),i1,nR)			
!			
!			nvar=8
!			write(fmt=*, unit= xy(nvar)), path_list(1:num_char,pos_ini_var(nvar))
!            xx=trim(output_folder)//'/'//adjustl(xy(nvar))
!			call write_FloatArrIJ(fileunit,xx,gam(i1,:),i1,nR)			

			nvar=4
			write(fmt=*, unit= xy(nvar)), path_list(1:num_char,pos_ini_var(nvar))
            xx=trim(output_folder)//'/'//adjustl(xy(nvar))
			call write_FloatArrIJ(fileunit,xx,v(i1,:),i1,nR)
			
			nvar=5
			write(fmt=*, unit= xy(nvar)), path_list(1:num_char,pos_ini_var(nvar))
            xx=trim(output_folder)//'/'//adjustl(xy(nvar))
			call write_FloatArrIJ(fileunit,xx,h(i1,:),i1,nR)			

			nvar=6
			write(fmt=*, unit= xy(nvar)), path_list(1:num_char,pos_ini_var(nvar))
            xx=trim(output_folder)//'/'//adjustl(xy(nvar))
			call write_FloatArrIJ(fileunit,xx,vf(i1,:),i1,nR)			
			
			nvar=7
			write(fmt=*, unit= xy(nvar)), path_list(1:num_char,pos_ini_var(nvar))
            xx=trim(output_folder)//'/'//adjustl(xy(nvar))
			call write_FloatArrIJ(fileunit,xx,w(i1,:),i1,nR)	
			
			nvar=8
			write(fmt=*, unit= xy(nvar)), path_list(1:num_char,pos_ini_var(nvar))
            xx=trim(output_folder)//'/'//adjustl(xy(nvar))
			call write_FloatArrIJ(fileunit,xx,wproy(i1,:),i1,nR)			

			nvar=9
			write(fmt=*, unit= xy(nvar)), path_list(1:num_char,pos_ini_var(nvar))
            xx=trim(output_folder)//'/'//adjustl(xy(nvar))
			call write_FloatArrIJ(fileunit,xx,c(i1,:),i1,nR)								

			nvar=10
			write(fmt=*, unit= xy(nvar)), path_list(1:num_char,pos_ini_var(nvar))
            xx=trim(output_folder)//'/'//adjustl(xy(nvar))
			call write_FloatArrIJ(fileunit,xx,h_s(i1,:),i1,nR)			
			
			nvar=11
			write(fmt=*, unit= xy(nvar)), path_list(1:num_char,pos_ini_var(nvar))
            xx=trim(output_folder)//'/'//adjustl(xy(nvar))
			call write_FloatArrIJ(fileunit,xx,h_sf(i1,:),i1,nR)			
			
			nvar=12
			write(fmt=*, unit= xy(nvar)), path_list(1:num_char,pos_ini_var(nvar))
            xx=trim(output_folder)//'/'//adjustl(xy(nvar))
			call write_FloatArrIJ(fileunit,xx,w_c(i1,:),i1,nR)			
			

			
!			nvar=21
!			write(fmt=*, unit= xy(nvar)), path_list(1:num_char,pos_ini_var(nvar))
!            xx=trim(output_folder)//'/'//adjustl(xy(nvar))
!			call write_FloatArrIJ(fileunit,xx,tsm_c(i1,:),i1,nR)						
			
			nvar=13
			write(fmt=*, unit= xy(nvar)), path_list(1:num_char,pos_ini_var(nvar))
            xx=trim(output_folder)//'/'//adjustl(xy(nvar))
			call write_FloatArrIJ(fileunit,xx,D50c(i1,:),i1,nR)			
			
			nvar=14
			write(fmt=*, unit= xy(nvar)), path_list(1:num_char,pos_ini_var(nvar))
            xx=trim(output_folder)//'/'//adjustl(xy(nvar))
			call write_FloatArrIJ(fileunit,xx,D84c(i1,:),i1,nR)			
			
			nvar=15
			write(fmt=*, unit= xy(nvar)), path_list(1:num_char,pos_ini_var(nvar))
            xx=trim(output_folder)//'/'//adjustl(xy(nvar))
			call write_FloatArrIJ(fileunit,xx,D50c_ss(i1,:),i1,nR)
			
			nvar=16
			write(fmt=*, unit= xy(nvar)), path_list(1:num_char,pos_ini_var(nvar))
            xx=trim(output_folder)//'/'//adjustl(xy(nvar))
			call write_FloatArrIJ(fileunit,xx,D84c_ss(i1,:),i1,nR)
			
			nvar=17
			write(fmt=*, unit= xy(nvar)), path_list(1:num_char,pos_ini_var(nvar))
            xx=trim(output_folder)//'/'//adjustl(xy(nvar))
			call write_FloatArrIJ(fileunit,xx,Qbf(i1,:),i1,nR)
			
			nvar=18
			write(fmt=*, unit= xy(nvar)), path_list(1:num_char,pos_ini_var(nvar))
            xx=trim(output_folder)//'/'//adjustl(xy(nvar))
			call write_FloatArrIJ(fileunit,xx,Qsbf(i1,:),i1,nR)						

			nvar=19
			write(fmt=*, unit= xy(nvar)), path_list(1:num_char,pos_ini_var(nvar))
            xx=trim(output_folder)//'/'//adjustl(xy(nvar))
			call write_FloatArrIJ(fileunit,xx,Q(i1,:),i1,nR)
			
			nvar=20
			write(fmt=*, unit= xy(nvar)), path_list(1:num_char,pos_ini_var(nvar))
            xx=trim(output_folder)//'/'//adjustl(xy(nvar))
			call write_FloatArrIJ(fileunit,xx,So(i1,:),i1,nR)				
			
			nvar=21
			write(fmt=*, unit= xy(nvar)), path_list(1:num_char,pos_ini_var(nvar))
            xx=trim(output_folder)//'/'//adjustl(xy(nvar))
			call write_FloatArrIJ(fileunit,xx,w_b(i1,:),i1,nR)			
			
			nvar=22
			write(fmt=*, unit= xy(nvar)), path_list(1:num_char,pos_ini_var(nvar))
            xx=trim(output_folder)//'/'//adjustl(xy(nvar))
			call write_FloatArrIJ(fileunit,xx,ts_c(i1,:),i1,nR)
			
			nvar=23
			write(fmt=*, unit= xy(nvar)), path_list(1:num_char,pos_ini_var(nvar))
            xx=trim(output_folder)//'/'//adjustl(xy(nvar))
			call write_FloatArrIJ(fileunit,xx,tc_c(i1,:),i1,nR)			

			nvar=24
			write(fmt=*, unit= xy(nvar)), path_list(1:num_char,pos_ini_var(nvar))
            xx=trim(output_folder)//'/'//adjustl(xy(nvar))
			call write_FloatArrIJ(fileunit,xx,Fr(i1,:),i1,nR)			

			nvar=25
			write(fmt=*, unit= xy(nvar)), path_list(1:num_char,pos_ini_var(nvar))
            xx=trim(output_folder)//'/'//adjustl(xy(nvar))
			call write_FloatArrIJ(fileunit,xx,aspect(i1,:),i1,nR)			
			
			nvar=26
			write(fmt=*, unit= xy(nvar)), path_list(1:num_char,pos_ini_var(nvar))
            xx=trim(output_folder)//'/'//adjustl(xy(nvar))
			call write_FloatArrIJ(fileunit,xx,subm(i1,:),i1,nR)			

			nvar=27
			write(fmt=*, unit= xy(nvar)), path_list(1:num_char,pos_ini_var(nvar))
            xx=trim(output_folder)//'/'//adjustl(xy(nvar))
			call write_FloatArrIJ(fileunit,xx,h_bf(i1,:),i1,nR)
			
			nvar=28
			write(fmt=*, unit= xy(nvar)), path_list(1:num_char,pos_ini_var(nvar))
            xx=trim(output_folder)//'/'//adjustl(xy(nvar))
			call write_FloatArrIJ(fileunit,xx,Is_pr_allD(i1,:),i1,nR)
			
			nvar=29
			write(fmt=*, unit= xy(nvar)), path_list(1:num_char,pos_ini_var(nvar))
            xx=trim(output_folder)//'/'//adjustl(xy(nvar))
			call write_FloatArrIJ(fileunit,xx,dvol(i1,:),i1,nR)			
			
			nvar=30
			write(fmt=*, unit= xy(nvar)), path_list(1:num_char,pos_ini_var(nvar))
            xx=trim(output_folder)//'/'//adjustl(xy(nvar))
			call write_FloatArrIJ(fileunit,xx,colm(i1,:),i1,nR)

!        	do k1=1,nD
!                nvar=27
!                write(fmt=*, unit= xy(nvar)), path_list(1:num_char,pos_ini_var(nvar)+k1-1)
!                xx=trim(output_folder)//'/'//adjustl(xy(nvar))
!        	    call write_FloatArrIJ(fileunit,xx,pi(i1,:,k1),i1,nR)
!			end do
        end do
        print*, '25apr20: firstly writing binaries, arrived to last line of subroutine model_sc'
	end subroutine		
end module
