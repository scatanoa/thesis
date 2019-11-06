module mod_model_sc
	implicit none	
	integer i,j,k,kk	!index for loops in 4D arrays: time, space, grainSize and widthPortion
	integer ii,jj,j3,j4,k3,k4	!auxiliar indexes in secondary subroutines
!	type string1D
!		character(len=:), allocatable :: str1D
!	end type string1D
	contains

	character(len=20) function str(k)
	!   "Convert an integer to string."
		integer, intent(in) :: k
		write (str, *) k	
		str = adjustl(str)
	end function str

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

	subroutine get_port(nPctil,port)	
		real, intent(out) :: port(11)
		integer nPctil
		do k3=1,nPctil
			port(k3)=(k3-1.)/10
		end do	
	end subroutine

	subroutine GSDkeyPctiles(P, D, Dp, D84, nD_in)
		!input: points of cumulative distribution
		!ouput: equally spaced grain size percentiles
		!percentil resolution: each .1
		!SUGGESTION: portion finer than Dmin of the sample?
		integer, intent(in) :: nD_in 
		real, intent(in) :: P(nD_in), D(nD_in)
		real, intent(out) :: D84, Dp(11)
		real port(11)

		!interpolation
		call get_port(11,port)
		do k=1,11
			jj=1
			do while (P(jj)<.99*port(k))
			    jj=jj+1
			end do
			if (jj>1) then
			    Dp(k)=log(D(jj-1))+(log(D(jj))-log(D(jj-1)))/(P(jj)-P(jj-1))*(port(k)-P(jj-1))
			    Dp(k)=exp(Dp(k))
			else
			    Dp(k)=(.062e-3*D(jj))**.5  !non cohesive sediment sizes
			endif
		end do
		!key pctiles, useful e.g. for friction formulas    
		D84=(Dp(9)+Dp(10))/2
	end subroutine
	
	subroutine backwater(Sf_min,w_alt,h_jdown_kk1,Head_down,Sf_down,So,Wflow,indW_max,s_elev,Iw_ini,Sw_ini,&
	dx,dt,roug,nW, Sf,Q,w_elev_new,v,h)
        integer, intent(in) :: nW
	    real, intent(in) :: Wflow(nW),s_elev(nW),Iw_ini,Sw_ini,dx,dt,roug(nW),h_jdown_kk1,So,Sf_min
	    real, intent(in) :: Head_down, Sf_down
	    !w_elev[m.a.ref.l]
	    integer, intent(in) :: indW_max
	    real, intent(out) :: Sf,Q,w_elev_new,v(nW),h(nW)
	    real, intent(inout) :: w_alt 
	    real cum, Asx, vm, hm, Fr_m, grad_h, Head_up
	    integer kk, cont_crit, indW
	    Fr_m=2
	    cont_crit=0
	    do while (Fr_m>1 .or. (cont_crit>1 .and. Fr_m<.9))   	    
	    !planned to enter only after any iteration become supercritical. Backwater is initialized with downstream   
	    !slope, which is assumed to be related to subcritical flow.	   
	        cont_crit=cont_crit+1
	        print *, 'cont_crit=', cont_crit
	        cum=0
            kk=1            
            do while (kk.lt.indW_max)
                if (w_alt.lt.s_elev(kk)) then
                    exit
                endif
                kk=kk+1
            end do
            indW=kk            
            Asx=0
            do kk=1,indW
                h(kk)=w_alt-s_elev(kk)
                Asx=Asx+h(kk)*Wflow(kk)
            end do            
            Q=(Iw_ini-(Asx*dx-Sw_ini))/dt            
            vm=Q/Asx
            hm=Asx/sum(Wflow(1:indW))
            Fr_m=vm/(10*hm)**.5
            if (cont_crit>1 .and. Fr_m<.9) then
                h(1)=h(1)/1.5
            else
                h(1)=2*h(1)
            endif
            w_alt=h(1)+s_elev(kk)    	    
	    end do
!        grad_h=(h_jdown_kk1-h(1))/dx    !based on central active channel of crosssection. Sf is same in Xsect
!        Sf=max(Sf_min,So-(1-Fr_m**2)*grad_h)
        print *, 'inside backwater subroutine, s_elev= ', s_elev
        print *, 'indW= ', indW
        print *, 'Q= ', Q
        print *, 'h_assum= ', h
        print *, 'vm= ', vm
        print *, 'hm= ', hm
        print *, 'Fr_m= ', Fr_m
        cum=0
        do kk=1,indW
            cum=cum+Wflow(kk)*h(kk)**2/sqrt(roug(kk))
        end do 
!        Sf=Q**2/(.38**2*10*cum**2)
        Sf=Sf_down
        print *, 'Sf= ', Sf
        do kk=1,indW
            v(kk)=.38*10**.5*h(kk)*Sf**.5
        end do
        print *, 'E= ', h+v**2/20
        print *, 'tao= ', 1e4*h*Sf
        h(1)=h_jdown_kk1-(So-Sf)/(1-Fr_m**2)*dx 
        print *, 'h_jdown= ', h_jdown_kk1
        print *, 'h(kk1)= ', h(1)
!        Head_up=Head_down+Sf*dx
        w_elev_new=h(1)+s_elev(1)
        print *, 'w_elev_new= ', w_elev_new
        do kk=2,indW    !simplification: sameLOW kinetic energy in cross-section, to avoid Xsect gradient in head or h
            h(kk)=w_elev_new-s_elev(kk)
        end do
        print *, 'h_backwatered= ', h
	end subroutine

	subroutine model_sc(iref_ini,iref_fin,jref_ini,jref_fin,&
	output_folder,path_list_py,pos_ini_var_py,flume_else_basin,dt_dQ_coarse,dt_minu_ev,dt_hours_postEv,dt_py,&
	wash_1sh_py,wash_2sh_py,slid_sh_py,Qhill_py, Wfeas_py, X_py, Wflow_py, Dcentr_py, Dpct_py, pctD_py, ini_alluv_py,&
	dum1, dum2, dum3, dum4,dum5, dum6,dum7,nt,n_coarse_dt, nW,nR,nD,nD_hillslope,num_paths,nvars,num_char)  
	!D_py, Qhill_py, dum,
		!inputs and outputs	
		integer, intent(in) :: n_coarse_dt, nt, nW, nR, nD, nD_hillslope 
		integer, intent(in) :: flume_else_basin, num_paths, nvars, num_char, dt_dQ_coarse 
		integer, intent(in) :: iref_ini, iref_fin, jref_ini, jref_fin , pos_ini_var_py(0:nvars-1)
		real, intent(in) :: dt_py(0:nt-1), dt_minu_ev, dt_hours_postEv
		real, intent(in) :: wash_1sh_py(0:n_coarse_dt-1,0:nR-1),wash_2sh_py(0:n_coarse_dt-1,0:nR-1)
        real, intent(in) :: slid_sh_py(0:n_coarse_dt-1,0:nR-1,0:nD_hillslope-1)
		!nD-2 is used given slide(,,nD) is triggered by creep due to incision in this morphodynamic model
		real, intent(in) :: Qhill_py(0:n_coarse_dt-1,0:nR-1),Wfeas_py(0:nR-1, 0:nW-1) 
		real, intent(in) :: X_py(0:nR-1, 0:8-1), Wflow_py(0:nR-1, 0:nW-1), pctD_py(0:nD-1), Dcentr_py(0:nD-1),Dpct_py(0:nD-1)
		real, intent(in) :: ini_alluv_py
!		character*10, intent(in) :: path_list_py(0:num_paths-1)
		character, intent(in), dimension(0:num_char-1, 0:num_paths-1) :: path_list_py
		real, intent(out) :: dum1(nt+1,nR), dum2(nt+1,nR), dum3(nt+1,nR), dum4(nt+1,nR)	!dummy output to debug code
		real, intent(out) :: dum5(nt+1,nR,nW), dum6(nt+1,nR,nW), dum7(nt+1,nR)
		!variables previous to timesteps
		character, dimension(num_char, num_paths) :: path_list
		character*255  :: output_folder
		integer :: pos_ini_var(nvars), nD_hill
		real X(nR,8), Wfeas(nR,nW), Wflow(nR,nW), Dcentr(nD), Dpct(nD), h_prev(nR, nW), iniXsect(nW), Dp(11), pctD(nD)
		real dt(nt)
		real wash_1sh(nt,nR), wash_2sh(nt,nR),slid_sh(nt,nR,nD),Qhill(nt,nR)
		real Sor(nR), dx(nR), zr(nR), thalwSeed, poros, SoMin, D84r, La(nt+1,nR,nW), f_sub_min, Fmin,D84a,D50a,D50ss,D84ss
		real thalwSeedIni(nW), slop_warm
		real So(nt+1,nR),c(nt+1,nR),ck(nt+1,nR,nD),v(nt+1,nR,nW),h(nt+1,nR,nW),vm(nt+1,nR),dm(nt+1,nR,nW),Fr(nt+1,nR,nW)
		real Head(nt+1,nR,nW),rel_submg(nt+1,nR,nW),W_h(nt+1,nR),Qs(nt+1,nR,nD,nW),Q(nt+1,nR)
		real Da84(nt+1,nR,nW),Da50(nt+1,nR,nW),Dss50(nt+1,nR,nW)
		real supl4(nt+1,nR),supl5(nt+1,nR),sumSedBed(nt+1,nR,nW),F(nt+1,nR,nD,nW),f_sub(nt+1,nR,nD,nW),Iw_prev(nR)
		real Is_prev(nR,nD),thratio84(nt+1,nR,nW)					
		real ref, betal, betah, gam0, gam1, gam2, fi0, taofm_taocref, beta	!GTM Recking16  	   
		integer k3,k3ant,kGrav
		integer indW_max(nR), indW_prev(nR), indW(nt+1,nR)
		!variables inside timesteps
		real Iw(nR),Is(nR,nD),jdown(nR),jup(nR),pa(nD,nW),pa_acu(nD,nW),pss(nD,nW),pss_acu(nD,nW), SedBalance(nW),tao(nW)
		real pIs(nD), qsbMx(nW),PblB_acu(nD,nW),pblb(nD,nW),Is_kk(nW), Pbl_acu(nD,nW),pbl(nD,nW),pbl_wc(nD,nW)
		real ff(11, nW),ff_nD(nD, nW), port(11),pa_i1(nD,nW),pblB_i1(nD,nW),pbl_i1(nD,nW),dzdt(nW),dLa_dt(nW),fI(nD,nW),Fs(nW)
		real qsbMx_wc(nW), F_interm(nD, nW), f_sub_interm(nD, nW)
		real sedIn(nW), sedOut(nW), zw  !, zw(nt+1,nR)
		real thm84a(nt,nR,nW),thm84B(nt,nR,nW), th50a(nt,nR,nW), th84a(nW)  
		real z_up, z_down, Sw_prev, fVegRoug(nW), roug(nW), yp, sum_y,sum_yz,sum_yz2,sum_w,sum_wz
		integer flood, i_h_err, j_h_err, kk_h_err, countDF(nR), cont,cont2, fixThalw, contBraidIJ, kk_braid, kk_steep
		real sum_Wext, sumWint, z_offset, dz_trasv, dz_adm, dz_exc, cum_thw, cum_area
		integer contOverLoadF, OverLoadF, contBF_artif,fileunit, contAdverseSlope, contBedrock, contRefreshText_LatReeq
		real aa,bb,cc,zw1,zw2,d1,d2,v1,v2,Fr1,Fr2,yy,zz, hm, WtotSup, thcr84(nW),s5,sumIs,volAct,sumIs_kk,Wstar,side_slp(nW)
		real volSubsup, D84blB, th84blB, M(nD), gam, pmx_mov, cum, zm, z1, z2, Qsb, fIs, Is_kk_k, fI_ok(nD), za
		real b, qsbMx_D(nD), initAlluv, alfaTEscobar, alfaActiveLayer,alfaWC,sedPlus, base_level, alfaGTM
		real t_supply1, t_supply2,t_supply3,t_supply4, Qpk, wash1pk, wash2pk, slidpk, t_pulse, grad_qsmx_y(nW)
		real ris_slop_Q, fall_slop_Q, ris_slop_wash1, fall_slop_wash1, ris_slop_wash2, fall_slop_wash2, Qhill_mean(nR)
		real sumFposit, sumF, ILay(nW), ssLay(nW), fcorrec, gravTpt, discr, So_DF, taorm(nW),taor(nD,nW),alfaVflume
!		character*255 path_Q, path_h_kk1, path_h_kk2,path_h_kk3, path_h_kk4, path_v_kk1, path_v_kk2
!		character*255 path_z_kk1, path_z_kk2,path_z_kk3, path_z_kk4
!		character*255 path_Da50_kk1, path_Da50_kk2,path_Da50_kk3, path_Da50_kk4
!		character*255 path_Dss50_kk1, path_Dss50_kk2,path_Dss50_kk3, path_Dss50_kk4
!		character*255 path_Da84_kk1, path_Da84_kk2,path_Da84_kk3, path_Da84_kk4
!		character*255 path_zw, path_So, path_indW,path_thm84a,path_thm84B,path_thratio84
!		character*255 path_c_k1,path_c_k2,path_c_k3,path_c_k4,path_Qhill,path_La
        character*255 path_c_k1
        character*255 :: path_c(nD) 						
        character*255 :: path_h(nW), path_v(nW), path_z(nW), path_Da50(nW), path_Dss50(nW), path_Da84(nW) 
        character*255 :: path_Q, path_zw, path_So, path_indW, path_thm84a, path_thm84B, path_thratio84,xx,xy(nvars) 
		integer flagSupply, flagDF, flagBackwat, flagWC, tb_hill,tpk_hill,t_calc(nR),flagCourant, trench, nvar
		integer niter_dune_mx, niter_dune_mn, flag_trench, refresh_texture(nW), kk_mx_tooDeep, kk2, sup_lim(nW)
		real tan_repos_ang, zmin, iniPort_D, zbw(nt,nR),latent(nt,nR),tbkwave, celer, t_cum, Hup, Hdown
		real w_elev_ref, Sf_asum, tol_dune, err_dune, tol_trench, dry_repos_ang, tan_stabl_ang(nW), sin_dry_ang
		real :: So_dune,Q_dune,v_dune(nW),h_dune(nW), Head_down, Sf_down, s_elev(nW), w_elev_dune, Sf, kdx
		real :: pulse_feed_rate(nD), feed0s(18), day_of_calc(nt), i_plot, jam_fact(nW), z_pre, ini_elev_repair
		real :: dz_lateral_reequil(nW), z_prev(nW), As_redistr(nD), hs_prev(nD), cum_ini_alluv, rel_err_ini_alluv
		real :: abs_dz_latReeq(nW), p_redistr_mix(nD), p_depos_up(nD), Asx_sed_pos_lat_reeq, Asx_sed_pre_lat_reeq
		real, allocatable :: w_elev(:)
	    integer :: indW_dune, flagDune, ramp(4), flagSo_smooth, flag_day_else_tStep, flag_printDetails, flag_pulse
	    integer :: flag_waterSupplyVar, flag_jamming, cont_debraid, contSupLim
!		character(len=255), dimension(nW) :: path_h, path_v
!		type(string1D) :: path_h(nW), path_v(nW)

		!Flags/parameters related to model type is 1st, do not move following line, so it can be visible:
		flag_waterSupplyVar=1   !flag_pulse=1: deprecated to avoid damage in dt structure.
		flagSupply=0; flagDF=0; flagWC=0; flagSo_smooth=1; 
		!output printing options:
		flag_printDetails=1; flag_day_else_tStep=0;
		flag_trench=1; flag_jamming=0
		
		if (flume_else_basin.eq.1) then
		    nD_hill=nD
		else
		    nD_hill=nD-1    !do not print boulder supply by creep in same binary file
		endif 
		!WC can be ommited when supply and alluvium exhaustion occur, to estimate fractional qsmx by Recking's GTM
		flagBackwat=0; tol_dune=.05; niter_dune_mx=20; niter_dune_mn=5     !use niter_dune_mn>3
        if (allocated(w_elev)) deallocate (w_elev)
        allocate(w_elev(niter_dune_mx))
		flagCourant=0 !1: hydraulic (related to v), 2: morphodynamic (related to vsb~qsb/h),0: binary dt
!        iref_ini=550
!        iref_fin=610
!		jref=10
		dry_repos_ang = 30*3.1416/180
		tan_repos_ang=tan(dry_repos_ang) !It can increase with time after event (coh,veg). f(g) after Kleinhans11
		sin_dry_ang = sin(dry_repos_ang)
        poros=.3	!cui07
        SoMin=5e-4  !calibrate to get partial mobility to coarsen bed
        D84r=.03	!Mar23, 2018: bedrock roughness. May27, 2018: 
        alfaTEscobar=.45 !% weight of active layer respect to bedload txt; it is .3 after Cui07, .45 after MullerHassan18
        alfaActiveLayer=2.    !1.5 was calib with ElguetaHassan2017 for ppt at RFG18. Value=2 was gotten by MullerHassan18.
        alfaWC=1. !1.8 was calib with ElguetaHassan2017 for ppt at RFG18
        alfaGTM = .36
        
        if (flume_else_basin.eq.1) then
            alfaVflume=2.     !2.5 was calib with ElguetaHassan2017 for ppt at RFG18: less disip  in flume
        else
            alfaVflume=1.
        endif
        
        tol_trench = .05    !admisible relative error
        slop_warm = 1.  !if its value is smaller than one, reduces n% the criticality of bank slope

        !parameters for fractional bedload transport 'GTM' Recking2016
        ref=84.	!Santi, take care, if you change this value, th84/th84c must be replaced, as pctil 84 does not hold anymore.
        betal=.05; gam0=.5; gam1=30.; gam2=.25; fi0 = .01   !10jul2019: vasquezTarrio18 leaves constant gam0 and gam1. fi0 can depend on nD.
        taofm_taocref=2.    ![2-4] according to pg4 ReckingGTM16
        
        
		!this new version with dynamic dt has new components, which are identified as 'ddt' (while I learn version code GitHub):
		!FOLLOWING 7 STEPS ARE DEPRECATED FROM APRIL 16TH, 2018
		!1.ddt. downscaling of input series 
		!2.ddt. creation of t_calc(nR) array to control number of timesteps calculated for each reach j. 
		!3.ddt. downscaling of declarated timeseries
		!4.ddt. filling of timeseries if dt>1min and 
		!5.ddt. actualize local t_calc(nR)
		!6.ddt. actualize downstream t_calc(nR) as MINIMUM t_calc from upstream 		
		!7.ddt. accelerator of code run by compressing data in recession. 		
		
		X = X_py
		Wfeas = Wfeas_py
		Wflow = Wflow_py
		Dcentr = Dcentr_py
		Dpct = Dpct_py
		do j=1, nR
    		Qhill_mean(j)=sum(Qhill_py(:,j-1))/nt
    		!Qhill_mean varies spatially to account for gradients in average rainfall
		end do
!		Qhill = Qhill_py
!        print *, 'in f90, dt_py= ', dt_py

		dt=3600*dt_py
		jam_fact=1
		
!		print *, 'in f90, dt= ', dt
	 	!list of paths is really a matrix built from python as chars/paths x nPaths (rows x cols). Element is character
        do i=1,num_char
            do j=1,num_paths
                path_list(i,j)=path_list_py(i-1,j-1)
            end do     
        end do			
		pos_ini_var=pos_ini_var_py
!		wash_1sh=wash_1sh_py 
!		wash_2sh=wash_2sh_py
!		slid_sh=slid_sh_py*0
		Sor= X(:,6)
		dx= X(:,7)
		zr= X(:,8)
!		print *, 'Sor= ', Sor		
!		print *, 'dx= ', dx		
!		print *, 'zr= ', zr
!		print *, 'in f90, dt=', dt
		if (flume_else_basin.eq.0) then
            initAlluv=200.  !.5
       	    thalwSeed=.05  !.05 used for basin, >>h if nW=2 to simulate ~nW=1
       	    !.01(flume SCatano's trench?) 	!Mar23/2018: vertical distance of low flow bed below the rest of valley	            
       	    thalwSeedIni=thalwSeed
        else
            thalwSeed=.001
            if (nW.gt.1) then   !trenched flume
                base_level = .2
                initAlluv= base_level + slop_warm*.5*Wfeas(1,nW)*tan_repos_ang  !UBC weir has h~.2m.
            else !fully 1D flume
                initAlluv=.1    !Marias setup, see her MsC14:
                !https://open.library.ubc.ca/cIRcle/collections/ubctheses/24/items/1.0165847
                !see fig 3.1 for longprofiles
            endif
            thalwSeedIni=Wflow(1,:)/2*tan_repos_ang*slop_warm  !90% of crit, so channel 'warms'
                !assuming Wflow uniformly discretized across valley
        endif    
	
		!1.ddt. downscaling of input series, from day to minute resolution 
		t_calc=1
!		Qhill_mean=sum(Qhill_py)/nt
		tb_hill=6 !default: 6 !hours. Must be related to hillslope area. 6 hours is assumed for 4km2	
		tpk_hill=tb_hill/2.7    !pg 39 Catano2015. SCS triangular Qt (http://www.bdigital.unal.edu.co/50573/1/1020427238_2015.pdf)
!		if (dt(1)<8e4) then !initialize counter for loop coming about downscaling supply when event occurs
!		    k3=24*60+1    !refine to minutes in day
!		else
!		    k3=1+1    !leave coarse daily resolution		    
!		endif
        k3=1    
!        print *, 'flume_else_basin= ', flume_else_basin
!        print *, 'Qhill(0,0)= ', Qhill_py(0,0)
        print *, 'dt_minu_ev= ', dt_minu_ev        
!        if (flume_else_basin.eq.0 .and. flag_waterSupplyVar.eq.0) then  
            !! DEPRECATED SAME DAY I PROPOSED IT: dt structure can not be changed
!            dt=86400
!        endif
        if (flume_else_basin.eq.0) then		   
!        print *, 'dt= ', dt
  		    k3ant=k3
            if (flag_day_else_tStep.eq.1 .and. i.ge.iref_ini-2 .and. i.le.iref_fin+2) then   
            !2 days before and after the ones to be zoomed		        
	            print *, 'alive to refine supply with dt, for i>1, in f90'
	            print *, 'in refining such, day i= ', 0
	            print *, 'k3= ', k3
	            print *, 'dt(k3)= ', dt(k3)
            endif  		    
		    if (dt(k3)<3e3) then !i=1 does not have antecedent to estimate pulse
                do j=1,nR    
		            Qhill(k3:k3+24*60/dt_minu_ev-1,j)=Qhill_py(0,j-1)			
		            wash_1sh(k3:k3+24*60/dt_minu_ev-1,j)=wash_1sh_py(0,j-1)
		            wash_2sh(k3:k3+24*60/dt_minu_ev-1,j)=wash_2sh_py(0,j-1)
		            do k=1,nD-1 !except boulders which are not coming from hillslope hidrology
            		    slid_sh(k3:k3+24*60/dt_minu_ev-1,j,k)=slid_sh_py(0,j-1,k-1)
            		end do		            
		        end do
		        k3=k3+24*60/dt_minu_ev

		    else !must be coarse time steps		        
		        do j=1,nR
	               	Qhill(k3,j)=Qhill_py(0,j-1)			
			        wash_1sh(k3,j)=wash_1sh_py(0,j-1)
			        wash_2sh(k3,j)=wash_2sh_py(0,j-1)
			        do k=1,nD-1
            			slid_sh(k3,j,k)=slid_sh_py(0,j-1,k-1)
			        end do			        
                end do
                k3=k3+1	        		    		    
		    endif
   		    day_of_calc(k3ant:k3-1)=0
		    do i=2,n_coarse_dt	!remember nt is days, that's why some 24 appear in these lines
		        k3ant=k3
                if (flag_day_else_tStep.eq.1 .and. i.ge.iref_ini-2 .and. i.le.iref_fin+2) then   
                !2 days before and after the ones to be zoomed		        
		            print *, 'alive to refine supply with dt, for i>1, in f90'
		            print *, 'in refining such, day i= ', i-1
		            print *, 'k3= ', k3
		            print *, 'dt(k3)= ', dt(k3)
                endif
    !			wash_1sh_py(0:i-1)
    !			wash_2sh_py(0:i-1),slid_1sh_py(0:i-1),slid_2sh_py(0:i-1)
    !			slid_3sh_py(0:i-1),slid_4sh_py(0:i-1),             
		        if (dt(k3)<3e3) then !.and. Qhill_py(i-1)>Qhill_py(i-2)  .and. i>1) then	
		        !intraday water pulse from rain event is assumed
		            t_supply1=k3+tpk_hill*60/dt_minu_ev    !(~arrives to t_peak)
		            t_supply2=k3+tb_hill*60/dt_minu_ev  !(~ends tb)
		            t_supply3=k3+24*60/dt_minu_ev-1		            
       		        do j=1,nR
                        if (i<n_coarse_dt &   !.and. flag_pulse.eq.1   !flag_waterSupplyVar.eq.1 .and.  
                        & .and. Qhill_py(i-2,j-1).lt.Qhill_py(i-1,j-1) .and. Qhill_py(i-1,j-1).gt.Qhill_mean(j)) then  
                            
                            !last condition avoids to build triangular hydrograph where there is not rain event			                    
		                    !Production of hydrograph differs from pg39_Catano2015, given here Q(t+1) is Qbase: suits better for simulQdaily 
		                    t_pulse=0   ![hours]
		                    ii=0
		                    Qpk=Qhill_py(i,j-1) + 48./tb_hill*(Qhill_py(i-1,j-1)-Qhill_py(i,j-1))
		                    wash1pk=wash_1sh_py(i,j-1) + 48./tb_hill*(wash_1sh_py(i-1,j-1)-wash_1sh_py(i,j-1))
		                    wash2pk=wash_2sh_py(i,j-1) + 48./tb_hill*(wash_2sh_py(i-1,j-1)-wash_2sh_py(i,j-1))
		                    ris_slop_Q=(Qpk-Qhill_py(i-2,j-1))/tpk_hill
		                    ris_slop_wash1=(wash1pk-wash_1sh_py(i-2,j-1))/tpk_hill
		                    ris_slop_wash2=(wash2pk-wash_2sh_py(i-2,j-1))/tpk_hill
		                    fall_slop_Q=(Qhill_py(i,j-1)-Qpk)/(tb_hill-tpk_hill)			                    
		                    fall_slop_wash1=(wash_1sh_py(i,j-1)-wash1pk)/(tb_hill-tpk_hill)
		                    fall_slop_wash2=(wash_2sh_py(i,j-1)-wash2pk)/(tb_hill-tpk_hill)
		                    do while (t_pulse.lt.tpk_hill)  !RISING limb in local (subbasin) hydrograph
		                        Qhill(k3+ii,j)=Qhill(k3+ii-1,j)+ris_slop_Q*dt(k3+ii)/3600
		                        wash_1sh(k3+ii,j)=wash_1sh(k3+ii-1,j)+ris_slop_wash1*dt(k3+ii)/3600
		                        wash_2sh(k3+ii,j)=wash_2sh(k3+ii-1,j)+ris_slop_wash2*dt(k3+ii)/3600
		                        ii=ii+1
		                        t_pulse=t_pulse+dt(k3+ii)/3600
		                    end do
		                    do while (t_pulse.lt.tb_hill)   !FALLING limb in local (subbasin) hydrograph
		                        Qhill(k3+ii,j)=Qhill(k3+ii-1,j)+fall_slop_Q*dt(k3+ii)/3600
		                        wash_1sh(k3+ii,j)=wash_1sh(k3+ii-1,j)+fall_slop_wash1*dt(k3+ii)/3600
		                        wash_2sh(k3+ii,j)=wash_2sh(k3+ii-1,j)+fall_slop_wash2*dt(k3+ii)/3600
		                        ii=ii+1
		                        t_pulse=t_pulse+dt(k3+ii)/3600
		                    end do			                    			                    
                            Qhill(k3+ii:t_supply3,j)= Qhill_py(i,j-1)
		                    !Sincronization at minute scale between water and sediment supply is assumed
		                    wash_1sh(k3+ii:t_supply3,j)=wash_1sh_py(i,j-1)
		                    wash_2sh(k3+ii:t_supply3,j)=wash_2sh_py(i,j-1)
		                    !slide is NOT supplied by triangular regime, but constant in day
		                    do k=1,nD-1 
               				    slid_sh(k3:t_supply3,j,k)=slid_sh_py(i-1,j-1,k-1)
            				end do
            		    else
               				Qhill(k3:t_supply3,j)=Qhill_py(i-1,j-1)			
		                    wash_1sh(k3:t_supply3,j)=wash_1sh_py(i-1,j-1)
		                    wash_2sh(k3:t_supply3,j)=wash_2sh_py(i-1,j-1)
		                    do k=1,nD-1
            				    slid_sh(k3:t_supply3,j,k)=slid_sh_py(i-1,j-1,k-1)
                            end do                                                		        
            		    end if				            	                                        				    
			        end do
			        k3=t_supply3+1 !jump to 1st timestep of next day. Next i will review if day has event.
		        elseif (dt(k3)>3e3 .and. dt(k3)<8e4) then !decreasing daily water supply or dt coarse		!flag_waterSupplyVar.eq.1 .and. 
		            t_supply4=k3+24/dt_hours_postEv-1
		            do j=1,nR
           				Qhill(k3:t_supply4,j)=Qhill_py(i-1,j-1)			
			            wash_1sh(k3:t_supply4,j)=wash_1sh_py(i-1,j-1)
			            wash_2sh(k3:t_supply4,j)=wash_2sh_py(i-1,j-1)
			            do k=1,nD-1
        				    slid_sh(k3:t_supply4,j,k)=slid_sh_py(i-1,j-1,k-1)
                        end do    				    
			        end do
    	            k3=t_supply4+1
		        else    !go here for the whole t, if flag_waterSupplyVar=0 even
		            do j=1,nR
           				Qhill(k3,j)=Qhill_py(i-1,j-1)			
			            wash_1sh(k3,j)=wash_1sh_py(i-1,j-1)
			            wash_2sh(k3,j)=wash_2sh_py(i-1,j-1)
			            do k=1,nD-1
        				    slid_sh(k3,j,k)=slid_sh_py(i-1,j-1,k-1)
			            end do			            
                    end do		
                    k3=k3+1		            
		        endif
		        day_of_calc(k3ant:k3-1)=i-1 
		    end do
		else    !flume. Experim duration is less than 1e4min so you don't need to coarsen inefficient time periods of process
		    do i=1,n_coarse_dt    	!dt_dQ_coarse is 24 in basin from daily lumped SHIA hydrologic model, or Maria's experiment
		        k3ant=k3
	            do j=1,nR	        
                    Qhill(k3:k3+dt_dQ_coarse*60/dt_minu_ev-1,j)=Qhill_py(i-1,j-1)
                    wash_1sh(k3:k3+dt_dQ_coarse*60/dt_minu_ev-1,j)=wash_1sh_py(i-1,j-1)			    		    
                    wash_2sh(k3:k3+dt_dQ_coarse*60/dt_minu_ev-1,j)=wash_2sh_py(i-1,j-1)			    		   
                    do k=1,nD
    				    slid_sh(k3:k3+dt_dQ_coarse*60/dt_minu_ev-1,j,k)=slid_sh_py(i-1,j-1,k-1)
	                end do			    		            		    
        		end do
!        		print *, 'Qhill(k3:k3+dt_dQ_coarse*60/dt_minu_ev-1,j)= ', Qhill(k3:k3+dt_dQ_coarse*60/dt_minu_ev-1,1)
        		k3=k3+dt_dQ_coarse*60/dt_minu_ev
        		day_of_calc(k3ant:k3-1)=i-1
    		end do		
	 	endif	 
	 	!Qramp starting flume to avoid exhaustion wave. 
	 	!(26MAY: it is not too relevant to avoid exhaustion)
	 	!Please make sure Qhill temporally downscaled has enough time steps to do ramp: 
	 	print *, 'alive after refining supply with dt, in f90'
	 	if (flume_else_basin.eq.1) then
	     	ramp=(/ 1,2,3,4 /)  ![minutes from start time of experiment]
	     	ramp=ramp/dt_minu_ev    !time steps where such minutes are located in timeseries    
	     	print *, 'ramp= ', ramp
            Qhill(1:ramp(1),1)=.1*Qhill(1:ramp(1),1)    !ramp is applied only to initial reach
            Qhill(ramp(1)+1:ramp(2),1)=.3*Qhill(ramp(1)+1:ramp(2),1)
            Qhill(ramp(2)+1:ramp(3),1)=.5*Qhill(ramp(2)+1:ramp(3),1)
            Qhill(ramp(3)+1:ramp(4),1)=.8*Qhill(ramp(3)+1:ramp(4),1)
        endif    
!        print *, 'in j=1 Qhill= ', Qhill(:,1)
        !Maria's feed pulses (ElguetaHassan2017):
        pctD=pctD_py; pctD=.01*pctD
        print *, 'in f90, pctD= ', pctD
        pulse_feed_rate(1)=83.3e-6/2.65*pctD(1)
        do k=2,nD
            pulse_feed_rate(k)=83.3e-6/2.65*(pctD(k)-pctD(k-1))
        end do
            ![hours from start time of experiment]:
        feed0s=(/0.,40., 81.,120., 120.25,130., 130.25,140., 140.25,150.,&
        150.25,160., 160.5,180., 180.5,200., 240.,280./)         
        feed0s=feed0s*60/dt_minu_ev  !time steps where such minutes are located in timeseries    
        if (flume_else_basin.eq.1 .and. flagSupply.eq.1) then !slid_sh comes with constant feed by default
            do i=1,9
                slid_sh(feed0s(2*i-1)+1:feed0s(2*i),1,:)=0
                do k=1,nD    
                    if (i.ge.2 .and. i.le.7) slid_sh(feed0s(2*i)+1:feed0s(2*i+1),1,k)=pulse_feed_rate(k)   
                end do
            end do
            do i=80*60/dt_minu_ev+1,feed0s(3)
                slid_sh(i,1,:)=pulse_feed_rate
            end do
       	 	!check supplied sediment [kg]:
       	 	cum=0
	     	do i=1,nt
	     	    cum=cum+sum(slid_sh(i,1,:))*dt(i)
	     	end do
            cum=cum*2650    ![kg]
            print *, 'feed [kg]= ', cum
        endif
        slid_sh=flagSupply*slid_sh
        if (flag_waterSupplyVar.eq.0) then
            do j=1,nR
                Qhill(:,j)=sum(Qhill_py(0:n_coarse_dt-1,j-1))/n_coarse_dt
            end do	 	
        endif            
		!PREVIOUS TO TIME STEPS:
!		path_La='doct_fictBasin/La'
!		path_Qhill='doct_fictBasin/Q_hill'
!		path_Q='doct_fictBasin/Q'
!		path_zw='doct_fictBasin/zw'
!		path_So='doct_fictBasin/So'
!		path_indW='doct_fictBasin/indW'
!		path_thm84a='doct_fictBasin/thm84a'
!		path_thm84B='doct_fictBasin/thm84B'
!		path_thratio84='doct_fictBasin/thratio84'
        
!        path_c_k1='doct_fictBasin/c_k1'
!        do k=1,nD
!!    		path_c(k)=path_list(pos_ini_var(13)+k-1)    !13th variable to plot as stated in dictionary in python
!!		path_c_k2='doct_fictBasin/c_k2'
!!		path_c_k3='doct_fictBasin/c_k3'
!!		path_c_k4='doct_fictBasin/c_k4'		
!        end do
!        
!        do kk=1,nW
!!    		path_v(kk)=path_list(pos_ini_var(3)+kk-1)
!!    		path_h(kk)= output_folder // '/' // path_h(kk)
!!		path_h_kk2='doct_fictBasin/h_kk2'
!!		path_h_kk3='doct_fictBasin/h_kk3'
!!		path_h_kk4='doct_fictBasin/h_kk4'
!!		path_v_kk1='doct_fictBasin/v_kk1'
!!		path_v_kk2='doct_fictBasin/v_kk2'
!!		path_Da84_kk1='doct_fictBasin/Da84_kk1'
!!		path_Da84_kk2='doct_fictBasin/Da84_kk2'
!!		path_Da84_kk3='doct_fictBasin/Da84_kk3'
!!		path_Da84_kk4='doct_fictBasin/Da84_kk4'
!!		path_Da50_kk1='doct_fictBasin/Da50_kk1'
!!		path_Da50_kk2='doct_fictBasin/Da50_kk2'
!!		path_Da50_kk3='doct_fictBasin/Da50_kk3'
!!		path_Da50_kk4='doct_fictBasin/Da50_kk4'
!!		path_Dss50_kk1='doct_fictBasin/Dss50_kk1'
!!		path_Dss50_kk2='doct_fictBasin/Dss50_kk2'
!!		path_Dss50_kk3='doct_fictBasin/Dss50_kk3'
!!		path_Dss50_kk4='doct_fictBasin/Dss50_kk4'		
!!		path_z_kk1='doct_fictBasin/z_kk1'
!!		path_z_kk2='doct_fictBasin/z_kk2'
!!		path_z_kk3='doct_fictBasin/z_kk3'
!!		path_z_kk4='doct_fictBasin/z_kk4'
!        end do 
		
		
		
!		do kk=1,nW
!			path_h(nW)="doct_fictBasin/h"//str(kk)
!!			path_h(nW)=trim(path_h(nW))
!			path_v(nW)="doct_fictBasin/v"//str(kk)
!!			path_v(nW)=trim(path_v(nW))
!		end do 
!		path_h='doct_fictBasin/h'
!		path_v='doct_fictBasin/v'
		contAdverseSlope=0
		contBF_artif=0
		contBraidIJ=0
		contBedrock = 0
		contRefreshText_LatReeq = 0
		contSupLim = 0
		contOverLoadF=0
		countDF=0
		indW_max=0
		indW_prev=0
		iniXsect=0
		h_prev=0
		So=0;c=0;vm=0;dm=0;W_h=0;Q=0;supl4=0;supl5=0;indW=0;	!i,j
		v=0;h=0;Fr=0;Head=0;rel_submg=0;Da84=0;Da50=0;Dss50=0;sumSedBed=0;thratio84=0	!i,j,kk
		Qs=0; F=0; f_sub=0	!i,j,kk,k
		ck=0 !i,j,k
		Iw_prev=0 
		Is_prev=0
		zw=0
		La=0
		zbw=0
		latent=0
		w_elev=0
		do j=1,nR
			kk=1
			do while (kk<=nW .and. Wflow(j,kk)>0)         
				kk=kk+1
				indW_max(j)=indW_max(j)+1       
			end do
		end do
		print *, 'indW_max= ', indW_max					
		do kk=1,nW
!			k3=nW-kk+1
			iniXsect(kk)=initAlluv-sum(thalwSeedIni(kk:nW))	!the more external width portion the thicker initial alluvium
		end do	
!		print *, 'iniXsect= ', iniXsect
		do j=1,nR
			sumSedBed(1,j,:)=iniXsect
!		    print *, 'j= ', j
!    		print *, 'initial sumSedBed= ', sumSedBed(1,j,:)			
		end do
		
		!check if initial volume of alluvium which forms initial trench is coherent (jul10, 2019)
		cum_ini_alluv = 0
		do j=1,nR
		    cum_ini_alluv = cum_ini_alluv + sum(sumSedBed(1,j,:)*Wflow(j,:))
		end do
		print *, 'thalwSeedIni, all kk: ', thalwSeedIni
		print *, 'iniXsect: ', iniXsect
        print *, 'cum_ini_alluv [m3]', cum_ini_alluv        
        print *, 'ini_alluv_py [m3]', ini_alluv_py
        rel_err_ini_alluv = (cum_ini_alluv - ini_alluv_py) / ini_alluv_py
        print *, 'rel_err cum_ini_alluv', rel_err_ini_alluv
        if (abs(rel_err_ini_alluv) .gt. .1) STOP 'abs(rel_err_ini_alluv) exceeds 10%'    
        ini_elev_repair = (ini_alluv_py - cum_ini_alluv) / (sum(dx)*sum(Wfeas(:, nW))/nR)
        print *, 'z ini thw, PRE py correx: ', sumSedBed(1,:,1)        
        sumSedBed(1,:,:) = sumSedBed(1,:,:) + ini_elev_repair
        print *, 'z ini thw, POS py correx: ', sumSedBed(1,:,1)
        base_level = base_level + ini_elev_repair 
        print *, 'elevation correction: ', ini_elev_repair
        print *, 'base_level: ', base_level                  
		cum_ini_alluv = 0
		do j=1,nR
		    cum_ini_alluv = cum_ini_alluv + sum(sumSedBed(1,j,:)*Wflow(j,:))
		end do
        print *, 'drainage L:', sum(dx)
        print *, 'mean W valley:', sum(Wfeas(:, nW))/nR     
        print *, 'repaired cum_ini_alluv [m3]', cum_ini_alluv
        rel_err_ini_alluv = (cum_ini_alluv - ini_alluv_py) / ini_alluv_py
        print *, 'repaired rel_err cum_ini_alluv', rel_err_ini_alluv      

		iniXsect=iniXsect-sum(iniXsect)/nW !rescale to leave it as template when flow attempt to create braided river
		!initial portions of each grain size class in alluvium (F: active layer, f: subsurface layer)
!		if (flume_else_basin.eq.0) then
!		    F(1,:,1:2,:)=1./(2*nD)	!half of uniform frequency is assumed to initialize a coarse alluvium 
!		    F(1,:,3:nD,:)=(1-(F(1,1,1,1)+F(1,1,2,1)))/(nD-3+1)
!		    f_sub(1,:,:,:)=1./nD
!        else
        iniPort_D=.01*pctD_py(0)    
        F(1,:,1,:)=iniPort_D
        f_sub(1,:,1,:)=iniPort_D
        do k=2,nD 
            iniPort_D=.01*(pctD_py(k-1)-pctD_py(k-1-1))
!                print *, 'iniPort_D', iniPort_D
            F(1,:,k,:)=iniPort_D
            f_sub(1,:,k,:)=iniPort_D  !initial alluvium is vertically mixed 
        end do    
!!      	    f_sub(1,:,1,:)=.1    !active layer starts well mixed with subsurface
!!       	    f_sub(1,:,4,:)=.1
!!       	    f_sub(1,:,2,:)=.4
!!       	    f_sub(1,:,3,:)=.4
!        endif
        print *, 'pctD_py= ', pctD_py 		    
!		La=20*Sor  !10*(Sor./.037/2).^(1/.68)
		i=0	!contSand
		do k=1,nD
			if (Dcentr(k)<=2e-3) then
				i=i+1
			endif
		end do
		kGrav=i+1
		print *, 'kGrav= ', kGrav 
		Fmin=.01
		f_sub_min=.25	!Skeleton, Takahashi
!		jup(1)=0
		jup=0
		do j=1,nR
		    if (X(j,4)>0) then  !if it is not basin outlet
    		    jup(int(X(j,4)))=j    
    		else
    		    jup(j)=nR-1
    		endif
		end do
		jdown=X(:,4)
		print *, 'jup= ', jup
		print *, 'jdown= ', jdown
!		day_of_calc=1
		!TIME STEPS *******************************************************************************************************
		do i=1,nt
!		    print *, 't_run [%]', int(float(i)/nt*100)			
			Iw=0; Is=0
!			jdown=0
			j=1
!			backWat=0	!declarate if you want to use it
!            print *, 'i= ', i
		    if (flag_day_else_tStep.eq.0) then
			    i_plot=i
			else
			    i_plot=day_of_calc(i)
		    endif            
		    
			do while (j .le. nR)				
!			    print *, 'i', i, 'j', j	
                if (j .eq. 1) print *, 'time step i: ', i
				!a) BED TEXTURE -------------------------------------------------------------------------------------
!				print *, 'j= ', j
				pa=0; pa_acu=0; pss=0; pss_acu=0
				do kk=1,indW_max(j)  !verify in each portion of valley             
				    if (sumSedBed(i,j,kk) .gt. 0) then  !is there alluvium?
				        do k=1,nD 
				            pa(k,kk)=F(i,j,k,kk)
				            pa_acu(k,kk)=sum(pa(1:k,kk))
				        end do
						call GSDkeyPctiles(pa_acu(:,kk), Dpct, Dp, D84a, nD)	!Dp requested to plot GSD to check
                		if (sum(pa_acu(:,kk) - pa_acu(:,kk)) .ne. 0) then
			                print *, 'err i= ', i
			                print *, 'j= ', j
			                print *, 'kk= ', kk                		    
                		    STOP 'NaN pa_acu'   !NaN identifier						
                		endif
						D50a=Dp(6)
						La(i,j,kk)=alfaActiveLayer*Dp(10)		!after discussion with Marwan. Apr 10, 2018			
				    else
				        D84a=D84r
				        D50a=D84r
				    endif
				    
!        		    if (j.le.jref_fin .and. j.ge.jref_ini .and. i_plot.le.iref_fin .and. i_plot.ge.iref_ini &
!        		    & .and. flag_printDetails.eq.1) then
!	                	print *, 'kk', kk
!	                	print *, 'f_sub(i,j,:,kk): ', f_sub(i,j,:,kk)
!	                endif    
				    
					Da84(i,j,kk)=D84a
					Da50(i,j,kk)=D50a
					if (sumSedBed(i,j,kk) .gt. La(i,j,kk)) then  !is there subsurface layer?
						do k=1,nD 
				            pss(k,kk)=f_sub(i,j,k,kk)
				            pss_acu(k,kk)=sum(pss(1:k,kk))	!plot GSD to check
				        end do
				        call GSDkeyPctiles(pss_acu(:,kk), Dpct, Dp, D84ss, nD)
                		if (sum(pss_acu(:,kk) - pss_acu(:,kk)) .ne. 0) then
			                print *, 'err i= ', i
			                print *, 'j= ', j
			                print *, 'kk= ', kk                		    
                		    STOP 'NaN pss_acu'   !NaN identifier						
                		endif				        
						D50ss=Dp(6)
					else
						D50ss=1e-10 !practically 0, but plotable in logscale                
					endif
					Dss50(i,j,kk)=D50ss
				end do			
	            
				if (i.eq.1 .and. j.eq.1) then
					pa_i1=pa
				endif
				!b) WATER DISCHARGE FROM BALANCE(velez01)------------------------------------------------------------------
				!cumulative Iw[m3] is coherent. X topology is ranked by w and upstr to downstream, given w
			    indW(i,j)=indW_prev(j)                 
			    Iw_prev(j)=Iw_prev(j)+Qhill(i,j)*dt(i)  !water supply from local hillslopes is added to what cAme from upstr
				!hypothesis: channel SLOPE determined by thalweg
!			    z_up=(sumSedBed(i,j,1)*Wflow(j,1)+sumSedBed(i,j,2)*Wflow(j,2))/(Wflow(j,1)+Wflow(j,2))                
                if (flagSo_smooth.eq.1) then
                    !spatially smoothed (central difference?) as Parker ebook, ch17, pg11
                    if (jdown(j).gt.0 .and. jup(j).gt.0) then  
                    !intermediate nodes, with neighbors both up and downstream    
                        z_up=sumSedBed(i,jup(j),1)
                        z_down=sumSedBed(i,jdown(j),1)    
    			        kdx=2   
    			    elseif (jdown(j).eq.0) then  !final node= basin outlet
    			        z_up=sumSedBed(i,j,1)     !jup(jup(j))
    			        z_down = base_level  !weir?
    			        kdx=1
    			    else   !initial node= headwater			        
			            z_up=sumSedBed(i,j,1)
			            z_down=sumSedBed(i,jdown(j),1)
			            kdx=1 !2
    			    endif
                else
                    z_up=sumSedBed(i,j,1)
                    kdx=1
                    if (j.lt.nR) then		                        
		                z_down = sumSedBed(i,jdown(j),1)
		            else
		                z_down = base_level     !initAlluv-(nW-1)*thalwSeed   !weir?
		            endif 
                endif 
			    So(i,j)=Sor(j)+(z_up-z_down)/(kdx*dx(j))			    
			    kk=1
			    Sw_prev=0			    
			    do while (kk<=indW(i,j))	!JIV mass balance equation must be solved lumped in complete wetted cross section
			        Sw_prev=Sw_prev+h_prev(j,kk)*Wflow(j,kk)*dx(j)
			        kk=kk+1
			    end do
    		    if (j.le.jref_fin .and. j.ge.jref_ini .and. i_plot.le.iref_fin .and. i_plot.ge.iref_ini &
    		    & .and. flag_printDetails.eq.1) then
    		        print *, 'i:', i
	            	print *, 'dt(i)= ', dt(i)	            	
	            	print *, 'days_calculated', day_of_calc(i)
                	print *, 'node j= ', j
   			        print *, 'zbw(i,j)= ', zbw(i,j)
   			        print *, 'latent(i,j)= ', latent(i,j)
   			        print *, 'zr(j)= ', zr(j)
                	print *, 'z_up= ', z_up
                	print *, 'z_down= ', z_down
                	print *, 'kdx: ', kdx
                	print *, 'So= ', So(i,j)
                	print *, 'La= ', La(i,j,:)  
	            endif
   			    flood=1 
			    indW(i,j)=1	!Gurnel12 dynamic model for vegetation is yet to be explored.
			    fVegRoug=0	!compatible with Paul's (Mr New Mexico) daughter?. She defended dissertation by March2018
			    roug=0
			    if (flagBackwat.eq.0 .and. So(i,j).le.SoMin) So(i,j)=SoMin
			    !May26, 2018. final factor 'alfaVflume' in yp traslates Recking11 regression. Flume implies less friction.
			    yp=1/.39*10**.5*So(i,j)**.5*alfaVflume
			    flagDune=0			    
			    !friction law of Recking2011, v~f(d/D), ~consistent dimensions
			    do while (flood .eq. 1) !water level over bank of active channel
			        do kk=1,indW(i,j)
			        	fVegRoug(kk)=1
			            if (sumSedBed(i,j,kk)<D84r) then	!some bedrock in roughness?
			                roug(kk)=fVegRoug(kk)*max(Da84(i,j,kk),D84r-sumSedBed(i,j,kk))   
			                !some fill but not smoother than texture of fill
			            else
			                roug(kk)=fVegRoug(kk)*Da84(i,j,kk)                    
			            endif
			        end do    
    			    if (flagBackwat.eq.1 .and. So(i,j).le.0) then
			            contAdverseSlope=contAdverseSlope+1
			            zw1=sumSedBed(i,j,kk)   !make hydraulic depth lowest non-negat, so backwater level wins
			            
			            !solve hydraulics because adverse slope might not implie backwaver in enough energy upstream
			            !consider Sf>0 <> So<0. Sf is assumed and iterated.			            
                         
			            Sf_asum=So(i-1,jdown(j))    !.5 !sum(So(i-1,:))/nR    
			            !downstream gradient suggested by Cui's code
			            k4=1
			            print *, 'ADVERSE SLOPE SO DUNE'
			            print *, 'Sf_asum= ', Sf_asum
			            print *, 'zr(jdown(j))= ', zr(jdown(j))
			            print *, 'zs_downstr= ', sumSedBed(i,jdown(j),1)
			            print *,  'decoupled -lagged- h_downstr', h(i-1,jdown(j),1)
			            !lacking set downstream control e.g. weir
			            Head_down=zr(jdown(j))+sumSedBed(i,jdown(j),1)+h(i-1,jdown(j),1)
			            w_elev(k4)=Head_down+Sf_asum*dx(j)			            
			            Head_down=Head_down+v(i-1,jdown(j),1)**2/20
			            s_elev=sumSedBed(i,j,:)+zr(j)
			            print *, 's_elev= ', s_elev
			            print *, '1stly printed w_elev= ', w_elev
			            print *, 'Head_down= ', Head_down
!			            print *, 'initial w_elev= ', w_elev(k4)
!			            indW_dune=indW(i,j)
			            err_dune=1.   !initialize high enough to enter to loop
			            do while (k4.le.niter_dune_mn .or. err_dune.ge.tol_dune) 
			            !.or. Sf_dune_converg(k4)>Sf_dune_converg(k4-1))
	                        print *, '.................next iter to improve h' 
                            k4=k4+1    			            
			                if (k4.gt.niter_dune_mx) then
			                    exit
			                endif
                            print *, 'k4= ', k4
                            print *, 'initial w_elev= ', w_elev(k4-1)
                            w_elev_ref=w_elev(k4-1)
                            call backwater(SoMin, w_elev_ref,h(i-1,jdown(j),1),Head_down,Sf_asum,So(i,j),&
                            Wflow(j,:),indW_max(j),s_elev,Iw_prev(j),Sw_prev,dx(j),dt(i),roug,nW, Sf,&
                            Q_dune,w_elev_dune,v_dune,h_dune)   
                            !given weird alternating convergence was founded, to Sf~1e-3 and Sf~1e-1:
!                            Sf_asum=So_dune
!                            Sf_dune_converg(k4)=Sf_asum
                            w_elev(k4)=w_elev_dune
       			            print *, 'new w_elev added= ', w_elev
                            if (k4>2) err_dune=abs(w_elev(k4)-w_elev(k4-2))/w_elev(k4)
                            print *, 'err_dune= ', err_dune                            
                        end do
                        Q(i,j)=Q_dune
                        v(i,j,:)=v_dune
                        h(i,j,:)=h_dune
                        So(i,j)=Sf
!                        So(i,j)=minval(Sf_dune_converg(k4-3:k4-1))
                        flagDune=1
			            exit    !avoid loop to get uniform flow (Sf=So)
    			    endif
			        sum_y=0; sum_yz=0; sum_yz2=0; sum_w=0; sum_wz=0; zw=0            
			        do kk=1,indW(i,j)
!	!		            roug(kk)=fVegRoug(kk)*Da84(i,j,kk);		            
			            yy=Wflow(j,kk)/roug(kk)**.5
			            zz=sumSedBed(i,j,kk)
			            sum_y=sum_y+yy
			            sum_yz=sum_yz+yy*zz
			            sum_yz2=sum_yz2+yy*zz**2
			            sum_w=Wfeas(j,kk)
			            sum_wz=sum_wz+Wflow(j,kk)*zz
			        end do
			        aa=dt(i)*yp*sum_y
			        bb=-2*dt(i)*yp*sum_yz+sum_w*dx(j)
			        cc=dt(i)*yp*sum_yz2-(Iw_prev(j)+Sw_prev)-sum_wz*dx(j)
			        discr=bb**2-4*aa*cc
			        zw1=(-bb+sqrt(discr))/(2*aa)
			        zw2=(-bb-sqrt(discr))/(2*aa)
			        d1=zw1-sumSedBed(i,j,1)
			        d2=zw2-sumSedBed(i,j,1)
			        v1=yp/roug(1)**.5*d1
			        v2=yp/roug(1)**.5*d2
			        Fr1=v1/sqrt(10*d1)
			        Fr2=v2/sqrt(10*d2)
			        if (discr.lt.0 .or. zw1.le.sumSedBed(i,j,indW(i,j))) then !water level in wider channel would be lower than bank
			        	contBF_artif=contBF_artif+1
			        	zw1=sumSedBed(i,j,indW(i,j))	!do not follow exactly friction law 			        	
			        	flood=0
			        elseif (indW(i,j).eq.nW .or. Wflow(j,indW(i,j)+1).eq.0 .or. zw1.le.sumSedBed(i,j,indW(i,j)+1)) then
			        	
			            flood=0
			        elseif (nW>1) then
			            indW(i,j)=indW(i,j)+1
			        endif
			    end do	            			    			    	            			    
!	!		    Q(i,j)=0;
!				zw(i,j)=zw1

!	            if (zw(i,j)>1e6) then
!	            	stop 'review why zw(i,j)>1e6'	           
!	            endif
!	            if (zw(i,j).eq.0 .and. i.le.nt) then
!	            	stop 'review why zw(i,j)=0'	           
!	            endif	            			    
                !no backwater profile is calculated. Binary process: uniform flow or baselevel control
                if (flagBackwat.eq.1 .and. int(latent(i,j))*zbw(i,j).gt.zr(j)+zw1) then   !baselevel wins
                	k4=1
                	err_dune=1.
                	do while (k4.le.niter_dune_mn .or. err_dune.ge.tol_dune) 
		            !.or. Sf_dune_converg(k4)>Sf_dune_converg(k4-1))
		                print *, 'solving h due to wave traveling upstream'
                        print *, '.................next iter to improve h' 
                        k4=k4+1    			            
		                if (k4.gt.niter_dune_mx) then
		                    exit
		                endif
                        print *, 'k4= ', k4
                        print *, 'initial w_elev= ', w_elev(k4-1)
                        w_elev_ref=w_elev(k4-1)
                        call backwater(SoMin, w_elev_ref,h(i-1,jdown(j),1),Head_down,SoMin,So(i,j),&
                        Wflow(j,:),indW_max(j),s_elev,Iw_prev(j),Sw_prev,dx(j),dt(i),roug,nW, Sf,&
                        Q_dune,w_elev_dune,v_dune,h_dune)   
                        !given weird alternating convergence was founded, to Sf~1e-3 and Sf~1e-1:
!                            Sf_asum=So_dune
!                            Sf_dune_converg(k4)=Sf_asum
                        w_elev(k4)=w_elev_dune
   			            print *, 'new w_elev added= ', w_elev
                        if (k4>2) err_dune=abs(w_elev(k4)-w_elev(k4-2))/w_elev(k4)
                        print *, 'err_dune= ', err_dune                            
                    end do                                    
                elseif (flagDune.eq.0) then    !otherwise h, v, S=Sf and Q have been already calculated                
			        do kk=1,indW(i,j)
			            h(i,j,kk)=zw1-sumSedBed(i,j,kk)			        
    			        v(i,j,kk)=yp/roug(kk)**.5*h(i,j,kk)
			            Q(i,j)=Q(i,j)+Wflow(j,kk)*h(i,j,kk)*v(i,j,kk)
			        end do
			    endif
			    indW_prev(j)=indW(i,j)
			    h_prev(j,:)=h(i,j,:)
			    if (jdown(j).gt.0 .and. i.lt.nt) then
  			        Iw(jdown(j))=Iw(jdown(j))+Q(i,j)*dt(i+1)  !add water to downstream reach. 
  			        !Iw enters there in water balance in next t, then volume must be cumulated with related dt(t+1)        
			    endif

			    WtotSup=sum(Wflow(j,1:indW(i,j)))
			    hm=0
			    do kk=1,indW(i,j)
			        Fr(i,j,kk)=v(i,j,kk)/sqrt(10*h(i,j,kk))	 !check1
			        rel_submg(i,j,kk)=h(i,j,kk)/roug(kk)    !check3
			        hm=hm+h(i,j,kk)*Wflow(j,kk)/WtotSup
			    end do        
			    W_h(i,j)=WtotSup/hm    !check4

			    !jun4, 2018: to make sure just enough width is used			    
                if (nW>1) then
    			    indW(i,j)=1                    
                    kk=1
			        do while (kk.lt.nW .and. h(i,j,kk+1).gt.1e-3)
			            kk=kk+1
			        end do
		            indW(i,j)=kk			        
			    endif
			    
			    Head(i,j,:)=zr(j)+sumSedBed(i,j,:)+h(i,j,:)+v(i,j,:)**2/20    !based on alluvium after work
			    if (j.le.jref_fin .and. j.ge.jref_ini .and. i_plot.le.iref_fin .and. i_plot.ge.iref_ini &
			    & .and. flag_printDetails.eq.1) then
	            	print *, 'indW= ', indW(i,j) 
	            	print *, 'roug= ', roug
	            	print *, 'Qhill_day [m3/s]= ', Qhill_py(day_of_calc(i),j-1) !python starts time index at 0 but day_of_calc also does
	            	print *, 'Qhill(i,j) [m3/s]= ', Qhill(i,j)
	            	print *, 'Qhill_sed(i,j,:) [g/s]= ', slid_sh(i,j,:)*2.65e6 
	            	print *, 'Iw_prev(j)= ', Iw_prev(j)
	            	print *, 'Sw_prev= ', Sw_prev
	            	print *, 'A*dx= ', (sum_w*zw1-sum_wz)*dx(j)	            	
	            	print *, 'aa= ', aa
	            	print *, 'bb= ', bb
	            	print *, 'cc= ', cc
                    print *, 'discr= ', discr	            		            		            		            	
	            	print *, 'zw1= ', zw1
	            	print *, 'zw2= ', zw2	            	
	            	print *, 'Sf= ', So(i,j)
	            	print *, 'Q= ', Q(i,j)
	            	print *, 'h= ', h(i,j,:)
	            	print *, 'v= ', v(i,j,:)	
	            	print *, 'E=', h(i,j,:)+v(i,j,:)**2/20
	            	print *, 'Before fluvial work in current t, Head=', Head(i,j,:)
	            	print *, 'Fr= ', v(i,j,:)/sqrt(10*h(i,j,:))             	
			    	do k=1,nD
			    	    print *, 'k= ', k
				    	print *, 'pa(k,:)= ', pa(k,:)
				    	print *, 'pss(k,:)= ', pss(k,:)		
			    	end do
	            endif
			    do kk=1,indW(i,j)
			        if (h(i,j,kk)<0) then
			            print *, 'err i= ', i
			            print *, 'j= ', j
			            print *, 'kk= ', kk
!			            i_h_err=i
!			            j_h_err=j
!			            kk_h_err=kk
!			            pause 1	!type 'go' to resume execution 
                        print *, 'h(i,j,kk) = ', h(i,j,kk)         
						STOP 'review where is failure such that h<0'
			        endif
			    end do
   			    Iw_prev(j)=Iw(j)   !input to this reach in next time comes from upstream in this time
				!c)SEDIMENT TRANSPORT CAPACITY --------------------------------------------------------------------   				
			    Is_prev(j,1)=Is_prev(j,1)+wash_1sh(i,j)*dt(i)*flagSupply
			    Is_prev(j,2)=Is_prev(j,2)+wash_2sh(i,j)*dt(i)*flagSupply
				do k=1,nD                
    			    Is_prev(j,k)=Is_prev(j,k)+(slid_sh(i,j,k))*dt(i)*flagSupply        			            
			    end do
				s5=0	!to be done: creep in valley bank pg~30 Integral notebook and Golly17 (SalettiFriend)
			    SedBalance=0
			    tao=0
			    dzdt=0; dLa_dt=0; fI_ok=0; fI=0; Fs=0; taorm=0; taor=0; qsbMx_wc=0; grad_qsmx_y=0
			    do kk=1,indW(i,j)
			    	tao(kk)=1e4*h(i,j,kk)*So(i,j)			        
			    end do			    
                if (kGrav>1) then			    
			        do kk=1,indW(i,j)			    
        			    Fs(kk)=sum(pa(1:kGrav-1,kk))                    
			        end do  
			    endif 
			    thcr84 = 0   
			    do kk=1,indW(i,j)
			        !positive exponent if S>.9: jamming?:			         
			         thm84a(i,j,kk)=alfaGTM*(5*So(i,j)+0.06)*(Da84(i,j,kk)/Da50(i,j,kk))**(4.4*sqrt(So(i,j))-1.5)	
    			     thcr84(kk) = alfaGTM*(1.32*So(i,j)+0.037)*(Da84(i,j,kk)/Da50(i,j,kk))**(-.93)   !eqn 2 Recking13, ASCE
!                     thm84a(i,j,kk) = .56*So(i,j) + 0.021 !GTMreking16, after recking09                                                                             
			         th50a(i,j,kk)=tao(kk)/(16500*Da50(i,j,kk))
			         th84a(kk) = tao(kk)/(16500*Da84(i,j,kk))
			         thratio84(i,j,kk) = th84a(kk) / thm84a(i,j,kk)  !check5: must be ~1 for Q=Qbf
			    end do			                 
!	!		    Is_prev(j,5)=Is_prev(j,5);%+s5;
			    sumIs=sum(Is_prev(j,:))			        
			    if (sumIs>0) then
			        pIs=Is_prev(j,:)/sumIs
			    else
			        pIs=0
			    endif
			    	!Recking equation
!			    thm84=1.5*So(i,j)**.75	!f(S) from PitonRecking16

			    qsbMx=0
			    PblB_acu=0; pblB=0
			    !Grain size portions in 'Budget'= supply + mobile alluvium
			    do kk=1,indW(i,j)
			        volAct=min(sumSedBed(i,j,kk),La(i,j,kk))*Wflow(j,kk)*dx(j)
			        sumIs_kk=Wflow(j,kk)/WtotSup*sumIs
			        if (sumIs_kk+volAct>0) then
			            do k=1,nD
			            	pblB(k,kk)=(pa(k,kk)*volAct+pIs(k)*sumIs_kk)/(sumIs_kk+volAct)
			            end do
			            if (thratio84(i,j,kk)>1) then	!SC: travelling bed load concept (Recking) is not being considered
			                 volSubsup=max(0.,sumSedBed(i,j,kk)-La(i,j,kk))*Wflow(j,kk)*dx(j)	!'max' to avoid consider 'negative' subsurface when zb<La
			                do k=1,nD
			                	pblB(k,kk)= (pblB(k,kk)*(sumIs+volAct)+pss(k,kk)*volSubsup)                 
							end do                  
							pblB(:,kk)= pblB(:,kk)/(sumIs+volAct+volSubsup)	!weighted average
			            endif
			            do k=1,nD
			            	PblB_acu(k,kk)=sum(pblB(1:k,kk))
			            end do
			            call GSDkeyPctiles(PblB_acu(:,kk), Dpct, Dp, D84blB, nD)
                		if (sum(PblB_acu(:,kk) - PblB_acu(:,kk)) .ne. 0) then
			                print *, 'err i= ', i
			                print *, 'j= ', j
			                print *, 'kk= ', kk                		    
                		    STOP 'NaN PblB_acu'   !NaN identifier						
                		endif
						thm84B(i,j,kk)=alfaGTM*(5*So(i,j)+0.06)*(D84blB/Dp(6))**(4.4*sqrt(So(i,j))-1.5)
			            th84blB=tao(kk)/(16500*D84blB)
			            if (j.le.jref_fin .and. j.ge.jref_ini .and. i_plot.le.iref_fin .and. i_plot.ge.iref_ini) then
			            	print *, 'kk= ', kk
			            	print *, 'tao(kk)= ', tao(kk)
			            	print *, 'D84blB= ', D84blB
			            	print *, 'D84/D50: ', D84blB/Dp(6)
			            	print *, 'thm84B(i,j,kk)= ', thm84B(i,j,kk)
			            endif    
			            qsbMx(kk)=sqrt(10*1.65*D84blB**3)*14*th84blB**2.5/(1+(thm84B(i,j,kk)/th84blB)**4)
			        endif
			    end do
			    !WilcockCrowe2003, surface based, alluvial
			    if (flagWC.eq.1) then
			        if (flag_jamming.eq.1) then
			            jam_fact(1)=max(1.,2.5-1.5/19*(Wflow(j,1)/Da84(i,j,1)-1)) !see 'ubc Experim cases for details. Based on Zimmermann2013'
			        endif
			    	do kk=1,indW(i,j)
			    		cum=1
                        if (sumSedBed(i,j,kk).gt.0) then			    		
			        		do k=1,nD
			        			cum=cum*Dcentr(k)**(pa(k,kk))
			        		end do
			        		dm(i,j,kk)=cum**1    !sum of exponents in 'cum' is 1=TotalProbability. 
!			        		dm(i,j,kk)=Da50(i,j,kk) 
			        	else
    			        	call GSDkeyPctiles(PblB_acu(:,kk), Dpct, Dp, D84blB, nD)
                    		if (sum(PblB_acu(:,kk) - PblB_acu(:,kk)) .ne. 0) then
			                    print *, 'err i= ', i
			                    print *, 'j= ', j
			                    print *, 'kk= ', kk                		    
                    		    STOP 'NaN PblB_acu'   !NaN identifier						
                    		endif    			        	
			        		do k=1,11
			        			cum=cum+Dp(k)   !Dp are percentiles in budget=supply(no alluvium) 
			        		end do
			        		dm(i,j,kk)=cum/11			        	
			        	endif        			    					        				    		
				    	taorm(kk)=jam_fact(kk)*alfaWC*(2.65-1)*1e3*10*dm(i,j,kk)*(.021+.015*exp(-20*Fs(kk)))
			    		do k=1,nD
			    			b=.67/(1+exp(1.5-Dcentr(k)/dm(i,j,kk)))
			    			taor(k,kk)=taorm(kk)*(Dcentr(k)/dm(i,j,kk))**b
			    			if (tao(kk)/taor(k,kk) .lt. 1.35) then
			    				Wstar=2e-3*(tao(kk)/taor(k,kk))**7.5
			    			else
			    				Wstar=14*(1-.894/sqrt(tao(kk)/taor(k,kk)))**4.5
			    			endif
			    			if (sumSedBed(i,j,kk).gt.0) then			    		
    			    			qsbMx_D(k)=(sqrt(tao(kk)/1000))**3*pa(k,kk)/((1.65-1)*10)*Wstar
   			    			else
   			    			    qsbMx_D(k)=(sqrt(tao(kk)/1000))**3*pblB(k,kk)/((1.65-1)*10)*Wstar
   			    			endif     			    			
!		                    if (j.le.jref_fin .and. j.ge.jref_ini .and. i_plot.le.iref_fin .and. i_plot.ge.iref_ini .and. &
!		                    & flag_printDetails.eq.1) then
!		                        print *, 'k= ', k
!		                    	print *, 'taor(k,kk)= ', taor(k,kk)
!		                    	print *, 'tao/taor (k)= ', tao(kk)/taor(k,kk)
!		                    	print *, 'Wstar (k)= ', Wstar
!		                    	print *, 'weighted by texture portion, qsbMx_D= ', qsbMx_D(k)		              
!		                    endif        			    		
			    		end do
			    		qsbMx_wc(kk)=sum(qsbMx_D)
!			    		if (qsbMx_wc(kk)>0) then
		    		    pbl_wc(:,kk)=qsbMx_D/qsbMx_wc(kk)				    		
!			    		else    !e.g. too low shear
!			    		    pbl_wc(1,kk)=1
!			    		    pbl_wc(2:nD,kk)=0			    		    
!			    		endif    
			    	end do
			    endif
			    
                if (j.le.jref_fin .and. j.ge.jref_ini .and. i_plot.le.iref_fin .and. i_plot.ge.iref_ini &
                & .and. flag_printDetails.eq.1) then                    
                	print *, 'across kk, dm= ', dm(i,j,:)
                endif
			    
   				if (i.eq.1 .and. j.eq.1) then
					pblB_i1=pblB    
				endif	
!			    !Grain size portions in bedload, by GTM proposed by Recking16		       
				call get_port(11,port)	!0 to 1, each .1
			    pbl=0; Pbl_acu=0; M=0; ff=0; ff_nD=0
			    do kk=1,indW(i,j)	!only in wetted zones of the cross section 		    			    
!			        if (qsbMx(kk) .gt. 0) then	!is there any available sediment (alluvial or from upstream)?
		            if (sumSedBed(i,j,kk).eq.0 .and. sumIs.eq.0) exit
		            !exit avoids NaNs. Existing values of pbl and qsbMx can remain. Anyway, Exner output is going to be 0
		            if (sumSedBed(i,j,kk) .gt. 0) then						                
		                if (thratio84(i,j,kk)>1) then
		                    betah = (13.4 - 2.91*log(ref)) / taofm_taocref
		                    beta=betah
		                else
		                    beta=betal        
		                endif		                
		                M(kk)=min(100.,ref*thratio84(i,j,kk)**beta)
		                gam=min(100.,gam0+gam1*thratio84(i,j,kk)**gam2)
		                pmx_mov=nint(M(kk)/10)+1  !resolution of grain pctiles each 10%
		                do ii=1,pmx_mov
							ff(ii, kk) = min(1.,fi0+(1-(100*port(ii)/M(kk))**gam))
		                end do
		                !average ff vs p/	M function to use with nD<11 grain sizes		                
		                k=1; jj=1 
		                cum=0; cont=0			                 			               
		                do while (PblB_acu(k,kk) .le. .1*(pmx_mov-1))
		                	cum=0; cont=0			                 
		                	do while (jj<=11 .and. port(jj)<=PblB_acu(k,kk) .and. ff(jj, kk)>fi0)
		                		cum=cum + ff(jj, kk)
		                		cont=cont+1
		                		jj=jj+1
		                	end do
							if (cont>0) then
			                	ff_nD(k, kk)=cum/cont
                            elseif (k .gt. 1) then    !two consecutive size PblB_acu values in same multiple of 10% in port, 
                            !e.g. PblB_acu = .05 .06 .16; 2nd value takes this 'else'
                                ff_nD(k, kk) = ff_nD(k-1, kk)
		                	endif
		                	k=k+1			                	
		                end do
!			                if (M(kk)/100<PblB_acu(1,kk)) then
                        if (flagWC.eq.0) then
!                            print *, 'eh, WHY ENTER HERE IF WC IS OFF?'
					        if (sum(ff_nD(:,kk)).eq.0) then
	                        	qsbMx(kk)=0 	!flow unable to move any grain of the classes modelled
	                        	do k=1,nD 	!only to avoid Nan fractional Qs 
	                        		if (k.eq.1) then
	                        			pbl(k,kk)=1
	                        		else
	                        			pbl(k,kk)=0
	                        		endif 
	                        	end do
	                        else
	                        	pbl(:,kk)=ff_nD(:,kk)*pblB(:,kk)
	                        	pbl(:,kk)=pbl(:,kk)/sum(pbl(:,kk)) !so sum=1, can be redundant
	                        endif               			                   			    
		                endif    
!			                if (i.eq.1819 .and. j.eq.jref) then
!			                	print *, 'port= ', port
!			                	print *, 'PblB_acu(:,kk)= ', PblB_acu(:,kk)
!			                	print *, 'thratio84(i,j,kk)= ', thratio84(i,j,kk)
!			                	print *, 'pmx_mov= ', pmx_mov
!			                	print *, 'kk= ', kk
!			                	print *, 'ff= ', ff
!			                	print *, 'cont= ', cont			
!			                	print *, 'ff_nD= ', ff_nD	
!			                endif               			                
                        if (j.le.jref_fin .and. j.ge.jref_ini .and. i_plot.le.iref_fin .and. i_plot.ge.iref_ini &
                        & .and. flag_printDetails.eq.1) then                    
                        	print *, 'kk:', kk, 'beta:', beta, 'gam: ', gam
                        	
                        endif
		            else
		                if (j.le.jref_fin .and. j.ge.jref_ini .and. i_plot.le.iref_fin .and. i_plot.ge.iref_ini) then
		                    print *, 'i= ', i 
		                    print *, 'eh, WHY ENTER HERE IF THERE IS ALLUVIUM IN THIS kk?'
		                endif    
		                pbl(:,kk)=pblB(:,kk)	!budget mixed in cross section if no alluvium
		            endif
	           		do k=1,nD
	                    Pbl_acu(k,kk)=sum(pbl(1:k,kk))
	                end do
			    end do
   				if (i.eq.1 .and. j.eq.1) then
					pbl_i1=pbl
				endif			    					   			
!				if (i.eq.nt .and. j.eq.jref) then
!					print *, 'from f90, i= ', i
!					print *, 'j= ', j
!					print *, 'j1, alluv_i1_bed portions= ', pa_i1
!					print *, 'j1, alluv_i1_budget portions= ', pblB_i1
!					print *, 'j1, alluv_i1_bedLoad portions= ', pbl_i1			    			   			    
!			    endif			       			  
!			end do 
!			do while (j .le. nR)    !to get every capacity to apply central difference scheme to solve qs gradients
			    !d) EXNER BULK MASS BALANCE ---------------------------------------------------------          
			    Is_kk=0; sedIn=0; sedOut=0; sup_lim = 0 
			    if (flagWC.eq.1 .and. sumSedBed(i,j,kk).gt.0) then
			    	qsbMx=qsbMx_wc
			    end if
			    grad_qsmx_y = Wflow(j,:)*tao**1.5
			    grad_qsmx_y = grad_qsmx_y/sum(grad_qsmx_y(1:indW(i,j)))    !portion (sum=1) of sediment capacity per zone of crosssection.
			    do kk=1,nW    !indW(i,j)
			        Is_kk = grad_qsmx_y*sumIs
			        sedIn(kk)=Is_kk(kk)/(Wflow(j,kk)*dx(j))	![m3/m2 = m]				
			        sedOut(kk)=qsbMx(kk)/dx(j)*dt(i)    ![m2/s/m*s = m]
!			        sedOut_jdown(kk)=
!			        qs_gradient=    !central difference, parker ebook ch17
			        SedBalance(kk)= sedIn(kk) - sedOut(kk)
			        sumSedBed(i+1,j,kk)=max(0.,sumSedBed(i,j,kk) + SedBalance(kk)/(1-poros))
			        if (sumSedBed(i+1,j,kk) .eq. 0.) then
			            sup_lim(kk) = 1
			            contSupLim = contSupLim + 1
                    endif			            
			    end do

			    
			    if (j.le.jref_fin .and. j.ge.jref_ini .and. i_plot.le.iref_fin .and. i_plot.ge.iref_ini .and. &
			    & flag_printDetails.eq.1) then
!			    	print *, 'qsbMx across valley: ', qsbMx(:)
			    	print *, 'along Xsection, jam_fact= ', jam_fact
			    	print *, 'IN (+) in balance [m]: ', sedIn/(1-poros)
			    	print *, 'grad_qsmx_y: ', grad_qsmx_y
			    	print *, 'SedBalance: ', SedBalance
			    	do kk = 1, nW
    			    	print *, 'kk:', kk, 'M(kk):', M(kk), 'ff_GTM_per_p:', ff(:, kk)
    			    	print *, 'PblB_acu(:,kk):', PblB_acu(:,kk)
    			    	print *, 'ff_nD:', ff_nD(:, kk)
                    end do    			    				    	
			    	print *, 'thm84a(i,j,:)', thm84a(i,j,:)
			    	print *, 'thratio84(i,j,:)', thratio84(i,j,:)
			    	print *, 'taorm= ', taorm			    	
			    	print *, 'tao= ', tao
			    	do k=1,nD
			    	    print *, 'pblB(k,:)', pblB(k,:)
!			    	    print *, 'pbl_wc(k,:)= ', pbl_wc(k,:)				    					    	
				    	print *, 'taor(k,:)= ', taor(k,:)
			    	end do	
			    	print *, 'qsbMx_wc: ', qsbMx_wc
			    	print *, 'OUT at capacity (-) in balance [m]: ', sedOut/(1-poros)	    	
			    	print *, 'sumSedBed(i,j,:)= ', sumSedBed(i,j,:)
			    	print *, 'AFTER simple Exner, BEFORE lateral reequil, sumSedBed(i+1,j,:)= ', sumSedBed(i+1,j,:)  
			    endif			    			    			    
			    
			    
			    
			    if (flagDF.eq.1) then
		            do kk=1,indW(i,j)
		            	So_DF=(sumSedBed(i+1,j,kk)-sumSedBed(i,jdown(j),kk))/dx(j)
				        if (So_DF .gt. .06) then	!Debris Flow mechanics in same time step, to reduce backwater in next time step
					        if (sumIs/dt(i)/Q(i,j) .gt. .05) then
						        countDF(j)=countDF(j)+1
						        qsbMx(kk)=v(i,j,k)*h(i,j,k)*min(.5,5.5*So_DF**2)	!assuming 8<h/D<25 (Mizuyama81; in Takahashi14, fig 2.30. See 2.51 too. See f3.30 for runout. f3.52 & 71 for scaling of lifetime and Qk of NatDam)
							        !eqn 6.8 assumed
						        print *, 'in DF, i= ', i
						        print *, 'in DF, j= ', j	        			        			      			        
						        print *, 'in DF, kk= ', kk
						        print *, 'in DF, cmx= ', min(.5,5.5*So_DF**2)
        !							stop
				            end if
				        endif			    
		            end do	
		        endif		    
			    if (flagDF.eq.1) then	!refresh alluvium balance
					do kk=1,indW(i,j)
					    sedOut(kk)=qsbMx(kk)/dx(j)*dt(i)
					    SedBalance(kk)=sumSedBed(i,j,kk)+(sedIn(kk)-sedOut(kk))/(1-poros)
					    sumSedBed(i+1,j,kk)=max(0.,SedBalance(kk))
					end do			    			    
			    endif			  			  			 			   			    			    			    
			    if (indW_max(j)>indW(i,j)) then
			        do kk=indW(i,j)+1,indW_max(j)
			            sumSedBed(i+1,j,kk)=sumSedBed(i,j,kk)
			        end do
			    endif

                Asx_sed_pre_lat_reeq = 0    ![m3]
                do kk=1, indW_max(j)
                    Asx_sed_pre_lat_reeq = Asx_sed_pre_lat_reeq + sumSedBed(i+1,j,kk)*Wflow(j,kk)
                end do
                
                                              
                    			   			   			   
			    !e) actual sediment output from reach, by grain sizes -------------------------------------------------------------
   			    do kk=1,indW(i,j)
       			    if (flagWC.eq.1 .and. sumSedBed(i,j,kk).gt.0) then
			        	pbl(:,kk)=pbl_wc(:,kk)
			        end if
			    end do    
			    do kk=1,indW(i,j)
			        if (sup_lim(kk) .eq. 0) then
			            Qsb=qsbMx(kk)*Wflow(j,kk)   !flow is capacity-limited
			        else
			            Qsb=(sumSedBed(i,j,kk)*Wflow(j,kk)*dx(j)+Is_kk(kk))/dt(i)   !flow is supply-limited
			        endif
			        do k=1,nD
			            Qs(i,j,k,kk)=Qsb*pbl(k,kk)
			        end do
			    end do
			    !!!what about %Qs_susp?
			    do k=1,nD
			        do kk=1,indW(i,j)
			            ck(i,j,k)=ck(i,j,k)+Qs(i,j,k,kk)            
			        end do
			        ck(i,j,k)=ck(i,j,k)/Q(i,j)
			    end do
!			    if (j.le.jref_fin .and. j.ge.jref_ini .and. i_plot.le.iref_fin .and. i_plot.ge.iref_ini) then
!			    	do k=1,nD
!!				    	print *, 'pblB(k,:)= ', pblB(k,:)				    				    	
!			    	end do
!			    endif			    
			    if (j.le.jref_fin .and. j.ge.jref_ini .and. i_plot.le.iref_fin .and. i_plot.ge.iref_ini .and. &
			    & flag_printDetails.eq.1) then
!			    	print *, 'pa(nD, nW)= ', pa 			    	
!			    	print *, 'pblB(nD, nW)= ', pblB
			    	do k=1,nD
!			    	    print *, 'pbl_wc(k,:)= ', pbl_wc(k,:)				    					    	
				    	print *, 'pbl(k,:)= ', pbl(k,:)				    					    	
			    	end do
			    	print *, 'ck(i,j,:)= ', ck(i,j,:)			    	
			    endif
			    do k=1,nD
			    	if (isnan(ck(i,j,k))) then 
			    		print *, 'i= ', i
			    		print *, 'j= ', j
			    		print *, 'k= ', k
			    		STOP 'ck is NaN'
			    	endif			    
				end do
			    if (jdown(j).gt.0 .and. i.lt.nt) then
			        do k=1,nD
			            Is(jdown(j),k)=Is(jdown(j),k)+ck(i,j,k)*Q(i,j)*dt(i+1)  !add sediment to downstream reach   
			        end do
			    endif
			    !f) actualize texture in surface and subsurface layer -------------------------------------------------------------				
			    !f.1) interface fractions, after Cui2007
			    do kk=1,indW(i,j)
			    	gravTpt=sum(pbl(kGrav:nD,kk))     
			        dzdt(kk) = (sumSedBed(i+1,j,kk) - sumSedBed(i,j,kk)) / dt(i)   
			            !as dzdt activates conditional codes only if it has no-null magnitude.
			        if (dzdt(kk).lt.0  .and. sumSedBed(i+1,j,kk)>La(i,j,kk)) then
			            fI(:,kk)=pss(:,kk)   !it does not matter if no real f when bedrock, cos used for ponder budget related to '0' subsurf volume             
			        elseif (dzdt(kk).gt.0 .and. sumSedBed(i+1,j,kk)>La(i,j,kk)) then
			            do k=1,kGrav-1 !sand
			                !not only 1 gravel size, as cui07, but kGrav-1. kGrav-1 DiffEqns to solve grain sizes
			                fI(k,kk)=Fs(kk)*pa(k,kk)/sum(pa(1:kGrav-1,kk)) !Fig2b (Cui2007). The more sand in active layer in last time step, the more infiltration to subsurface (?) !min(1.,Fs(kk)+.1)
			            end do                
			            fIs=sum(fI(1:kGrav-1,kk))                    
			            do k=kGrav,nD 	!gravel
			                if (pa(k,kk).gt.0 .and. gravTpt.gt.0) then	!toroEscobar ponder (cui07)
		                        fI(k,kk)=(1-fIs)*(alfaTEscobar*pa(k,kk)/sum(pa(kGrav:nD,kk))+&
		                        &(1-alfaTEscobar)*pbl(k,kk)/sum(pbl(kGrav:nD,kk)))             
		                    elseif (pa(k,kk).gt.0) then !alluvium but no gravel output                                
		                        fI(k,kk)=(1-fIs)*pa(k,kk)/sum(pa(kGrav:nD,kk))
			                elseif (gravTpt .gt. 0) then 	!no alluvium but gravel output
		                        fI(k,kk)=(1-fIs)*pbl(k,kk)/sum(pbl(kGrav:nD,kk))
			                else !only gravel could come from gravel stored supply in this time step
		                        fI(k,kk)=(1-fIs)/(nD-kGrav+1)                   
			                endif
			            end do
!			        else !erosion or deposition, but NOT existing subsurface alluvium in next t
!			        	fI(:,kk)=0
!!			        	do k=2,nD	!fI=1 for k1 only to get sum(fI)=1			        		
!!			        		fI(k,kk)=0	
!!			        	end do
!						contDeposNoSubsurf=contDeposNoSubsurf+1
			        endif
			    end do
			    
			    if (j.le.jref_fin .and. j.ge.jref_ini .and. i_plot.le.iref_fin .and. i_plot.ge.iref_ini .and. &
			    & flag_printDetails.eq.1) then
                    print *, 'dzdt(all kk):', dzdt
!                    do kk=1,indW(i,j)
!                        print *, 'kk', kk, 'fI:', fI(:,kk)
!                    end do
			    endif				    
			    
!   			    if (dzdt>0) then
!   			    	print *, 'i= ', i
!   			    	print *, 'j= ', j
!   			    	print *, 'fI= ', fI
!   			    	stop
!   			    endif
 			   				    			    			    			    			  
			    do kk=1,indW(i,j)
			        if (abs(dzdt(kk))>0 .and. sumSedBed(i+1,j,kk)>La(i,j,kk) .and. abs(1-sum(fI(:,kk)))>.05) then	!only if subsurf exists and change was generated
!	!		            i_fI_err=i
!	!		            j_fI_err=j
!	!		            kk_fI_err=kk
         
			        endif
			    end do			   
			    !f2) GSD in active layer (dLa/dt~0), La=active layer thickness
			    do kk=1,indW(i,j)
			        if (i.gt.1) then
    			        dLa_dt(kk)=(La(i,j,kk)-La(i-1,j,kk))/dt(i)  
    			        !if z(i-1) was bedrock then La(i-1)=0 (not assigned value), so actLayer is born
			        else !parker ebook, ch17, pg15.
!			            dLa_dt(kk)=La(i,j,kk)/dt(i)     !active layer is 'born'
                        dLa_dt(kk) = 0
			        endif
			    end do
			    overLoadF=0; F_interm = 0
			    do k=1,nD
			        do kk=1,indW(i,j)
			            Is_kk_k=Is_prev(j,k)*grad_qsmx_y(kk)
			        end do
                end do			    
		        do kk=1,indW(i,j)
	                if (sumSedBed(i+1,j,kk)>0) then
!			                Is_kk_k=Is_prev(j,k)*Wflow(j,kk)/WtotSup
	                    if (sumSedBed(i+1,j,kk)>La(i,j,kk)) then
	                        !if there is subsurface layer, then consider flux to/from subsurface (zb-La):
	                        fI_ok=fI(:,kk)  !i.e. it is not 0 anymore
	                        za=La(i,j,kk)
	                    else 
	                        za=sumSedBed(i+1,j,kk)
	                    endif
        			    do k=1,nD		                			                
	                        F_interm(k,kk)=F(i,j,k,kk)+&
	                        dt(i)/za*dLa_dt(kk)*(F(i,j,k,kk)-fI_ok(k))+dt(i)/za/(1-poros)*&
	                        &((Is_kk_k/dt(i)-Qs(i,j,k,kk))/(dx(j)*Wflow(j,kk))-(1-poros)*fI_ok(k)*dzdt(kk))
	                        F_interm(k,kk)=max(Fmin,F_interm(k,kk))            
	                        if (F_interm(k,kk) > 1) then
	                            overLoadF=overLoadF+1
!                                print *, 'i', i, 'j', j, 'k', k, 'kk', kk
!                                STOP 'overLoadF'
						        contOverLoadF=contOverLoadF+1
	                        endif
	                    end do
	                else
	                	F_interm(k,kk)=0
!			                if (k>1) then
!			                    F(i+1,j,k,kk)=0
!			                else
!			                    F(i+1,j,k,kk)=1 !to determine little D84 in next time step to ensure bedrock roughness
!			                endif
	                endif
		        end do
			    do k=1,nD
			        Is_prev(j,k)=Is(j,k)
			    end do 
			    
			    !F bug 1: exhausted alluvium?
!			    do kk=1,indW(i,j)
!			        if (sum(F(i+1,j,:,kk)).eq.0) then 
!			            do k=1,nD
!			                if (k>1) then 
!			                    F(i+1,j,k,kk)=0
!			                else
!			                    F(i+1,j,k,kk)=1 !to determine little D84 in next time step to ensure bedrock roughness
!			                endif
!			            end do            
!			        endif
!			    end do
			    !F bug 2: overload		        
			    if (overLoadF>0) then 
			        do kk=1,indW(i,j)
			            sumFposit=0
			            do k=1,nD
			                if (F_interm(k,kk)>0) then 
			                    sumFposit=sumFposit+F_interm(k,kk)
			                else
			                    F_interm(k,kk)=0	!absent size. This line do not seem to be necessary
			                endif
			            end do
			            F_interm(:,kk)=F_interm(:,kk)/sumFposit !!!!!!!!!!!!!!!questionable            
			        end do
			    endif
			    !F bug 3: rigurous balance
			    do kk=1,indW(i,j)
			        sumF=sum(F_interm(:,kk))
			        if (abs(1-sumF)>.0001) then 
			        	F_interm(:,kk)=F_interm(:,kk)/sumF
			        endif
			    end do
			    !F bug 4: unmodified size fractions in 0flow zones of crosssection
			    if (indW_max(j)>indW(i,j)) then 
			        do kk=indW(i,j)+1,indW_max(j)
			            F_interm(:,kk)=F(i,j,:,kk)
			        end do
			    endif
			    !F bug 5: does balance work?
			    do kk=1,indW(i,j)
			        if (sumSedBed(i+1,j,kk)>0 .and. abs(1-sum(F_interm(:,kk)))>.05) then 
!	!		            i_F_err=i
!	!		            j_F_err=j
!	!		            kk_F_err=kk
!	!		            pause;
			            print *, 'err i= ', i
			            print *, 'j= ', j
			            print *, 'kk= ', kk
			            print *, 'F_interm(:,kk)', F_interm(:,kk)
						STOP 'review where is failure in F_interm'						         
			        endif
			    end do
!			    !f3) fractions in subsurface layer
				ILay=0; ssLay=0; p_depos_up = 0; f_sub_interm = 0
			    do kk=1,indW(i,j)		            
		            if (sumSedBed(i+1,j,kk)>La(i,j,kk) .and. dzdt(kk).lt.0) then
	                    f_sub_interm(:,kk)=pss(:,kk)
	                elseif (sumSedBed(i+1,j,kk)>La(i,j,kk)) then	!deposition on alluvium	                    
	                    ssLay(kk)=max(0.,sumSedBed(i,j,kk)-La(i,j,kk))	                                        	                    
                        p_depos_up = fI(:,kk)
                        ILay(kk)=dzdt(kk)*dt(i)
!           			    if (j.le.jref_fin .and. j.ge.jref_ini .and. i_plot.le.iref_fin .and. i_plot.ge.iref_ini .and. &
!           			    & flag_printDetails.eq.1) then
!       			    		print *, 'kk', kk, 'p_redistr_mix:', p_redistr_mix
!       			    		print *, 'p_depos_up: ', p_depos_up
!       			    		print *, 'ILay(kk): ', ILay(kk)
!       			    		print *, 'ssLay(kk): ', ssLay(kk)
!			            endif				            	                    
	                    do k=1,nD
!		                	print *, 'i=', i
!		                	print *,'enter to D loop: OK'
	                        f_sub_interm(k,kk) = (ILay(kk)*p_depos_up(k) + ssLay(kk)*pss(k,kk)) / (ILay(kk) + ssLay(kk))
	                    end do
		            else	!alluvium is thinner than active layer
!						f_sub(i+1,j,:,kk)=0			        commented on June1, 2018
                        f_sub_interm(:,kk)=pa(:,kk) !10jul2019: it matters if lateral redistribution occur.
!		                do k=1,nD
!		                    if (k<nD) then 
!		                        f_sub_interm(k,kk)=0
!		                    else
!		                        f_sub_interm(k,kk)=1 !just to sum 1, cos no effect when ponders in next step cos vol subsurf=0
!		                    endif           
!		                end do
		            endif
			    end do
			    !fsub bug 1: no coarse structure with voids 
!			    do kk=1,indW(i,j)            
!			        if (sum(f_sub(i+1,j,1:nD-1,kk))<f_sub_min) then 
!			            fcorrec=fmin-sum(f_sub(i+1,j,1:k,kk))
!			            f_sub(i+1,j,1,kk)=f_sub(i+1,j,1,kk)+fcorrec    !!!!'comes' new fine to fill skeleton (from susp load?)..hypoth: enough fine supply so no organal
!			            f_sub(i+1,j,nD,kk)=f_sub(i+1,j,nD,kk)-fcorrec
!			        endif            
!			    end do
				!fsub bug 2: unmodified size fractions in 0flow zones of crosssection
			    if (indW_max(j)>indW(i,j)) then 
			        do kk=indW(i,j)+1,indW_max(j)
			            f_sub_interm(:,kk)=f_sub(i,j,:,kk)
			        end do
			    endif
			    
			    




			    
			    
			    dz_lateral_reequil=0
                z_prev=sumSedBed(i+1,j,:)
                
       			!if trenching exceeds repose angle, then redistribute volume in crosssection
       			trench=1
       			if (flag_trench.eq.1) then
           			do while (trench .eq. 1)       			
    !			        do kk=2,indW_max(j)    !Wflow(j,kk)+1    !,indW_max(j)
    !			            if (sumSedBed(i+1,j,kk)-sumSedBed(i+1,j,kk-1) .gt. .25*(Wflow(j,kk)+Wflow(j,kk-1))*tan_repos_ang) then
    !!			                sumSedBed(i+1,j,kk) .lt. sumSedBed(i,j,kk) .and.
    !			            
    !			                !dz>dz_max?
    !			                trench=1			            
    !			                kk_steep=kk
    !			            endif
    !			        end do
                        trench = 0; tan_stabl_ang = 0			        
		                do kk=2,indW_max(j) !from center to wall
		                    if (thratio84(i,j,kk) .gt. 1) then
		                        tan_stabl_ang(kk) = 0
		                    elseif (h(i,j,kk) .gt. 0) then
		                        tan_stabl_ang(kk) = tan(asin(sin_dry_ang*sqrt(1-(th84a(kk)/thcr84(kk))**2)))  !adapted from millar1993
		                    else
		                        tan_stabl_ang(kk) = tan_repos_ang
		                    endif
		                    dz_trasv = sumSedBed(i+1,j,kk) - sumSedBed(i+1,j,kk-1)
		                    dz_adm = .25*(Wflow(j,kk) + Wflow(j,kk-1))*tan_stabl_ang(kk)
		                    if (dz_trasv .gt. dz_adm) then
                                cum_area = Wflow(j,kk-1)*sumSedBed(i+1,j,kk-1) + Wflow(j,kk)*sumSedBed(i+1,j,kk)
    !			                sumSedBed(i+1,j,kk) = sumSedBed(i+1,j,kk) - dz_exc
		                        sumSedBed(i+1,j,kk-1) = (cum_area - Wflow(j,kk)*dz_adm) / (Wflow(j,kk) + Wflow(j,kk-1))
		                        sumSedBed(i+1,j,kk) = sumSedBed(i+1,j,kk-1) + dz_adm
                            endif
		                end do
		                do kk2 = 2,indW_max(j) !from wall to center
		                    kk = indW_max(j) - kk2 + 2
		                    dz_trasv = sumSedBed(i+1,j,kk) - sumSedBed(i+1,j,kk-1)
		                    dz_adm = .25*(Wflow(j,kk) + Wflow(j,kk-1))*tan_repos_ang
		                    if (dz_trasv .gt. dz_adm) then
                                cum_area = Wflow(j,kk-1)*sumSedBed(i+1,j,kk-1) + Wflow(j,kk)*sumSedBed(i+1,j,kk)
		                        sumSedBed(i+1,j,kk-1) = (cum_area - Wflow(j,kk)*dz_adm) / (Wflow(j,kk) + Wflow(j,kk-1))
		                        z_pre = sumSedBed(i+1,j,kk) 
		                        sumSedBed(i+1,j,kk) = sumSedBed(i+1,j,kk-1) + dz_adm
		                        if (abs(sumSedBed(i+1,j,kk)-z_pre)/sumSedBed(i+1,j,kk) .gt. tol_trench) then
		                            trench = 1
!		                            print *, 'inside return trenching'
		                        endif  
                            endif
		                end do			        
    !			        sumSedBed(i+1,j,1) = sumSedBed(i+1,j,1) + cum_area/Wflow(j,1)
                    end do
		        endif

!			    if (trench.eq.1) then			        
!!	                do while (trench.eq.1)		        
!       			    	cum=0
!			        	do kk=1,kk_steep
!			        	    cum=cum+sumSedBed(i+1,j,kk)*Wflow(j,kk)
!			        	end do
!			            zm=cum/Wfeas(j,kk_steep)
!			            zmin=zm-.25*Wfeas(j,kk_steep)*tan_repos_ang !small base of trapecium
!			            sumSedBed(i+1,j,1)=zmin+.25*Wflow(j,1)*tan_repos_ang
!			            do kk=2,kk_steep
!			                sumSedBed(i+1,j,kk)=sumSedBed(i+1,j,kk-1)+.25*(Wflow(j,kk-1)+Wflow(j,kk))*tan_repos_ang    
!			            end do
!			            trench=0
!!			            do kk=2,indW_max(j)    !Wflow(j,kk)+1    !,indW_max(j)
!!			                if (sumSedBed(i+1,j,kk)-sumSedBed(i+1,j,kk-1) .gt. .5*(Wflow(j,kk)+Wflow(j,kk-1))*tan_repos_ang) then
!!			                    !dz>dz_max?
!!			                    trench=1			            
!!			                    kk_steep=kk
!!			                endif
!!			            end do			            
!!			        end do                    			        			        
!			    endif
			    
			    if (j.le.jref_fin .and. j.ge.jref_ini .and. i_plot.le.iref_fin .and. i_plot.ge.iref_ini .and. &
			    & flag_printDetails.eq.1) then
			        print *, 'tan_stabl_ang(all kk)', tan_stabl_ang
                    print *, 'AFTER trench, BEFORE de-braid, sumSedBed(i+1,j,:)= ', sumSedBed(i+1,j,:)  
			    endif			    

			    
			    fixThalw=0; kk_braid = 0; kk = 1; cont_debraid = 0
				!not sure but...: get unique channel in cross section, preserving its sediment volume
			    
			    do while ((kk_braid .eq. 0) .and. (kk .lt. indW_max(j)))
    			    kk = kk + 1
			        if (sumSedBed(i+1,j,kk) .lt. sumSedBed(i+1,j,kk-1)) then
			            fixThalw=1
			            kk_braid = kk-1   !'this is the most central zone of crosssection where island emerges'
                        cont_debraid = cont_debraid + 1			            
		                if (j.le.jref_fin .and. j.ge.jref_ini .and. i_plot.le.iref_fin .and. i_plot.ge.iref_ini .and. &
		                & flag_printDetails.eq.1) then
                            print *, 'kk_braid, before big while: ', kk_braid
                        endif                            
			        endif
			    end do
!			    zm=sum(sumSedBed(i+1,j,1:indW(i,j)))/indW_max(j)
			    if (fixThalw.eq.1) then			        
    			    contBraidIJ=contBraidIJ+1
	                do while (fixThalw.eq.1)
	                    cum = 0			        	
			        	!assume island can be smoothed as zm in its external remaining crosssection. 
			        	!set thalweg seed after smoothing			        	
!			        	sum_seed=0                                               
!                        do while (sumSedBed(i+1,j,kk_tooDeep).lt.sumSedBed(i+1,j,kk_braid) .and. kk_tooDeep.le.indW_max(j))
                        do kk = kk_braid + 1 , indW_max(j)
                            if (sumSedBed(i+1,j,kk).lt.sumSedBed(i+1,j,kk_braid)) then  !find gap closest to valley wall
                                kk_mx_tooDeep = kk                            
                            endif
			                if (j.le.jref_fin .and. j.ge.jref_ini .and. i_plot.le.iref_fin .and. i_plot.ge.iref_ini .and. &
			                & flag_printDetails.eq.1) then
			                    print *, 'inside while of de-braid'
			                endif
                        end do
		                if (j.le.jref_fin .and. j.ge.jref_ini .and. i_plot.le.iref_fin .and. i_plot.ge.iref_ini .and. &
		                & flag_printDetails.eq.1) then
                            print *, 'kk_braid= ', kk_braid
                            print *, 'kk_mx_tooDeep= ', kk_mx_tooDeep
                        end if                        
                        do kk=kk_braid,kk_mx_tooDeep
                            cum = cum + sumSedBed(i+1,j,kk)*Wflow(j,kk)                            
                        end do
                        sum_Wext = sum(Wflow(j,kk_braid:kk_mx_tooDeep))                        
			            zm=cum/sum_Wext
			            cum_thw = 0
			            do kk=kk_braid,kk_mx_tooDeep
                            sumSedBed(i+1,j,kk) = zm - (kk_mx_tooDeep-kk)*thalwSeed        
                            cum_thw = cum_thw + sumSedBed(i+1,j,kk)*Wflow(j,kk)
			        	end do
			        	z_offset = (cum - cum_thw)/sum_Wext !to restore mass balance
		                if (j.le.jref_fin .and. j.ge.jref_ini .and. i_plot.le.iref_fin .and. i_plot.ge.iref_ini .and. &
		                & flag_printDetails.eq.1) then
                            print *, 'zm: ', zm
                            print *, 'z_offset: ', z_offset                                                                
                        end if			        	
			        	do kk=kk_braid,kk_mx_tooDeep
			        	    sumSedBed(i+1,j,kk) = sumSedBed(i+1,j,kk) + z_offset
			        	end do

			            if (j.le.jref_fin .and. j.ge.jref_ini .and. i_plot.le.iref_fin .and. i_plot.ge.iref_ini .and. &
			            & flag_printDetails.eq.1) then			                
			                print *, 'zm: ', zm
			            endif			    	
			        	

    !			        z1=max(0.,zm-thalwSeed)
    !			        z2=(zm*Wfeas(j,indW_max(j))-z1*Wflow(j,1))/(Wfeas(j,indW_max(j))-Wflow(j,1))
    !			        sumSedBed(i+1,j,1)=z1
    !			        do kk=2,indW_max(j)
    !			            sumSedBed(i+1,j,kk)=z2
    !			        end do
                        fixThalw=0
			            do while ((kk_braid .eq. 0) .and. (kk .lt. indW_max(j)))
            			    kk = kk + 1
			                if (sumSedBed(i+1,j,kk) .lt. sumSedBed(i+1,j,kk-1)) then
			                    fixThalw=1
			                    kk_braid = kk-1   !'this is the most central zone of crosssection where island emerges'
			                    cont_debraid = cont_debraid + 1
		                        if (j.le.jref_fin .and. j.ge.jref_ini .and. i_plot.le.iref_fin .and. i_plot.ge.iref_ini .and. &
		                        & flag_printDetails.eq.1) then
                                    print *, 'kk_braid, after debraiding: ', kk_braid
                                    print *, 'it is coming debraiding #', cont_debraid
                                endif                            
			                endif
			            end do                      
                    end do
			    endif
			    
			    if (j.le.jref_fin .and. j.ge.jref_ini .and. i_plot.le.iref_fin .and. i_plot.ge.iref_ini .and. &
			    & flag_printDetails.eq.1) then
                    print *, 'AFTER lateral reequil, sumSedBed(i+1,j,:)= ', sumSedBed(i+1,j,:)
                    do kk=1,indW_max(j)
                        if (sumSedBed(i+1,j,kk) .ne. sumSedBed(i+1,j,kk)) then   !it is a way to identify NaN values
                            print *, 'err i= ', i
                            print *, 'j= ', j
                            print *, 'kk= ', kk                        
                            print *, 'sumSedBed(i+1,j,kk)', sumSedBed(i+1,j,kk)
                            STOP 'NaN bed elevation'
                        endif                                                   
                    end do
			    endif			    
			    
			    
	            !if negative storage results applying iniXsect rule, then redistribute in cross section
	            !june7, 2018: very low likely cos I started to use too thick alluvium
	            if (sumSedBed(i+1,j,1)<0) then
	                contBedrock=contBedrock+1
	                sedPlus=-sumSedBed(i+1,j,1)*Wflow(j,1)
	                cont=1	
	                !sum material to redistribute		        
	                do kk=2,indW_max(j)
	                    contBedrock=contBedrock+1
	                    if (sumSedBed(i+1,j,kk)<0) then
			                sedPlus=sedPlus-sumSedBed(i+1,j,kk)*Wflow(j,kk)
			                cont=cont+1
			            endif
	                end do 
	                !redistribute in no-bedrock portions of crosssection
                    do kk=1,indW_max(j)
                        if (kk>cont) then
    		                sumSedBed(i+1,j,kk) = sumSedBed(i+1,j,kk) + &
    		                & (sedPlus*Wflow(j,kk)/sum(Wflow(j,cont+1:indW_max(j)))) / Wflow(j,kk)
                        else
                            sumSedBed(i+1,j,kk)=0
                        endif
                    end do			            
	            endif			    
			    
			    dz_lateral_reequil=sumSedBed(i+1,j,:)-z_prev    !collects changes (+/-) in elevation along crosssection by lateral reequilibrium
			    
			    !if any lateral reequilibrium then turn ON flag to refresh bed texture along crosssection:
			    refresh_texture = 0 !flag per zone: break evolution where lateral redistribution occurs.
			    do kk=1,indW_max(j)
			        if (abs(dz_lateral_reequil(kk))>1e-3) then
                        refresh_texture(kk)=1
			        endif
                    if (abs(dz_lateral_reequil(kk)) .gt. sumSedBed(i+1,j,kk)) then  !negative z_prev??
                        print *, 'err i= ', i
                        print *, 'j= ', j
                        print *, 'kk= ', kk 
                        print *, 'dz_lateral_reequil: ', dz_lateral_reequil                           
                        STOP 'there is not enough material to erode in this lateral redistribution (& negative z_prev??)'
                    end if			        
			    end do
			    
                if (sum(refresh_texture) .gt. 0) contRefreshText_LatReeq = contRefreshText_LatReeq + 1
			    
			    if (j.le.jref_fin .and. j.ge.jref_ini .and. i_plot.le.iref_fin .and. i_plot.ge.iref_ini .and. &
			    & flag_printDetails.eq.1) then			        
			        print *, 'dz_lateral_reequil(all kk): ', dz_lateral_reequil
                    print *, 'refresh_texture(all kk): ', refresh_texture
!                    if (sum(refresh_texture) .eq. 0) then
!                        print *, 'i= ', i
!                        print *, 'j= ', j
!                        print *, 'kk= ', kk                    
!                        STOP 'look: there is case when simple texture evolution happen cos no lateral redistri'
!                    endif
			    endif				   		             
                
                Asx_sed_pos_lat_reeq = 0    ![m3]
                do kk=1,indW_max(j)
                    Asx_sed_pos_lat_reeq = Asx_sed_pos_lat_reeq + sumSedBed(i+1,j,kk)*Wflow(j,kk)
                end do
                
			    if (j.le.jref_fin .and. j.ge.jref_ini .and. i_plot.le.iref_fin .and. i_plot.ge.iref_ini .and. &
			    & flag_printDetails.eq.1) then			        
			        print *, 'Asx_sed_pre_lat_reeq (i+1): ', Asx_sed_pre_lat_reeq
			        print *, 'Asx_sed_pos_lat_reeq (i+1): ', Asx_sed_pos_lat_reeq
			        if (abs(Asx_sed_pos_lat_reeq - Asx_sed_pre_lat_reeq) / Asx_sed_pos_lat_reeq .gt. .01) then
			            print *, 'i', i, 'j', j, 'rel_err >1% in Asx_sed due to lat_reeq' 
                        STOP 
			        endif
			    endif			    
			    
			    
			    


			    
			    
			    
			    
                if (sum(refresh_texture).ge.1) then     !is there any lateral reequilibrium?                
                !texture refresh by lateral reequilibrium (due to both 'DE'braiding and/or trenching):
                !refreshing is (1) mixing all eroded zones to (2) distribute it among depositing ones.
                !note no change in surface texture occurs at kk such that dz_lateral_reequil(kk)=0.
                    As_redistr=0    !solid area of removed material, cumulated in eroding zones of crosssection, per grain size.
                    abs_dz_latReeq = abs(dz_lateral_reequil)
                    do kk=1,indW_max(j)
                    !(1) mixing: cumulate erosion amount by size classes, and reduce active layer thickness:
                        if (dz_lateral_reequil(kk).lt.0) then   !erosion                              
                            if (-dz_lateral_reequil(kk).gt.La(i,j,kk)) then  !erosion deeper than active layer
                                F(i+1,j,:,kk)=f_sub_interm(:,kk)     !subsurface layer gets exposed to alluvium surface

                                As_redistr = As_redistr + (1-poros)*Wflow(j,kk)*(F_interm(:,kk)*La(i,j,kk)+ &
                                & f_sub_interm(:,kk)*(abs_dz_latReeq(kk)-La(i,j,kk)))
                            else 
                                F(i+1,j,:,kk) = (f_sub_interm(:,kk)*(La(i,j,kk)-abs_dz_latReeq(kk)) + &
                                & F_interm(:,kk)*abs_dz_latReeq(kk))/La(i,j,kk) !texture of active layer becomes a mixture with subsurface layer
                                As_redistr = As_redistr + (1-poros)*F_interm(:,kk)*abs_dz_latReeq(kk)*Wflow(j,kk)
                            endif
                        endif                 
                    end do
                    p_redistr_mix = As_redistr / sum(As_redistr)
                    do kk=1,indW_max(j)    
                    !(2) redistribution: increase La (hypothesis) and refresh texture in depositing zones of crosssection:
                        if (dz_lateral_reequil(kk).gt.0) then   !only where fill occurs
                            if (dz_lateral_reequil(kk).gt.La(i,j,kk)) then  !deposition deeper than active layer
                                F(i+1,j,:,kk) = p_redistr_mix
                            else
                                F(i+1,j,:,kk) = (p_redistr_mix*abs_dz_latReeq(kk) + &
                                & F_interm(:,kk)*(La(i,j,kk)-abs_dz_latReeq(kk))) / La(i,j,kk)
                            endif                                                                                                                    
                        endif
		                if (j.le.jref_fin .and. j.ge.jref_ini .and. i_plot.le.iref_fin .and. i_plot.ge.iref_ini .and. &
		                & flag_printDetails.eq.1) then                        
		                    if (abs(1-sum(F_interm(:,kk)))>.05) then
		                        print *, 'err i= ', i
		                        print *, 'j= ', j
		                        print *, 'kk= ', kk         
		                        print *, 'all kk, dz_lateral_reequil= ', dz_lateral_reequil
		                        print *, 'all k, F_interm(:,kk)= ', F_interm(:,kk)
					            STOP 'failure in F inside refreshing texture after lateral reequilibrium'
		                    endif
		                endif
                    end do		            
                endif
                
                do kk=1,indW_max(j) !pass F_interm to next time step for kk_zones which did not suffer lat_reeq (dz_lateral_reequil(kk)=0)
                    if (refresh_texture(kk) .eq. 0) F(i+1,j,:,kk) = F_interm(:,kk)
                end do
                
			    do kk=1,indW(i,j)
			        if (sumSedBed(i+1,j,kk)>0 .and. abs(1-sum(F(i+1,j,:,kk)))>.05) then 
			            print *, 'err i= ', i
			            print *, 'j= ', j
			            print *, 'kk= ', kk         
						STOP 'review where is failure in F'						         
			        endif
			    end do                
                
!			    !fractions in subsurface layer, for lateral redistribution
				ILay=0; ssLay=0; p_depos_up = 0
			    do kk=1,indW_max(j)		            
		            if (sumSedBed(i+1,j,kk)>La(i,j,kk) .and. dz_lateral_reequil(kk).lt.0) then  !erosion in lat_reequil
	                    f_sub(i+1,j,:,kk)=f_sub_interm(:,kk)
	                elseif (sumSedBed(i+1,j,kk)>La(i,j,kk)) then	!deposition on alluvium	                    
	                    ssLay(kk)=max(0.,sumSedBed(i,j,kk)-La(i,j,kk))	                                        	                    
	                    if (refresh_texture(kk) .eq. 1) then    !this 'if' might be avoided if f_sub ponderation were set as subroutine.
	                        p_depos_up = p_redistr_mix
                            ILay(kk) = dz_lateral_reequil(kk)
	                    endif	                    
!           			    if (j.le.jref_fin .and. j.ge.jref_ini .and. i_plot.le.iref_fin .and. i_plot.ge.iref_ini .and. &
!           			    & flag_printDetails.eq.1) then
!       			    		print *, 'kk', kk, 'p_redistr_mix:', p_redistr_mix
!       			    		print *, 'p_depos_up: ', p_depos_up
!       			    		print *, 'ILay(kk): ', ILay(kk)
!       			    		print *, 'ssLay(kk): ', ssLay(kk)
!			            endif				            	                    
	                    do k=1,nD
!		                	print *, 'i=', i
!		                	print *,'enter to D loop: OK'
	                        f_sub(i+1,j,k,kk) = (ILay(kk)*p_depos_up(k) + ssLay(kk)*f_sub_interm(k,kk)) / (ILay(kk) + ssLay(kk))
	                    end do
		            else	!deposition on bedrock or thin alluvium with no bedrock
!						f_sub(i+1,j,:,kk)=0			        commented on June1, 2018
		                do k=1,nD
		                    if (k<nD) then 
		                        f_sub(i+1,j,k,kk)=0
		                    else
		                        f_sub(i+1,j,k,kk)=1 !just to sum 1, cos no effect when ponders in next step cos vol subsurf=0
		                    endif                
		                end do
		            endif
			    end do
!				!fsub bug 2: unmodified size fractions in 0flow zones of crosssection
!			    if (indW_max(j)>indW(i,j)) then 
!			        do kk=indW(i,j)+1,indW_max(j)
!			            if (f_sub(i+1,j,k,kk) .eq. 0) f_sub(i+1,j,:,kk) = f_sub_interm(:,kk)
!			        end do
!			    endif
                do kk=1,indW_max(j) !even in dry channel
                !pass f_sub_interm to next time step for kk_zones which did not suffer lat_reeq (dz_lateral_reequil(kk)=0)
                    if (refresh_texture(kk) .eq. 0) f_sub(i+1,j,:,kk) = f_sub_interm(:,kk)
                end do
         
                
                
                
                
                
                
                
                
                
                
                
                
                
                			    
			    
!			    if (i.gt.295 .and. i.lt.306 .and. j.eq.18) then
!		            print *, 'i= ', i
!		            print *, 'j= ', j
!		            print *, 'dzdt(kk:)=', dzdt						
!		            print *, 'Fs(kk=1)=', Fs(1)
!		            print *, 'fIs=', fIs
!		            print *, 'pss(:,kk=1)=', pss(:,1)
!	            	print *, 'fI(:,kk=1)=', fI(:,1)
!	            	print *, 'ILay(kk=1)=', ILay(1)
!	            	print *, 'ssLay(kk=1)=', ssLay(1)
!					print *, 'f_sub(i+1,j,:,kk=1)', f_sub(i+1,j,:,1)
!			    	print *, '----------------another i coming----------------------------'
!			    	print *, '----------------another i coming----------------------------'						            			            
!		            if (i .eq. 305) then
!			            STOP 'review why deposit does not change subsurface texture'
!		            endif
!			    endif

   			    if (j.le.jref_fin .and. j.ge.jref_ini .and. i_plot.le.iref_fin .and. i_plot.ge.iref_ini .and. &
   			    & flag_printDetails.eq.1) then
   			    	do k=1,nD
   			    		print *, 'k', k, 'fI(k,:)= ', fI(k,:)
   			    		print *, 'F_interm(k,:)= ', F_interm(k,:)
   			    		print *, 'F(i+1,j,k,:)= ', F(i+1,j,k,:)
   			    		print *, 'f_sub_interm(k,:)= ', f_sub_interm(k,:)
   			    		print *, 'f_sub(i+1,j,k,:)= ', f_sub(i+1,j,k,:)
   			    	end do
			    endif
			    
			    do kk=1,indW_max(j)
   			    	do k=1,nD
                        if (F(i+1,j,k,kk) .lt. 0) then
                            print *, 'i', i, 'j', j, 'k', k, 'kk', kk
                            STOP 'F(i+1,j,k,kk) < 0'
                        endif   			    	
                        if (f_sub(i+1,j,k,kk) .lt. 0) then
                            print *, 'i', i, 'j', j, 'k', k, 'kk', kk
                            STOP 'f_sub(i+1,j,k,kk) < 0'
                        endif                         
   			    	end do                    			    
			    end do
			    
			     
			    !does balance work?
			    do kk=1,indW(i,j)
			        if (sumSedBed(i+1,j,kk)>La(i,j,kk) .and. abs(1-sum(f_sub(i+1,j,:,kk)))>.05) then 
!	!		            i_f_err=i
!	!		            j_f_err=j
!	!		            kk_f_err=kk
!	!		            pause;
			            print *, 'err i= ', i
			            print *, 'j= ', j
			            print *, 'kk= ', kk
			            print *, 'La(j)= ', La(i,j,kk)
			            print *, 'f_sub(i+1,j,:,kk)', f_sub(i+1,j,:,kk)       
						STOP 'review where is failure in f_sub'          
			        endif
			    end do
			    !can I consider Cui2007's entrainment eqns, to actualize f_sub to couple with Lisle02 (armour but no scour)?
!	!		    g) adverse Sf? backwater?. Assume conditions in active (inner) channel as representative of crosssection -----------------------------
                Head(i,j,:)=zr(j)+sumSedBed(i+1,j,:)+h(i,j,:)+v(i,j,:)**2/20    !based on alluvium after work
                if (j.eq.nR .and. flagBackwat.eq.1) then 
                    !set new backwater pulses where energy obstacles appear
                    do jj=1,nR-1      			        
      			        !alluvium refreshed to foresee backwater
!                        Hup=E(i,jj,1)+sumSedBed(i,jj,1)+zr(jj)
!                        Hdown=E(i,jdown(jj),1)+sumSedBed(i,jdown(jj),1)+zr(jdown(jj))                        
                        if (Head(i,jdown(jj),1)>Head(i,j,1)) then !send signal upstream
                            j3=jdown(jj)
                            zbw(i,j3)=zr(j3)+sumSedBed(i,j3,1)+h(i,j3,1) !initial pulse
                            latent(i,j3)=1. !pulse is present now. It is not waiting to arrive
                            celer=sqrt(10*h(i,j3,1)) !wave celerity c~sqrt(gh) if h/L<20                            
                            if (celer*dt(i).gt.dx(j3)) then
                                !signal travels upstream faster than dx/dt                                
                                t_cum=dx(j3)/celer
                                do while (dt(i)>t_cum) !available time to travel
                                    if ((jup(j3).eq.0).or.(zbw(i,j3).lt.sumSedBed(i,jup(j3),1)+zr(jup(j3)))) then
                                        exit    !alluvium or headwater were reached by projected backwater
                                    endif
                                    latent(i+1,j3)=1.
                                    zbw(i+1,jup(j3))=max(zbw(i+1,jup(j3)),int(latent(i,j3))*zbw(i,j3))    
                                    !max signal wins, latent is activated only if travel was completed
                                    j3=jup(j3)  !go upstream for next position
                                    celer=sqrt(10*h(i,j3,1))
                                    t_cum = t_cum + dx(j3)/celer               
                                end do                                
                            else  
                                zbw(i+1,jup(j3))=zbw(i,j3)
                                latent(i+1,jup(j3))=min(1.,latent(i,jup(j3))+celer*dt(i)/dx(j3))
                                !memory of previous partial travels to get next dx upstream
                            endif                        
                        endif    
                    end do
                endif
                if (j.le.jref_fin .and. j.ge.jref_ini .and. i_plot.le.iref_fin .and. i_plot.ge.iref_ini &
                & .and. flag_printDetails.eq.1) then
                    print *, 'after refreshing alluvium, Head=', Head(i,j,:)
                    print *, '----------------another j coming----------------------------'
                endif    	    	
				j=j+1
			end do						
			fileunit=12
			
			nvar=1
			write(fmt=*, unit= xy(nvar)), path_list(1:num_char,pos_ini_var(nvar))
!			write(fmt=*, unit= xx), output_folder//'/'//xy
            xx=trim(output_folder)//'/'//adjustl(xy(nvar))
			call write_FloatArrIJ(fileunit,xx,Q(i,:),i,nR)			
			nvar=8
			write(fmt=*, unit= xy(nvar)), path_list(1:num_char,pos_ini_var(nvar))
			xx=trim(output_folder)//'/'//adjustl(xy(nvar))
			call write_FloatArrIJ(fileunit,xx,So(i,:),i,nR)
			nvar=9
			write(fmt=*, unit= xy(nvar)), path_list(1:num_char,pos_ini_var(nvar))
			xx=trim(output_folder)//'/'//adjustl(xy(nvar))
			call write_FloatArrIJ(fileunit,xx,real(indW(i,:)),i,nR)
			nvar=10
			write(fmt=*, unit= xy(nvar)), path_list(1:num_char,pos_ini_var(nvar))
			xx=trim(output_folder)//'/'//adjustl(xy(nvar))
			call write_FloatArrIJ(fileunit,xx,th50a(i,:,1),i,nR)
			nvar=11
			write(fmt=*, unit= xy(nvar)), path_list(1:num_char,pos_ini_var(nvar))
			xx=trim(output_folder)//'/'//adjustl(xy(nvar))
			call write_FloatArrIJ(fileunit,xx,thm84B(i,:,1),i,nR)
			nvar=12
			write(fmt=*, unit= xy(nvar)), path_list(1:num_char,pos_ini_var(nvar))
			xx=trim(output_folder)//'/'//adjustl(xy(nvar))
			call write_FloatArrIJ(fileunit,xx,thratio84(i,:,1),i,nR)
			nvar=14
		    write(fmt=*, unit= xy(nvar)), path_list(1:num_char,pos_ini_var(nvar))               
		    xx=trim(output_folder)//'/'//adjustl(xy(nvar))
            call write_FloatArrIJ(fileunit,xx,Qhill(i,:),i,nR)			
												                
!			call write_FloatArrIJ(fileunit,path_Q,Q(i,:),i,nR)  !record = time step, recordLenght = #reaches
!        	call write_FloatArrIJ(fileunit,path_So,So(i,:),i,nR)
!			call write_FloatArrIJ(fileunit,path_indW,real(indW(i,:)),i,nR)
!			call write_FloatArrIJ(fileunit,path_thm84a,thm84a(i,:,1),i,nR)
!			call write_FloatArrIJ(fileunit,path_thm84B,thm84B(i,:,1),i,nR)
!			call write_FloatArrIJ(fileunit,path_thratio84,thratio84(i,:,1),i,nR)			
			
!			print *, 'path_list= ', path_list
			do kk=1,nW	
			    nvar=2
			    write(fmt=*, unit= xy(nvar)), path_list(1:num_char,pos_ini_var(nvar)+kk-1)                
			    xx=trim(output_folder)//'/'//adjustl(xy(nvar))
			    call write_FloatArrIJ(fileunit,xx,h(i,:,kk),i,nR)
!			    print *, 'kk= ', kk
!			    print *, 'path_h(kk)=', path_h(kk)
                nvar=3
			    write(fmt=*, unit= xy(nvar)), path_list(1:num_char,pos_ini_var(nvar)+kk-1)                
			    xx=trim(output_folder)//'/'//adjustl(xy(nvar))
!                print *, 'given kk, path_v= ', xy(nvar)			    
			    call write_FloatArrIJ(fileunit,xx,v(i,:,kk),i,nR)
			    nvar=4
			    write(fmt=*, unit= xy(nvar)), path_list(1:num_char,pos_ini_var(nvar)+kk-1)               
			    xx=trim(output_folder)//'/'//adjustl(xy(nvar))
			    call write_FloatArrIJ(fileunit,xx,sumSedBed(i,:,kk),i,nR)
			    nvar=5
			    write(fmt=*, unit= xy(nvar)), path_list(1:num_char,pos_ini_var(nvar)+kk-1)                
			    xx=trim(output_folder)//'/'//adjustl(xy(nvar))
			    call write_FloatArrIJ(fileunit,xx,Da84(i,:,kk),i,nR)
			    nvar=6
			    write(fmt=*, unit= xy(nvar)), path_list(1:num_char,pos_ini_var(nvar)+kk-1)            
			    xx=trim(output_folder)//'/'//adjustl(xy(nvar))
			    call write_FloatArrIJ(fileunit,xx,Da50(i,:,kk),i,nR)			    
			    nvar=7
			    write(fmt=*, unit= xy(nvar)), path_list(1:num_char,pos_ini_var(nvar)+kk-1)                
			    xx=trim(output_folder)//'/'//adjustl(xy(nvar))
			    call write_FloatArrIJ(fileunit,xx,Dss50(i,:,kk),i,nR)			     			    
			    			    
!			    call write_FloatArrIJ(fileunit,path_h_kk2,h(i,:,2),i,nR)
!			    call write_FloatArrIJ(fileunit,path_h_kk3,h(i,:,3),i,nR)
!			    call write_FloatArrIJ(fileunit,path_h_kk4,h(i,:,4),i,nR)			
!			    call write_FloatArrIJ(fileunit,path_v_kk1,v(i,:,1),i,nR)
!			    call write_FloatArrIJ(fileunit,path_v_kk2,v(i,:,2),i,nR)
!			    call write_FloatArrIJ(fileunit,path_z_kk1,sumSedBed(i,:,1),i,nR)
!			    call write_FloatArrIJ(fileunit,path_z_kk2,sumSedBed(i,:,2),i,nR)
!			    call write_FloatArrIJ(fileunit,path_z_kk3,sumSedBed(i,:,3),i,nR)
!			    call write_FloatArrIJ(fileunit,path_z_kk4,sumSedBed(i,:,4),i,nR)
!    !			call write_FloatArrIJ(fileunit,path_zw,zw(i,:),i,nR)						
!			    call write_FloatArrIJ(fileunit,path_Da84_kk1,Da84(i,:,1),i,nR)
!			    call write_FloatArrIJ(fileunit,path_Da84_kk2,Da84(i,:,2),i,nR)
!			    call write_FloatArrIJ(fileunit,path_Da84_kk3,Da84(i,:,3),i,nR)
!			    call write_FloatArrIJ(fileunit,path_Da84_kk4,Da84(i,:,4),i,nR)			
!			    call write_FloatArrIJ(fileunit,path_Da50_kk1,Da50(i,:,1),i,nR)
!			    call write_FloatArrIJ(fileunit,path_Da50_kk2,Da50(i,:,2),i,nR)
!			    call write_FloatArrIJ(fileunit,path_Da50_kk3,Da50(i,:,3),i,nR)
!			    call write_FloatArrIJ(fileunit,path_Da50_kk4,Da50(i,:,4),i,nR)			
!			    call write_FloatArrIJ(fileunit,path_Dss50_kk1,Dss50(i,:,1),i,nR)
!			    call write_FloatArrIJ(fileunit,path_Dss50_kk2,Dss50(i,:,2),i,nR)
!			    call write_FloatArrIJ(fileunit,path_Dss50_kk3,Dss50(i,:,3),i,nR)
!			    call write_FloatArrIJ(fileunit,path_Dss50_kk4,Dss50(i,:,4),i,nR)									
			end do
        	do k=1,nD
!        	    xy=path_c(k)
!                print *, 'path_c(k)=', path_c(k)
!                print *, 'path_c_k1=', path_c_k1
!                print *, 'k= ', k
                nvar=13
                write(fmt=*, unit= xy(nvar)), path_list(1:num_char,pos_ini_var(nvar)+k-1)                
!                print *, 'given k, path_c= ', xy(nvar)
                xx=trim(output_folder)//'/'//adjustl(xy(nvar))
        	    call write_FloatArrIJ(fileunit,xx,ck(i,:,k),i,nR)       	            	            	    
			end do
        	do k=1,nD_hill
                nvar=15
                write(fmt=*, unit= xy(nvar)), path_list(1:num_char,pos_ini_var(nvar)+k-1)                
                xx=trim(output_folder)//'/'//adjustl(xy(nvar))
        	    call write_FloatArrIJ(fileunit,xx,slid_sh(i,:,k),i,nR)        	            	            	    
			end do			
						
!				fileunit=fileunit+1
!				call write_FloatArrIJ(fileunit,path_h(kk),h(i,:,kk),i,nR)
!				fileunit=fileunit+1
!				call write_FloatArrIJ(fileunit,path_v(kk),v(i,:,kk),i,nR)						
            
!            day_of_calc=day_of_calc+dt(i)/86400.            
            
            if (j.le.jref_fin .and. j.ge.jref_ini .and. i_plot.le.iref_fin .and. i_plot.ge.iref_ini &
            & .and. flag_printDetails.eq.1) then
                print *, '-------------------ANOTHER i COMING-----------------------------------------------'
            endif    	    	            
            
		end do
		
		print *, 'contAdverseSlope= ', contAdverseSlope
		print *, 'contBraidIJ= ', contBraidIJ
		print *, 'contBedrock= ', contBedrock
		print *, 'contRefreshText_LatReeq', contRefreshText_LatReeq	
		print *, 'contSupLim', contSupLim	
		print *, 'contOverLoadF= ', contOverLoadF
		print *, 'contBF_artif= ', contBF_artif
		print *, 'iniXsect= ', iniXsect
!		do j=1,nR
!			print *, 'j= ', j
!			print *, 'zw_mean= ', sum(zw(:,j))/nt			
!		end do 				

!		print *, 'path_c= ', path_c
		
!		dum1=h(:,:,1); dum2=v(:,:,1); dum3=Fr(:,:,1); dum4=Q; dum5=E; dum6=rel_submg; dum7=W_h;		
	end subroutine
	
	
end module
