module ddt_class_fld
    implicit none
    integer i,j,k
    
    contains
    
    function temporal_averager_2d(nx, nt, array)
        integer, intent(in) :: nx, nt
        real, intent(in) :: array(nt, nx)
        real temporal_averager_2d(nx)
		do j=1, nx !accounts for gradients in average rainfall
    		temporal_averager_2d(j) = sum(array(:,j))/nt
		end do
    end function
    
    function FloodWhenWhere(n_coarse_dt, nR, Qhill_mean, Qhill_in)
        integer, intent(in) :: n_coarse_dt, nR
        real, intent(in) :: Qhill_mean(nR), Qhill_in(n_coarse_dt, nR)
        logical :: FloodWhenWhere(n_coarse_dt, nR)
        do j=1, nR
            FloodWhenWhere(:, j) = Qhill_in(:, j) .gt. Qhill_mean(j)
        end do
    end function
        
    subroutine printDdtBasin(range2zoomDdt, coarse_step, nt, dt, k3)
        integer, intent(in) :: range2zoomDdt(3), coarse_step, nt
        real, intent(in) :: dt(nt)
        integer, intent(inout) :: k3
        integer flag_day_else_tStep, iref_ini, iref_fin
        
        flag_day_else_tStep = range2zoomDdt(1); iref_ini = range2zoomDdt(2); iref_fin = range2zoomDdt(3)
        if (flag_day_else_tStep.eq.1 .and. coarse_step.ge.iref_ini-2 .and. coarse_step.le.iref_fin+2) then   
        !2 days before and after the ones to be zoomed		                    
!            print *, 'k3= ', k3
!            print *, 'dt(k3)= ', dt(k3)
        endif        
    end subroutine
         
    subroutine applyDdt_flume(nt, nR, n_coarse_dt, ddt_args, dt_dQ_crs_arr, k3, array4ddt, refinedArray, day_of_calc)
        integer, intent(in) :: nt, nR, n_coarse_dt, ddt_args(6)
        real, intent(in) :: array4ddt(n_coarse_dt, nR), dt_dQ_crs_arr(n_coarse_dt)
        integer, intent(inout) :: k3
        real, intent(out) :: refinedArray(nt, nR)
        integer, intent(out) :: day_of_calc(nt)
        integer k3ant, ind_t
        integer :: dt_dQ_coarse, dt_ev
        
!        print *, '6mar20, vectorcrudo j1: ', array4ddt(:,1)
!        print *, '6mar20, vectorcrudo j2: ', array4ddt(:,2)
        dt_ev = ddt_args(1)     !        dt_dQ_coarse = ddt_args(7);  !9apr2020: variable dt_dQ_coarse for flood flume        
        
	    do i=1, n_coarse_dt    	!dt_dQ_coarse is 24 in basin from daily lumped SHIA hydrologic model, or Maria's experiment
!	        print *, '---------new coarse t, i=', i, 'k3 = ', k3
	        k3ant = k3
            ind_t = int(dt_dQ_crs_arr(i)/dt_ev)
!            print*,'-------------------------------------------------------------'
!            print *, '9apr2020, in applyDdt_flume, at i=', i,',  ind_t =', ind_t
            do j=1, nR
                refinedArray(k3 : k3+ind_t-1, j) = array4ddt(i, j)
!                print *, '9apr2020, in applyDdt_flume, at j=', j,',  array4ddt(i, j) =', array4ddt(i, j)
!                print *, '9apr2020, in applyDdt_flume, at j=', j,',  refinedArray(k3 : k3+ind_t-1, j) =', &
!                    refinedArray(k3 : k3+ind_t-1, j)
    		end do
    		k3 = k3 + ind_t
    		day_of_calc(k3ant:k3-1) = i
		end do
!        print *, 'vectorRefined j1: ', refinedArray(1:65, 1)
!		print *, 'sizes of crude and refined vector: ', size(array4ddt(:,1)), size(refinedArray(:,1))
    end subroutine
    
    subroutine applyDdt_basin_iniDay(fpd, dt_day1, nt, nR, steps_event_day, k3, array4ddt_iniDay, refinedArray)
        integer, intent(in) :: fpd, nt, nR, steps_event_day
        real, intent(in) :: array4ddt_iniDay(nR), dt_day1
        integer, intent(inout) :: k3
        real, intent(out) :: refinedArray(nt, nR)
        integer :: advanc_ind = 0
!        print*,'6may20; in applyDdt_basin_iniDay; array4ddt_iniDay: ', array4ddt_iniDay
        refinedArray = 0.0
	    if (dt_day1<3e3) advanc_ind = steps_event_day  !i=1 does not have antecedent to estimate pulse.
!	    if (fpd.eq.1) print*, 'dt_day1 = ', dt_day1
!	    if (fpd.eq.1) print*, 'advanc_ind= ', advanc_ind
	    do i=k3, k3 + advanc_ind
            refinedArray(i, :) = array4ddt_iniDay
        end do            
        k3 = k3 + advanc_ind + 1
    end subroutine
    
    
    subroutine applyDdt_basin_pulse_slopes(fpd, tb_hill, tpk_hill, array4ddt_xSt_2t, ris_slop, fall_slop)
        integer, intent(in) :: fpd
        real, intent(in) :: tb_hill, tpk_hill, array4ddt_xSt_2t(2)
        real, intent(out) :: ris_slop, fall_slop
        real Qpeak, Qbase, Qmean_today
        integer :: thpd = 48*3600  !twice_hperday_secs

        Qbase = array4ddt_xSt_2t(1);  Qmean_today = array4ddt_xSt_2t(2)
!        if (fpd.eq.1) print*, 'Qbase', Qbase, 'Qmean_today', Qmean_today, 'thpd/tb_hill= ', thpd, '/', tb_hill
        Qpeak = Qbase + thpd/tb_hill*(Qmean_today - Qbase)
!        if (fpd.eq.1) print*, 'in applyDdt_basin_pulse_slopes, Qpeak= ', Qpeak        
        !this hydrograph differs from pg39_Catano2015, given here Q(t+1) is Qbase: it suits better for simulQdaily.
        ris_slop = (Qpeak - Qbase) / tpk_hill
        fall_slop = (Qbase - Qpeak) / (tb_hill - tpk_hill)
!        if (fpd.eq.1) print*, 'ris_slop/fall_slop= ', ris_slop, '/', fall_slop
    end subroutine
        
        
    subroutine applyDdt_basin_pulse_trend(fpd, nt, k3, slope, dt_today, tRef, ii, t_pulse, refinedArray_xSt, Qbase)
        integer, intent(in) :: fpd, nt, k3
        real, intent(in) :: slope, tRef, Qbase, dt_today
        integer, intent(inout) :: ii
        real, intent(inout) :: t_pulse, refinedArray_xSt(nt)
        real dy, prev_val
        if (ii.eq.0) then
!            if (fpd.eq.1) print*, 'just inside applyDdt_basin_pulse_trend by 1st time: rise'
            prev_val = Qbase
        else
!            if (fpd.eq.1) print*, 'just inside applyDdt_basin_pulse_trend by 2nd time: fall'
            prev_val = refinedArray_xSt( k3 + ii - 1)
        endif            
!        if (fpd.eq.1) print*, 'max lim for t_pulse [hours]: tRef= ', tRef
!        if (fpd.eq.1) print*, 'dt_today = ', dt_today
        dy = slope*dt_today
        do while (t_pulse.lt.tRef)  !RISING limb in local (subbasin) hydrograph
!            if (fpd.eq.1) print*, 'ii= ', ii
!            if (fpd.eq.1) print*, 't_pulse[hour]= ', t_pulse
            refinedArray_xSt( k3 + ii ) = max(Qbase, prev_val  + dy )
            ii = ii + 1
            prev_val = refinedArray_xSt( k3 + ii - 1)
            t_pulse = t_pulse + dt_today
        end do
    end subroutine
    
    subroutine applyDdt_basin_pulse_here(fpd, nt, k3, t_coarse, steps_event_day, tb_hill, tpk_hill, dt, &
        array4ddt_xSt, refinedArray_xSt)
        integer, intent(in) :: fpd, nt, k3, t_coarse, steps_event_day
        real, intent(in) :: tb_hill, tpk_hill, array4ddt_xSt(nt), dt(nt)
        real, intent(inout) :: refinedArray_xSt(nt)
        real :: ris_slop, fall_slop, array4ddt_xSt_2t(2), t_pulse, Qbase
        integer :: ii
        ii = 0;  t_pulse = 0.0
!        if (fpd.eq.1) print*, 'just started applyDdt_basin_pulse_here so ii==0, so ii= ', ii
!        if (fpd.eq.1) print*, 'checking other initializations, t_pulse==0.0, so it = ', t_pulse        
        array4ddt_xSt_2t = array4ddt_xSt(t_coarse-1 : t_coarse)   !yesterday & today.
        Qbase = array4ddt_xSt_2t(1)
        call applyDdt_basin_pulse_slopes(fpd, tb_hill, tpk_hill, array4ddt_xSt_2t, ris_slop, fall_slop)
        call applyDdt_basin_pulse_trend(fpd, nt, k3, ris_slop, dt(k3), tpk_hill, ii, t_pulse, refinedArray_xSt, Qbase)
!        if (fpd.eq.1) print*, 'just between limbs ii>0 is expected, so ii= ', ii        
        call applyDdt_basin_pulse_trend(fpd, nt, k3, fall_slop, dt(k3), tb_hill, ii, t_pulse, refinedArray_xSt, Qbase)
        refinedArray_xSt( k3+ii : k3 + steps_event_day) = Qbase   !must be filled with Qbase = yesterday
    end subroutine
    
    subroutine applyDdt_basin_pulse(fpd, needPulse, ic, nt, nR, n_coarse_dt, steps_event_day, k3, tb_hill, tpk_hill, &
        range2zoomDdt, dt, FloodTodayHere, array4ddt, refinedArray)   !dt smaller than 1 hour.
        !intraday water pulse from rain event is assumed
        integer, intent(in) :: fpd, ic, nt, nR, steps_event_day, n_coarse_dt, range2zoomDdt(3)
        logical, intent(in) :: needPulse, FloodTodayHere(n_coarse_dt, nR)
        real, intent(in) :: tb_hill, tpk_hill, array4ddt(n_coarse_dt, nR), dt(nt)
        integer, intent(inout) :: k3
        real, intent(inout) :: refinedArray(nt, nR)
        real :: Qtoday, Qyesterday
        logical fth
          
        if (fpd.eq.1) call printDdtBasin(range2zoomDdt, ic, nt, dt, k3)   !ic: i_coarse, e.g. day.
        if (needPulse) then
            do j=1, nR
                fth = FloodTodayHere(ic, j)
                Qtoday = array4ddt(ic, j);  Qyesterday = array4ddt(ic-1, j)
!                if (fpd.eq.1) print*,  'in reach j= ', j, &
!                    ', pulseConds: FloodTodayHere/Qy>Qt: ', fth, Qyesterday, '/', Qtoday
                if (ic.lt.n_coarse_dt .and. Qyesterday.lt.Qtoday) then  !fth .and. 
                    !build triangular hydrograph only where (j) it rained during day ic, which implies Qmean today > yesterday
!                    if (fpd.eq.1) print*, '-------------pulse in reach j= ', j
                    call applyDdt_basin_pulse_here(fpd, nt, k3, ic, steps_event_day, tb_hill, tpk_hill, dt, &
                        array4ddt(:, j), refinedArray(:, j))
		        else
                    refinedArray(k3 : k3 + steps_event_day, j) = array4ddt(ic, j)
		        end if				            	                                        				    
            end do        
        else
            do j=1, nR        
!                if (array4ddt(ic, j) .gt. 1e-9) then
!                    print*, 'landslide supply[m3/s]=', array4ddt(ic, j), ', toDay', ic, 'here/nR', j, '/', nR, &
!                        '; run will stop'; stop
!                endif
                refinedArray(k3 : k3 + steps_event_day, j) = array4ddt(ic, j)
            end do                                
        endif
!        print*, '6may20; in last line of applyDdt_basin_pulse()'
    end subroutine   
        
    subroutine applyDdt_basin(fpd, wv, needPulse, nt, nR, n_coarse_dt, ddt_args, k3, tb_hill, &
        range2zoomDdt, dt, FloodTodayHere, array4ddt, refinedArray, day_of_calc)
	    integer, intent(in) :: fpd, wv, nt, nR, n_coarse_dt, range2zoomDdt(3), ddt_args(6)
	    logical, intent(in) :: needPulse, FloodTodayHere(n_coarse_dt, nR)	    
	    real, intent(in) :: tb_hill, array4ddt(n_coarse_dt, nR), dt(nt)
	    real, intent(out) :: refinedArray(nt, nR)
        integer, intent(out) :: day_of_calc(nt)	    
        integer, intent(inout) :: k3
        integer steps_event_day
        real final_step_in_day, tpk_hill, bal   !tpk_hill to get SCS triangular Qt as pg 39 in master thesis of SCatano2015.
        integer k3ant, jp
        character*10 :: which_var
!        print*, '6may20: just IN applyDdt_basin() of ddt f90 module.'
        jp = 1  !reach ID to print example of mass balance check.
        k3ant = k3        
        tpk_hill = tb_hill/2.7;  steps_event_day = sum(ddt_args(4:6)) - 1
        if (fpd.eq.1) call printDdtBasin(range2zoomDdt, 1, nt, dt, k3)
!        if (fpd.eq.1) print*, 'Santi, welcome to f90 (ddt_class). steps_event_day= ', steps_event_day
!        if (fpd.eq.1) print*, 'en applyDdt_basin, justo despues de salir de printDdtBasin'
        if (wv.eq.1) then
            which_var = 'Qhill'                            
        elseif (wv.eq.2) then
            which_var = 'wash'
        else
            which_var = 'slid'
        endif
!        if (fpd.eq.1) print*, 'for var ', which_var, 'in 1st step of day ', 1, ', dt will be ', dt(1), 's'
        call applyDdt_basin_iniDay(fpd, dt(1), nt, nR, steps_event_day, k3, array4ddt(1, :), refinedArray)
	    day_of_calc(k3ant: k3-1) = 1
!        print*, 'for reach', jp, ', must be 1 for balance: Qday_input / Qday_asFineMean = ', &
!            array4ddt(1, 1), ' / ', sum(refinedArray(k3ant: k3-1, 1) * dt(k3ant: k3-1)) / sum(dt(k3ant: k3-1))
            !ok: reviewed jan10, 2020; 8:02 GMT-5.
	    do i=2, n_coarse_dt
!            if (fpd.eq.1) print*, '--------------------------------------------------------------'
!            if (fpd.eq.1) print*, 'for var ', which_var, 'in 1st step of day ', i, ', dt will be ', dt(k3), 's, and k3=', k3
!	        print*, '6may20; in applyDdt_basin(), before choose if event day for i=', i, 'n_coarse_dt=', n_coarse_dt
	        if (dt(k3)<3e3) then   !As k3 is 1st step in day, it contains an event. Remember that [dt] = seg.
!                print*, '6may20; in applyDdt_basin(); event day before go into pulse'
                call applyDdt_basin_pulse(fpd, needPulse, i, nt, nR, n_coarse_dt, steps_event_day, k3, tb_hill, tpk_hill, &
                    range2zoomDdt, dt, FloodTodayHere, array4ddt, refinedArray)
                final_step_in_day = k3 + steps_event_day
	        else
    			refinedArray(k3, :) = array4ddt(i, :)
                final_step_in_day = k3
	        endif
	        k3ant = k3        
	        k3 = final_step_in_day + 1
	        day_of_calc(k3ant:k3-1) = i
!            if (wv.eq.1) then   !ok: reviewed jan10, 2020; 8:02 GMT-5.
!                print*, 'for reach', jp, ', must be 1 for balance: Qday_input / Qday_asFineMean = ', &
!                    array4ddt(i, 1), ' / ', sum(refinedArray(k3ant: k3-1, 1) * dt(k3ant: k3-1)) / sum(dt(k3ant: k3-1))
!            endif
            if (wv.eq.1) then   !ok: reviewed jan10, 2020; 8:02 GMT-5.
                bal =  array4ddt(i, 1) / ( sum( refinedArray(k3ant: k3-1, 1) * dt(k3ant: k3-1) ) / sum(dt(k3ant: k3-1)) )
                if (abs(1-bal)>.1) then
!                    print*, 'will STOP for reach', jp, '. Ratio must be 1 for balance: Qday_input / Qday_asFineMean : ', &
!                        array4ddt(i, 1), ' / ', sum(refinedArray(k3ant: k3-1, 1) * dt(k3ant: k3-1)) / sum(dt(k3ant: k3-1))
                    print*, 'in applyDdt_basin(); for DEBUG; i_coarse=', i ,'Q abs(err)=', abs(1-bal)
                    stop 'err balance >10%'
                endif                        
            endif
	    end do
!        if (fpd.eq.1) print*, 'var ', which_var, ', reach', jp, ', ratio must be 1 for balance: Qday_input / Qday_asFineMean = ', &
!            sum(array4ddt(:, 1))/size(array4ddt(:, 1)), &
!            ' / ', sum(refinedArray(:, 1) * dt(:)) / sum(dt)
    end subroutine
    
	subroutine applyDdt(fpd, wv, needPulse, flume_else_basin, nt, nR, ddt_args, n_coarse_dt, &
        tb_hill, dt_dQ_crs_arr, dt, range2zoomDdt, FloodTodayHere, array4ddt, refinedArray, day_of_calc)
	    integer, intent(in) :: fpd, wv, flume_else_basin, nt, nR, ddt_args(6)
        integer, intent(in) :: n_coarse_dt, range2zoomDdt(3)
	    logical, intent(in) :: needPulse, FloodTodayHere(n_coarse_dt, nR)
	    real, intent(in) :: tb_hill, array4ddt(n_coarse_dt, nR), dt_dQ_crs_arr(n_coarse_dt), dt(nt)
	    real, intent(out) :: refinedArray(nt, nR)
        integer, intent(out) :: day_of_calc(nt)	    
		integer k3  !k3 serves for ddt subroutine, called from hgy (Q) & from sed (Qhill).
        k3=1
!        print*, '6may20: just IN applyDdt() of ddt f90 module.'
        if (flume_else_basin .eq. 1) then
            !flume. Experim duration is less than 1e4min so you don't need to coarsen inefficient time periods of process.
            call applyDdt_flume(nt, nR, n_coarse_dt, ddt_args, dt_dQ_crs_arr, k3, array4ddt, refinedArray, day_of_calc)
		else
!		    if (fpd.eq.1) print*, 'en applyDdt, sí lee que está en cuenca'
            call applyDdt_basin(fpd, wv, needPulse, nt, nR, n_coarse_dt, ddt_args, k3, tb_hill, &
                range2zoomDdt, dt, FloodTodayHere, array4ddt, refinedArray, day_of_calc)
	 	endif
	 	!Qramp starting flume to avoid exhaustion wave. (26MAY 2018 or 2019: it is not too relevant to avoid exhaustion)
        !here there was module for lab sediment feed adapted to Maria's experiment; see texturalModel.        
    end subroutine 
end module
