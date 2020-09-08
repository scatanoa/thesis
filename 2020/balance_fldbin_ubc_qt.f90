module balance_fldbin_ubc_qt
    implicit none
    integer i,j,k, k2
    
    contains
    
    real function fv(qu, coef_fu, sq, sq3)
        real, intent(in) :: qu, coef_fu(12), sq, sq3
        real :: b_1, b_3, b_2, a_1, a_2 !, ta    !, w_b !D_84, S, g, 
!        S = coef_fu(3); g = coef_fu(4) !D_84 = coef_fu(2); 
        b_1 = coef_fu(6); b_3 = coef_fu(7); b_2 = coef_fu(8); a_1 = coef_fu(9); a_2 = coef_fu(10)
        fv = a_1 * (qu/sq3)**b_1 * sq * ((qu/(a_2*sq3))**b_2 + 1)**b_3
    end function    

    real function fv2(qu, coefrix, sq, sq3)
        real, intent(in) :: qu, coefrix(5), sq, sq3
        real :: b_1, b_3, b_2, a_1, a_2
        b_1 = coefrix(1); b_3 = coefrix(2); b_2 = coefrix(3); a_1 = coefrix(4); a_2 = coefrix(5)
        fv2 = a_1 * (qu/sq3)**b_1 * sq * (  ( qu / (a_2*sq3) )**b_2 + 1  )**b_3
    end function
    
    real function f_dcy(a1_2, a2_2, a12_2, y) !see dedux at p206d
        real, intent(in) :: y, a12_2, a1_2, a2_2
        f_dcy = 8*y**2 * ( a1_2 + a2_2/y**1.7 ) / a12_2
    end function    
    
    real function fv_dcy(sg8, f, R)
        real, intent(in) :: sg8, f, R
        fv_dcy = sqrt(sg8*R/f)
    end function
    
    real function fQ(cp3, cp4, hfld)
        real, intent(in) :: cp3, cp4, hfld
        fQ = cp3 - cp4*hfld
    end function    
    
    real function fbar(Drp, Q, qf, wf, b_1, wcb1, qcs_a2, b_2, a2_sq3f, b_3) !deduct p
        real, intent(in) :: Drp, Q, qf, wf, b_1, wcb1, qcs_a2, b_2, a2_sq3f, b_3
        fbar = Drp * (Q/qf - wf)**b_1/wcb1 * ( (qcs_a2**b_2 + 1) / ((qf/a2_sq3f)**b_2 + 1) )**b_3
    end function
    real function fej( x, coef) !wiki ejemplo: https://es.wikipedia.org/wiki/M%C3%A9todo_de_Newton
        real, intent(in) :: x, coef
        fej = coef*cos(x) - x**3
    end function
        
    subroutine find_root_fej( xinit, coef, tol, maxiter, result, success )  !f, 
!        interface   !17mar20: deprecated as f2py seems unable to catch this interface.
!            real function f(x)
!                real, intent(in) :: x
!            end function f
!        end interface
        real, intent(in)     :: xinit, coef, tol
        integer, intent(in)  :: maxiter
        real, intent(out)    :: result
        logical, intent(out) :: success
        real                 :: eps = 1e-2 !default: 1.0e-4
        real                 :: fx1, fx2, fprime, x, xnew
        integer              :: i

        result  = 0.0
        success = .false.
        x = xinit
        do i = 1,max(1,maxiter)
            fx1    = fej(x, coef)
            fx2    = fej(x+eps, coef)
            write(*,*) i, fx1, fx2, eps
            fprime = (fx2 - fx1) / eps
            xnew = x - fx1 / fprime
            if ( abs(xnew-x) <= tol ) then
                success = .true.
                result  = xnew
                exit
            endif
            x = xnew
         enddo
    end subroutine    
    
    real function fu(w, coefs, coef_fu)
        real, intent(in) :: w, coefs(4), coef_fu(12)
        real :: sq, sq3, qs
        real :: Q, D_84, S, g, w_b, b_1, b_3, b_2, a_1, a_2, ta, a1_sq, Q_sq3, ta4, wb2
        Q = coef_fu(1);
        b_1 = coef_fu(6); b_3 = coef_fu(7); b_2 = coef_fu(8); a_1 = coef_fu(9); a_2 = coef_fu(10)
        Q_sq3 = coefs(1);         wb2 = coefs(2);         ta4 = coefs(3);         a1_sq = coefs(4); 
        qs = Q_sq3/w
        fu =  Q   -   a1_sq*(qs)**b_1 * (w**2 - wb2) * ((qs/a_2)**b_2 + 1)**b_3 * ta4
    end function
    
    real function fu2(compos, q_f, coefs2, coef_fu2, fx_1_or_2, print_cnd)
    !fu2 q_f from sympy_v_law_contin.py/floodplain_ok
        logical, intent(in) :: print_cnd, compos
        integer, intent(in) :: fx_1_or_2
        real, intent(in) :: q_f, coefs2(10), coef_fu2(16)
        real :: sq3c, sq3f, Acd, qcs, qcs_a2, t1, t2, t3, qfs, sqc, sqf
        real :: Q, b_1, b_3, b_2, a_1, a_2, cw, cQ, a2_sq3f, Acd_a1_sqc, Drp, wc1b1, wf
        Q = coef_fu2(1); b_1 = coef_fu2(8); b_3 = coef_fu2(9); b_2 = coef_fu2(10); 
        a_1 = coef_fu2(11); a_2 = coef_fu2(12)
        Acd_a1_sqc = coefs2(1); wc1b1 = coefs2(2); Drp = coefs2(3); 
        a2_sq3f = coefs2(4); cQ = coefs2(5); cw = coefs2(6); wf = coefs2(7); sq3f = coefs2(8)
        sqf = coefs2(9); sqc = coefs2(10)
        qcs = cQ - cw*q_f
        qcs_a2 = qcs/a_2
        if (.not. compos) then    !vc>vf
            t1 = Acd_a1_sqc * qcs**b_1 * (qcs_a2**b_2 + 1)**b_3
            t2 = q_f* wc1b1 * Drp * (Q/q_f - wf)**b_1 * &
                ((qcs_a2**b_2 + 1) / ((q_f/a2_sq3f)**b_2 + 1))**b_3
            t3 = q_f*wf
            if (print_cnd) print*, 'in fu2; a2_sq3f: ', a2_sq3f        
            if (print_cnd) print*, 'in fu2; fx_1_or_2,qcs,t1,t2,t3: ', fx_1_or_2,qcs,t1,t2,t3
            fu2 = Q-t1-t2-t3
        else !vc>vf
            t1 = sqc * qcs**b_1 * (qcs_a2**b_2 + 1)**b_3
            qfs = q_f/sq3f
            t2 = sqf * qfs**b_1 * ( (qfs/a_2)**b_2 + 1 )**b_3
            fu2 = t1-t2
        end if            
    end function
    
    real function fu3(q, coefs3, coef_fu3)
        real, intent(in) :: q, coefs3(8), coef_fu3(14)
        real :: ta2, ta22, wb2, cp, dx, dt, qv, v
        ta2 = coefs3(1); ta22 = coefs3(2); wb2 = coefs3(3); cp = coefs3(4); dx = coefs3(5); dt = coefs3(6)
        v = fv2(q, coef_fu3(8:12), coefs3(7), coefs3(8))
        qv = q/v
        fu3 =  q*( ta2*qv + sqrt( ta22*qv*qv + wb2 ) ) - cp*v / (dx + v*dt)
    end function
    
    real function fu4(qf, coefs4, coef_fu4, fx_1_or_2, print_cnd)
    !7may20; see eq (V) phdSCnotes p194c.
        logical, intent(in) :: print_cnd
        integer, intent(in) :: fx_1_or_2
        real, intent(in) :: qf, coefs4(14), coef_fu4(17)
        real :: sq3c, sq3f, Acd, qcs, qcs_a2, qfs, sqc, sqf, coefrix(5)
        real :: Q, b_1, b_3, b_2, a_2, cw, cQ, a2_sq3f, Acd_a1_sqc, Drp, wcb1, wf, wc, vf, vc, vc_vf
        real :: vc_a1, vf_a1, cp3, cp4, wc_sq3c, sum_vA, Qc_hg, Qc_hc, Qc_hcd, Qc_hcu, Qf_hg
        coefrix = coef_fu4(1:5)
        b_1 = coefrix(1); b_3 = coefrix(2); b_2 = coefrix(3); a_2 = coefrix(5)
        cp3 = coefs4(1); cp4 = coefs4(2); 
        Drp = coefs4(3); wf = coefs4(4); wcb1 = coefs4(5); cw = coefs4(6); a2_sq3f = coefs4(7);
        sqc = coefs4(8); sq3c = coefs4(9); sqf = coefs4(10); sq3f = coefs4(11); wc = coefs4(12)
        wc_sq3c = coefs4(13); Acd = coefs4(14)
        vf = fv2(qf, coefrix, sqf, sq3f)
        Q = fQ(cp3, cp4, qf/vf)
        cQ = Q/wc_sq3c
        qcs = cQ - cw*qf
        qcs_a2 = qcs/a_2
        
!        if (.not. compos) then    !vc>vf
        Qf_hg = qf*wf
        Qc_hg = Q-Qf_hg
        vc = fv2(Qc_hg/wc, coefrix, sqc, sq3c)
        vc_vf = fbar(Drp, Q, qf, wf, b_1, wcb1, qcs_a2, b_2, a2_sq3f, b_3) 
            !4jun20; for fbar theory see p144d, today was checked to result same number as vf/vc in f4sln.
            !7may20; 10:40: fbar is too complex and dont work; why OK for flumeFlood?
!            vc_vf = vc/vf
        Qc_hcu = qf*vc_vf*wc
        Qc_hcd = vc*Acd
        Qc_hc = Qc_hcu + Qc_hcd
        sum_vA = Qc_hc + Qf_hg
!            sum_vA = vc*Acd + qf*( vc_vf*wc + wf )      
        fu4 = Q - sum_vA
!        print*, '--'
!        if (print_cnd) print*, 'in fu4 NoCompos; : sqc,sqf,Q,Qc_hcu,Qc_hcd,Qf_hg:', &
!            sqc,sqf,Q,Qc_hcu,Qc_hcd,Qf_hg
        if (print_cnd) print*, 'in fu4 NoCompos; : vc,vf,vc_vf_these,vc_vf_fbar:', vc,vf,vc/vf,vc_vf
        if (print_cnd) print*, 'in fu4 NoCompos; fx_1_or_2,Q,Qc_hcd,Qc_hcu,Qf_hg,sum_vA: ', &
            fx_1_or_2,Q,Qc_hcd,Qc_hcu,Qf_hg,sum_vA
!        else !vc=vf
!            vc_a1 = sqc * qcs**b_1 * (qcs_a2**b_2 + 1)**b_3
!            qfs = qf/sq3f
!            vf_a1 = sqf * qfs**b_1 * ( (qfs/a_2)**b_2 + 1 )**b_3
!            if (print_cnd) print*, 'in fu4 Compos; Q,qf,qcs,vc_a1,vf_a1: ', Q,qf,qcs,vc_a1,vf_a1
!            if (print_cnd) print*, 'only if finer fld, sq3c>sq3f, then vf can exceed vc: ', sq3c,sq3f
!            fu4 = vc_a1 - vf_a1
!        end if
    end function
    
    real function fu5(hfld, coefs5, print_cnd, rockw) 
        !ferguson07 VPE friction law for composite channel (equiv Darcy f). See p2019_207u.
        logical, intent(in) :: print_cnd, rockw
        real, intent(in) :: hfld, coefs5(13)
        real hmcd, a1_2, a2_2, a12_2, wf, Pc, Acd, wv, sg8, cpp, dx_dt, Df, Dc
        real hmc, ff, fc, Pf, P, fe, A, R, v_hc, v_hg
        hmcd = coefs5(1); a1_2 = coefs5(2); a2_2 = coefs5(3); a12_2 = coefs5(4); wf = coefs5(5)
        Pc = coefs5(6); Acd = coefs5(7); wv = coefs5(8); sg8 = coefs5(9); cpp = coefs5(10); 
        dx_dt = coefs5(11); Df = coefs5(12); Dc = coefs5(13)
        hmc = hfld + hmcd
        if (rockw) then
            ff = .3 !21may20; p228. Ferguson17. Bedrock wall friction.
        else
            ff = f_dcy(a1_2, a2_2, a12_2, Df/hfld)
        endif
        fc = f_dcy(a1_2, a2_2, a12_2, Dc/hmc)
        Pf = wf + hfld
        P = Pc + Pf
        fe = (Pc*fc + Pf*ff) / P !dedux p155, analogous to horton law n_eq ~ 2/3 powers ponder.
        A = Acd + wv*hfld
        R = A/P
        v_hc = fv_dcy(sg8, fe, R)
        v_hg = cpp/A - dx_dt
        fu5 = v_hg - v_hc
        print*, 'in fu5 Compos; : hfld,fe,R,v_hc,v_hg,fu5', hfld,fe,R,v_hc,v_hg,fu5 !if (print_cnd) 
    end function

    subroutine find_root_fu(print_cnd, xinit, coefs, coef_fu, tol, maxiter, result, success )  !f, 
        logical, intent(in) :: print_cnd    
        real, intent(in)     :: xinit, coef_fu(12), coefs(4), tol
        integer, intent(in)  :: maxiter
        real, intent(out)    :: result
        logical, intent(out) :: success
        real                 :: eps = 1e-2 !default: 1.0e-4
        real                 :: fx1, fx2, fprime, x, xnew
        integer              :: i
        result  = 0.0
        success = .false.
        x = xinit
        if (print_cnd) print*, 'in find_root_fu, to start Newton iters'
        do i = 1,max(1,maxiter)
            fx1    = fu(x, coefs, coef_fu)
            fx2    = fu(x+eps, coefs, coef_fu)
            if (print_cnd) write(*,*) i, x, fx1, fx2, eps
            fprime = (fx2 - fx1) / eps
            xnew   = x - fx1 / fprime
            if ( abs(xnew-x) <= tol ) then
                success = .true.
                result  = xnew
                exit
            endif
            x = xnew
         enddo
    end subroutine

    subroutine find_root_fu2(compos, print_cnd, xinit, coefs2, coef_fu2, tol, maxiter, result, success, & 
        xnew_neg)  !f, 
        logical, intent(in) :: print_cnd, compos
        real, intent(in)     :: xinit, coefs2(10), coef_fu2(16), tol
        integer, intent(in)  :: maxiter
        real, intent(out)    :: result
        logical, intent(out) :: success, xnew_neg
        real                 :: eps_fact = 1e-5 !units, if no composite channel: 
            ![f2]=m3/s, [x] = m2/s; so [fprime=df/eps]=m; 
            !so eqs=d(qf) & [eps]=m2/s
        real                 :: eps  ! = 1e-2 !default: 1.0e-4
        real                 :: fx1, fx2, fprime, x, xnew, dx
        integer              :: i
        result  = 0.0
        success = .false.
        xnew_neg = .false.
        x = xinit
        if (print_cnd) print*, 'in find_root_fu2, to start Newton iters'
        do i = 1,max(1,maxiter)
            eps = eps_fact*x
            fx1    = fu2(compos, x, coefs2, coef_fu2, 1, print_cnd)
            fx2    = fu2(compos, x+eps, coefs2, coef_fu2, 2, print_cnd)
            fprime = (fx2 - fx1) / eps
            dx = - fx1 / fprime
            xnew = x + dx
            if (xnew < 0) then            
!                print*, '29ago20; in xnew_neg for f2 flume fld'
!                stop                
                xnew_neg = .true.
                success = .true.
                if (print_cnd) print*, '29ago20; success. Flood was solved by setting min qf, as f2 root xnew_neg.'
                result  = 1e-6
                exit
            endif                
            xnew = max(1e-6, xnew)  !qf > 1e-2 l/s/m !!
                !29ago2020: to be consistent with seed (min asumed) qf_try0 value for fu2 before called by 1st time.
            if (print_cnd) print*, 'i,x,fx1,fx2,eps,fprime,dx,xnew', i,x,fx1,fx2,eps,fprime,dx,xnew
            if (print_cnd) print*, '---next iter'            
            if ( abs(dx)/x <= tol ) then
                success = .true.
                if (print_cnd) print*, 'success; it means flood was solved'
                result  = xnew
                exit
            endif
            x = xnew
         enddo
    end subroutine
    
    subroutine find_root_fu3(print_cnd, xinit, coefs3, coef_fu3, tol, maxiter, result, success )  !f, 
        logical, intent(in) :: print_cnd    
        real, intent(in)     :: xinit, coef_fu3(14), coefs3(8), tol
        integer, intent(in)  :: maxiter
        real, intent(out)    :: result
        logical, intent(out) :: success
        real                 :: eps = 1e-2 !default: 1.0e-4
        real                 :: fx1, fx2, fprime, x, xnew
        integer              :: i
        result  = 0.0
        success = .false.
        x = xinit
        if (print_cnd) print*, 'in find_root_fu3, to start Newton iters with xinit: ', x
        do i = 1,max(1,maxiter)
            fx1    = fu3(x, coefs3, coef_fu3)
            fx2    = fu3(x+eps, coefs3, coef_fu3)
            if (print_cnd) write(*,*) i, x, fx1, fx2, eps
            fprime = (fx2 - fx1) / eps
            xnew   = x - fx1 / fprime
            if ( abs(xnew-x) <= tol ) then
                success = .true.
                result  = xnew
                exit
            endif
            x = xnew
         enddo
    end subroutine
    
    subroutine find_root_fu4(qfmin, print_cnd, xinit, coefs4, coef_fu4, errQadm_j, &
        tol, maxiter, result, &
        success, fff4)
        logical, intent(in) :: print_cnd
        real, intent(in)     :: xinit, coefs4(14), coef_fu4(17), tol, errQadm_j
        real, intent(inout)  :: qfmin
        integer, intent(in)  :: maxiter
        real, intent(out)    :: result, fff4
        logical, intent(out) :: success
!        real                 :: eps_fact = 1e-5 !portion of current x to set dx for calc derivative fprime=df/dx.
        real                 :: eps  ! = 1e-2 !default: 1.0e-4
        real                 :: fx1, fx2, fprime, x, xnew, dx, fx1_prev, x_prev, a, b, c, fc
        integer              :: i, j
        logical              :: bisec, succ_cnd, cnd1, cnd2
        result  = 0.0; fff4=0; bisec = .false.
        success = .false.
        x = xinit
        eps = qfmin
        if (print_cnd) print*, 'in find_root_fu2, to start Newton iters'
        fx1_prev = 0
        do i = 1,max(1,maxiter)
!            eps = eps_fact*x
            fx1    = fu4(x, coefs4, coef_fu4, 1, print_cnd)
            fx2    = fu4(x+eps, coefs4, coef_fu4, 2, print_cnd)
            if (fx1 .eq. fx2) then !6jun20; curve lacks slope, p243cl
                eps=10*eps !coarsen dx. Assume not many iters are needed from here.
                fx1    = fu4(x, coefs4, coef_fu4, 1, print_cnd)
                fx2    = fu4(x+eps, coefs4, coef_fu4, 2, print_cnd)                
            endif
            if( (x .eq. qfmin) .and. (fx1 .ne. fx1) )  qfmin=.1*qfmin !21may20; p226u(a).
            cnd1 = (fx1_prev .gt. 0) .and. ( (fx1 .ne. fx1) .or. (fx2 .ne. fx2) )
            cnd2 = (fx1 .gt. 0) .and. (fx2 .ne. fx2)  !9jun20; p253c; I know df4/dx<0; NoFld&geom indics too low qf;
                !so f4(x1=qfmin)>0, f4(x2=qfmin+eps)=NaN
            if ( cnd1 .or. cnd2 ) then 
                !if f+ goes to NaN without f- i.e. no RootCut shown.
                bisec = .true.
                if (print_cnd) print*,'  need BISEC'
                exit !go to bisection, knowing f4 decreases with x.
            endif
            if (x .eq. qfmin .and. fx1 .lt. 0) then
                !10may20; wvmx=1km and wc near 0 carry only 1l/s if qfmin=1e-6. 
                !f4 always decrease in range (qfmin,qf_root)
                result = x;  fff4 = fx1;  exit
            endif
            fprime = (fx2 - fx1) / eps
            dx = - fx1 / fprime
            xnew   = x + dx
!            xnew = max(qfmin*(1+i*1./maxiter), xnew)  !7may20; variable xmin to avoid circular iteration.
!            if ((fx1_prev .gt. 0) .and. (fx1 /= fx1)) then !2nd conditions is NaN check. See p204u.
                !make rougher secant as tangent approx, to reduce chance of projection out of valid domain.
                !I know f4 falls with qf and seems that proyects to f4 NaN due to f4 downward concave.
!                stop '12may20; coarse secant'
!                eps = eps*10
!                xnew = x
!            else
            xnew = max(qfmin, xnew)  !dont bother circular iteration; I need to get when f4(qfmin)<0.            
!                fx1_prev = fx1 !if coarsening secant avoid update fx1_prev to preserve noNaN positive val.
!            endif
            if (print_cnd) print*, 'i,x,fx1,fx2,eps,fprime,dx,xnew', i,x,fx1,fx2,eps,fprime,dx,xnew
            if (print_cnd) print*, '---next iter'
            succ_cnd = ( ( abs(dx)/x .lt. tol ) .and. ( abs(fx1) .lt. errQadm_j ) )  .or. (fx1*fx2 .lt. 0) 
                !2nd cond; 6jun20; p244ur; solution between x and x + eps. 
            if (succ_cnd) then !4jun20; and instead or.
                !12may20; vertical check is set for cases when even with fine dqf next qf to small +f4 is NaN.                
!                print*, '4jun20; SUCCESS COMES; tol,abs(dx)/x,abs(fx1),errQadm_j', &
!                    tol,abs(dx)/x,abs(fx1),errQadm_j                
                success = .true.
                if (print_cnd) print*, 'success; it means flood was solved'
                result  = x !p226dl; deprec: xnew
                exit
            endif
            fx1_prev = fx1
            x_prev = x
            x = xnew
        enddo
        if (bisec) then
            if (cnd2) then
                a = x ; b = x+eps
            else
                a = x_prev; b = x
            endif

            do i = 1, maxiter
                c = (a+b)/2
                fc = fu4(c, coefs4, coef_fu4, 0, print_cnd)
                if (print_cnd) print*, '  f4 BISEC; i,a,b,c,fc:', i,a,b,c,fc !15may20
                if (print_cnd) print*, '---next iter'
                if ( abs(fc) .lt. errQadm_j ) then
                    success = .true.
                    if (print_cnd) print*, ' f4 BISEC success: flood was solved'
                    result  = c
                    exit                    
                endif
                if ((fc .lt. 0) .or. (fc .ne. fc)) then !choose which half (interval) for next iter.
                    b = c !left.
                else
                    a = c !right.
                endif
            enddo
        endif
    end subroutine
    
    subroutine find_root_fu5(print_cnd, rockw, xinit, coefs5, tol, maxiter, result, success )
        logical, intent(in) :: print_cnd, rockw
        real, intent(in)     :: xinit, coefs5(13), tol
        integer, intent(in)  :: maxiter
        real, intent(out)    :: result
        logical, intent(out) :: success
        real                 :: eps = 1e-2 !default: 1.0e-4
        real                 :: fx1, fx2, fprime, x, xnew, dx, fx1_prev, x_prev, a, b, c, fc
        integer              :: i
        logical              :: bisec        
        result  = 0.0
        success = .false.
        x = xinit
        if (print_cnd) print*, 'in find_root_fu3, to start Newton iters with xinit: ', x
        
        print*, '  solving COMPOS F5; i,x,fx1,fx2,eps,fprime,dx,xnew', i,x,fx1,fx2,eps,fprime,dx,xnew
        print*, '---next iter'
        do i = 1,max(1,maxiter)
            fx1    = fu5(x, coefs5, print_cnd, rockw)
            fx2    = fu5(x+eps, coefs5, print_cnd, rockw)
            if (print_cnd) write(*,*) i, x, fx1, fx2, eps
            
            if ( (fx1_prev .gt. 0) .and. ( (fx1 .ne. fx1) .or. (fx2 .ne. fx2) ) ) then 
                !if f+ goes to NaN without f- i.e. no RootCut shown.
                bisec = .true.
                if (print_cnd) print*,'  f5_compos needs BISEC'                   
                exit !go to bisection, knowing f4 decreases with x.
            endif            
            
            fprime = (fx2 - fx1) / eps
            dx = - fx1 / fprime
            xnew   = x + dx
            if ( abs(dx/x) <= tol ) then
                success = .true.
                result  = xnew
                exit
            endif
            x = xnew
            fx1_prev = fx1
            x_prev = x            
        enddo
        if (bisec) then
            a = x_prev; b = x
            do i = 1, maxiter
                c = (a+b)/2
                fc = fu5(c, coefs5, print_cnd, rockw)
                if (print_cnd) print*, '  f5 BISEC; i,a,b,c,fc:', i,a,b,c,fc !15may20
                if (print_cnd) print*, '---next iter'
                if ( abs(fc) .lt. .001 ) then !tol: 1mm/s; so low as some floods are reaaally shallow.
                    success = .true.
                    if (print_cnd) print*, ' f5 BISEC success: flood was solved'
                    result  = c
                    exit                    
                endif
                if ((fc .lt. 0) .or. (fc .ne. fc)) then !choose which half (interval) for next iter.
                    b = c !left.
                else
                    a = c !right.
                endif
            enddo
        endif                 
    end subroutine
    
    subroutine solve_flow_flume(i1,j1,print_cnd, hprev, wc, h_bf, coef_fu, w, h, v, Fr, aspect, sqc, sq3c, isFlood)
        logical, intent(in) :: print_cnd
        integer, intent(in) :: i1,j1
        real, intent(in) :: hprev, wc, h_bf, coef_fu(12)
        real, intent(out) :: w, h, v, Fr, aspect, sqc, sq3c
        logical, intent(out) :: isFlood
        real           :: coefs(4), w_try0, w_try, wb, root, spiral_ampl, Q, g, ta, retry_fact, hh, A
        real           :: a_1, D_84, Q_sq3, S, wb2, ta4, a1_sq
        integer :: it, itmx, spiral_dir
        integer        :: maxiter = 20  !for internal loop of newton solver, given initial value.
        logical        :: success, condFr

        D_84 = coef_fu(2); S = coef_fu(3); g = coef_fu(4)
        wb = coef_fu(5); Q = coef_fu(1); a_1 = coef_fu(9); ta = coef_fu(11); ta4 = coef_fu(12)
        sqc = sqrt(D_84*S*g);  sq3c = sqrt(D_84**3*S*g);  wb2 = wb*wb
        Q_sq3 = Q/sq3c;  a1_sq = a_1*sqc
        coefs = (/Q_sq3, wb2, ta4, a1_sq/)
        w_try0 = min(wb + 2*hprev/ta, wc)
        it = 0 !; spiral_ampl = .3
        retry_fact = .2
        w_try = w_try0
        condFr = .false.
        if (print_cnd) print*, 'IN solve_flow_flume, condFr=', condFr, '; wc,w_try0=', wc,w_try0
        itmx = maxiter !max # of iters changing initial values spirally. 
            !This is an external loop, which repeats if num sln fail condFr.
        success = .false.            
        do while (.not.condFr)
            it = it + 1
            if (it>itmx) then
                print*, 'take this t,x to edit print_range for DEBUG: i,j:', i1,j1
                print*, 'itmx=', itmx
                stop 'In NOflood, tried many w_try, but now it>itmx'
            endif
            call find_root_fu(print_cnd, w_try, coefs, coef_fu, 1.0e-3, maxiter, w, success)
            v = fv(Q/w, coef_fu, sqc, sq3c)
            A = Q/v
            hh = A/w !hydraulic depth
            Fr = Q/A/sqrt(g*hh)
            aspect = w/hh
            condFr = (Fr.lt.3 .and. v.gt.0)
!            if (mod(it, 2).eq.0) then
!                spiral_dir = 1
!            else
!                spiral_dir = -1
!            endif
            if (.not.success .or. .not.condFr) then
                w_try = (1 + retry_fact)*w_try
!                retry_fact = maxin(1, 1 + spiral_dir*spiral_ampl*it) 
!                w_try = retry_fact * w_try0  !make spiral around initial guess.
            endif                    
        end do
        h = .5*(w - wb)*ta
        isFlood = h .gt. 1.05*h_bf
        if (h .gt. h_bf .and. h .lt. 1.05*h_bf) then            !velocity update if 'threshold' flow
            !29ago20; where 1.03 repaired to 1.05.
            h = .99*h_bf
            v = fv(Q/wc, coef_fu, sqc, sq3c)
            A = Q/v
            hh = A/wc   !hydraulic depth
            Fr = Q/A/sqrt(g*hh)
            aspect = wc/hh
            condFr = (Fr.lt.3 .and. v.gt.0)                                    
            if (print_cnd) print*, '    in NOflood, entered to update v just below hbf'
        endif
        if (print_cnd) print*, '    final NOflood sln for this x,t; h=', h, 'w=', w, 'v=', v, 'Fr=', Fr
        if (print_cnd .and. isFlood) print*, '....but isFlood with really 5%tol, as h,h_bf: ',  h,h_bf
    end subroutine    
    
    subroutine solve_flood(which_flume, ta_inv_sqd, i1,j1, nR, print_cnd, coef_fu, coef_fu2, h, vc, vf, Frc, aspect, w, sq3f, Q_f)
        logical, intent(in) :: print_cnd
        integer, intent(in) :: i1,j1, nR, which_flume
        real, intent(in) :: coef_fu(12), coef_fu2(16), ta_inv_sqd
        real, intent(out) :: h, vc, vf, Frc, aspect, w, sq3f, Q_f
        logical :: condFr, success, compos, xnew_neg
        integer :: it, itmx     !for internal loop of newton solver, given initial value.        
        real :: qf_try0, coefs2(10), b_1, a_1, a_2, a2_sq3f, Ac, Acd, Acd_a1_sqc, Acu, Af, qfmx
        real :: cw, cQ, D84f, D84c, Drp, Frf, g, hbf, Q, Qc, qf, qf_try, Sc, Sf, dqf
        real :: sq3c, sqf, sqc, ta, wc, wc1b1, wf, wv, h_fld
!        real           :: w_try0, w_try, wb, root, spiral_ampl, Q, g, ta, retry_fact, hh, A
        Q = coef_fu2(1); D84f = coef_fu2(2); hbf = coef_fu2(4); 
        g = coef_fu2(5); wc = coef_fu2(6); wv = coef_fu2(7); a_2 = coef_fu2(12); b_1 = coef_fu2(8)
        a_1 = coef_fu2(11); ta = coef_fu2(13); sqc = coef_fu2(14); sq3c = coef_fu2(15)
        D84c = coef_fu2(16); Sc = coef_fu2(3)
        sq3f = sqrt(D84f**3*Sc*g);         
        if (j1 .lt. nR) then
            Acd = hbf*(wc-hbf/ta)
            sqf = sqrt(D84f*Sc*g)    !proust09wrr says typically S is shared for both ch and fl zones.
            wf = wv-wc
            cQ = Q/(wc*sq3c); cw = wf/(wc*sq3c)
            Drp = (D84f/D84c)**(1.5*b_1 - 0.5); wc1b1 = wc**(1-b_1); !Dratio grouped by 4jun20.
            Acd_a1_sqc = Acd*a_1*sqc; a2_sq3f = a_2*sq3f
            coefs2 = (/Acd_a1_sqc, wc1b1, Drp, a2_sq3f, cQ, cw, wf, sq3f, sqf, sqc/)
            
    !        vcp_vfp = 3 !assumed ratio of velocities in channel and floodplain.
    !        vcp = wc/wv*Q/Acd !no numerator in eqn for qf_try0 is greater than 0; the lower wf the more water in channel hbf.
    !        qf_try0 = (Q - vcp*Acd)/(vcp_vfp*wc + wf) !see deduction in phdSC notes p152.       
            itmx = 20  !max # of iters changing initial values for Newton numSolver. 
            qf_try0 = 1e-6 !min
!            if (wf/wv<.1) then !this cond is for basin
            qfmx = Q/wf
!            else !8jun20; pp249. Narrow valley can convey Q>Qnf; nf=noFlood (up-projected trapeze crosssection.)
!                qfmx = v2
!            endif
            dqf = (qfmx - qf_try0) / itmx             
            qf_try = qf_try0
            condFr = .false.
            xnew_neg = .false.
            compos = .false. !see phdSCnotes p156: when no more momentum is transferable form ch to fl as vf>vc.
            it = 0        
            if (print_cnd) print*, 'IN solve_flood, qf_try0,qfmx,itmx,dqf=', qf_try0,qfmx,itmx,dqf
                !This is an external loop, which repeats if num sln fail condFr.
            do while (.not.condFr .and. .not.xnew_neg)
                it = it + 1
                if (it>itmx) then
                    print*, 'take this t,x to edit print_range for DEBUG: i,j:', i1,j1
                    print*, 'itmx=', itmx
                    stop 'In flood, tried many qf_try, but now it>itmx'
                endif
                if (print_cnd) print*, 'in solve_flood; prev to find_root_fu2, it,qf_try:', it,qf_try
                call find_root_fu2(compos, print_cnd, qf_try, coefs2, coef_fu2, 1e-2, itmx, qf, success, &
                    xnew_neg)
!                if (xnew_neg) then
!                    print*, 'qf=', qf
!                    stop !29ago20; to see if qf is min.                    
!                endif
                Q_f = qf*wf
                Qc = Q - Q_f
                vc = fv(Qc/wc, coef_fu, sqc, sq3c)
                vf = fv(qf, coef_fu, sqf, sq3f)
                h_fld = qf/vf
                Acu = wc*h_fld
                Ac = Acd + Acu
                Frc = vc/sqrt(g*Ac/wc)
                Af = wf*h_fld
                Frf = vf/sqrt(g*Af/wf)
                if (.not.compos  .and.  vf .gt. vc) compos = .true. 
                    !13may20; note basin flood has vpe approach which is the valid one. 
                    !This compos code is deprecated, as calibrated pitlick13flood never shows vf>vc
                if (print_cnd) print*, 'in solve_flood (29ago20), after found root_fu2; Frc,Frf,vc,vf,compos: ', &
                    Frc,Frf,vc,vf,compos
                condFr = (Frc.lt.2 .and. Frf.lt.2 .and. vf.gt.0 .and. vf .le. 1.01*vc) 
                    !vf.le.vc as if compos then they are equal. 5% tol due to error solving composite flow.
                if (.not.condFr) then
                    qf_try = qf_try + dqf
                    if (qf_try .gt. qfmx) then
                        print*, 'take this t,x to edit print_range for dqf DEBUG: i,j:', i1,j1
                        print*, &
                            '15apr20; STOP, as all qf domain was explored in steps dqf=', dqf
                        stop
                    endif                    
                endif                    
            end do
            
            Af = wf*h_fld
            aspect = wv**2/(Af + Ac)
            h = hbf + h_fld
            w = wv
            if (print_cnd) print*, '  FLOODsln_j<nR: h=', h, ', h_fld=', h_fld, &
                ', vc=', vc, ', vf=', vf, ', Frc=', Frc, ', Frf=', Frf, ', w=', w, 'Qc,Q_f:', Qc,Q_f
        else    !critical flow at outlet with rigid walls.
            if (which_flume .eq. 2) then !27ago2020; pitlick13flood experiment; with rectangle crosssection.
                !see Pitlick email to me, about boundary condition
                w=wv; h=(Q/(w*sqrt(g)))**.67
            else !ubc trench; triangular outlet
                h=( 2*Q**2 / (g*ta_inv_sqd) )**.2
                w= 2*h/ta  
            endif
            vc=sqrt(g*h); vf=0; Frc=1; Frf=0; aspect=w/h; Q_f=0; h_fld =0
            if (print_cnd) print*, '  FLOODsln_j=nR: h=', h, ', h_fld=', h_fld, &
                ', vc=', vc, ', vf=', vf, ', Frc=', Frc, ', Frf=', Frf, ', w=', w, 'Q=Qc,Q_f:', Q,Q_f
        endif
    end subroutine
    
    subroutine solve_flow_basin(minifld, i1,j1,print_cnd, ta,qtry0, wc, h_bf, coef_fu3, w, h, v, Fr, aspect, sqc, &
        sq3c, &
        Q, Qf_try0, A, isFlood, prt_cd)
        logical, intent(in) :: print_cnd
        integer, intent(in) :: i1,j1
        real, intent(in) :: ta, qtry0, wc, h_bf, coef_fu3(14)
        integer, intent(inout) :: minifld !if it is 'on'=1, it is turned here to 'off'=0 as it's unable to go to fld.
        integer, intent(out) :: prt_cd
        real, intent(out) :: w, h, v, Fr, aspect, sqc, sq3c, Q, Qf_try0, A
        logical, intent(out) :: isFlood
        real           :: coefs3(8), qtry, wb, root, spiral_ampl, g, retry_fact, hh, qv, dq
        real           :: D_84, S, wb2, ta2, cp, dt, dx, qmx, ta22
        integer :: it, itmx, spiral_dir
        integer        :: maxiter = 20  !for internal loop of newton solver, given initial value.
        logical        :: success, condFr, vpos, fld
        Q = 0
        prt_cd = 0 
        cp = coef_fu3(1);  dx = coef_fu3(2);  dt = coef_fu3(3)
        D_84 = coef_fu3(4); S = coef_fu3(5); g = coef_fu3(6); wb = coef_fu3(7); 
        ta2 = coef_fu3(13); ta22 = coef_fu3(14)
        sqc = sqrt(D_84*S*g);  sq3c = sqrt(D_84**3*S*g);  wb2 = wb*wb
        coefs3 = (/ta2, ta22, wb2, cp, dx, dt, sqc, sq3c/)
!        print*, '7may20; coefs3: ', coefs3
        it = 0 !; spiral_ampl = .3
        retry_fact = .2
        qtry = qtry0
        qmx = 10 ![m3/s/m]
        itmx = maxiter !max # of iters changing initial values spirally. 
        dq = (qmx - qtry0) / itmx         
        condFr = .false.
        if (print_cnd) print*, 'IN solve_flow_basin, h_bf,condFr:', h_bf,condFr, '; qtry0=', qtry0
            !This is an external loop, which repeats if num sln fail condFr.
        success = .false.
        fld = .false.; vpos = .false.
        do while (.not.condFr)
            it = it + 1
            if (it>itmx) then
                print*, 'take this t,x to edit solve_flow_basin DEBUG: i,j:', i1,j1
                print*, 'itmx=', itmx
                prt_cd = 1
                exit                
!                stop 'In NOflood, tried many w_try, but now it>itmx'
            endif
            call find_root_fu3(print_cnd, qtry, coefs3, coef_fu3, 1.0e-3, maxiter, q, success)           
            v = fv2(q, coef_fu3(8:12), sqc, sq3c)
            qv = q/v
            w = ta2*qv + sqrt( ta22*qv*qv + wb2 ) !phdSCnotes p194u.
            Q = q*w
            A = Q/v
            hh = A/w !hydraulic depth
            Fr = Q/A/sqrt(g*hh)
            h = .5*(w - wb)*ta
            vpos = v .gt. 0
            if (h .gt. 1.03*h_bf .and. vpos) then
                fld = .true.
                exit
            endif
            condFr = (Fr.lt.5 .and. vpos)            
            if (.not.success .or. .not.condFr) then
                qtry = qtry + dq
            endif                    
        end do
        aspect = w/hh
        isFlood = ( fld .and. (minifld .lt. 1) .and. (prt_cd .lt. 1) )
            !dont recalc flood, as qf is negligible: minifld. 
            !Go out from this subrout and avoid floodCalcs
                !where prt_cd can be edited, to be able to find the bug which set prt_cd=1 in this subrout.
        if (minifld .eq. 1) minifld = 0
        if (isFlood) Qf_try0 = (h-h_bf)*v*.5*(w+wc)  !p197ur
        if (h .gt. h_bf .and. h .lt. 1.03*h_bf) then            !velocity update if 'threshold' flow
            !21may20; water balance error remains to be evaluated here. 
                !That's  why I prefer to lower qftry_min in solving f4.
            h = .99*h_bf
            v = fv2(Q/wc, coef_fu3(8:12), sqc, sq3c)
            A = Q/v
            hh = A/wc   !hydraulic depth
            Fr = Q/A/sqrt(g*hh)
            aspect = wc/hh
            condFr = (Fr.lt.5 .and. v.gt.0)                                    
            if (print_cnd) print*, '    in NOflood, entered to update v just below hbf'
        endif
        if (print_cnd) print*, '    final NOflood sln for this x,t; Q,h=', Q,h, 'w,wc,wb=', w,wc,wb, 'v=', v, 'Fr=', Fr
!        if (minifld .eq. 1) then
!            print*, 'i1,j1:', i1,j1, 'coge minifld pa evitar q diagnostiq flood q ya dio qf<<1e-6~0'
!            print*, 'Q,h,h_bf:', Q,h,h_bf
!            stop 'stop to celebr!'
!        endif
        if (print_cnd .and. isFlood) print*, '....but isFlood with 3%tol, as h,h_bf,Qf_try0: ',  h,h_bf,Qf_try0
    end subroutine          

    subroutine solve_flood_basin(i1,j1, nR, print_cnd, hnf, wnf, vnf, sg8, Qup, Qf_try0nf, Acd, dx_dt, cpp, sa, &
        a1_2, a2_2, a12_2, coef_fu4, errQadm_j, &
        h, vc, vf, Frc, aspect, w, &
        sq3f, Q_f, Q, Aw, prt_cd, minifld, repeater)
        logical, intent(in) :: print_cnd
        integer, intent(in) :: i1, nR, j1
        real, intent(in) :: Qup, Qf_try0nf, Acd, dx_dt, cpp, sa, coef_fu4(17), errQadm_j
        real, intent(in) :: a1_2, a2_2, a12_2, sg8, hnf, wnf, vnf
        integer, intent(inout) :: repeater
        real, intent(out) :: h, vc, vf, Frc, aspect, w, sq3f, Q_f, Q, Aw
        integer, intent(out) ::prt_cd, minifld
        logical :: condFr, success, compos, rockw, persist
        integer :: it, itmx     !for internal loop of newton solver, given initial value.        
        real :: coefrix(5), coefs4(14), coefs5(13)
        real :: b_1, b_3, a_2, a2_sq3f, Ac, Acu, Qcd, Qcu, Af, qf_try_ini, qf_tryNF, qfmin, fff4
        real :: cw, cQ, D84f, D84c, Drp, Frf, g, hbf, Qc, q_c, qf, qf_try, qft, Sc, Sf, dqf, qfmx
        real :: sq3c, sqf, sqc, ta, wc, wc1b1, wf, wv, h_fld, cp3, cp4, wcb1, wc_sq3c, ff4, m, relerr
        real :: vm, hfld_try, hfmin, hfmx, hft, ff5, wb, dhf, vcomp, hadd, Avirt, v2, Qmx
        persist = .true.
        prt_cd = 0; fff4=0
        coefrix = coef_fu4(1:5)
        b_1 = coef_fu4(1); b_3 = coef_fu4(2); a_2 = coef_fu4(5)
        g = coef_fu4(6); wc = coef_fu4(7); wv = coef_fu4(8);  D84c = coef_fu4(9); D84f = coef_fu4(10);
        ta = coef_fu4(11); hbf = coef_fu4(12); Sc = coef_fu4(13); cp3 = coef_fu4(14); cp4 =coef_fu4(15);
        sqc = coef_fu4(16); sq3c = coef_fu4(17)
        sq3f = sqrt(D84f**3*Sc*g);
        minifld = 0
        rockw = .false.              
        if (j1 .lt. nR) then
            sqf = sqrt(D84f*Sc*g)    !proust09wrr says typically S is shared for both ch and fl zones.
            wf = wv-wc
            if (wf .eq. 0)  rockw = .true.
            wc_sq3c = wc*sq3c            
            cw = wf/wc_sq3c !cQ = Q/(wc*sq3c); #7may20; cQ must be calculated for each qf tried.
            Drp = (D84f/D84c)**(1.5*b_1 - 0.5); wcb1 = wc**b_1; 
            a2_sq3f = a_2*sq3f
            coefs4 = (/ cp3, cp4, Drp, wf, wcb1, cw, a2_sq3f, sqc, sq3c, sqf, sq3f, wc, wc_sq3c, Acd/)
            itmx = 20  !max # of iters changing initial values for Newton numSolver. 
            qf_tryNF = Qf_try0nf/wv
            qfmin = 1e-6
            qf_try = .1*qf_tryNF !9may20; assume vf~10%vc.
            qf_try_ini = qf_try            
            if ( (.not. rockw) .and. (wnf<wv) ) then
                qfmx = Q/wf !not more than ALL flow might travel via floodplain. Q is NoFloodSolution, by f90 runMemory.
            else !8jun20; pp249. Narrow valley can convey Q>Qnf; nf=noFlood (up-projected trapeze crosssection.)
                Avirt = .25*(wnf+wv)*(wnf-wv) !pair of triangles out of valley walls.
                hadd = Avirt/wv !hadd may be overstimated; if v2>vnf then wet area falls, though Q may increase by w-.
                v2 = vnf * ( hnf+hadd ) / hnf !local hydraulic geometry for v and h scales similar with Q, so v/h~const.
                Qmx = v2 * ( wv*hnf + Avirt ) !neglect alluvial banks, i.e. assumes rectangular crosssection.
                qfmx = Qmx/Wv
                print*, '8jun20; choosing qfmx pre iter in narrow valley; wnf,wv,vnf,v2,hnf,hadd:', &
                    wnf,wv,vnf,v2,hnf,hadd
                print*, ' 2nd must be larger; Qnf,Qmx,qfmx:', Q,Qmx,qfmx
                print*, ' wf/wv:', wf/wv
                print*, ' likely 2nd is larger, as 1st grows a lot (low wf/wv); qfmx(wide_wv),qfmx:', Q/wf,qfmx
!                stop
            endif            
            dqf = (qfmx - qf_try) / itmx
            condFr = .false.
            compos = .false. !see phdSCnotes p156: when no more momentum is transferable form ch to fl as vf>vc.
            it = 0
            if (print_cnd) print*, 'coefs4:', coefs4
            if (print_cnd) print*, 'IN solve_flood; wv,wc,wf: ', wv,wc,wf
            if (print_cnd) print*, 'IN solve_flood; preWhile; qf_tryNF,qf_try,qfmx,Q_NoFld,itmx,dqf: ', &
                qf_tryNF,qf_try,qfmx,Q,itmx,dqf
            if (.not.rockw) then
                do while (.not.condFr)
                    !This is an external loop, which repeats if num sln fail condFr.            
                    it = it + 1

                    if (it>itmx) then
                        print*, 'take this t,x to edit print_range for DEBUG: i,j:', i1,j1
                        print*, 'itmx=', itmx
                        print*, 'VIENE f4 for all qf: '
                        prt_cd = 1
                        persist = .false.
    !                    
    !                    !8may20: see f4 shape.
    !!                    if (i1 .eq. 64 .and. j1 .eq. 14) then
!                        m = (1e1 - 1e-2) / 1000
                        m = (qfmx - qf_try_ini) / 1000
                        do k=1,1000
                            qft = qfmin + (k-1)*m
                            ff4 = fu4(qft, coefs4, coef_fu4, 1, print_cnd)
    !                        if (print_cnd) then
                            print*, '**********************************************************'                    
                            print*, 'k,qft,f4', k,qft,ff4
    !                        endif                        
                        end do
    !!                    endif                    
    !                    
    !                    stop 'In flood, tried many qf_try, but now it>itmx'
                        exit
                    endif
                    
                    if (print_cnd) then
                        print*, 'in solve_flood; prev to find_root_fu4, it,qf_try:', it,qf_try
                        if (compos) print*, 'to seek vf reduction to vc; qfmx: ', qfmx
                    endif
                    if (persist) then
                        call find_root_fu4(qfmin, print_cnd, qf_try, coefs4, coef_fu4, errQadm_j, 1e-2, 10, &
                            qf, success, fff4)
                        vf = fv2(qf, coefrix, sqf, sq3f)                            
                    else !10jun20; maybe qf is too low to not be NaN in x range for f4(x).
                        qf=0; vf = 0
                        print*,'10jun20, NO MAS PERSIST CON qf' 
                    endif                        
                    if (.not.compos) then                    
                        Q_f = qf*wf
                        Q = fQ(cp3, cp4, qf/vf)
        !                vc_vf = fbar(Drp, fQ(cp3, cp4, qf/vf), qf, wf, b_1, wcb1, (cQ-cw*q_f)/a_2, b_2, &
        !                    a2_sq3f, b_3)
        !                Q = wc*Acd + vc_vf*qf*wc + Qf
                        Qc = Q - Q_f
                        vc = fv2(Qc/wc, coefrix, sqc, sq3c)
                    else
                                                                
                    endif                
                    h_fld = qf/vf
                    Acu = wc*h_fld
                    Ac = Acd + Acu
                    Frc = vc/sqrt(g*Ac/wc)
                    Af = wf*h_fld
                    Frf = vf/sqrt(g*Af/wf)
                    q_c = Qc/wc
                    Qcd = vc*Acd
                    Qcu = vc*Acu
                    if (.not.compos  .and.  vf .gt. vc) then 
                        !as Sf=Sc=S, h=hc>hfld and v(q,S,D): so compos (vf>vc) only if Df<Dc; see p197cl.
                        compos = .true.
    !                    vm = (vf + vc) / 2                    
                        hfmx = Q/vc-Acd !21may20; v_compos, i.e. vf=vc relates with no compos v as: >vc=vmin and <vf.
    !                    hfld_try = max(.01,(Q/vm - Acd)/wv) !max added by 21may20, to avoid neg hf deprec 'if' below.
                        if (hfmx .lt. 0) then !UNDEPREC SAME DAY. deprec 21may20; vm is too simple averaging.
                            print*, 'COMPOS DEBUG:'
                            print*, 'COMPOS: wc,wv,vf,vc,h_fld: ', wc,wv,vf,vc,h_fld
                            print*, 'COMPOS: Q, Qcd, Qcu, Qc, Q_f: ', Q, Qcd, Qcu, Qc, Q_f
                            print*, 'compare above: Qcd+Qcu vs Qc, also Q=sumQi'
                            print*, 'debe ser q no era sln; success,fff4:: ', success,fff4
                            print*, 'COMPOS: vc,Q/vc: ', vc,Q/vc, '; while Acd: ', Acd                            
                            stop ' then hfld_try<0'
                        endif
                        exit !leads to VPE composite channel solution; see phdSCnotes p207.
                    endif                    
                    if (print_cnd) print*, 'in solve_flood, after found root_fu4; qc,qf,Qc,Q_f,Q: ', &
                        q_c,qf,Qc,Q_f,Q
                    if (qf .le. 1.01*qfmin)  then !all f4 negative?, so f4 nearest to 0 is for qfmin as f4~1/qf.
                        !(fff4 .lt. 0) can occurs with any qf if newton passed root and needed to come back.
                        relerr = fff4/Q

    !                    print*, 'NOTE qf<qfmin i.e. its TOO LOW: f4= errQ= Qhg-Qhc:', fff4, &
    !                    'relerr:', relerr, 'as Q,qf:', Q,qf

                        if (abs(relerr) .gt. .1) then
    !                        prt_cd = 1
                            minifld = 1 !fld is so small than bankfull flow will be assumed to avoid Qpeak understim.
                            print*, 'take this t,x to edit print_range for this qf<qfmin DEBUG: i,j:', i1,j1
    !                        stop 'assuming qf=qfmin instead real qf<qfmin gets relerr>10%' !stop in repeater in hg.f90
                            exit
                        endif
                    endif                                        
                    if (print_cnd) print*, 'in solve_flood, after found root_fu4; Frc,Frf,vc,vf,compos: ', &
                        Frc,Frf,vc,vf,compos

                    !14may20; I think this is now solved by compositeVPE.
    !                if (q_c .lt. qf .and. sqf .gt. 1.03*sqc) then
    !                    print*, 'in solve_flood_basin; after Newton; edit print_range for DEBUG with i,j:', i1,j1
    !                    print*, 'q_c,qf,sqc,sqf:', q_c,qf,sqc,sqf
    !                    stop 'why slower in channel if its bed is not coarser?'
    !                endif
    
                    condFr = (Frc.lt.5 .and. Frf.lt.5 .and. vf.gt.0 .and. vf .le. 1.01*vc) 
                        !vf.le.vc as if compos then they are equal. 5% tol due to error solving composite flow.
                    if (.not.condFr) then
                        if (print_cnd) print*, 'in solve_flood_basin; condFr failed as Frc,Frf,vc,vf:',Frc,Frf,vc,vf
                        qf_try = qf_try_ini + it*dqf
                        
                        if (qf_try .gt. qfmx) then
                            prt_cd = 1
                            m = (1e-1 - 1e-4) / 1000                            
    !                        do k=1,1000
    !                            qft = qfmin + (k-1)*m
    !                            ff4 = fu4(.false., qft, coefs4, coef_fu4, 1, print_cnd)
    !    !                        if (print_cnd) then
    !                            print*, '**********************************************************'
    !                            print*, '**********************************************************'                    
    !                            print*, 'k,qft,f4', k,qft,ff4
    !    !                        endif                        
    !                        end do                                            
                            print*, &
                                '15apr20; STOP, as all qf domain was explored in steps dqf=', dqf
!                        print*, ''8jun20 6:00, prt verbose here?, prt f4 chorrera?..see added qfmx calcs                                
                            print*, 'take this t,x to edit print_range for DEBUG: i,j:', i1,j1
                            stop
                        endif                    
                    endif                    
                end do
            else
                compos = .true.
                hfmx = .0125/Sc !v>.1, f_dcy<100 (ferguson17bedrock), p228cl
            endif
            Af = wf*h_fld
            Aw = Af + Ac
            aspect = wv**2/Aw
            h = hbf + h_fld
            w = wv; wb = wc-2*hbf/ta
            
!            sg8 = 80*Sc !used also to get Qbf for all ij|j<nR.
            
            if (compos) then
!                if (j1 .lt. nR) stop '4jun20; in compos, j<nR, no hf<0 bug'
                hfmin = 1e-4! ; hfmx = 1e0  21may20; known that v>vmin
                condFr = .false.
                it = 0                
                coefs5 = (/ Acd/wc, a1_2, a2_2, a12_2, wf, wb + 2*hbf/sa, Acd, wv, sg8, &
                    cpp, dx_dt, D84f, D84c /)
                
!                m = (hfmx - hfmin) / 1000  
!                do k=1,1000
!                    hft = hfmin + (k-1)*m
!                    ff5 = fu5(hft, coefs5, print_cnd)
!!                        if (print_cnd) then
!                    print*, '**********************************************************'
!                    print*, '**********************************************************'                    
!                    print*, 'k,hft,f5', k,hft,ff5
!!                        endif                        
!                end do
                hfld_try = hfmin
                do while ( .not. condFr)
                    it = it + 1
                    dhf = (hfmx - hfmin) / itmx                    
                    if (it>itmx) then
                        print*, 'take this t,x to edit print_range for compos DEBUG: i,j,rockw:', i1,j1,rockw
                        print*, 'itmx=', itmx
                        prt_cd = 1
!                        stop 'In flood, tried many hfld_try, but now it>itmx'
                        exit
                    endif                    
                    if (print_cnd) then
                        print*, 'in solve_flood; COMPOSITE channel, prev to find_root_fu5, it,hfld_try:', &
                            it,hfld_try
                    endif
                    call find_root_fu5(print_cnd, rockw, hfld_try, coefs5, 1e-2, 10, &
                        h_fld, success)
                    Aw = Acd + wv*h_fld
!                    print*, '15may20; in COMPOS, i1,j1,h_fld,Aw: ', i1,j1,h_fld,Aw
                    vcomp = cpp/Aw - dx_dt; vc = vcomp; vf = vcomp
                    Af = wf*h_fld                    
                    Ac = Aw - Af
                    Frc = vcomp/sqrt(g*Ac/wc)
                    Frf = vcomp/sqrt(g*h_fld)
                    if (rockw) then
                        condFr = Frc.lt.5 !4jun20; as wf=0.
                    else
                        condFr = (Frc.lt.5 .and. Frf.lt.5)
                    endif
                    print*, ' 21may20; FLD COMPOS, i1,j1,it,Frc,Frf,hfld_try,h_fld_result:', &
                        i1,j1,it,Frc,Frf,hfld_try,h_fld                        
                    if (.not.condFr) then
                        if (print_cnd) print*, '  COMPOSITE channel; condFr failed as Frc,Frf:',Frc,Frf
                        hfld_try = hfmin + it*dhf
                    endif
                end do
                if (repeater .eq. 1 .and. Aw .eq. Aw) then
                    repeater = 0 
                    print*, '8jun20; fld sln was repited and now Aw is ok, so repeater is set to 0 for +running'
                endif
                    !if came from Aw_NaN and Aw_ok in reiter,  then continue run.
                Q = vcomp*Aw; Q_f = vcomp*Af; Qc = Q-Q_f; q_c= Qc/wc; qf = Q_f/wf
            endif
            aspect = wv**2 / Aw
            h = hbf + h_fld
            w = wv
            
            if (print_cnd) print*, 'in solve_flood_basin; 15may20; re-calc to debug Aw; Aw:', Aw
            if (Aw .ne. Aw) then !15may20; Aw NaN check
                print*, 'take this t,x to edit print_range for Aw DEBUG: i,j,Aw:', i1,j1,Aw
                prt_cd = 1 !to repeat i1j1 calcs to print details and debug, 
                    !as bugs dont emerge in same i1j1 then whole model run is done 
                    !(f90 numNoise Lorenzian chaos?).
!                stop
            endif
        else    !critical flow at outlet with rectangle crosssection and rigid walls.
            Q = Qup
            w=wv; h=(Q/(w*sqrt(g)))**.67; vc=sqrt(g*h); vf=0; Frc=1; Frf=0; aspect=w/h; Q_f=0; h_fld =0 
            if (print_cnd) print*, '  FLOODsln_j=nR: h=', h, ', h_fld=', h_fld, & 
                ', vc=', vc, ', vf=', vf, ', Frc=', Frc, ', Frf=', Frf, ', w=', w, 'Q=Qc,Q_f:', Q,Q_f
!            print*, '15may20; as j1:', j1, ', code must come here for critical section: Fr=1.'
        endif
!        print*, 'mijito, see this vars at end of solve_flood_basin at i1,j1,compos,minifld,qfmin?fff4:', &
!            i1,j1,compos,minifld,fff4
!        print*, '  h,h_fld,wb,wc,wv: ', h,h_fld,wb,wc,wv        
!        print*, '  sqc,sqf,Q: ', sqc,sqf,Q
!        print*, '  Qc,Q_f,q_c,qf,vc,vf: ', Qc,Q_f,q_c,qf,vc,vf
    end subroutine
    
    real function tsm84(calib, Dr, S)
        real, intent(in) :: calib, Dr, S
        if (S .lt. 3e-2) then !eq 11 recking13asce, for So<.03 see eq5 reckParker15. 
            !For steeper, eq 4 (Lamb type) of reckParker15.
            tsm84 = (5*S + 0.06)*Dr**(4.4*sqrt(S) - 1.5)
        else
            tsm84 = 1.5*S**.75
        endif
        tsm84 = calib*tsm84
    end function
    
    subroutine update_La(print_cnd, dt, kLa, D90, La, La_fut, dLa_dt)
        logical :: print_cnd
        real, intent(in) :: dt 
        real, intent(in) :: kLa, D90, La
        real, intent(out) :: La_fut, dLa_dt
        La_fut = kLa*D90
        dLa_dt = (La_fut - La) / dt
        if (print_cnd) print*, '  in update_La; La,kLa,D90,dLa_dt', La,kLa,D90,dLa_dt
    end subroutine
    
	function get_port(np) result(port)
		integer np
		real :: port(np), den
		den = np-1
		do k=1,np
			port(k)=(k-1.)/den
		end do	
	end function
	
	subroutine get_GSDpctiles(print_cnd, nD, np, port, P, D, D50, D84, D90, indDp)
		!input: points of cumulative distribution
		!produces equally spaced grain size percentiles, each dp=.1.
		logical print_cnd
		integer, intent(in) :: nD, np, indDp(3)
		real, intent(in) :: port(np), P(nD), D(nD)
		real, intent(out) :: D50, D84, D90
		real :: Pmx, dD_dp, dp, Dpct(np) !to get pctiles each 1/(np-1)*100%         , D80
		do k=1,np
			k2=1
			Pmx = .99*port(k)
			do while (P(k2)<Pmx)
			    k2=k2+1
			end do
			if (k2>1) then
			    dD_dp = (log(D(k2))-log(D(k2-1))) / (P(k2)-P(k2-1))
			    dp = (port(k)-P(k2-1))
			    Dpct(k)=log(D(k2-1))+ dD_dp*dp
			    Dpct(k)=exp(Dpct(k))
			else
			    Dpct(k)=(.062e-3*D(k2))**.5  !non cohesive sediment sizes
			endif
		end do
		!key pctiles, useful e.g. for friction formulas.
		D50 = Dpct(indDp(1))
        D90 = Dpct(indDp(3)) !D80 = Dpct(indDp(2));
		D84 = Dpct(indDp(2))
!		print*, '28mar, in get_GSDpctiles; D50,D84,D90 = ',D50,D84,D90
	end subroutine	
    
    subroutine update_text(i1, j1, print_cnd, dt, nD, psolid, h_s, La, wbed, dhs_dt, dLa_dt, dwbed_dt, dx, &
        fa, fs, Qsi_k, Qso_k,  &
        fa_fut, fs_fut)
        logical print_cnd
        integer i1, j1
        integer, intent(in) :: nD
        real, intent(in) :: psolid, h_s, La, wbed, dhs_dt, dLa_dt, dwbed_dt, dx, dt
        real, intent(in), dimension(nD) :: fa, fs, Qsi_k, Qso_k
        real, intent(out), dimension(nD) :: fa_fut, fs_fut
        real, dimension(nD) :: fi, div_qs,dfa_dt, dfs_dt
        real fa_sum_abs, fs_sum_abs, hss, dhss_dt, dLar, dhssr1, dhssr2, dwr, fak, fsk
        real fmin
        if (wbed>0.01) then        
            if (print_cnd) print*, 'in update_text; wbed>.01m: there is width over which bed texture evolution can happen'
            if (dhs_dt.gt.0) then   !simple texture interchange criteria of gasparini04
                fi=fa !deposition
            else
                fi=fs   !erosion
            endif
            hss = h_s - La
            dhss_dt = dhs_dt - dLa_dt
            !following are vectorial equations (nD):
            div_qs = (Qso_k - Qsi_k) / (wbed*dx)
            dLar = dLa_dt/La
            dhssr1 = dhss_dt/La
            dhssr2 = dhss_dt/hss
            dwr = dwbed_dt/wbed
            dfa_dt = -div_qs/(psolid*La)*(1 + dwr) - dLar*fa - dhssr1*fi - dwr*(fa + fi*(h_s/La-1))
                !19jun20; deduced about 25mar20 at p126
            dfs_dt = -(dhssr2 + dwr)*(fi + fs)
            fa_fut = fa + dfa_dt*dt
            fs_fut = fs + dfs_dt*dt
        else
            if (print_cnd) print*, 'in update_text; wbed<.01m: too narrow to have credible bed texture evolution'
            fa_fut = fs     !mixed alluvium 
            fs_fut = fs
        endif

        if (print_cnd) then
            print*, '  eros(-)/depos(+): dhs_dt=', dhs_dt, ', so fI=', fi
            print*, '  div_qs=', div_qs
            print*, '  dwr,dhssr1,dhssr2,dLar=', dwr,dhssr1,dhssr2,dLar
            print*, '  La,hss,wbed=', La,hss,wbed
            print*, '  fa=', fa, ', sum(fa)=', sum(fa)
            print*, '  dfa_dt=', dfa_dt            
            print*, '  fa_fut=', fa_fut
            print*, '  fs=', fs, ', sum(fs)=', sum(fs)
            print*, '  dfs_dt=', dfs_dt            
            print*, '  fs_fut=', fs_fut
        endif
        fmin = minval(fa_fut)
        if (fmin .lt. 0) then !2nd type of remormalization, via traslation and compression as described in pen notes p128. 
                                    !The more negative the less abundant.
            fa_fut = (fa_fut-fmin) / (maxval(fa_fut) - fmin)
        endif
        fmin = minval(fs_fut)
        if (fmin .lt. 0) then !2nd type of remormalization, via traslation and compression as described in pen notes p128. 
                                    !The more negative the less abundant.
            fs_fut = (fs_fut-fmin) / (maxval(fs_fut) - fmin)
        endif

!        do k=1,nD
!            fak = fa_fut(k)
!            fsk = fs_fut(k)
!            if (fak.lt.0 .or. fak.gt.1) then
!                print*, 'take this t,x to edit print_range for DEBUG: i,j:', i1,j1            
!                print*, 'checking active txt for k=', k, 'fak=', fak 
!                stop 'for this k, active layer texture evolved out of (0,1) range'                
!            endif                           
!            if (fsk.lt.0 .or. fsk.gt.1) then
!                print*, 'take this t,x to edit print_range for DEBUG: i,j:', i1,j1            
!                print*, 'checking subsurf txt for k=', k, 'fsk=', fsk
!                stop 'for this k, subsurface layer texture evolved out of (0,1) range'                
!            endif                            
!        end do
        fa_fut = fa_fut/sum(fa_fut)
        fs_fut = fs_fut/sum(fs_fut)
!        fa_sum_abs = abs(sum(fa_fut))
!        fs_sum_abs = abs(sum(fs_fut))
        if (print_cnd) then
            print*, 'in textCheck after output remormalization:'            
            print*, '  fa_fut=', fa_fut            
            print*, '  fs_fut=', fs_fut
        endif
!        if (fa_sum_abs .lt. .95 .or. fa_sum_abs .gt. 1.05) then
!            print*, 'take this t,x to edit print_range for DEBUG: i,j:', i1,j1        
!            print*, 'fa_fut=', fa_fut
!            print*, 'far from 1 ACTIVE texture sum:', fa_sum_abs
!            stop
!        endif
!        if (fs_sum_abs .lt. .95 .or. fs_sum_abs .gt. 1.05) then
!            print*, 'take this t,x to edit print_range for DEBUG: i,j:', i1,j1        
!            print*, 'fs=', fs, ', sum(fs)=', sum(fs)
!            print*, 'dfs_dt=', dfs_dt
!            print*, 'fs_fut=', fs_fut
!            print*, 'far from 1 SUBSURF texture sum:', fs_sum_abs
!            stop
!        endif        
    end subroutine   
    
    subroutine gtm(print_cnd, nD, np, xnp, ynp, port, fl_gtm, tr, pblB, PblB_acu, pbl, M, gam)
        logical :: print_cnd
        integer, intent(in) :: nD, np
        real, intent(in) :: fl_gtm(7), tr, pblB(nD), PblB_acu(nD), port(np), xnp, ynp
        real, intent(out) :: pbl(nD), M, gam
        integer :: cont, pmx_mov
        real :: ff(np), ff_nD(nD), cum, betah, beta
        real :: ref, betal, gam0, gam1, gam2, fi0, taofm_taocref, fic
        !parameters for fractional bedload transport 'GTM' Recking2016
        ref=fl_gtm(1); 
        betal=fl_gtm(2); gam0=fl_gtm(3); gam1=fl_gtm(4); gam2=fl_gtm(5); fi0 = fl_gtm(6); 
        taofm_taocref=fl_gtm(7) 
        pbl = 0; M=0; ff=0; ff_nD=0
        if (tr>1) then
            if (print_cnd) print*, '23abr20; taofm_taocref: ', taofm_taocref
            if (taofm_taocref .gt. 1e-2) then !1% with num tol. Remember gtm recking16 method requires 50<ref<90, 2<tfm/tcr<4.
                betah = (13.4 - 2.91*log(ref)) / taofm_taocref
            else
                betah = betal !as ref grain is moved almost all other move too.
                if (print_cnd) print*, '23abr20; debiera entrar aqui'
            endif                
            beta=betah
        else
            beta=betal
        endif
        M = min(100.,ref*tr**beta)
        gam=min(100.,gam0+gam1*tr**gam2)
        pmx_mov = nint( xnp*M/100 ) + 1  !28mar20: resolution of grain pctiles each 1/(np-1)*100% 
            !***for pmx_mov deduction, see notes p130ul.
        if (print_cnd) print*,'in gtm(), tr,M,pmx_mov,gam,betah: ',tr,M,pmx_mov,gam,betah
        do k=1,pmx_mov
!            print*,'ants d usar port'
            fic = 1 - (100*port(k)/M)**gam
!            print*,'tras usar port'
		    ff(k) = min(1., fi0 + max(0., fic))
        end do
        k=1; k2=1 
        cum=0; cont=0
        do while ( (k .le. nD) .and. (PblB_acu(k) .le. ynp*(pmx_mov-1)) ) !see comment ***.
        	cum=0; cont=0
        	do while (k2 <= np .and. port(k2)<=PblB_acu(k) .and. ff(k2)>fi0)
        		cum=cum + ff(k2)
        		cont=cont+1
        		k2=k2+1
        	end do
			if (cont>0) then
            	ff_nD(k)=cum/cont
            elseif (k .gt. 1) then    !two consecutive size PblB_acu values in same multiple of 10% in port, 
            !e.g. PblB_acu = .05 .06 .16; 2nd value takes this 'else'
                ff_nD(k) = ff_nD(k-1)
        	endif
        	k=k+1
!        	if (print_cnd) print*, '12apr20: in gtm 2xwhile, k=',k
        end do
        if (print_cnd) print*,'  ff',ff
        if (print_cnd) print*,'  ff_nD',ff_nD

        pbl = ff_nD*pblB
        if (sum(pbl) .eq. 0) pbl(1)=1 !only to avoid Nan fractional Qs 
        if (print_cnd) print*,'  pblB=',pblB
        if (print_cnd) print*,'  raw (no normalized) pbl=',pbl
!            pbl(nD) = 0
    	pbl = pbl/sum(pbl) !so sum=1, can be redundant
    	if (print_cnd) print*, '  normalized pbl: ', pbl
!        if (print_cnd) print*, '12apr20: passed norm pbl if'
!        if (print_cnd) print*, '12apr20: pbl(nD), size(pbl), M, gam:', pbl(nD), size(pbl), M, gam
    end subroutine

    subroutine get_Qs_capac(print_cnd, coefs_wqsy, ss_trg, calib, Q, Qf, v, w_sf, w_v, w_c, w_b, w, h, h_bf, sa, ta, D50c, D50f, &
        D84c, D84f, Sc, i1,j1, &
        a_qs, b1_qs, b2_qs, flume_else_basin, Acd, a1_2, a2_2, a12_2, sg8, &
        fl_gtm, nD, np, xnp, ynp, sq3c, sq3f, port, fac, fac_acu, faf, faf_acu, &
        Qs_c, Qs_fl, pbl, pblf, ts50h, tc50, tsm, subm, Mgtm, gam, qsbx, fic50, fif50, fibc50, cmx, Qbf, Qsbf)
        logical, intent(in) :: print_cnd
        integer, intent(in) :: flume_else_basin, nD, np, i1, j1
        real, intent(in) :: coefs_wqsy(8), ss_trg(3), Qf, Acd, a1_2, a2_2, a12_2, sg8
        real, intent(in) :: calib, Q, v, w_sf, w_v, w_c, w_b, w, h, h_bf, sa, ta, D50c, D50f, D84c, D84f, Sc, a_qs, b1_qs, b2_qs
        real, intent(in) :: fl_gtm(7), port(11), fac(nD), fac_acu(nD), faf(nD), faf_acu(nD), xnp, ynp, sq3c, sq3f
        real, intent(out) :: Qs_c, Qs_fl, pbl(nD), pblf(nD), ts50h, tc50, tsm, subm, Mgtm, gam, qsbx, fibc50, cmx
        real, intent(out) :: fic50, fif50, Qbf, Qsbf
        real :: Awet_c, Awet_f, Pc, Pf, Rc, hfld, Rf, w_eff_ch, c1, eps, fim, fibm, tbt, sq50_11, Rp, asp_b_inv, Pcbf
        real :: ts50hf, ctot, cr, hmc, fc, vbf, ts50bf, fibf50, cts, cbf
        !same critical stress for FL and CH, as both share same energy gradient Sc, even if FL/terrace may not be parallel to CH.
        cmx = 4e-2 !p243; maybe this supplies too many sed to jdown leading to errSedBal due to high dhs.
            !calhoun17espl says limit c above which flood mech is affected: 4% by vol...guarde doi!!!
        tc50 = calib*(1.32*Sc + .037) !eq 23 recking09 (D50)
!        tc = calib*(.56*Sc + .021) !eq 12 recking16gtm, from eq 24 recking09 (D84)
        tbt = coefs_wqsy(2)  !tbt is the same lower case fi, that is ratio tao_bank/tao; this name tbt is used instead
              !to avoid confusion with fi=tao/taop, p= c or m                      
        Pcbf = w_b + 2*h_bf*sqrt(1 + 1/ta**2) !note interface with floodplain flow in not accounted.              
        if (h < h_bf) then
            Pc = w_b + 2*h*sqrt(1 + 1/ta**2)
            Awet_c = Q/v
            Pf = 0; Awet_f = 0; Qs_fl = 0; pblf = 0
            w_eff_ch = w
            asp_b_inv = h/w_b
            Rp = h*(1 + asp_b_inv/ta) / (1 + 2*tbt/sa*asp_b_inv)     !like Cantelli07 c2 coef, but for whole instead half channel; 
                                !after Qunif momentum eqn. (phdSCnotes p 167d).
            if (print_cnd) print*,'1may20; in get_Qs_capac; insumos de Rp h,asp_b_inv,ta,sa,tbt,Rp:', h,asp_b_inv,ta,sa,tbt,Rp
        else !20mar20: pendiente pues friction law must come form different eqn.
            if (print_cnd) print*, 'In get_Qs_capac, FL ::'
            hfld = h - h_bf
            Pc = Pcbf
!            Pf = w_sf + 2*hfld
            Awet_c = .5*(w_c + w_b)*h_bf + w_c*hfld
!            Awet_f = hfld*w_sf
!            Rf = Awet_f/Pf
            
!            ts = Rf*Sc/(1.65*D84f)
            tsm = tsm84(calib, D84f/D50f, Sc)
!            Qs_fl = w_sf * sq3f * a_qs*ts**b1_qs/(1+(tsm/ts)**b2_qs)

            asp_b_inv = h_bf/w_b
            Rp = ( h_bf*(1 + asp_b_inv/ta) + w_v/w_b*hfld ) / ( 1 + 2*tbt*( hfld/h*(hfld+w_sf)/w_b + asp_b_inv/sa) )
                !phdSCnotes p167d.
            ts50h = Rp*Sc/(1.65*D50c)
            ts50hf = ts50h*hfld/h*D50c/D50f
            sq50_11 = sqrt(16.5*D50f**3)*11.2 !for silica r*g is 16.5; r=s-1; s = gamma_s/gamma_w.
            fif50 = ts50hf/tc50
            if (fif50>1) then
                Qs_fl = w_sf*sq50_11*ts50hf**1.5*(1- 1/fif50)**4.5
            else !22apr20: to avoid negative base of power that is NaN
                Qs_fl = 0                                
            endif

            if (Qs_fl .ne. Qs_fl) then !10jun20; to solve c= NaN bug, p255ul.
                print*, 'take this t,x to edit print_range for Qs_f NaN, in get_Qs_capac(), DEBUG: i,j:', &
                    i1,j1
                print*,'w_eff_ch,sq50_11,ts50h,1- 1/fic50:',w_sf,sq50_11,ts50hf,1- 1/fif50
                stop
            endif
            
            w_eff_ch = w_c

            if (print_cnd) print*, 'In get_Qs_capac, FL: Sc,D50f,D84f,hfld,ts50hf: ', Sc,D50f,D84f,hfld,ts50hf
            if (print_cnd) print*, '  14may20; before go in gtm pblB=faf: ', faf
            if (print_cnd) print*, '  shields_f,tc50,Qs_fl,c: ', ts50hf,tc50,Qs_fl,Qs_fl/Qf
            call gtm(print_cnd, nD, np, xnp, ynp, port, fl_gtm, ts50hf/tc50, faf, faf_acu, pblf, Mgtm, gam)
            if (sum(pblf).eq.0) Qs_fl=0    !no flow competence.

        end if
        if (print_cnd) print*, 'In get_Qs_capac, CH ::'
        Rc = Awet_c/Pc
!        ts = Rc*Sc/(1.65*D84c)
        tsm = tsm84(calib, D84c/D50c, Sc)
!        fim = ts/tsm
        c1 = sq3c * a_qs
        
!        Qs_c = w_eff_ch * c1 * ts**b1_qs / (1 + (1/fim)**b2_qs)
        !22apr20: parker79 try, as more peaky c in initial new Q is requiered for pitlick13flood calibration.
        !based on cantelli form of parker eqn (coincides with pitlick13; see phdSCnotes p169)
        cts = Sc/(1.65*D50c)
        ts50h = Rp*cts
        sq50_11 = sqrt(16.5*D50c**3)*11.2
        fic50 = ts50h/tc50
        if (fic50>1) then        
            Qs_c = w_eff_ch*sq50_11*ts50h**1.5*(1- 1/fic50)**4.5
        else
            Qs_c = 0
        endif            
        
        if (Qs_c .ne. Qs_c) then !10jun20; to solve c= NaN bug, p255ul.
            print*, 'take this t,x to edit print_range for Qs_c NaN, in get_Qs_capac(), DEBUG: i,j:', &
                i1,j1
            print*,'w_eff_ch,sq50_11,ts50h,1- 1/fic50:',w_eff_ch,sq50_11,ts50h,1- 1/fic50
            stop
        endif
        
        !longitudinal transport on the sloping bank. See phdSCnotes p164
        eps = f_eps(coefs_wqsy, ss_trg) !taoc_sloping/taoc_flat (eq. 2.69 garcia08asceCh2), with intermediate side slope.
!        fibm = fim*tbt/eps !assume eps traslates equally from channel to bank the vars: tsm and tc.
            !Subscript b means in bank, then tao fluid and resistant is affected by fi and eps.
            
!        qsbx = c1 * (tbt*ts)**b1_qs / (1 + (1/fibm)**b2_qs) !using c1 and ts implies D84c is also in sloping banks.        
        !parker79; see reasons about 5-10 lines above:        
!        print*, '21may20; tbt,eps,tbt/eps: ', tbt,eps,tbt/eps
!        stop       
        fibc50 = tbt/eps*fic50
        if (fibc50 > 1) then
            qsbx = sq50_11*(tbt*ts50h)**1.5*(1 - 1/fibc50)**4.5
        else
            qsbx = 0
        endif
        
            !Subscript x indicates x (longitudinal) direction, as qsbx will affect qsby when estimating widening.
        if (print_cnd) print*, 'in get_Qs_capac, CH; lateralVars; eps,fi,fibc50,qsbx', eps,tbt,fibc50,qsbx
            
        call gtm(print_cnd, nD, np, xnp, ynp, port, fl_gtm, ts50h/tc50, fac, fac_acu, pbl, Mgtm, gam)
!        if (print_cnd) print*, '12apr20; In get_Qs_capac, CH, passed gtm'
        if (sum(pbl).eq.0) Qs_c=0    !no flow competence.        
        subm = Rc/D84c
        if (print_cnd) print*, 'In get_Qs_capac, CH, So_c,,D50c,D84c,Rp: ', Sc,D50c,D84c,Rp
        if (print_cnd) print*, '  Shields,tc50,Qs_c: ', ts50h,tc50,Qs_c
            !shields is eq4 pitlick13flood.
        ctot = (Qs_c + Qs_fl) / Q
        cr = ctot/cmx
        if (cr .gt. 1) then        
            Qs_c = Qs_c/cr;  Qs_fl = Qs_fl/cr
            if (print_cnd) print*, ' 19may20; rescaled Qs to be below cmx=', cmx
        endif
        
        !26may20;  Qbf (bankfull) to validate with HG parker07.        
        if (flume_else_basin .eq. 0) then
            hmc = Acd/w_c
            Rc = Acd/Pcbf
            fc = f_dcy(a1_2, a2_2, a12_2, D84c/hmc)
            vbf = fv_dcy(sg8, fc, Rc)
            Qbf = vbf*Acd
            if ( Qbf .eq. 0) then
                print*, 'sg8, fc, Rc,D84c,w_c: ', sg8, fc, Rc,D84c,w_c
                print*, 'Acd, vbf: ', Acd, vbf
                stop '29may20; only way Qbf null is if Qcd does. Why this xt gives Qbf=0 ?'
            endif
            ts50bf = Rc*cts
            fibf50 = ts50bf/tc50
            if (fibf50>1) then
                Qsbf = w_c*sq50_11*ts50bf**1.5*(1- 1/fibf50)**4.5
            else
                Qsbf = 0                                    
            endif
            cbf = Qsbf/Qbf
    !        print*, 'In get_Qs_capac; 27may20; i1,j1,Qbf: ', i1,j1,Qbf
    !        print*, 'In get_Qs_capac; 27may20; fc,vbf,hmc,ts50bf,cbf: ', fc,vbf,hmc,ts50bf,cbf
        endif
    end subroutine
    
    subroutine get_wlog(Qs_c, Q, w, dt, par_wlog, wlog)
        !empirical logistic grouth of flow top width. See trench_flume_2020.xlsx
        real, intent(in) :: Qs_c, Q, w, par_wlog(6), dt
        real, intent(out) :: wlog
        real :: tsat, wsat, dt_nd, w_nd, w_nd_next, dw_dt_nd, r
        real :: a_tsat, a_wsat, a_r_wlog, b_tsat, b_wsat, b_r_wlog
        a_tsat = par_wlog(1); a_wsat = par_wlog(2); a_r_wlog = par_wlog(3)
        b_tsat = par_wlog(4); b_wsat = par_wlog(5); b_r_wlog = par_wlog(6)
        tsat = a_tsat*Q**b_tsat
        wsat = a_wsat*Q**b_wsat
        dt_nd = dt/tsat
        w_nd = w/wsat
        r = a_r_wlog*(Qs_c/Q)**b_r_wlog
        dw_dt_nd = r*w_nd*(1-w_nd)
        w_nd_next = w_nd + dw_dt_nd*dt_nd
        wlog = w_nd_next*wsat
    end subroutine

    real function f_eps(coefs_wqsy, ss_trg)
        !deduction: notes p148d. Variable ensapsulation before x,t calcs did not include bank slope 'a',
        !as it will be dynamic in future due to cohesion/vegetation.
        real, intent(in) :: coefs_wqsy(8), ss_trg(3)
        real sa, ca, a, c_det, c_det2, cb_ng, miu2, inv_2a, b_ng, ca2, trg1, trg2, det
        sa = ss_trg(1); ca = ss_trg(2) !trigonometry updated for intermediate slope in trapeze corner to get qsy to widen.
                                        !see phdSCnotes p165dr.
        a = coefs_wqsy(4); c_det = coefs_wqsy(5); c_det2 = coefs_wqsy(6); cb_ng = coefs_wqsy(7); miu2 = coefs_wqsy(8)
        inv_2a = 1/(2*a)
        b_ng = cb_ng*ca
        ca2 = ca**2
        trg1 = c_det2*ca2
        trg2 = sa**2/miu2 - ca2
        det = c_det * (trg1 - 4*a*trg2)
        f_eps = inv_2a * ( b_ng + sqrt(det) ) !I chose + as b_ng<0.
        !seminaraParkerSolari02 have more general formulae for high longitudinal slope; I studied their neat statics.
    end function
    
    subroutine get_w_qsy(print_cnd, flume_else_basin, isFlood, coefs_wqsy, ss_trg, qssx, tpt_stg, psolid, hbf, wc, &
        dt, w, pWreveg, wc_qsy)
        !transverse -y- sediment transport qsy in small slopes. qsy = Qsy per unit 'width' that, for y direction, is dx.
        !see notes p146 and pitlick13, cantelli07 and ferrer14. 
        !This is a linear theory for gentle x&y bed gradients, and small secondary currents i.e. tao_y/tao_x.
        logical, intent(in) :: print_cnd, isFlood
        integer, intent(in) :: flume_else_basin       
        real, intent(in) :: coefs_wqsy(8), ss_trg(3), qssx, tpt_stg, psolid, hbf, wc, dt, w, pWreveg
        real, intent(out) :: wc_qsy
        real ta, b, fi, qssy_qssx, qssy, dwc_dt,w_expos
        ta = ss_trg(3) !intermediate slope smaller than repose, so feasible qsy so grow wc.
        b = coefs_wqsy(1)
            !physics vary if sources of variable ta like **: wetness, cohesion or vegetation.
        qssy_qssx = b * sqrt(1/tpt_stg) * ta    !tpt_stg tao/taoc comes already traslated to bank via eps and fi coefs.
        qssy = qssy_qssx * qssx
        dwc_dt = qssy/(psolid*hbf)   !. For chances of narrowing, see pitlick13 discussion. 
        wc_qsy = wc + dwc_dt*dt
        if (print_cnd) print*, 'in get_w_qsy; ta_eff,b: ', ta,b
        if (print_cnd) print*, 'in get_w_qsy; hbf,qssy_qssx,dwc_dt(cm/min): ', hbf,qssy_qssx,dwc_dt*6e3
        if (print_cnd) print*, 'in get_w_qsy; wc,wc_qsy: ', wc,wc_qsy
            !note Ferrer, Cantelli and Croissant incision narrows!. Do they talk about wc or w (flow)?        

!        print*, '25may20; qssx:', qssx
        w_expos = wc - w !25may20; grows with incision. pWreveg factor involves %t with h over plants colmating &not eroding'em.
        if ( (flume_else_basin .eq. 0) .and. (qssx .eq. 0) .and. (w_expos .gt. 0) ) then !26may20; depr: .not.isFlood (w<wc).
            !w- via reveget. p231.
            wc_qsy = wc_qsy - w_expos*pWreveg
!            !Modified from DavidsonEaton2018.
        endif
        if (print_cnd) print*, 'in get_w_qsy; after reveget; w_expos,pWreveg,wc,wc_qsyReveg: ', w_expos,pWreveg,wc,wc_qsy
    end subroutine    
    
    subroutine get_alluv_derivs(i1,j1,print_cnd, h_bf, Q, w_c, w_v, h_s, h_sf, w_b, colm_prev, ta, dx, dvol_dwc, &
        dvol_dhsf, dvol_dhs, w_sf)
        logical, intent(in) :: print_cnd
        integer, intent(in) :: i1,j1
        real, intent(in) :: Q, w_v, h_s, h_sf, w_b, colm_prev, ta, dx, w_c, h_bf
        real, intent(out) :: dvol_dwc, dvol_dhsf, dvol_dhs, w_sf
        real :: w_ss
!        h_bf = h_sf - h_s
        w_ss = 2*h_bf/ta
!        w_c = w_b + w_ss !7jun20; deprec (wc in instead), 
            !to avoid initial numf90 bug in colm:-2 at get_hs. (wc>wv)
        if ( abs(w_c-w_ss-w_b)/w_c .gt. .05 ) then !8jun20; tol 1%  !10jun20; see analogous check in get_hs().
            print*, 'take this t,x to edit print_range for wc vs wb DEBUG: i,j:', i1,j1
            print*, 'w_c,h_bf,w_ss,w_b, err= w_c-w_ss-w_b:', w_c,h_bf,w_ss,w_b, w_c-w_ss-w_b
            stop 'check geom, as wb is not coherent with wc'
        endif            
        w_sf = w_v - w_c
        if (w_sf<-.05) then
            print*, 'take this t,x to edit print_range for DEBUG: i,j:', i1,j1
            print*, 'Q=', Q, 'wv=', w_v, ' and wc=', w_c
            print*, 'colm_prev:', colm_prev !p244d
            stop 'wsf=wv-wc<0'
        endif
        dvol_dwc = -dx*h_bf
        dvol_dhsf = dx*(w_sf + w_ss)
        dvol_dhs = dx*w_b !6jun20; p224d(w_c - w_ss) equiv
        if (print_cnd) print*, 'in get_alluv_derivs; wv,wc, dvodwc,dvodhsf,dvodhs', &
            w_v,w_c, dvol_dwc, dvol_dhsf, dvol_dhs
    end subroutine
    
    subroutine get_wc_hsf(print_cnd, wc_prj, Q, Qf, v, ta, wv, h, h_bf, wb, wc, h_s, h_sf, &
        w_sf, psolid, Is_prev, Qs_fl, dt, dx, wb_fut, wc_fut, h_sf_fut, dwc_dt, dhsf_dt, pf, dhsf_type, wc_ambit)
        logical, intent(in) :: print_cnd
        real, intent(in) :: dt
        real, intent(in) :: wc_prj, Q, Qf, v, ta, wv, h, h_bf, wb, wc, h_s, h_sf, w_sf, psolid, Is_prev, Qs_fl, dx
        real, intent(out) :: wb_fut, wc_fut, h_sf_fut, dwc_dt, dhsf_dt, pf, wc_ambit
        integer, intent(out) :: dhsf_type
        real :: det, h_sf_change, Is_pr_allD_f
        dhsf_type = 0; wc_ambit = 0

!        if (Qf.eq.0) then  !deprecated wlog (which was flow w) as widening mech was changed to pitlick13 qsy.
!            det = wc_prj**2 - 4*Q/(v*ta) !see deduction in phdSCnotes p120d.
!            if (det.gt.0) then
!                if (print_cnd) print*, 'in get_wc_hsf, choosing wb_fut: Qf.eq.0, det.gt.0'
!                wb_fut = max(sqrt(det), wb)
!                    !max is used because det (f(w_log)) produced smaller wc_fut so channel has not known way to get narrow,
!                    !like vegetation encroachment (EatonDavidson18?).                    
!            else
!                if (print_cnd) print*, 'in get_wc_hsf, choosing wb_fut: Qf.eq.0, det.lt.0'
!                wb_fut = wb
!            endif
!        else
!            if (print_cnd) print*, 'in get_wc_hsf, choosing wb_fut: Qf.gt.0'
!            wb_fut = wb + (wc_prj-wc) !bank retreats paralell i.e. with same slope angle.
!                !though ubc flume experiments at 2018 did not included floodplain, pitlick13 did as saw dwcdt>0 (see qsy).
!        endif
        wb_fut = wb + (wc_prj-wc)
!        wc_fut = wb_fut + 2*h_bf/ta  !19may20; deprec: to avoid num errors when checking no widening.
        wc_fut = wc_prj
        if (print_cnd) print*, 'in get_wc_hsf, wc_fut(virt)=',wc_fut
        if (print_cnd) print*, 'in get_wc_hsf, Qf=',Qf
        pf = Qf/Q
        Is_pr_allD_f = pf*Is_prev
        if (wc_fut .gt. wv) then
            h_sf_fut = h_sf !h_sf - .5*(wc_fut - wv)*ta !deprec 6jun20, p245. Leave it to be changed in get_hs, 
                !so dhsf_dt = 0 from this subroutine .
                !Though this nonull deprec formulae was valid for ubcTrenchExp (wf=0), is not for basin wf>=0.
            dhsf_type = 1 !touch wall leads to sediment fall, but without alluvWallExhaustion. 
                !Observed in UBCtrench (wf=0=const).
            wc_ambit = wc_fut !6jun20; ambitioned even with valley constraint.
!            if (h_sf_fut .lt. h_s) then !deprec by 6jun20
!                h_sf_fut = h_s !6jun20; only to make code enter to colm1 for thalwegization.
!                dhsf_type = 2 !touch wall and ExhaustWall
!                if (print_cnd) print*, '  15may20; EXHAUSTwall: not all sed demanded via widening is viable, otherwise hsf<hs.'
!            endif
            if (print_cnd) print*, 'in get_wc_hsf, choosing h_s_fut: wc_fut.gt.wv. wv,wc_fut(virt),h_sf,h_sf_fut=',  &
                wv,wc_fut,h_sf,h_sf_fut
        else
            if (pf .lt. 1e-3) then
                if (print_cnd) print*, 'in get_wc_hsf, choosing h_s_fut: wc_fut.le.wv, Qf/Q.lt.1e-3'
                h_sf_fut = h_sf
                dhsf_type = 3 !no flood
            else
                if (print_cnd) print*, 'in get_wc_hsf, choosing h_s_fut: wc_fut.le.wv, Qf/Q.ge.1e-3'
                    !assumes vertical concentration profile, i.e. constant in depth, despite c being bedload.
                        !Yang singapur vertical currents of 2nd  flow verticalize profile.
                h_sf_change = (Is_pr_allD_f - Qs_fl*dt)/(psolid*w_sf*dx)
                h_sf_fut = h_sf + h_sf_change
                dhsf_type = 4 !erosion by flow shear on floodplain
            endif            
        endif
        wc_fut = min(wv, wc_fut)
        dhsf_dt = (h_sf_fut - h_sf) / dt !6jun20; no null only if dhsf_type=4 (fld flow shear)
        if (pf < 1e-3) dhsf_dt = 0 !wall fall due to wc reaching w like in ubc trench comes in get_hs().
        dwc_dt = (wc_fut - wc) / dt
        if (print_cnd) print*, 'in get_wc_hsf, w_b,w_b_fut,h_sf,h_sf_fut', wb,wb_fut,h_sf,h_sf_fut
    end subroutine

    real function zc(i1,j1,wv, ta, hbf, wb_min)
        !see deduction in notes, p2019_121u. Assumed wb_min=0.01m, to preserve trapezial shape.
        integer, intent(in) :: i1,j1
        real, intent(in) :: wv, ta, hbf, wb_min
        real :: mb, c, term2, zc1, zc2
        mb = wv*ta + 2*hbf !mb: -b, i.e. 'minus'b.
        c = hbf*ta*(wv - wb_min)
        term2 = sqrt(mb**2 - 4*c)
!        zc1 = .5*(mb + term2)
        zc2 = .5*(mb - term2) !10may20: in an example run I found + root gives zc>hbf: nonsense!.
!        print*, 'take this t,x to edit print_range for DEBUG: i,j:', i1,j1
!        print*, 'solving zc for thalweg, zc1=', zc1, 'zc2=', zc2
!        stop 'check which zc to choose'        
        zc = zc2
    end function
    
    real function dw_rw(print_cnd, b, c)
        !see deduction, and selection of + instead -, in notes, p2019_123d.
        logical, intent(in) :: print_cnd
        real, intent(in) :: b, c !a=1 to reduce 2nd order polinomia to solve.
        real :: mb, term2
        mb = -b
        term2 = sqrt(mb**2 - 4*c)
        dw_rw = .5*(mb + term2)
        if (print_cnd) print*, '      in fn dw_rw; b,c,term2:', b,c,term2
    end function    
    
    subroutine triangle_incis(i1, j1, print_cnd, w_c, h_sf, h_s, dt, dvol_dwc, dvol_dhsf, dvol_dhs, dwc_dt, dhsf_dt, ni, dx, &
        sa, ca, ta, wv, h_sf_fut2,z,dwc)
        logical, intent(in) :: print_cnd
        integer, intent(in) :: i1, j1
        real, intent(in) :: dt
        real, intent(in) :: w_c, h_sf, h_s, dvol_dwc, dvol_dhsf, dvol_dhs, dwc_dt, dhsf_dt, ni, dx, sa, ca, ta, wv
        real, intent(out) :: h_sf_fut2, z, dwc !z is used instead of dhs to write eqns faster. 
        real :: z_adm, niadm,nirw !admisible incision [m], and net and via-rewiden input [m3sed/s].
        real :: Ai, x1, x2, dwrw
        z_adm = h_sf - h_s - .5*w_c*ta     !which reduces wb of trapeze to 0.
        niadm = dvol_dhs*z_adm/dt + dvol_dhsf*dhsf_dt + dvol_dwc*dwc_dt
        nirw = ni - niadm
        Ai = -nirw*dt/dx !Ai must be positive for geometry calcs, then minus is used as I know ni, niadm and nirw are <0.
        dwrw = dw_rw(print_cnd, 4*w_c, -4*Ai/ta)
        if (print_cnd) print*, '  in triangle_incis; z_adm,ni,niadm,nirw,Ai,dwrw,w_c,wv',z_adm,ni,niadm,nirw,Ai,dwrw,w_c,wv
        if (w_c + dwrw .lt. wv) then
            if (print_cnd) print*, '    this incision preserves hsf i.e. elevation of alluvium touching the valley wall'            

            z = dwrw*ta
            h_sf_fut2 = h_sf
            dwc = dwrw
            
!            print*, '!!take this t,x to edit print_range for DEBUG: i,j:', i1,j1
!            print*, '  2may20; dhs,dwc=dw_rw: ', z,dwc
!            stop 'why incision preserves hsf in trenchExp??'
        else
            if (print_cnd) print*, '    this incision decreases hsf i.e. elevation of alluvium touching the valley wall'
            
!            print*, '!!take this t,x to edit print_range for DEBUG: i,j:', i1,j1
!            stop 'OK; incision changes hsf in trenchExp'
            
            x1 = (wv - w_c)/2
            x2 = dwrw/2 - x1
            z = - (Ai + (x1**2 + x2**2)*ta) / wv !intented negative output as incision geometry (dhs<0) was assumed for geom.
            h_sf_fut2 = h_sf - x2*ta
            dwc = 2*x1
        endif
    end subroutine

    subroutine get_hs(i1,j1,nR,print_cnd, h_bf, diag_un, ta, ca, sa, wv, dx, psolid, h_s, h_sf, w_b, w_c, &
        fic50, fif50, fibc50, Is_prev, Q, Qs, dt, wb_thw, hbf_thw, cmx, dhsf_type, wc_ambit, &
        dvol_dwc, dvol_dhsf, dvol_dhs, dwc_dt, dhsf_dt, w_c_fut, h_sf_fut, w_b_fut, h_s_fut, net_Qs_flow, &
        dhs_dt, dhsf2_dt, dwb_dt, dwf_dt, colm, Qspp)
        logical, intent(in) :: print_cnd
        integer, intent(in) :: i1,j1, nR, dhsf_type
        real, intent(in) :: dt, wb_thw, hbf_thw, cmx, wc_ambit
        real, intent(in) :: h_bf, diag_un, ta, ca, sa, wv, dx, psolid, h_s, Is_prev, Q, Qs, dvol_dwc, dvol_dhsf
        real, intent(in) :: dwc_dt, dhsf_dt, w_b, w_c, h_sf, dvol_dhs, fibc50, fic50,  fif50
        real, intent(inout) :: w_c_fut, h_sf_fut
        real, intent(out) :: w_b_fut, h_s_fut, net_Qs_flow, dhs_dt, dhsf2_dt, dwb_dt, dwf_dt, colm, Qspp
        real :: hs_fld, h_bf2, A_trap2, dhs, zch, wb_min, h_sf_fut2, dhs2, dwc, dz, lss, wcgeom, wcgeomFut, x2, A3
        real :: dhsf_dt2, dwc_dt2, x1, dwb, dhs_loc, dhs_ups, hs_fut_loc, hs_fut_ups, Qsu, dhs_adm, c_add, Qsmx
        real :: Aoff, c_cap, relerr_ups, hbf_fut, dhsf, dhsf_shear, dhs_adm_dt, dhsf_tw, dvol_tw, wc_oov, ddhs, ddwc
        real :: dvol_E, k11, k10, x4, wc_adm, k12, A_wbRecov, dhsf_orig, dwc_adm, Arw, x6, Qsp2, Qsp3, dQs_hsAdm
        real:: relerrGeom, wss, hsf_fut_precolm, net_Qs_flow_colm2, dwc_colm2, dwc_exn_colm2, Qs_hsf_orig, ta_2
        integer :: flg_wb_min
        logical :: tooMuchIncis, cnd1, cnd2, cond_hbf, cond_wb, cond_hsfhs, cond_colm_gt2, cond_inc
        logical :: condNan
        tooMuchIncis = .false.
            !sufix '2' means 'virtual' while channel bed is higher than floodplain.
        Qspp = 0            
        colm=0            
        wb_min = wb_thw     !21may20; to undeprec wb_thw.      !vs Dmx
        Qsu = Is_prev/dt !from upstream
        if (Qsu < 1e-8) Qsu=0 !means c~1e-5 for Q=1l/s
        net_Qs_flow = (Qsu - Qs) / psolid
!        h_bf2 = h_sf_fut - h_s
        dhs_ups =  net_Qs_flow*dt / dvol_dhs
        Qs_hsf_orig = dvol_dhsf*dhsf_dt         
        dhs_loc = -(Qs_hsf_orig + dvol_dwc*dwc_dt)*dt / dvol_dhs
        dhs = dhs_ups + dhs_loc
!        dhs = (net_Qs_flow - dvol_dhsf*dhsf_dt - dvol_dwc*dwc_dt)*dt / dvol_dhs  !exner, try 1/3: wb>min
            !17may20; if dhs fills channel (to hbf shallower than thw_seed) then colmat exner is needed.
            !  If colmat is due to upstream load: thalwegize. Otherwise causes are local:
            !  If local trend to colmat (wc+, hsf-; p215dl) expulse exceess with gained w. 
            !  Expulsion due to solariParker vectorial qs.
            !wv>=wc was garanteed from get_wc_hsf().
        h_s_fut = h_s + dhs
        hs_fut_loc = h_s + dhs_loc
        hs_fut_ups = h_s + dhs_ups
        w_b_fut = w_c_fut - 2*(h_sf_fut - h_s_fut)/ta
        cnd1 = (w_b_fut.lt. wb_min) !6jun20; .99*
        if (cnd1 .and. dwc_dt .lt. 0) then
            !9jun20; p252c, reveget is not able to incise to rewiden --> solve again colm0|reveg=0
            colm = -5
            dhs_loc = -Qs_hsf_orig*dt / dvol_dhs !null dwc
            w_c_fut = w_c
            !rest of if block is pure update..:
            dhs = dhs_ups + dhs_loc
            h_s_fut = h_s + dhs
            hs_fut_loc = h_s + dhs_loc
            hs_fut_ups = h_s + dhs_ups
            w_b_fut = w_c_fut - 2*(h_sf_fut - h_s_fut)/ta
            cnd1 = (w_b_fut.lt. wb_min)
        endif        
        dhsf_orig = dhsf_dt*dt 
!        print*, 'pre check wb_min; w_b_fut:', w_b_fut
        print*, 'in get_hs:' !only dhsf_type=4 brings nonull dhsf_dt.
        print*, ' at i1,j1:', i1,j1
        wss = 2*h_bf/ta
        relerrGeom = (w_c-wss - w_b)/w_c !10jun20; p254dr; denominador wc no wb pa mas flex.
        if ( abs(relerrGeom) .gt. .05 ) then !9jun20; tol 1%; 10jun20; tol5% (1% fail at ~50yr)
            print*, ' geom DEBUG at i1,j1:', i1,j1
            print*, ' relerrGeom= (w_c-wss - w_b)/w_b:', relerrGeom
            print*, '  h_bf,hbf_thw: ', h_bf,hbf_thw
            print*, '  w_b,wss,w_c,wv: ', w_b,wss,w_c,wv
            stop
        endif
        if ((fibc50 .lt. 0) .and. (dwc_dt .gt. 0)) then !p218d
            print*, 'in get_hs; w+_motor DEBUG at i1,j1:', i1,j1
            print*, 'Qs, fibc50, fic50, fif50', Qs, fibc50, fic50, fif50 
                !as fibc50 = tbt/eps*fic50 (tbt~.85, eps~.5 for reposAng~30deg),
                    !side longit qssx may occur (fibc50>1) with c_cap=0 (fic50<1).
            print*, 'w_c, w_c_fut, dwc_dt', w_c, w_c_fut, dwc_dt
            stop 'how widens if no qssx (longit qs near bank) was triggered?'
        endif        
        print*, ' colm= 0 or -5; colm,Is_prev, Q, c, net_VolSed_InFlow:', colm,Is_prev, Q, Qs/Q, net_Qs_flow*dt
        print*, '  dt, dhs, dhsf_dt, dwc_dt(pre colm:-5):', dt, dhs, dhsf_dt, dwc_dt
        print*, '  dvol_dwc,dvol_dhs,wv:', dvol_dwc,dvol_dhs,wv
        print*, '  dhs_ups, dhs_loc:', dhs_ups, dhs_loc        
        print*, '  w_b, w_c, h_s, h_sf:', w_b, w_c, h_s, h_sf
        print*, '  w_b_fut, w_c_fut(byGet_wc_hsf), h_s_fut, h_sf_fut:', &
            w_b_fut, w_c_fut, h_s_fut, h_sf_fut
        flg_wb_min = 0
        cnd2 = (dhsf_type .eq. 1) !7jun20; if w ambitioned by qsy is reached, alluvium falls at valley wall.
        
        !INCISION processes in narrow (wb|wv) setting: 1: wbmin, 2: wcReq>wv.
            !hsf_fut came from own floodplain exner balance.
            !6jun20; added condition 2 after the OR cnd on dhsf_type to include WidthAmbitiousFlows. 
                !E.g. ubcTrench must be colm-2 to destroy wall as real. In fact cnd2 is set only to go to colm-2.
        if (cnd1 .and. .not. cnd2) then 
            !cnd2 implies diferent geometry of falling alluvial wall  to calculate vol due to dhsf.
            colm=-1
            flg_wb_min = 1
            print*, ' colm=-1; pre-solve it; i1,j1,w_b_fut,wb_min:', i1,j1,w_b_fut,wb_min
            w_b_fut = wb_min
!            call triangle_incis(i1, j1, print_cnd, w_c_fut, h_sf_fut, h_s, dt, dvol_dwc, dvol_dhsf, dvol_dhs, dwc_dt, &
!                dhsf_dt, net_Qs_flow, dx, sa, ca, ta, wv, h_sf_fut2, dhs2, dwc)
!                !note partial derivatives are biased to the past, i.e. to the left in time axis.
!                !2may20: note 
!            h_sf_fut = h_sf_fut2;  dhs = dhs2

            !wb_min trapeze incision as Cui06DREAMS. Incision narrows until wb reaches wbmin. wv>wc. Eqn deduced at p187.
            dwb = w_b_fut - w_b !7jun20; p246ul; if wbmin then no reveget alone can drive wb- via hs- (2ndary flow).
!            dhs = (net_Qs_flow - dvol_dwc*dwb/dt - Qs_hsf_orig)*dt / (dvol_dhs  - 2*dvol_dwc/ta)
                !exner, try 2/3: 
                !note Qs convergence, base narrowing and floodplain lowering mitigates bed deposition.
!                if (print_cnd) print*, '!!take this t,x to edit print_range for DEBUG: i,j:', i1,j1
            ta_2 = 2/ta
            dhs = (net_Qs_flow - dvol_dwc*(dwb/dt + ta_2*dhsf_dt) - Qs_hsf_orig )*dt / &
                (dvol_dhs  - ta_2*dvol_dwc) 
                !10jun20; updated p187d, eq10.

            h_s_fut = h_s + dhs !hs and hsf are lowered parallel.
            h_sf_fut = h_sf + dhsf_orig  !in wide valley, bank retreats without lowering floodplain elevation.
            dwc = ta_2*(dhsf_orig - dhs) + dwb !updated 10jun20; see p187d eq 9.   
                !widens if dhs incises more than dhsf be negative, though negative dwb mitigates w gain.
            w_c_fut = w_c + dwc  !10jun20; 3rd term added. 
                 !overwrites widening calculated from transverse transport qsy.            
                !6jun20; altough well fld may erode by shear
            print*, '  dhs,dhsf_orig,dwc,dwb,wv:', dhs,dhsf_orig,dwc,dwb,wv
            print*, '  w_b, w_c, h_s, h_sf:', w_b, w_c, h_s, h_sf
            print*, '  w_b_fut, w_c_fut(vs wv?), h_s_fut, h_sf_fut:', &
                w_b_fut, w_c_fut, h_s_fut, h_sf_fut
!            if (dhs>0) stop 'why at colm=-1 dhs>0?' !7jun20; cnd1 and p187 deduction implies dhs<0; 
                !8jun20; unless numErrors like 2nd-order corners dydz and derivs in t-1.
        endif
        if ((w_c_fut .gt. wv) .or. (cnd2)) then !6jun20
            !10jun20; as all dhsf_type=1 come here, then here come all chances to ambition wcfut>wv.
                !Therefore, colm3 must not produce wc>wv.
            !widening desire of channel base in narrow valley leads to alluvium decrease 
            !if cnd2 w_c_fut came already limited by wv instead of the ambitioned by qsy. 
                !So how to get ambition to set x2?
            !Condition w_c_fut .gt. wv can be true if resulting for the last colm-1 if (exner2/3).
            !in valley wall, as trench ubc experiment SCatano&Hassan 2018.
            if (cnd2) then
                wc_oov = wc_ambit
!                wc_oov = min(w_c_fut, wc_ambit) !8jun20; deprec. 7jun20; supply (as said by colm=-1 exner) vs
                    !capac as said by qsy. oov=out of valley
                    !6jun20; though qsy wants high w, it is restricted by latIncisSupply.
            else
                wc_oov = w_c_fut
            endif
            colm=-2
            x2 = .5*(wc_oov-wv)
            x1 = .5*(wv - w_c)
            if ((x1<0) .or. (x2<0)) then
                print*, '!!take this t,x to edit print_range for x1|x2 neg colm=-2 DEBUG: i,j:', i1,j1
                print*, 'x1,x2:', x1,x2
                stop
            endif
            dhsf_tw = -x2*ta !6jun20; tw: touch wall, implying fall of alluvium top via slide due to dw+.
!            dhsf_dt2 = dhsf_dt + dhsf_tw/dt    !add wall decrease x2*ta (p188u) to floodplain eros/depos.
            dwc_dt2 = 2*x1/dt !replace qsy power by desire by what is possible: wv. 2 channel sides.
!                dhs = (net_Qs_flow - dvol_dhsf*dhsf_dt2 - dvol_dwc*dwc_dt2)*dt / dvol_dhs !deprec by 6jun20. p245cr.
!                dvol_tw = 2*dx*(x2*abs(dhsf_tw)/2) !6jun20; 2 channel sides.
                !6jun20 WRONG; deduc p245cr triangle area of virtual volume to add to balance.
            dvol_tw = 2*dx*x2*h_bf !7jun20; deduc at p246c. Added out of vol partialDerivs, as tw changes volEqn.
                !8jun20; 2 banks.
            k10 = dvol_dhsf*dhsf_dt + dvol_dwc*dwc_dt2
            dhs = ( (net_Qs_flow - k10)*dt + dvol_tw )/ dvol_dhs
                !7jun20; fig at p245c and p246cl; add paralelogram vol due to wall fall. 
                !6jun20 RIGHT: deprec virtual volume. Note dwc was already decreased due to wv lim.
                !note only real erosion adds sediment as dhs; instead, dhsf_tw only changes geometry.

                !exner, try 3/3: wb=min, wv<wc for w_c_fut that resulted from exner2.
                !note if wc cte (trench) so hsf fall (-) deposits in bed (+).                
            h_s_fut = h_s + dhs
            dhsf = dhsf_orig + dhsf_tw
            h_sf_fut = h_sf + dhsf !see y4 in p188 draw. As flood vanishes and 
                !w_c_fut goes to wv (x2=0) hsf tends to unchanged .
            
            !3may20; provisional checks:
            w_c_fut = wv
            hbf_fut = h_sf_fut - h_s_fut
            w_b_fut = w_c_fut - 2*hbf_fut/ta
            print*, ' colm=-2; i1,j1,dhsf_type,wc_oov,wc_ambit,wv:', i1,j1,dhsf_type,wc_oov,wc_ambit,wv
            print*, '  x1,x2,dhs,dhsf_orig,dhsf_tw,dhsf:', x1,x2,dhs,dhsf_orig,dhsf_tw,dhsf
            print*, '  Q,net_Qs_flow:', Q,net_Qs_flow            
            print*, '  w_b, w_c, h_s, h_sf:', w_b, w_c, h_s, h_sf
            print*, '  w_b_fut, w_c_fut, h_s_fut, h_sf_fut:', &
                w_b_fut, w_c_fut, h_s_fut, h_sf_fut
            if (h_s_fut > h_sf_fut) tooMuchIncis = .true. !8jun20; bypass to colm=1. 
                !9jun20; must be true if hbf_fut<thweg?
            if ((w_b_fut<wb_min) .and. .not. tooMuchIncis) then !excessExpulsion as colm=2, to gain width to get wbmin.
                !8jun20; colm-3 and -4 are motivated by trigger process of golly17 buttress effect.
                colm=-3 !p247dl, note both banks are sumed.
                x4 = .5*(wb_min - w_b_fut) !>0 given the preceding if.
                A_wbRecov = 2*x4*hbf_fut
                Qspp =  psolid*dx*A_wbRecov/dt
                net_Qs_flow = (Qsu - Qs - Qspp) / psolid
                c_add = Qspp/Q
                c_cap = Qs/Q
                !updating geometry: hs remains, hsf (wall) falls as wc=wv, so wc remains and wb=min
                w_b_fut = wb_min
                h_sf_fut = h_sf_fut - x4*ta !8jun20; as x4 included both banks.
                print*, ' colm=-3: i1,j1,x4=wb_min-w_b_fut, dhsf, hbf_fut:', i1,j1,x4, h_sf_fut-h_sf, hbf_fut
                print*, '  A_wbRecov,c_add,c_cap,net_Qs_flow:', A_wbRecov,c_add,c_cap,net_Qs_flow
                print*, '  w_b, w_c, h_s, h_sf:', w_b, w_c, h_s, h_sf
                print*, '  w_b_fut, w_c_fut, h_s_fut, h_sf_fut:', &
                    w_b_fut, w_c_fut, h_s_fut, h_sf_fut
                if (c_add + c_cap > cmx) then !near debris flow. 
                    colm=-4 !p248; resemble Golly17 butress effect. 
                    Qsmx = cmx*Q !19may20; as cmx limit was set from qs calc at get_capac(), deprec: max(Qs, cmx*Q)
                    Qspp = Qsmx-Qs
                    net_Qs_flow = (Qsu - Qsmx) / psolid
                    dvol_E = ( net_Qs_flow - k10 )*dt
                    k12 =  .5*ta + h_bf*dx/dvol_dhs
                    k11 = h_sf + dhsf_orig + wv*k12 - h_s - dvol_E/dvol_dhs
                    wc_adm = (k11 - (wv-wb_min)*.5*ta) / k12
                    x2 = .5*(wc_adm - wv)
                    dhsf_tw = -x2*ta
                    dvol_tw = 2*dx*x2*h_bf !8jun20; 2 banks.
                    dhs = (dvol_E + dvol_tw) / dvol_dhs
                    dhsf = dhsf_orig + dhsf_tw                    
                    !updating geometry: wc remains wv, wb remains min, hsf (wall) falls as wc=wv and hs varies given dhs.
                    h_s_fut = h_s + dhs
                    h_sf_fut = h_sf + dhsf
                    print*, 'colm=-4; please check if OK: less incis than hs(colm-2) and hsf(colm-3)'
                    print*, '  >0: i1,j1,k12,x2', i1,j1,k12,x2
                    print*, '  =0 if in initially touch wall and no fld; k10:', k10                    
                    print*, '  ok wc_adm>wv: wc_adm,wv', wc_adm,wv
                    print*, '  dhs,dhsf:', dhs,dhsf
                    print*, '  net_Qs_flow:', net_Qs_flow
                    print*, '  w_b, w_c, h_s, h_sf:', w_b, w_c, h_s, h_sf
                    print*, '  w_b_fut, w_c_fut, h_s_fut, h_sf_fut:', &
                        w_b_fut, w_c_fut, h_s_fut, h_sf_fut
                    if (wc_adm<wv) then
                        print*, '!!take this t,x to edit print_range for colm=-4 wc_adm DEBUG: i,j:', i1,j1
                        print*, 'wc_adm,wv:', wc_adm,wv
                        stop 'deduce p248theory for NoWallFall (x2=0)'
                    endif
                endif                    
            endif                
        endif

        
        !COLMATATION processes: 1.thalwegization, 2.excessExpulsion, 3.cmx_dwcMitig or 4.cmx_dhsf(dwc0).
        if (j1.lt.nR) then
            !thalwegization (if so, then no chance of colmat excess expulsion: morphology will be reset to thw):
!            relerr_ups = ( hbf_thw - (h_sf_fut - hs_fut_ups))/hbf_thw !'dam' colmat if hbfreq<thw (relerr_ups>0)            
            hbf_fut = h_sf_fut-h_s_fut
            cond_inc = (dhs .lt. 0)
            cond_hbf = hbf_fut .lt. hbf_thw
            cond_colm_gt2 = ( (colm .eq. 5) .or. (colm .eq. 0) ) .and. cond_hbf
                !9jun20; colm=0: eqns deduced diff/ for colmat and incision procs.
            cond_hsfhs = cond_inc .and. cond_colm_gt2   !9jun20; p252d: hbf_fut<thw as dhsf<dhs<0.
                !as local supply is not widening but floodplain erosion, FL-CH interchange > inertia2channel.
            dhs_adm = (h_sf_fut-hbf_thw) - h_s   !new: 9jun20; p251. 
                !9jun20; p250d if cond_hbf then too low hbf_fut so enlarge it.
!            cond_dhsAdm = cond_colm_gt2 .and. (dhs_adm .ge. 0)    !9jun20; p253u.
            if (cond_hbf) then
                print*, 'as cond_hbf=T; check if goes to colm1 or 2; h_sf_fut,hbf_thw,h_s,dhs_adm:', &
                    h_sf_fut,hbf_thw,h_s,dhs_adm
                print*, 'cond_colm_gt2,cond_inc:',cond_colm_gt2,cond_inc
            endif            
            if (  ( cond_colm_gt2 .and. (dhs_adm .lt. 0) ) .or. cond_hsfhs .or. (h_s .eq. h_sf_fut) .or. &
                (h_sf_fut - hs_fut_ups .lt. hbf_thw) .or. tooMuchIncis  ) then 
                !7jun20;  deprec: .or. tooMuchIncis. 8jun20: dedeprec
                !20may20; ~0 instead hbf_thw in 2nd cond, p224dl. 
                !19may20; 1st cnd(wallExhaust): see p220c.
                !why hs_ups&loc and why 'dam': p220ur.
                !note no variable in this cnd was altered by prev 'excess' if.
                !see notes p121 for thalwegization and 
                    !p215dl & p216 for widening and expulsion of exceess above hbf_thw.
                !deprec: ...h_s_fut.gt.h_sf_fut....15may20:  .and. w_c_fut.lt.wv
                !14may20; thalwegization if too shallow flow, i.e. if deposition leads to hbf<hbf_thw
                !whole channel colmatation: crosssection is flattened across section, so it seeds artificial thalweg
                colm=1
                if (print_cnd) print*, 'in get_hs; h_s_fut.gt.h_sf_fut, so whole channel colmatation'
                A_trap2 = .5*(w_c + w_b)*h_bf !10may20; due to errBalance: *h_bf2
                hs_fld = h_sf + (net_Qs_flow*dt/dx - A_trap2) / wv !!10may20; due to errBalance: +hsf_fut (p202cr)
                    !15may20; hs_fld is calculated from reset, it did not feel widening EXHAUSTwall correction (hsf>=h_s).
                    !Instead, this correction must have affected dhs by lowering dhsf_dt so maybe avoiding colmatation.
                if (print_cnd) print*, '  w_c,w_b,h_bf,A_trap2', w_c,w_b,h_bf,A_trap2            
                if (print_cnd) print*, '  h_sf,net_Qs_flow,dt,dx,wv', h_sf,net_Qs_flow,dt,dx,wv
                !flat section:
                h_s_fut = hs_fld !h_s + h_bf + hs_fld !10may20; due to errBalance: +h_bf2
                h_sf_fut = hs_fld !h_sf_fut + hs_fld
                !artificial thalweg:rai
    !            wb_min = .01
                zch = zc(i1,j1,wv, ta, hbf_thw, wb_min)
                h_s_fut = h_s_fut - zch
                h_sf_fut = h_sf_fut + (hbf_thw - zch)
                if (h_sf_fut-h_s_fut .lt. .0001) then
                    print*, ' in thwColmat: why h_sf_fut=h_s_fut?; h_sf_fut,h_s_fut:', &
                        h_sf_fut,h_s_fut
                    stop                        
                endif                 
!                if (colm>1) then
!                    print*, '!!take this t,x for doubleColmat DEBUG: i1,j1:', i1,j1
!                    print*, 'h_sf_fut,hs_fut_ups,hbfReq,hbf_thw:',h_sf_fut,hs_fut_ups,h_sf_fut-hs_fut_ups,hbf_thw
!                    stop 'eh, after total colmat (also with localSuppl), also colm due to Qsu?'
!                endif
                w_b_fut = wb_min
                w_c_fut = w_b_fut + 2*hbf_thw/ta
                print*, ' colm=1; hs_fld,hbf_thw,zch:', hs_fld,hbf_thw,zch
                print*, '  w_b, w_c, h_s, h_sf:', w_b, w_c, h_s, h_sf
                print*, '  w_b_fut, w_c_fut, h_s_fut, h_sf_fut:', &
                    w_b_fut, w_c_fut, h_s_fut, h_sf_fut
            endif
            cond_colm_gt2 = (colm .eq. 0) .and. cond_hbf !9jun20; colm=0: eqns deduced diff/ for colmat and incision procs
                !evaluated again here as maybe colm changed to 1.
            !excess expulsion, if sed concentration is not debris flow.
!            cond_wb = w_b_fut .lt. wb_min !this happens only if colm=-1 or -2.
            if ( cond_colm_gt2 .and. (dhs_adm .ge. 0) .and. .not. cond_inc ) then
                !9jun20; colm0 avoids to come here if previously colm1.
                !7jun20; deprec: (colm .eq. 0) .and. (revived as colm-3 will solve incis.)
                !7jun20; p246d; if too small (hztal|vertical) channel evacuate sed excess as its source is local (inertia).
                    !sed sources if colm<0: bankPlanarFail in paralel incision with wbmin (colm=-1), 
                    !or wall fall at colm=-2.
                !20may20; ~0 instead hbf_thw in 2nd cond, p224dl...DISCARDED: Aw_NaN
                if (colm .eq. 1) then
                    print*, 'pre stop: cond_colm_gt2,dhs,cond_inc:  ', cond_colm_gt2,dhs,cond_inc
                    stop 'why are you here at colm2; see cnd to enter!'
                endif
                !if Qsu+widening&floodplain feeds so much, leads2fill channel:
                !deduction: p216u.
!                colm_pre2 = colm !depr 7jun20; go to colm-3.
!                if (dhs_adm<0 .and. dhs_adm/hbf_thw>-.03) dhs_adm=0 !not too negative (5% tol)
                    !19may20; this conditional is to avoid ever decresing wb, which goes in t too left from wbmin.
                ddhs = 0; ddwc = 0
!                if (cond_wb) then !deprec by 7jun20; go for colm=-3.
!                    dwc_adm = wb_min - w_b 
!                    if (dwc_adm>0) stop 'dwc_adm<0'
!                    dwc = w_c_fut - w_c
!                    ddwc = dwc - dwc_adm                
!                    if (ddwc<0) stop 'ddwc<0'
!                endif
                colm=2; x6 = 0; Arw = 0; Qsp3 = 0
                    !20may20: to avoid having hbf at colmat limit, deprec: dhs_adm = h_bf - hbf_thw
!                    dhs_adm = h_bf - hbf_thw !7jun20; recovered.                 
!                    dhs = h_s_fut - h_s !7jun20; no need, as dhs comes clean from 1st exner due to colm0 condition.
                ddhs = dhs - dhs_adm
                Qsp2 = psolid * dvol_dhs * ddhs / dt  !p216u.
!                h_s_fut = h_sf_fut-hbf_thw !7jun20; -h_bf.    !20may20, note 'updtHc p224ul'; deprec: hbf_thw
                    !7jun20; note h_sf_fut was updated by any dhsf_type if colm<0; e.g by fld tao and/or wall fall.
                    !deprec (19may20) as uses hsf (no fut) and larger error propag assumed: h_s + dhs_adm. See p221cr.
                w_b_fut = w_c_fut-2*hbf_thw/ta !7jun20; -2*h_bf/ta. !20may20, note 'updtHc p224ul'; deprec: - 2*hbf_thw/ta
                h_s_fut = h_s + dhs_adm !9jun20; 16:00, p254. 
                if (w_b_fut<wb_min) then !9jun20; p251d
                    colm=5 !rewiden (rw) to wb min and evacuate the consequent new excess. Work on hs set by colm2 (thwg)
                    x6 = wb_min-w_b_fut
                    Arw = 2*x6*hbf_thw
                    Qsp3 = Arw*dx/dt
                    w_b_fut = wb_min
                    w_c_fut = w_b_fut + 2*hbf_thw/ta
                endif
                Qspp = Qsp2 + Qsp3  !9jun20; p251; add vertical and horizontal sed excesses.
                net_Qs_flow = (Qsu - Qs - Qspp) / psolid
                    !as new hs is higher: wb>min as wb grows when up in trapeze.
                c_add = Qspp/Q !19may20; in hg.f90, there is another var called cpp.
                c_cap = Qs/Q !!19may20; due to capacity eqn.                
                print*, ' colm=2 or 5; colm, dt, Q:', colm,dt,Q
                print*, '  dhs,dhs_adm,dhsf_orig(tao):', dhs,dhs_adm,dhsf_orig
                print*, '  if colm5 x6>0 and Qsp3>0; x6,Qsp2,Qsp3:', x6,Qsp2,Qsp3
                print*, '  c_cap,c_add:', c_cap,c_add
                print*, '  h_bf,hbf_fut(colm0),hbf_thw:', h_bf,hbf_fut,hbf_thw
                print*, '  w_b, w_c, h_s, h_sf:', w_b, w_c, h_s, h_sf
                print*, '  w_b_fut, w_c_fut, h_s_fut, h_sf_fut:', &
                    w_b_fut, w_c_fut, h_s_fut, h_sf_fut
                if (h_s_fut .eq. h_sf_fut) then
                    print*, '!!take this t,x to edit print_range for colmatQs++ DEBUG: i,j:', i1,j1
                    print*, 'why are the same h_s_fut,h_sf_fut:?', h_s_fut,h_sf_fut
                    stop
                endif
                if (w_b_fut .le. .99*wb_min) then !19may20; 3%tol..7jun20:0% !8jun20; tol1%
                    print*, '!!take this t,x to edit print_range for colmatQs++ DEBUG: i,j:', i1,j1
                    print*, 'why w_b_fut < wb_min:?', w_b_fut,wb_min
                    stop
!                elseif (w_b_fut .le. wb_min) then !19may20; eh, why dhs_adm<0 leading to this?.Depr:SeekDefiSln!
!                    w_b_fut=wb_min !if small error (<3% below wbmn) then reset to avoid cumul(t) bug.
                endif   
                if (c_add + c_cap > cmx) then !near debris flow. 
                            !p216d; occur for almost any sed tpt qsx, as qsy/qsx~1, as dvoldhs high as 
                            !dx/w=1000 instead 1 (as flume calib).
                    colm=3 !p216cr; see note on colm3 at beginning of colm-2; on p254cl 'latent bug'.
                    !see p217/8: real widening is localized; this code widens few as c allows, cos it's mean reach w.
                    Qsmx = cmx*Q !19may20; as cmx limit was set from qs calc at get_capac(), deprec: max(Qs, cmx*Q)
                    Qspp = Qsmx-Qs
                    if (Qspp/Qsmx .le. -1e-2) then !too negative (1%tol)
                        print*, '!!take this t,x to edit print_range for colmatQs++ DEBUG: i,j:', i1,j1
                        print*, 'why Qsmx<Qs; Qsmx,Qs:?', Qsmx,Qs
                        stop
                    endif
                    dhs_adm_dt = dhs_adm/dt
                    net_Qs_flow = (Qsu - Qsmx) / psolid
                    dQs_hsAdm = dvol_dhs*dhs_adm_dt
                    dwc_colm2 = w_c_fut - w_c
                    dwc = ( net_Qs_flow - Qs_hsf_orig - dQs_hsAdm )*dt/dvol_dwc
                        !cmx mitigates dwc. P218.                    
                    w_c_fut = w_c + dwc
                    w_b_fut = w_c_fut - 2*hbf_thw/ta !20may20, note 'updtHc p224ul'; deprec: - 2*hbf_thw/ta
                        !as new hs is higher: wb>min as wb grows when up in trapeze.
                    print*, ' colm=3; Qsu,Qs,Qsmx,Qspp:', Qsu,Qs,Qsmx,Qspp
                    print*, '  dwc terms; dwc_exner,net_Qs_flow,Qs_hsf_orig,dQs_hsAdm,dvol_dwc:', &
                        dwc,net_Qs_flow,Qs_hsf_orig,dQs_hsAdm,dvol_dwc
                    net_Qs_flow_colm2 = (Qsu - (c_add + c_cap)*Q) / psolid
                    dwc_exn_colm2 = ( net_Qs_flow_colm2 - Qs_hsf_orig - dQs_hsAdm )*dt/dvol_dwc
                    print*, '  exner check, if c>mx (colm2) dwcExner,dwc_colm2:', &
                        dwc_exn_colm2, dwc_colm2
                    print*, '  colm2 dwc terms; dwc_exn_colm2,net_Qs_flow,Qs_hsf_orig,dQs_hsAdm,dvol_dwc:', &
                        dwc_exn_colm2,net_Qs_flow_colm2,Qs_hsf_orig,dQs_hsAdm,dvol_dwc
                    print*, '  w_b, w_c, h_s, h_sf:', w_b, w_c, h_s, h_sf
                    print*, '  w_b_fut, w_c_fut, h_s_fut, h_sf_fut:', &
                        w_b_fut, w_c_fut, h_s_fut, h_sf_fut
                    if (dwc<0) then !stop 'dwc<0'
                        print*, '   19may20; as cmx implies dwc<0, lets restrain FL erosion dhsf here in colm4:'
                        colm=4 !p221d.
                        dhsf = ( net_Qs_flow - dQs_hsAdm )*dt/dvol_dhsf !dwc=0
                            !cmx short to gain dwc, so cmx also mitigates dhsf
                        h_sf_fut = h_sf + dhsf
                        w_c_fut = w_c
                        w_b_fut = w_c_fut - 2*hbf_thw/ta !20may20, note 'updtHc p224ul'; deprec: - 2*hbf_thw/ta
                        h_s_fut = h_sf_fut - hbf_thw  !9jun20; hbf_thw instead h_bf. !6jun20; as h_s_fut was updated
                        dhsf_shear= dhsf_dt*dt
                        print*, ' colm=3; dhsf_orig,dhsfPossible:',dhsf_shear,dhsf
                        print*, '  w_b, w_c, h_s, h_sf:', w_b, w_c, h_s, h_sf
                        print*, '  w_b_fut, w_c_fut, h_s_fut, h_sf_fut:', &
                            w_b_fut, w_c_fut, h_s_fut, h_sf_fut
                        if (dhsf_orig>dhsf) stop 'no sense: colm4 must make less negat dhsf from that due to shear'
                        if (h_sf_fut<0) stop 'hsf_fut<0?'
!                        print*,'check this colm4 pa i1,j1:', i1,j1
!                        stop  
                    endif
                endif
            endif
        endif
        condNan = (w_b_fut .ne. w_b_fut) .or. (w_c_fut .ne. w_c_fut) .or. &
            (h_s_fut .ne. h_s_fut) .or. (h_sf_fut .ne. h_sf_fut)   
        if (condNan) stop 'a future geom variable in NaN'
        dhs_dt = (h_s_fut - h_s) / dt
        dhsf2_dt = (h_sf_fut - h_sf) / dt
        dwb_dt = (w_b_fut - w_b) / dt
        dwf_dt = (w_c - w_c_fut) / dt   !',t' being t deriv; wf,t=-wc,t as valley wv is constant and wf=wv-wc.
        hbf_fut = h_sf_fut - h_s_fut
!        if (hbf_fut .lt. .99*hbf_thw) then !20may20; p224d so deprec.
!            print*, '!!take this t,x to edit print_range for hbf_fut DEBUG: i,j:', i1,j1
!            print*, ' why <thw (1%tol)?: hbf_fut,hbf_thw,colm:', hbf_fut,hbf_thw,colm
!            stop    
!        elseif (hbf_fut .lt. hbf_thw) then 
!            !if err <1% below, to avoid cumulative wb lefting bug leading to <wbmin. p220d.
!            h_sf_fut = h_s_fut + hbf_thw
!        endif
        if (print_cnd) print*, 'in get_hs; after colmat check; w_b,w_c,h_s,h_sf:', &
            w_b,w_c,h_s,h_sf
        if (print_cnd) print*, 'in get_hs; after colmat check; w_b_fut,w_c_fut,h_s_fut,h_sf_fut:', &
            w_b_fut,w_c_fut,h_s_fut,h_sf_fut
        if (print_cnd) print*, 'compare lines above against boundaries wv,hbf_thw:', wv,hbf_thw
        if (h_sf_fut<h_s_fut) then
            print*, '!!take this t,x to edit print_range for DEBUG: i,j:', i1,j1
            print*, 'h_sf_fut=', h_sf_fut, '<h_s_fut=', h_s_fut
            print*, '15may20; due to need to thalwegization with wc=wv; exhaustWall'
            print*, '8jun20; in colm-2 due to excesive dhs+'            
            stop 'why?'
        endif            
        if (w_b<.99*wb_min) then
            print*, '!!take this t,x to edit print_range for DEBUG: i,j:', i1,j1
            print*, 'wb=', w_b, '<wb_min=', wb_min
            stop 'why?'
        endif
!        if (w_c<.98) then
!            print*, '!!take this t,x to edit print_range for DEBUG: i,j:', i1,j1
!            print*, 'w_c_fut: ', w_c_fut 
!            stop 'why w_c_fut<wv=1 in trenchExp??' !5may20: deprecated stop, so code works for any flume; 
!                !Anyway, remain aware in wc(t) trench plots at locations/abscissas.
!        endif
    end subroutine
    
    real function As(w_v, h_sf, w_c, w_b, h_s)
        real, intent(in) :: w_v, h_sf, w_c, w_b, h_s
        real :: h_bf
        h_bf = h_sf - h_s
        As = w_v*h_sf - .5*(w_c + w_b)*h_bf
    end function

    subroutine check_wat_balance(i1,j1, q_dl, tcum,cumWat_m3,cumWatErr_m3,print_cnd, isFlood, Qhc, hc, mp, dt, &
        dx_dt, cpp, Aw, err_m3)
        logical, intent(in) :: print_cnd, isFlood !, lowQ_CH_narrower
        integer, intent(in) :: i1,j1  
        real, intent(in) :: Qhc, dt, dx_dt, cpp, Aw, tcum, mp(9), hc(7),cumWat_m3,cumWatErr_m3, q_dl
        real, intent(out) :: err_m3
        real Qhg, relerr, smooth_rat, vm, tol, dQ
        if (q_dl<.1) then    !out of range for recking11 friction law, so v often < 1cm/s unrealistic.
            tol = .5 !50% tol!, just for calc to continue. Check rating v vs Q if the too low v are rare.
        else
            if (tcum .lt. 3e7) then !6jun20; first month channel hydraulics is 'warming'. 
                !width~1e1m, dx~1e3, h~1e0, Qh~1e-2 takes ~10 days for channel to fill itself. 
                !Often h ~1e-2m! (too shallow). 10jun20; lets give 1yr~3e7s
                tol=.3
            else
                tol=.2
            endif
        endif
        vm = Qhc/Aw
        smooth_rat = vm / (vm + dx_dt)
        Qhg = cpp*smooth_rat
        dQ = Qhc - Qhg
        relerr = dQ / Qhg
        err_m3 = dQ*dt
        if (print_cnd) print*, 'smooth_rat, watBalance relerr=', smooth_rat,relerr
        if (abs(relerr) .gt. tol) then
            print*, 'take this t,x to edit print_range for watBal DEBUG: i,j:', i1,j1
            print*, 'hydraulics now: S,w,h,hbf,vc,vf,Frc:', hc(1),hc(2),hc(3),hc(4),hc(5),hc(6),hc(7)
            print*, 'morphol now: colm, w_b, w_c, h_s, h_sf:', mp(1),mp(2),mp(3),mp(4),mp(5)
            print*, 'morphol future: w_b_fut, w_c_fut, h_s_fut, h_sf_fut:', &
                mp(6),mp(7),mp(8),mp(9)
            print*, 'tcum,tcum1month(highTol),tdt, dx_dt: ', tcum,2.592e6,dt,dx_dt
            print*, 'vm, smooth_rat, Q=Qhc, Qhg, relerr: ', vm, smooth_rat, Qhc, Qhg, relerr
!            print*, 'lowQ_CH_narrower: ', lowQ_CH_narrower
            print*, '%err in cumul water balance; cumWatErr_m3/cumWat_m3:', cumWatErr_m3/cumWat_m3
            stop 'please check water balance'
        endif
    end subroutine
        
    subroutine check_sed_balance(i1,j1,print_cnd, colm, wv, ta, dx, net_sedVolFlow, &
        dt, Q, psolid, &
        h_sf, w_c, w_b, h_s, &
        h_sf_fut, &
        w_c_fut, w_b_fut, h_s_fut, dvol)
        logical, intent(in) :: print_cnd 
        integer, intent(in) :: i1,j1           
        real, intent(in) :: wv, ta, dx, net_sedVolFlow, h_sf, w_c, w_b, h_s, psolid
        real, intent(in) :: h_sf_fut, w_c_fut, w_b_fut, h_s_fut, colm, dt, Q
        real, intent(out) :: dvol
        real vol_prev, vol_fut, relerr
        vol_prev = dx*As(wv, h_sf, w_c, w_b, h_s)
        vol_fut = dx*As(wv, h_sf_fut, w_c_fut, w_b_fut, h_s_fut)
        dvol = vol_fut - vol_prev
        relerr = (dvol - net_sedVolFlow) / vol_prev
        if (print_cnd) print*, 'dvol,sedbalance relerr=', dvol,relerr
        if (abs(relerr) .gt. .1) then
            print*, 'take this t,x to edit print_range for sedBal DEBUG: i,j:', i1,j1
            print*, 'geom: h_s, h_sf, w_b, w_c.   wv, colm: ', wv, colm
            print*, 'geomNow; ', h_s, h_sf, w_b, w_c
            print*, 'geomFut; ', h_s_fut, h_sf_fut, w_b_fut, w_c_fut
            print*, 'vol_fut,vol_prev,net_sedVolFlow,relerr', vol_fut,vol_prev,net_sedVolFlow,relerr
            print*, 'Q,dt', Q,dt
            stop 'please check sediment balance'
        endif
        dvol = dvol*psolid !2jun20; export to plot only alluvium store change of solid sediment, with no voids.
    end subroutine
    
!    subroutine calc_pondBackwat() !steps: phdSCnotes p158.
!    end subroutine
end module
