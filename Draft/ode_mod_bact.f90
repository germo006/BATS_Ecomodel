
    !-----------------------------------------------------------------------
    !      Bacterial Processes
    !-----------------------------------------------------------------------
    ! 1. Gross Grow
    ! Maximum possible C amount for bacterial use
    ALC = y(iLDOMc)
    !temp = MIN(c1, exp(b_SDONlabi * (y(iSDOMn)/y(iSDOMc)/q_BA_n - c1)) )
    !temp = MIN(temp, exp(b_SDOPlabi * (y(iSDOMp)/y(iSDOMc)/q_BA_p - c1)) )
    ASC = y(iSDOMc) * r_SDOM
    ! Carbon Usage
    Nfunc_ba_n = y(iBAn)/y(iBAc)/q_BA_n
    Nfunc_ba_p = y(iBAp)/y(iBAc)/q_BA_p
    temp = min(Nfunc_ba_n, Nfunc_ba_p)
    temp = min(temp, c1)
    growBAldoc = mu_BA * y(iBAc) * temp * ALC &
                          / (ALC+ k_DOM + ASC) 
    growBAsdoc = mu_BA * y(iBAc) * temp * ASC &
                          / (ASC+ k_DOM + ALC)
    ! DON and DOP usage
    growBAldon = growBAldoc/y(iLDOMc)*y(iLDOMn) ! available labile N
    growBAldop = growBAldoc/y(iLDOMc)*y(iLDOMp) ! available labile P
    growBAsdon = growBAsdoc * min(q_BA_n, &
               (y(iSDOMn)/y(iSDOMc) + f_BAslct/Nfunc_ba_n*(q_BA_n-y(iSDOMn)/y(iSDOMc))))
    growBAsdop = growBAsdoc * min(q_BA_p, &
               (y(iSDOMp)/y(iSDOMc) + f_BAslct/Nfunc_ba_p*(q_BA_p-y(iSDOMp)/y(iSDOMc))))

    ! inorganic nutrients uptake
    growBAnh4 = growBAldon / y(iLDOMn) * y(iNH4) * min(c1, 1/Nfunc_ba_n)
    if (Nfunc_ba_n.lt.c1) then
        growBAno3 = min(0.1 * (growBAldon + growBAsdon) / (y(iLDOMn) + y(iSDOMn)) &
                            * y(iNO3) * min(c1, 1/Nfunc_ba_n), &
			    (growBAldon + growBAsdon) / (y(iLDOMn) + y(iSDOMn)) &
                            * (y(iNO3)  + y(iNH4)) - growBAnh4)
        growBAno3 = max(c0, growBAno3)
    else
        growBAno3 = c0
    end if
    growBApo4 = growBAldop / y(iLDOMp) * y(iPO4) * min(c1, 1/Nfunc_ba_p)
    ! Bacteria gross growth
    growBAc = growBAldoc + growBAsdoc
    growBAn = growBAldon + growBAsdon + growBAnh4 + growBAno3
    growBAp = growBAldop + growBAsdop + growBApo4
    ! 2. respiration
    respBA = r_BAresp_1 * y(iBAc) + zeta * growBAno3 + &
        (r_BAresp_min + (r_BAresp_max - r_BAresp_min)*EXP(-b_BAresp*growBAc) ) &
        * growBAc
    ! 3. excreting refractory DOM
    refrBAc = r_BArefr * y(iBAc)
    refrBAn = q_refrDOM_n * refrBAc
    refrBAp = q_refrDOM_p * refrBAc
    ! 4. excreting semi-labile DOM and regenerating DIN
    IF ( (y(iBAc) < y(iBAn)/q_BA_n) .AND. &
         (y(iBAc) < y(iBAp)/q_BA_p) ) THEN  !Cabon in short
         excrBAc = c0
         excrBAn = c0
         excrBAp = c0
         remiBAn = r_BAremi * (y(iBAn) - y(iBAc) * q_BA_n)
         remiBAp = r_BAremi * (y(iBAp) - y(iBAc) * q_BA_p)
    ELSE IF ( (y(iBAc) > y(iBAn)/q_BA_n) .AND. &
         (y(iBAp)/q_BA_p > y(iBAn)/q_BA_n) ) THEN !Nitrogen in short
         excrBAc = r_BAadju * (y(iBAc) - y(iBAn)/q_BA_n)
         excrBAn = c0
         excrBAp = r_BAadju * (y(iBAp) - y(iBAn)/q_BA_n * q_BA_p)
         remiBAn = c0
         remiBAp = c0
    ELSE !Phosphorus in short
         excrBAc = r_BAadju * (y(iBAc) - y(iBAp)/q_BA_p)
         excrBAn = r_BAadju * (y(iBAn) - y(iBAp)/q_BA_p * q_BA_n)
         excrBAp = c0
         remiBAn = c0
         remiBAp = c0
    END IF
    !6. removal by grazing
    grazBAc = mu_PRT * y(iPRTc) * y(iBAc) * y(iBAc) &
        / (y(iBAc) * y(iBAc) + g_ba * g_ba + &
           y(iSPc)*y(iSPc)/g_sp/g_sp*g_ba*g_ba + &
           y(iUNc)*y(iUNc)/g_un/g_un*g_ba*g_ba)
    grazBAn = grazBAc / y(iBAc) * y(iBAn)
    grazBAp = grazBAc / y(iBAc) * y(iBAp)
    !6b. Mortality due to viruses
    mortBAc = r_BAmort * y(iBAc)
    mortBAn = r_BAmort * y(iBAn)
    mortBAp = r_BAmort * y(iBAp)
    !7. BA Derivs
    dydtt(iBAc) = (growBAc - refrBAc - excrBAc                   - grazBAc &
                   - respBA - mortBAc)/ SecPerDay
    dydtt(iBAn) = (growBAn - refrBAn - excrBAn - remiBAn - grazBAn - mortBAn) &
                   / SecPerDay
    dydtt(iBAp) = (growBAp - refrBAp - excrBAp - remiBAp - grazBAp - mortBAp) &
                   / SecPerDay
    !8. Flux of inorganic nutrients through bacteria
    fluxBAnh4 = growBAnh4 - remiBAn
    fluxBAno3 = growBAno3
    fluxBApo4 = growBApo4 - remiBAp
