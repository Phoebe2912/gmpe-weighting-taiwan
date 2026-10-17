gmpe_models <- function(tectonic) {
  if (tectonic == "cru") list(
    ASK14 = function(r, k, tp)
      ASK14(r$M[k], r$Rrup[k], r$Rjb[k], tp, r$Dip[k], r$ftype[k], r$width[k], r$Vs30[k], r$hwflag[k],
            r$Ztor[k], r$Vs30_class[k], r$Z1.0[k], r$Rx[k], r$Ry[k], regionflag = 1, msasflag = 0),
    BSSA14 = function(r, k, tp)
      BSSA14(r$M[k], r$Rjb[k], tp, r$Vs30[k], r$ftype[k], r$Z1.0[k], regionflag = 0, r$basinflag[k]),
    CB14 = function(r, k, tp)
      CB14(r$M[k], r$Rrup[k], r$Rjb[k], r$ftype[k], tp, r$Vs30[k], r$Ztor[k], r$Z2.5[k], r$Dip[k],
           r$depth[k], hwflag = 0, r$Rx[k], r$width[k], regionflag = 0),
    CY14 = function(r, k, tp)
      CY14(r$M[k], r$Rrup[k], r$Rjb[k], tp, r$Vs30[k], r$Dip[k], r$Ztor[k], r$ftype[k], r$Z1.0[k],
           r$Vs30_class[k], hwflag = 0, r$Rx[k], regionflag = 0),
    I14 = function(r, k, tp)
      I14(r$M[k], r$Rjb[k], tp, r$ftype[k], r$Vs30[k]),
    PhungCru20 = function(r, k, tp)
      PhungCru18(r$M[k], r$Rrup[k], r$Rjb[k], r$Rx[k], tp, r$Vs30[k], r$Dip[k], r$Ztor[k], r$ftype[k],
                 r$Z1.0[k], r$Vs30_class[k], regionflag = 1, r$hwflag[k]),
    Chao20 = function(r, k, tp)
      Chao2018(r$M[k], r$Rrup[k], tp, r$ftype[k], r$Vs30[k], r$Vs30_class[k], r$Ztor[k], r$Z1.0[k],
               r$sourcetype[k], r$msasflag[k]),
    Lin09 = function(r, k, tp)
      Lin09(r$M[k], r$Rrup[k], tp, r$Vs30[k], r$ftype[k]),
    Lin11 = function(r, k, tp)
      Lin11(r$M[k], r$Rrup[k], tp, r$hwflag[k], r$Vs30[k])
  ) else list(
    AG2020 = function(r, k, tp)
      AG2020(r$M[k], r$Rrup[k], tp, r$ftype[k], r$Vs30[k], r$Z2.5[k], r$Ztor[k], 7, 1),
    KBCG2020 = function(r, k, tp)
      KBCG2019(r$M[k], r$Rrup[k], tp, r$ftype[k], r$Vs30[k], r$Z2.5[k], mbInter = 7.1, mbSlab = 7.7,
               7, r$Ztor[k], r$Z1.0[k], 0),
    PSHAB2022 = function(r, k, tp)
      PSHAB2019(r$M[k], r$Rrup[k], tp, r$ftype[k], r$Vs30[k], r$Z2.5[k], r$depth[k],
                mbInter = 7.1, mbSlab = 7.7, 9, 0),
    SI2022 = function(r, k, tp)
      SI2020(r$M[k], r$Rrup[k], tp, r$ftype[k], r$Vs30[k], r$Z2.5[k], r$depth[k], r$Dmoho[k]),
    Chao20 = function(r, k, tp)
      Chao2018(r$M[k], r$Rrup[k], tp, r$ftype[k], r$Vs30[k], r$Vs30_class[k], r$Ztor[k], r$Z1.0[k],
               r$sourcetype[k], r$msasflag[k]),
    PhungSub20 = function(r, k, tp)
      PhungSub18(r$M[k], r$Rrup[k], tp, r$ftype[k], r$Vs30[k], r$Z1.0[k], r$Ztor[k], 1),
    LL08 = function(r, k, tp)
      LL08(r$M[k], r$Rrup[k], r$depth[k], r$ftype[k], tp, r$Vs30[k])
  )
}

gmpe_distance <- c(BSSA14 = "Rjb", I14 = "Rjb")
