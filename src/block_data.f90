!  block data:  initiallizes some physical constants, and also the
!  constants used by the conductive opacity regimes.  the numbers
!  have been taken from itoh et al (1983; 1984). 
!  Also added Itoh and company Bremsstrahlung neutrino coefficients.
!  See Itoh etal 1989, ApJ, 339, 354. + erratum 360, 741.
!***********************************************************************

module misc_const

  real*8, parameter :: an0 = 6.022169d23
  real*8, parameter :: bk_const = 1.380622d-16
  real*8, parameter :: ck = 1.5d0
  real*8, parameter :: pi = 3.14159265358979d0
  real*8, parameter :: fnat = 4.342945d-1
  real*8, parameter :: amsun = 1.989d+33    !Solar mass in g
  real*8, parameter :: alsun = 3.90d+33
  real*8, parameter :: g = 6.67259d-8 
  real*8, parameter :: mp = 1.67d-24
  real*8, parameter :: elec = 4.8d-10
  real*8, parameter :: c = 2.997925d10
  real*8, parameter :: a = 7.56566d-15 

end module misc_const

!***********************************************************************
! common /hydatl/ alhy(4),blhy(3),clhy(3),dlhy(4),elhy(3),flhy(3) 

module hydatl
  
  real*8, dimension(4) :: alhy, dlhy
  real*8, dimension(3) :: blhy,clhy,elhy,flhy
  
end module hydatl

!***********************************************************************
! common /hedatl/ alhe(4),blhe(3),clhe(3),dlhe(4),elhe(3),flhe(3) 

module hedatl

  real*8, dimension(4) :: alhe, dlhe
  real*8, dimension(3) :: blhe,clhe,elhe,flhe

end module hedatl

!***********************************************************************
! common/hedats/ sahe(5),sbhe(3),sehe(5),sfhe(3),sihe(5),sjhe(3), 
!                sphe(5), sqhe(3), alphahe(4), betahe(4), sche,
!                sdhe, sghe, shhe, skhe, slhe, srhe, sshe
 
module hedats
  
  real*8, dimension(5) :: sahe, sehe, sihe, sphe
  real*8, dimension(4) :: alphahe, betahe
  real*8, dimension(3) :: sbhe, sfhe, sjhe, sqhe
  real*8 :: sche, sdhe, sghe, shhe, skhe, slhe, srhe, sshe

end module hedats

!************************************************************************
! the following are the fitting coefficients quoted by itoh 
! "cdats" ==> carbon, solid state
! "cdatl" ==> carbon, liquid metal state
! "odats" ==> oxygen, solid state
! "odatl" ==> oxygen, liquid metal state
!  and so on. 
!************************************************************************

! common/cdats/ sac(5), sbc(3), sec(5), sfc(3), sic(5), sjc(3),
!               spc(5), sqc(3), alphac(4), betac(4), scc,
!               sdc, sgc, shc, skc, slc, src, ssc 

module cdats

  real*8, dimension(5) :: sac, sec, sic, spc
  real*8, dimension(4) :: alphac, betac
  real*8, dimension(3) :: sbc, sfc, sjc, sqc
  real*8 ::  al, scc, sdc, sgc, shc, skc, slc, src, ssc 

end module cdats

! common /cdatl/ alc(4), blc(3), clc(3), dlc(4), elc(3), flc(3)

module cdatl

  real*8, dimension(4) :: alc, dlc
  real*8, dimension(3) :: blc, clc, elc, flc

end module cdatl

! common/odats/ sao(5), sbo(3), seo(5), sfo(3), sio(5), sjo(3),
!               spo(5), sqo(3), alphao(4), betao(4), sco,
!               sdo, sgo, sho, sko, slo, sro, sso 

module odats
  
  real*8, dimension(5) :: sao, seo, sio, spo
  real*8, dimension(4) :: alphao, betao
  real*8, dimension(3) :: sbo, sfo, sjo, sqo
  real*8 :: sco, sdo, sgo, sho, sko, slo, sro, sso 

end module odats

! common /odatl/ alo(4), blo(3), clo(3), dlo(4), elo(3), flo(3)

module odatl
  
  real*8, dimension(4) :: alo, dlo
  real*8, dimension(3) :: blo, clo, elo, flo

end module odatl
   
!**********************************************************************
! the following are the Bremsstrahlung neutrino rate fitting 
! coefficients quoted by Itoh. Common'd blocks are:
! Helium: soldath, phondath, liqdath, quandath
! Carbon: soldatc, phondatc, liqdatc, quandatc
! Oxygen: soldato, phondato, liqdato, quandato
!   Iron: soldatf, phondatf, liqdatf, quandatf 
!**********************************************************************

! common/soldath/ snah(5), snbh(3), sneh(5), snfh(3), snih(5), 
!                 snjh(3), snph(5), snqh(3), alphasnh(4), 
!                 betasnh(4), snch, sndh, sngh, snhh, snkh, 
!                 snlh, snrh, snsh

module soldath
  
  real*8, dimension(5) :: snah, sneh, snih, snph
  real*8, dimension(4) :: alphasnh, betasnh
  real*8, dimension(3) :: snbh, snfh, snjh, snqh
  real*8 :: snch, sndh, sngh, snhh, snkh, snlh, snrh, snsh

end module soldath

! common/phondath/ pnah(5), pnbh(3), pnih(5), pnjh(3), 
!                  alphapnh(4), betapnh(4), pnch, pndh, pnkh, pnlh

module phondath

  real*8, dimension(5) :: pnah, pnih
  real*8, dimension(4) :: alphapnh, betapnh
  real*8, dimension(3) :: pnbh, pnjh
  real*8 :: pnch, pndh, pnkh, pnlh

end module phondath

! common/liqdath/ lnah(6), lnbh(4), lneh(6), lnfh(4), lnih(6), 
!                 lnjh(4), lnph(6), lnqh(4), alphalnh(4), 
!                 betalnh(4), lnch, lndh, lngh, lnhh, lnkh, 
!                 lnlh, lnrh, lnsh
! common/quandath/ qnah(6), qnbh(4), qneh(6), qnfh(4), 
!                  alphaqnh(5), qnch, qndh, qngh, qnhh
! For these two common blocks, created a data block called bremliqh_data. Code 
! takes the variables directly from there.

! common/liqdatc/ lnac(6), lnbc(4), lnec(6), lnfc(4), lnic(6), 
!                 lnjc(4), lnpc(6), lnqc(4), alphalnc(4), 
!                 betalnc(4), lncc, lndc, lngc, lnhc, lnkc, 
!                 lnlc, lnrc, lnsc 
! common/quandatc/ qnac(6), qnbc(4), qnec(6), qnfc(4), 
!                  alphaqnc(5), qncc, qndc, qngc, qnhc
! For these two common blocks, created a data block called bremliqc_data. Code 
! takes the variables directly from there.

! common/liqdato/ lnao(6), lnbo(4), lneo(6), lnfo(4), lnio(6), 
!                 lnjo(4), lnpo(6), lnqo(4), alphalno(4), 
!                 betalno(4), lnco, lndo, lngo, lnho, lnko, 
!                 lnlo, lnro, lnso
! common/quandato/ qnao(6), qnbo(4), qneo(6), qnfo(4), 
!                  alphaqno(5), qnco, qndo, qngo, qnho
! For these two common blocks, created a data block called bremliqo_data. Code 
! takes the variables directly from there.

! common/soldatc/ snac(5), snbc(3), snec(5), snfc(3), snic(5), 
!                 snjc(3), snpc(5), snqc(3), alphasnc(4), 
!                 betasnc(4), sncc, sndc, sngc, snhc, snkc, 
!                 snlc, snrc, snsc 

module soldatc
  
  real*8, target, dimension(5) :: cb_snac, cb_snec, cb_snic, cb_snpc
  real*8, target, dimension(4) :: cb_alphasnc, cb_betasnc
  real*8, target, dimension(3) :: cb_snbc, cb_snfc, cb_snjc, cb_snqc
  real*8, target :: cb_sncc, cb_sndc, cb_sngc, cb_snhc, cb_snkc, cb_snlc, &
       cb_snrc, cb_snsc 

end module soldatc

! common/phondatc/ pnac(5), pnbc(3), pnic(5), pnjc(3), 
!                  alphapnc(4), betapnc(4), pncc, pndc, pnkc, pnlc

module phondatc

  real*8, target, dimension(5) :: cb_pnac, cb_pnic
  real*8, target, dimension(4) :: cb_alphapnc, cb_betapnc
  real*8, target, dimension(3) :: cb_pnbc, cb_pnjc
  real*8, target :: cb_pncc, cb_pndc, cb_pnkc, cb_pnlc

end module phondatc




! common/soldato/ snao(5), snbo(3), sneo(5), snfo(3), snio(5), 
!                 snjo(3), snpo(5), snqo(3), alphasno(4), 
!                 betasno(4), snco, sndo, sngo, snho, snko, 
!                 snlo, snro, snso

module soldato

  real*8, target, dimension(5) :: cb_snao, cb_sneo, cb_snio, cb_snpo
  real*8, target, dimension(4) :: cb_alphasno, cb_betasno
  real*8, target, dimension(3) :: cb_snbo, cb_snfo, cb_snjo, cb_snqo
  real*8, target :: cb_snco, cb_sndo, cb_sngo, cb_snho, cb_snko, cb_snlo, &
       cb_snro, cb_snso

end module soldato

! common/phondato/ pnao(5), pnbo(3), pnio(5), pnjo(3), 
!                  alphapno(4), betapno(4), pnco, pndo, pnko, pnlo

module phondato

  real*8, target, dimension(5) :: cb_pnao, cb_pnio
  real*8, target, dimension(4) :: cb_alphapno, cb_betapno
  real*8, target, dimension(3) :: cb_pnbo, cb_pnjo
  real*8, target :: cb_pnco, cb_pndo, cb_pnko, cb_pnlo

end module phondato

! common/sfodatf/ snaf(5), snbf(3), snef(5), snff(3), snif(5), 
!                 snjf(3), snpf(5), snqf(3), alphasnf(4), 
!                 betasnf(4), sncf, sndf, sngf, snhf, snkf, 
!                 snlf, snrf, snsf

module soldatf

  real*8, dimension(5) :: snaf, snef, snif, snpf
  real*8, dimension(4) :: alphasnf, betasnf
  real*8, dimension(3) :: snbf, snff, snjf, snqf
  real*8 :: sncf, sndf, sngf, snhf, snkf, snlf, snrf, snsf
  
end module soldatf

! common/phondatf/ pnaf(5), pnbf(3), pnif(5), pnjf(3), 
!                  alphapnf(4), betapnf(4), pncf, pndf, pnkf, pnlf

module phondatf

  real*8, dimension(5) :: pnaf, pnif
  real*8, dimension(4) :: alphapnf, betapnf
  real*8, dimension(3) :: pnbf, pnjf
  real*8 :: pncf, pndf, pnkf, pnlf

end module phondatf

! common/liqdatf/ lnaf(6), lnbf(4), lnef(6), lnff(4), lnif(6), 
!                 lnjf(4), lnpf(6), lnqf(4), alphalnf(4), 
!                 betalnf(4), lncf, lndf, lngf, lnhf, lnkf, 
!                 lnlf, lnrf, lnsf

module liqdatf

  real*8, dimension(6) :: lnaf, lnef, lnif, lnpf
  real*8, dimension(4) :: lnbf, lnff, lnjf, lnqf, alphalnf, betalnf
  real*8 :: lncf, lndf, lngf, lnhf, lnkf, lnlf, lnrf, lnsf

end module liqdatf

! common/quandatf/ qnaf(6), qnbf(4), qnef(6), qnff(4), 
!                  alphaqnf(5), qncf, qndf, qngf, qnhf

module quandatf  

  real*8, dimension(6) :: qnaf, qnef
  real*8, dimension(5) :: alphaqnf
  real*8, dimension(4) :: qnbf, qnff
  real*8 :: qncf, qndf, qngf, qnhf

end module quandatf
 
!**********************************************************************
!  variables used in cndlhy
!**********************************************************************

module cndlhy_data
  
  real*8, dimension(4) :: a, d
  real*8, dimension(3) :: b, c, e, f

  data a / 0.34542, -.28157, 0.09184, -.03734 /     !alhy
  data b / -.61919, 0.40004, -.16585 /              !blhy
  data c / 0.35742, -.41151, 0.21552 /              !clhy
  data d / 0.21512, -.10843, -.00596, -.00950 /     !dlhy
  data e / -.36667, 0.14040, -.04588 /              !elhy
  data f / 0.10493, -.09537, 0.04682 /              !flhy

end module cndlhy_data

!**********************************************************************
!  variables used in cndlhe
!**********************************************************************

module cndlhe_data

  real*8, dimension(4) :: a, d
  real*8, dimension(3) :: b, c, e, f

  data a / 0.62199, -.16110, 0.15574, -.02893 /  !alhe
  data b / -.65222, 0.48601, -.18266 /           !blhe
  data c / 0.36580, -.52176, 0.26240 /           !clhe
  data d / 0.36090, 0.02576, 0.05061, 0.00015 /  !dlhe
  data e / -.40559, 0.15316, -.04058 /           !elhe
  data f / 0.12140, -.11621, 0.04939 /           !flhe

end module cndlhe_data

!**********************************************************************
!  variables used in cndshe
!**********************************************************************

! NOTE: the values stored in a(1), e(1), i(1) and p(1) are 1/2 the
!  values listed in the paper.  this is because they are used in the
!  formulae in the paper as a(1)/2, etc.  (same for carbon and oxygen).

module cndshe_data
  
  real*8, dimension(5) :: sahe, sehe, sihe, sphe, a, e, i, p
  real*8, dimension(4) :: alphahe, betahe, alpha, beta
  real*8, dimension(3) :: sbhe, sfhe, sjhe, sqhe, b, f, j, q
  real*8, parameter :: c   = -.16902, d = 1.4525
  real*8, parameter :: ggg = -.20406, h = 1.36576
  real*8, parameter :: k   = -.22355, l = 1.70151
  real*8, parameter :: r   = -.25532, s = 1.90219

  data sahe       / -.01006, 0.05975, -.04290, -.00341, -.00338 / !sahe
  data sbhe       / 0.11689, 0.01560, -.00713 /                   !sbhe
  data sehe       / -.05144, 0.10844, -.04380, -.00663, -.00657 / !sehe
  data sfhe       / 0.12937, 0.02888, -.00413 /                   !sfhe
  data sihe       / 0.12429, -.04498, -.06418, -.00931, -.00582 / !sihe
  data sjhe       / 0.11189, 0.00562, -.01176 /                   !sjhe
  data sphe       / 0.09218, -.00644, -.06552, -.01190, -.00832 / !sphe
  data sqhe       / 0.12575, 0.01598, -.00950 /                   !sqhe
  data alphahe    / 1.5611, -20.6104, 247.0830, -1003.2600 /      !alphahe
  data betahe     / 1.5874, -21.4772, 255.7040, -1028.9000 /      !betahe
  
end module cndshe_data

!**********************************************************************
!  variables used in cndlc
!**********************************************************************

module cndlc_data

  real*8, dimension(4) :: a, d
  real*8, dimension(3) :: b, c, e, f
  
  data a/  0.9896, -0.1851,  0.1019, -0.0360/  !alc
  data b/ -0.8825,  0.6657, -0.3798/           !blc
  data c/ -0.0915, -1.5848,  1.1882/           !clc
  data d/  0.4406, -0.0161, -0.0093, -0.0028/  !dlc
  data e/ -0.4821,  0.0826, -0.0557/           !elc
  data f/ -0.5193, -0.0830,  0.0147/           !flc

end module cndlc_data

!**********************************************************************
!  variables used in cndsc
!**********************************************************************

module cndsc_data

  real*8, dimension(5) :: a, e, i, p
  real*8, dimension(4) :: alpha, beta
  real*8, dimension(3) :: b, f, j, q
  real*8, parameter :: c   = -0.09378, d = 0.94180
  real*8, parameter :: ggg = -0.19961, h = 1.63031
  real*8, parameter :: k   = -0.10127, l = 2.10721
  real*8, parameter :: r   = -0.20326, s = 2.80479

  data a /   0.03387,  0.01795, -0.04324, -0.00437, -0.00422/ !sac
  data b /   0.04865, -0.00473, -0.00865/                     !sbc
  data e /  -0.00329,  0.08497, -0.06909, -0.00582, -0.00677/ !sec
  data f /   0.12695,  0.00854, -0.01197/                     !sfc
  data i /   0.177955,-0.08588, -0.07995, -0.00957, -0.00256/ !sic
  data j /   0.03417, -0.03395, -0.02519/                     !sjc
  data p /   0.18606, -0.06547, -0.10440, -0.01172, -0.00448/ !spc
  data q /   0.10304, -0.02601, -0.02922/                     !sqc
  data alpha/0.4786,  16.4194, -139.1920, 184.7580/            !alphac
  data beta/ 0.6146, 11.7296, -89.8319, 32.0399/               !betac

end module cndsc_data

!**********************************************************************
! variables used in cndlo
!**********************************************************************

module cndlo_data

  real*8, dimension(4) :: a, d
  real*8, dimension(3) :: b, c, e, f

  data a/ 1.0779, -0.1838, 0.1059, -.0290/  !alo
  data b/ -.9743,  0.6955, -.3966/          !blo
  data c/ -.1040, -1.7692, 1.3546 /         !clo
  data d/ 0.4486, -0.0160, -.0014, 0.0039 / !dlo
  data e/ -.5193,  0.0822, -.0467 /         !elo
  data f/ -.5403, -0.1022, 0.0416 /         !flo
  
end module cndlo_data

!**********************************************************************
!  variables used in cndso
!**********************************************************************

module cndso_data

  real*8, dimension(5) :: a, e, i, p
  real*8, dimension(4) :: alpha, beta
  real*8, dimension(3) :: b, f, j, q
  real*8, parameter :: c   = -.07178, d = 0.85073
  real*8, parameter :: ggg = -.19226, h = 1.66185
  real*8, parameter :: k   = -.06064, l = 2.14532
  real*8, parameter :: r   = -.17568, s = 2.99352

  data a /     0.02312, 0.02087, -.03684, -.00366, -.00348 /  !sao
  data b /     0.04177, -.00712, -.00865 /                    !sbo
  data e /     -.001185, 0.08552, -.07206, -.00554, -.00673 / !seo
  data f /     0.12586, 0.00371, -.01360 /                    !sfo
  data i /      .160905, -.08292, -.07377, -.00577, 0.00154/  !sio
  data j /     0.03816, -.04206, -.02870 /                    !sjo
  data p /      .19600, -.07598, -.10917, -.00942, -.00142 /  !spo
  data q /     0.10403, -.03774, -.03483 /                    !sqo
  data alpha/    .4682, 18.0197, -171.0090, 313.8890 /         !alphao
  data beta/    0.6085, 12.9705, -115.3180, 136.3590 /         !betao

end module cndso_data

!**********************************************************************
!  variables used in hlcnc
!**********************************************************************

module hlcnc_data

  real*8, dimension(25) :: a

  data a/-5.18321e0,2.88323e-2,-1.68700e0,-3.80388e-1,-1.42808e0, &
       1.42606e0,-6.77026e-2,2.99633e-1,3.02198e-1,-1.94700e-1,   &
       -2.72927e-3,4.24382e-2,-9.77295e-2,-1.50302e-2,8.16693e-3, &
       4.80447e-4,1.04512e-3,-8.95778e-3,1.27411e-2,-3.37754e-4,  &
       4.17955e-5,-7.52181e-5,-1.38498e-4,5.92069e-4,-5.46719e-4/

end module hlcnc_data

!**********************************************************************
!  variables used in hlcnhe
!**********************************************************************

module hlcnhe_data

  real*8, dimension(5,5) :: a

  data a /9.246777e+0,    -2.097036,    -1.267003,  5.490678e-2,  2.463979e-3,&
         -1.265387e+1,  1.716563e-1,  1.088986e+0, -4.576479e-2, -6.966255e-3,&
          4.539412e+0, -1.427574e-1, -3.499339e-1,  1.589662e-2,  3.421386e-3,&
         -5.842722e-1,  3.930197e-2,  4.702301e-2, -2.529960e-3, -5.814372e-4,&
          2.617848e-2, -2.791279e-3, -2.238884e-3,  1.442103e-4,  3.234929e-5/

end module hlcnhe_data

!**********************************************************************
!  variables used in hlcnh
!**********************************************************************

module hlcnh_data

  real*8, dimension(5,5) :: a

  data a/   12.40928,    -2.545420,    -3.019045,  6.746942e-2,  4.686838e-2,&
           -15.29732,     .8154491,     2.387849, -6.939705e-2, -4.129124e-2,&
            5.162517, -3.839103e-1, -6.899530e-1,  2.568522e-2,  1.281765e-2,&
           -.6347533,  7.203151e-2,  8.405360e-2, -3.965926e-3, -1.652045e-3,&
         2.699066e-2, -4.243979e-3, -3.648879e-3,  2.127303e-4,  7.509063e-5/

end module hlcnh_data

!**********************************************************************
! variables used in bremsolh
!**********************************************************************

!!$      data snah       / -.01148, 0.01601, -.00433, 0.00015, -.00034 / 
!!$      data snbh       / 0.01558, 0.00191, -.00055 /
!!$      data snch, sndh / -.01694, 0.10649 /
!!$      data sneh       / -.01827, 0.02395, -.00448, -.00033, -.00088 / 
!!$      data snfh       / 0.01730, 0.00402, -.00005 /
!!$      data sngh, snhh / -.02222, 0.13969 /
!!$      data snih       / -.00323, 0.00440, -.00110, 0.00001, -.00007 / 
!!$      data snjh       / 0.00294, 0.00059, -.00018 /
!!$      data snkh, snlh / -.00337, 0.02116 /
!!$      data snph       / -.00469, 0.00610, -.00114, -.00010, -.00018 / 
!!$      data snqh       / 0.00320, 0.00107, -.00008 /
!!$      data snrh, snsh / -.00442, 0.02775 /
!!$      data alphasnh   / 1.6449, -23.2588, 272.1670, -1074.7000 /
!!$      data betasnh    / 1.6443, -23.2414, 272.0080, -1074.2500 /
!!$      data pnah       / -.00687, 0.00957, -.00204, -.00005, -.00003 / 
!!$      data pnbh       / 0.00661, 0.00135, -.00035 /
!!$      data pnch, pndh / -.00811, 0.05098 /
!!$      data pnih       / -.00169, 0.00231, -.00047, -.00003, 0.00000 / 
!!$      data pnjh       / 0.00111, 0.00042, -.00010 /
!!$      data pnkh, pnlh / -.00161, 0.01013 /
!!$      data alphapnh   / -.1394, 7.0680, -115.5940, 619.9170 /
!!$      data betapnh    / -.1394, 7.0664, -115.5800, 619.8790 /

!**********************************************************************
! variables used in bremsolc
!**********************************************************************

module bremsolc_data

  real*8, dimension(5) :: snac, snec, snic, snpc, pnac, pnic
  real*8, dimension(4) :: alphasnc, betasnc, alphapnc, betapnc
  real*8, dimension(3) :: snbc, snfc, snjc, snqc, pnbc, pnjc
  real*8, parameter :: sncc = -.01093, sndc = 0.12431
  real*8, parameter :: sngc = -.02259, snhc = 0.20343
  real*8, parameter :: snkc = -.00398, snlc = 0.02499
  real*8, parameter :: snrc = -.00650, snsc = 0.04087
  real*8, parameter :: pncc = -.00729, pndc = 0.06630
  real*8, parameter :: pnkc = -.00212, pnlc = 0.01332

  data snac       / 0.01838, -.01066, -.00458, -.00177, -.00138 / 
  data snbc       / -.00244, -.00206, -.00037 /
  data snec       / 0.02360, -.01353, -.00619, -.00211, -.00176 / 
  data snfc       / 0.00456, -.00174, -.00031 /
  data snic       / 0.00053, -.00048, -.00022, 0.00019, -.00001 / 
  data snjc       / 0.00658, -.00180, 0.00036 /
  data snpc       / -.00024, 0.00063, -.00064, 0.00030, -.00006 / 
  data snqc       / 0.01013, -.00247, 0.00052 /
  data alphasnc   / 0.6252, 10.6819, -70.6879, -44.3349 /
  data betasnc    / 0.6307, 10.4966, -68.7973, -50.0581 /
  data pnac       / 0.01116, -.00589, -.00279, -.00073, -.00043 / 
  data pnbc       / -.00095, -.00059, 0.00002 /
  data pnic       / 0.00012, 0.00018, -.00028, 0.00012, -.00004 / 
  data pnjc       / 0.00339, -.00082, 0.00015 /
  data alphapnc   / 0.5481, -20.4731, 223.9220, -534.9400 /
  data betapnc    / 0.5413, -20.2069, 220.7060, -524.1240 /

end module bremsolc_data

!**********************************************************************
! variables used in bremsolo
!**********************************************************************

module bremsolo_data

  real*8, dimension(5) :: snao, sneo, snio, snpo, pnao, pnio
  real*8, dimension(4) :: alphasno, betasno, alphapno, betapno
  real*8, dimension(3) :: snbo, snfo, snjo, snqo, pnbo, pnjo
  real*8, parameter :: snco = -.00791, sndo = 0.13980
  real*8, parameter :: sngo = -.02610, snho = 0.26993
  real*8, parameter :: snko = -.00447, snlo = 0.02811
  real*8, parameter :: snro = -.00861, snso = 0.05414
  real*8, parameter :: pnco = -.00776, pndo = 0.08995
  real*8, parameter :: pnko = -.00287, pnlo = 0.01803
  
  data snao       / 0.01616, -.00874, -.00413, -.00190, -.00139 / 
  data snbo       / -.00344, -.00261, -.00070 /
  data sneo       / 0.02210, -.00883, -.00857, -.00257, -.00214 / 
  data snfo       / 0.00629, -.00210, -.00099 /
  data snio       / 0.00100, -.00112, -.00003, 0.00014, 0.00001 / 
  data snjo       / 0.00745, -.00209, 0.00044 /
  data snpo       / -.00056, 0.00110, -.00074, 0.00024, -.00004 / 
  data snqo       / 0.01286, -.00281, 0.00057 /
  data alphasno   / 0.4889, 16.1962, -138.4860, 185.7060 /
  data betasno    / 0.5111, 15.4195, -130.1540, 159.6050 /
  data pnao       / 0.00800, -.00191, -.00330, -.00075, -.00047 / 
  data pnbo       / 0.00088, -.00098, -.00036 /
  data pnio       / -.00009, 0.00055, -.00038, 0.00011, -.00003 / 
  data pnjo       / 0.00429, -.00088, 0.00014 /
  data alphapno   / 0.3173, -14.4048, 186.9100, -476.8100 /
  data betapno    / 0.3073, -13.7973, 176.9940, -438.7520 /

end module bremsolo_data

!**********************************************************************
! variables used in bremsolf
!**********************************************************************

!!$      data snaf       / 0.02096, -.01768, -.00007, -.00241, -.00080 / 
!!$      data snbf       / -.01705, -.00268, -.00141 /
!!$      data sncf, sndf / 0.00818, 0.13629 /
!!$      data snef       / 0.04628, -.03290, -.00523, -.00539, -.00276 / 
!!$      data snff       / -.02574, -.00630, -.00285 /
!!$      data sngf, snhf / 0.00022, 0.31871 /
!!$      data snif       / 0.00489, -.00653, 0.00171, -.00024, 0.00017 / 
!!$      data snjf       / 0.00869, -.00323, 0.00075 /
!!$      data snkf, snlf / -.00439, 0.02766 /
!!$      data snpf       / 0.00732, -.00957, 0.00222, -.00022, 0.00025 / 
!!$      data snqf       / 0.01867, -.00615, 0.00133 /
!!$      data snrf, snsf / -.01023, 0.06442 /
!!$      data alphasnf   / 0.6798, 12.7527, -140.1800, 268.8290 /
!!$      data betasnf    / 0.7783, 10.2315, -124.2640, 241.3060 /
!!$      data pnaf       / 0.01292, -.00894, -.00097, -.00125, -.00045 / 
!!$      data pnbf       / -.00739, -.00190, -.00062 /
!!$      data pncf, pndf / 0.00167, 0.09950 /
!!$      data pnif       / 0.00240, -.00271, 0.00047, 0.00001, -.00001 / 
!!$      data pnjf       / 0.00604, -.00197, 0.00044 /
!!$      data pnkf, pnlf / -.00320, 0.02012 /
!!$      data alphapnf   / 0.2847, -13.0828, 192.5030, -543.1480 /
!!$      data betapnf    / 0.3221, -13.8640, 201.0700, -573.1160 /

!**********************************************************************
! variables used in bremliqh
!**********************************************************************

module bremliqh_data

  real*8, dimension(6) :: lnah, lneh, lnih, lnph, qnah, qneh
  real*8, dimension(5) :: alphaqnh
  real*8, dimension(4) :: lnbh, lnfh, lnjh, lnqh, alphalnh, betalnh, qnbh, qnfh
  real*8, parameter :: lnch = 0.00671, lndh = 0.28130
  real*8, parameter :: lngh = -.02199, lnhh = 0.17300
  real*8, parameter :: lnkh = -.01021, lnlh = 0.06417
  real*8, parameter :: lnrh = -.00561, lnsh = 0.03522
  real*8, parameter :: qnch = -.02408, qndh = 0.28763
  real*8, parameter :: qngh = -.00764, qnhh = 0.53346
  
  data lnah       / 0.04518, -.03009, -.00564, -.00544, -.00290, -.00224 /
  data lnbh       / -.02148, -.00817, -.00300, -.00170 /
  data lneh       / -.01003, 0.01790, -.00783, -.00021, 0.00024, -.00014 / 
  data lnfh       / 0.00538, -.00175, -.00346, -.00031 /
  data lnih       / 0.00096, -.00301, -.00073, 0.00182, 0.00037, 0.00116 / 
  data lnjh       / 0.01706, -.00753, 0.00066, -.00060 /
  data lnph       / -.00556, 0.00603, -.00149, 0.00047, 0.00040, 0.00028 / 
  data lnqh       / 0.00422, -.00009, -.00066, -.00003 /
  data alphalnh   / -.07980, 0.17057, 1.51980, -.61058 /
  data betalnh    / -.05881, 0.00165, 1.82700, -.76993 /
  data qnah       / -.02260, 0.02634, -.00424, 0.00055, 0.00017, -.00042 / 
  data qnbh       / 0.00550, 0.00105, -.00179, 0.00028 /
  data qneh       / -.00539, 0.00617, -.00205, 0.00089, 0.00022, 0.00032 / 
  data qnfh       / 0.00615, -.00149, -.00055, -.00004 /
  data alphaqnh   / 2.5344, -7.7188, -9.5605, -5.6433, 1.2671 /

end module bremliqh_data

!**********************************************************************
! variables used in bremliqc
!**********************************************************************

module bremliqc_data

  real*8, dimension(6) :: lnac, lnec, lnic, lnpc, qnac, qnec
  real*8, dimension(5) :: alphaqnc
  real*8, dimension(4) :: lnbc, lnfc, lnjc, lnqc, alphalnc, betalnc, qnbc, qnfc
  real*8, parameter :: lncc = 0.00945, lndc = 0.34529
  real*8, parameter :: lngc = -.02342, lnhc = 0.24819
  real*8, parameter :: lnkc = -.01259, lnlc = 0.07917
  real*8, parameter :: lnrc = -.00829, lnsc = 0.05211
  real*8, parameter :: qncc = -.02458, qndc = 0.25843
  real*8, parameter :: qngc = -.01194, qnhc = 0.34437

  data lnac       / 0.08973, -.05821, -.01089, -.01147, -.00656, -.00519 / 
  data lnbc       / -.04969, -.01584, -.00504, -.00281 /
  data lnec       / 0.03390, -.00944, -.01289, -.00589, -.00404, -.00330 / 
  data lnfc       / -.02213, -.01136, -.00467, -.00131 /
  data lnic       / 0.00383, -.00710, -.00028, 0.00232, 0.00044, 0.00158 / 
  data lnjc       / 0.02300, -.01078, 0.00118, -.00089 /
  data lnpc       / -.00384, 0.00356, -.00184, 0.00146, 0.00031, 0.00069 / 
  data lnqc       / 0.01052, -.00354, -.00014, -.00018 /
  data alphalnc   / -.05483, -.01946, 1.86310, -.78873 /
  data betalnc    / -.06711, 0.06859, 1.74360, -.74498 /
  data qnac       / -.01900, 0.02446, -.00419, -.00030, -.00047, -.00099 / 
  data qnbc       / 0.00093, 0.00060, -.00176, 0.00018 /
  data qnec       / -.00690, 0.01046, -.00351, 0.00010, -.00003, -.00026 / 
  data qnfc       / 0.00289, -.00082, -.00155, -.00008 /
  data alphaqnc   / 0.9859, 2.2468, -10.6808, 11.4582, -4.0100 /

end module bremliqc_data

!**********************************************************************
! variables used in bremliqo
!**********************************************************************

module bremliqo_data

  real*8, dimension(6) :: lnao, lneo, lnio, lnpo, qnao, qneo
  real*8, dimension(5) :: alphaqno
  real*8, dimension(4) :: lnbo, lnfo, lnjo, lnqo, alphalno, betalno, qnbo, qnfo
  real*8, parameter :: lnco = 0.00952, lndo = 0.36029 
  real*8, parameter :: lngo = -.02513, lnho = 0.27480
  real*8, parameter :: lnko = -.01314, lnlo = 0.08263
  real*8, parameter :: lnro = -.00921, lnso = 0.05786
  real*8, parameter :: qnco = -.02491, qndo = 0.25207
  real*8, parameter :: qngo = -.01541, qnho = 0.32800

  data lnao       / 0.10466, -.06740, -.01293, -.01352, -.00776, -.00613 / 
  data lnbo       / -.05950, -.01837, -.00567, -.00310 /
  data lneo       / 0.04652, -.01656, -.01489, -.00778, -.00520, -.00418 / 
  data lnfo       / -.03076, -.01390, -.00522, -.00161 /
  data lnio       / 0.00476, -.00838, -.00011, 0.00244, 0.00046, 0.00168 / 
  data lnjo       / 0.02455, -.01167, 0.00132, -.00097 /
  data lnpo       / -.00350, 0.00295, -.00184, 0.00166, 0.00032, 0.00082 / 
  data lnqo       / 0.01231, -.00445, 0.00002, -.00026 /
  data alphalno   / -.06597, 0.06048, 1.74860, -.74308 /
  data betalno    / -.07356, 0.10865, 1.70150, -.73653 /
  data qnao       / -.01792, 0.02395, -.00425, -.00055, -.00065, -.00115 / 
  data qnbo       / -.00039, 0.00043, -.00177, 0.00015 /
  data qneo       / -.00866, 0.01344, -.00448, -.00002, -.00009, -.00041 /
  data qnfo       / 0.00291, -.00091, -.00200, -.00010 /
  data alphaqno   / 0.8978, 2.5443, -10.3549, 10.2665, -3.3537 /

end module bremliqo_data

!**********************************************************************
! variables used in bremliqf
!**********************************************************************

!!$      data lnaf       / 0.17444, -.11076, -.02349, -.02283, -.01250,
!!$     1                  -.00971 / 
!!$      data lnbf       / -.10661, -.02860, -.00785, -.00385 /
!!$      data lncf, lndf / 0.00766, 0.40991 /
!!$      data lnef       / 0.11580, -.05891, -.02531, -.01747, -.01021,
!!$     1                  -.00778 / 
!!$      data lnff       / -.07767, -.02516, -.00711, -.00260 /
!!$      data lngf, lnhf / -.03076, 0.36908 /
!!$      data lnif       / 0.00946, -.01493, -.00125, 0.00262, 0.00055,
!!$     1                  0.00209 / 
!!$      data lnjf       / 0.03034, -.01519, 0.00204, -.00135 /
!!$      data lnkf, lnlf / -.01494, 0.09395 /
!!$      data lnpf       / -.00006, -.00222, -.00104, 0.00225, 0.00037,
!!$     1                  0.00140 / 
!!$      data lnqf       / 0.02013, -.00883, 0.00090, -.00071 /
!!$      data lnrf, lnsf / -.01249, 0.07850 /
!!$      data alphalnf   / -.07608, 0.11559, 1.75730, -.79677 /
!!$      data betalnf    / -.08034, 0.14368, 1.73140, -.79467 /
!!$      data qnaf       / -.01326, 0.02179, -.00484, -.00146, -.00135, 
!!$     1                  -.00175 / 
!!$      data qnbf       / -.00567, -.00040, -.00174, 0.00012 /
!!$      data qncf, qndf / -.02710, 0.23276 /
!!$      data qnef       / -.00430, 0.01443, -.00677, -.00142, -.00119,
!!$     1                  -.00151 /
!!$      data qnff       / -.00396, -.00245, -.00278, -.00023 /
!!$      data qngf, qnhf / -.02288, 0.26787 /
!!$      data alphaqnf   / 1.0379, 1.0195, -5.1242, 3.9369, -0.8701 /

      
