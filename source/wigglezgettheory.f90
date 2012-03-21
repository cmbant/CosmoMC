!!! moved this routine from modules.f90 in CAMB here by preference of Antony.
!!! needs to be separate from Pktheory.f90 because that file depends on CMB_Cls somehow

! adapated version of Reid et al code for computing no-wigglez non-linear power spectra
! note: the outputs are very different

module wigglezgettheory
use precision
use Transfer
implicit none

real(dl), dimension(4) :: transferscalefid  
integer, dimension(4)  :: izwigglez
 !! in h/Mpc.  they are needed for normalizing nowiggs power spectrum. 
logical :: use_wigz10 = .false. 
contains

        subroutine Transfer_GetMatterPowerNWandNL(MTrans,outpower, itf, in, minkh, dlnkh, &
                npoints, outpowernwlin, outpowernwhf )

          !Allows for non-smooth priordial spectra
          !if CP%Nonlinear/ = NonLinear_none includes non-linear evolution
          !Get total matter power spectrum at logarithmically equal intervals dlnkh of k/h starting at minkh
          !in units of (h Mpc^{-1})^3.
          !Here there definition is < Delta^2(x) > = 1/(2 pi)^3 int d^3k P_k(k)
          !We are assuming that Cls are generated so any baryonic wiggles are well sampled and that matter power
          !sepctrum is generated to beyond the CMB k_max
          Type(MatterTransferData), intent(in) :: MTrans
          Type(MatterPowerData) :: PKnw

          integer, intent(in) :: itf, in, npoints 
          real, intent(out) :: outpower(npoints) ! the linear power spectra
          real, intent(out) :: outpowernwlin(npoints) ! the no-wigglez linear Pk
          real, intent(out) :: outpowernwhf(npoints) ! the no-wiiglez halofit
          real, intent(in) :: minkh, dlnkh
          real(dl), parameter :: cllo=1.e30_dl,clhi=1.e30_dl
          integer ik, llo,il,lhi,lastix
          real(dl) matpower(MTrans%num_q_trans), kh, kvals(MTrans%num_q_trans), ddmat(MTrans%num_q_trans)
          real(dl) atransfer,xi, a0, b0, ho, logmink,k, h, fbaryon,omegam
          real(dl) matpowernwlin(MTrans%num_q_trans), matpowernwhf(MTrans%num_q_trans), &
          & atransfernwlin, atransfernwhf, &
          &ddmatnwlin(MTrans%num_q_trans), ddmatnwhf(MTrans%num_q_trans)

          !!added for splining.
          real(dl) :: mykvals(MTrans%num_q_trans),mylnpklinear(MTrans%num_q_trans),mylnpksmooth(MTrans%num_q_trans)

          integer :: nwi,tempi, setabs
          Type(MatterTransferData) :: MTransnw

          MTransnw%num_q_trans = MTrans%num_q_trans
          allocate(MTransnw%q_trans(MTransnw%num_q_trans))
          allocate(MTransnw%TransferData(Transfer_max,MTransnw%num_q_trans,CP%Transfer%num_redshifts))
          allocate(MTransnw%sigma_8(CP%Transfer%num_redshifts, CP%InitPower%nn))

          h = CP%H0/100
          do nwi = 1, MTransnw%num_q_trans
            MTransnw%q_trans(nwi) = MTrans%q_trans(nwi)  !! not ever referenced.
            MTransnw%TransferData(Transfer_kh,nwi,1) = MTrans%TransferData(Transfer_kh,nwi,1)
            kh = MTrans%TransferData(Transfer_kh,nwi,1)
            k = kh*h
            do tempi=2,Transfer_tot
              MTransnw%TransferData(tempi,nwi,1) = 0.0d0
            end do
            mykvals(nwi) = k
            atransfer=MTrans%TransferData(transfer_power_var,nwi,itf)
            mylnpklinear(nwi) = log(atransfer**2*k*pi*twopi*h**3*ScalarPower(k,in))
          end do


          ! produce linear no-wigglez power spectrum
#ifdef WIGZ
          call dopksmoothbspline(mykvals,mylnpklinear,mylnpksmooth, MTrans%num_q_trans)
#else
         call MpiStop('mpk: edit makefile to have "EXTDATA = WIGZ" to inlude LRGs')
#endif

          do nwi = 1, MTransnw%num_q_trans
            kh = MTrans%TransferData(Transfer_kh,nwi,1)
            k = kh*h
            MTransnw%TransferData(transfer_power_var,nwi,1) = sqrt(exp(mylnpksmooth(nwi))/(k*pi*twopi*h**3*ScalarPower(k,in)))
          end do
          if (npoints < 2) stop 'Need at least 2 points in Transfer_GetMatterPower'

          !! get nonlinear on Pnw
          call Transfer_GetMatterPowerData(MTransnw,PKnw, in, 1)
          Pknw%redshifts(1) = CP%Transfer%Redshifts(itf)
          call NonLinear_GetRatios(Pknw)

          h = CP%H0/100
          logmink = log(minkh)
          do ik=1,MTrans%num_q_trans
             kh = MTrans%TransferData(Transfer_kh,ik,itf)
             k = kh*h
             kvals(ik) = log(kh)
             atransfer=MTrans%TransferData(transfer_power_var,ik,itf)
             atransfernwlin=MTransnw%TransferData(transfer_power_var,ik,1)
             atransfernwhf=MTransnw%TransferData(transfer_power_var,ik,1)
             atransfernwhf = atransfernwhf * PKnw%nonlin_ratio(ik,1)
             matpower(ik) = log(atransfer**2*k*pi*twopi*h**3)
                 !Put in power spectrum later: transfer functions should be smooth, initial power may not be

             matpowernwlin(ik) = log(atransfernwlin**2*k*pi*twopi*h**3)
             matpowernwhf(ik) = log(atransfernwhf**2*k*pi*twopi*h**3)
          end do
          call spline(kvals,matpower,MTrans%num_q_trans,cllo,clhi,ddmat)
          call spline(kvals,matpowernwlin,MTrans%num_q_trans,cllo,clhi,ddmatnwlin)
          call spline(kvals,matpowernwhf,MTrans%num_q_trans,cllo,clhi,ddmatnwhf)


            llo=1
            lastix = npoints + 1
            do il=1, npoints
               xi=logmink + dlnkh*(il-1)
               if (xi < kvals(1)) then
                 outpower(il)=-30.
                 outpowernwlin(il)=-30.
                 outpowernwhf(il)=-30.
                 cycle
               end if
               do while ((xi > kvals(llo+1)).and.(llo < MTrans%num_q_trans))
                  llo=llo+1
                  if (llo >= MTrans%num_q_trans) exit
               end do
               if (llo == MTrans%num_q_trans) then
                   lastix = il
                   exit
               end if
               lhi=llo+1
               ho=kvals(lhi)-kvals(llo)
               a0=(kvals(lhi)-xi)/ho
               b0=(xi-kvals(llo))/ho

               outpower(il) = a0*matpower(llo)+ b0*matpower(lhi)+((a0**3-a0)* ddmat(llo) &
                       +(b0**3-b0)*ddmat(lhi))*ho**2/6
               outpowernwlin(il) = a0*matpowernwlin(llo)+ b0*matpowernwlin(lhi)+((a0**3-a0)* ddmatnwlin(llo) &
                       +(b0**3-b0)*ddmatnwlin(lhi))*ho**2/6
               outpowernwhf(il)  = a0*matpowernwhf(llo) + b0*matpowernwhf(lhi) +((a0**3-a0)* ddmatnwhf(llo)  &
                       +(b0**3-b0)*ddmatnwhf(lhi))*ho**2/6

            end do

            do while (lastix <= npoints)
               !Do linear extrapolation in the log
               !Obviouly inaccurate, non-linear etc, but OK if only using in tails of window functions
               outpower(lastix) = 2*outpower(lastix-1) - outpower(lastix-2)
               outpowernwlin(lastix) = 2*outpowernwlin(lastix-1) - outpowernwlin(lastix-2)
               outpowernwhf(lastix) = 2*outpowernwhf(lastix-1) - outpowernwhf(lastix-2)
               lastix = lastix+1
            end do

            outpower = exp(max(-30.,outpower))
            outpowernwlin = exp(max(-30.,outpowernwlin))
            outpowernwhf = exp(max(-30.,outpowernwhf))


            do il = 1, npoints
               k = exp(logmink + dlnkh*(il-1))*h
               outpower(il)      = outpower(il) * ScalarPower(k,in)
               outpowernwhf(il)  = outpowernwhf(il) * ScalarPower(k,in)
               outpowernwlin(il) = outpowernwlin(il) * ScalarPower(k,in)
            end do

          call MatterPowerdata_Free(PKnw)
          deallocate(MTransnw%q_trans)
          deallocate(MTransnw%TransferData)
          deallocate(MTransnw%sigma_8)

        end subroutine Transfer_GetMatterPowerNWandNL


end module
