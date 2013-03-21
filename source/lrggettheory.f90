!!! moved this routine from modules.f90 in CAMB here by preference of Antony.
!!! needs to be separate from Pktheory.f90 because that file depends on CMB_Cls somehow

module lrggettheory
use precision
use Transfer
implicit none

real(dl), parameter ::  aNEAR = 0.809717d0, aMID = 0.745156d0, aFAR = 0.70373d0
real(dl), parameter ::  z0 = 0.0d0, zNEAR = 0.235d0, zMID = 0.342d0, zFAR = 0.421d0
real(dl), parameter ::  sigma2BAONEAR = 86.9988, sigma2BAOMID = 85.1374, sigma2BAOFAR = 84.5958
real(dl), parameter :: zeffDR7 = 0.312782  !! effective redshift of the LRG sample
real(dl), dimension(4) :: transferscalefid  !! this is set in LRGinfo_init
real(dl), dimension(4) :: powerscaletoz0
 !! this is to scale the amplitude of the redshift slices power spectra to the z=0 amplitude;
 !this is the assumption of the model.
real(dl), parameter :: kmindata = 0.02
 !! in h/Mpc.  they are needed for normalizing nowiggs power spectrum.
! Hard coded for the SDSS DR7 values.
integer :: iz0lrg, izNEARlrg, izMIDlrg, izFARlrg
logical :: use_dr7lrg = .false.

contains

        subroutine Transfer_GetMatterPowerAndNW(MTrans,outpower, itf, in, minkh, dlnkh, &
                npoints, kmindata, getabstransferscale, outpowernw, outpowerrationwhalofit)

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
          real, intent(out) :: outpower(npoints)
          real, intent(out) :: outpowernw(npoints), outpowerrationwhalofit(npoints)
          real, intent(in) :: minkh, dlnkh
          real(dl), intent(in) :: kmindata
          real(dl), intent(out) :: getabstransferscale
          real(dl), parameter :: cllo=1.e30_dl,clhi=1.e30_dl
          integer ik, llo,il,lhi,lastix
          real(dl) matpower(MTrans%num_q_trans), kh, kvals(MTrans%num_q_trans), ddmat(MTrans%num_q_trans)
          real(dl) atransfer,xi, a0, b0, ho, logmink,k, h, fbaryon,omegam
          real(dl) matpowernw(MTrans%num_q_trans), matpowernwhalofit(MTrans%num_q_trans), &
          & atransfernw, atransfernwhalofit, &
          &ddmatnw(MTrans%num_q_trans), ddmatnwhalofit(MTrans%num_q_trans)

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

#ifdef DR71RG
          call dopksmoothbspline(mykvals,mylnpklinear,mylnpksmooth, MTrans%num_q_trans)
#else
         call MpiStop('mpk: edit makefile to have "EXTDATA = LRG" to inlude LRGs')
#endif
          setabs = 0
          do nwi = 1, MTransnw%num_q_trans
            kh = MTrans%TransferData(Transfer_kh,nwi,1)
            if(kh > kmindata .and. setabs == 0) then
              getabstransferscale = sqrt(exp(mylnpklinear(nwi)))
              setabs = 1
           end if
            k = kh*h
            MTransnw%TransferData(transfer_power_var,nwi,1) = sqrt(exp(mylnpksmooth(nwi))/(k*pi*twopi*h**3*ScalarPower(k,in)))
          end do
          if (npoints < 2) stop 'Need at least 2 points in Transfer_GetMatterPower'
          if (minkh*exp((npoints-1)*dlnkh) > MTrans%TransferData(Transfer_kh,MTrans%num_q_trans,itf) &
                .and. FeedbackLevel > 0 ) &
                    write(*,*) 'Warning: extrapolating matter power in Transfer_GetMatterPower'

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
             atransfernw=MTransnw%TransferData(transfer_power_var,ik,1)
             atransfernwhalofit=MTransnw%TransferData(transfer_power_var,ik,1)
             atransfernwhalofit = atransfernwhalofit * PKnw%nonlin_ratio(ik,1)
             matpower(ik) = log(atransfer**2*k*pi*twopi*h**3)
                 !Put in power spectrum later: transfer functions should be smooth, initial power may not be

             matpowernw(ik) = log(atransfernw**2*k*pi*twopi*h**3)
             matpowernwhalofit(ik) = log(atransfernwhalofit**2*k*pi*twopi*h**3)
          end do
          call spline(kvals,matpower,MTrans%num_q_trans,cllo,clhi,ddmat)
          call spline(kvals,matpowernw,MTrans%num_q_trans,cllo,clhi,ddmatnw)
          call spline(kvals,matpowernwhalofit,MTrans%num_q_trans,cllo,clhi,ddmatnwhalofit)


            llo=1
            lastix = npoints + 1
            do il=1, npoints
               xi=logmink + dlnkh*(il-1)
               if (xi < kvals(1)) then
                 outpower(il)=-30.
                 outpowernw(il)=-30.
                 outpowerrationwhalofit(il)=-30.
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
               outpowernw(il) = a0*matpowernw(llo)+ b0*matpowernw(lhi)+((a0**3-a0)* ddmatnw(llo) &
                       +(b0**3-b0)*ddmatnw(lhi))*ho**2/6
               outpowerrationwhalofit(il) = a0*matpowernwhalofit(llo)+ b0*matpowernwhalofit(lhi)+((a0**3-a0)* ddmatnwhalofit(llo) &
                       +(b0**3-b0)*ddmatnwhalofit(lhi))*ho**2/6

            end do

            do while (lastix <= npoints)
               !Do linear extrapolation in the log
               !Obviouly inaccurate, non-linear etc, but OK if only using in tails of window functions
               outpower(lastix) = 2*outpower(lastix-1) - outpower(lastix-2)
               outpowernw(lastix) = 2*outpowernw(lastix-1) - outpowernw(lastix-2)
               outpowerrationwhalofit(lastix) = 2*outpowerrationwhalofit(lastix-1)&
                         - outpowerrationwhalofit(lastix-2)
               lastix = lastix+1
            end do

            outpower = exp(max(-30.,outpower))
            outpowernw = exp(max(-30.,outpowernw))
            outpowerrationwhalofit = exp(max(-30.,outpowerrationwhalofit))


            do il = 1, npoints
               k = exp(logmink + dlnkh*(il-1))*h
               outpower(il) = outpower(il) * ScalarPower(k,in)
               outpowerrationwhalofit(il) = outpowerrationwhalofit(il)/outpowernw(il)
                !! do this first because the ScalarPower calls cancel.
               outpowernw(il) = outpowernw(il) * ScalarPower(k,in)
               !print *,k/h,outpower(il),outpowernw(il)
               !print *,k/h,outpowerrationwhalofit(il)
            end do

          call MatterPowerdata_Free(PKnw)
          deallocate(MTransnw%q_trans)
          deallocate(MTransnw%TransferData)
          deallocate(MTransnw%sigma_8)

        end subroutine Transfer_GetMatterPowerAndNW


end module
