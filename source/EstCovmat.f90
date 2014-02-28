!Estimate the covariance matrix by making steps on a grid in nd
! slb 5 aug 2004
! AML 2012, updates to have chance of working with bounded and unconstrained parameters

! Overwrites whatever is in propose_matrix with the newly estimated matrix
! (In future could add/replace some entries with any better info on the covmat
! which is already known.)
! Uses the steps from the params.ini file to choose grid

module EstCovmatModule
 use Random
 use CalcLike
 use ParamPointSet
 use settings
 use Matrixutils
 use BaseParameters
 implicit none

 integer, allocatable :: step_positions(:,:) !step, parameter
 !step is -1,0,1 for normal case, but one-sided e.g. 0,1,2 in case of hard prior boundary nearby
 logical, allocatable  :: is_weakly_constrained(:)
contains


 function EstCovmat(Params, bestdll_in, sucess)
   integer :: ntries = 10 ! Maximum number of attempts to get a covmat with all +ve evals
   integer itries ! Number of tries so far
   integer, intent(out) :: sucess
   real(mcp) parsteps(num_params_used), bestdll
   real(mcp), intent(in) :: bestdll_in  ! the dloglike value used to determine grid spacing eg. try 4
   Type(ParamSet) Params, CenterParams
   !double precision Gaussian1
   real(mcp) Hess(num_params_used,num_params_used), U(num_params_used,num_params_used)
   real(mcp) D(num_params_used),EstCovmat(num_params_used,num_params_used)
   real(mcp) CenterLike

   itries=0
   sucess=0
   bestdll = bestdll_in


   if (Feedback>1) write(*,*)
   if (Feedback>1) write(*,*) 'Estimating the covariance matrix by finding likelihood values on a grid'

   ! Set up the grid position
   CenterParams = Params
   allocate(step_positions(3,num_params_used))
   allocate(is_weakly_constrained(num_params_used))
   parsteps=BaseParams%PWidth(params_used)
   CenterLike = PlaceGrid(CenterParams,parsteps,bestdll)

   ! Loop over attempts to find the covariance matrix
   do while (sucess.eq.0)
      itries = itries + 1
      call GetHess(CenterParams,CenterLike,parsteps,Hess)

      U=Hess
      call Matrix_Diagonalize(U, D, num_params_used)
      if (any(D <= 0)) then
        write(*,*)
        write(*,*) '!!! -ve evals, so spreading the grid out by a factor of ~2'
        bestdll=bestdll*2 ! try using a more spread out grid
        CenterLike= PlaceGrid(CenterParams,parsteps,bestdll)
      else
             sucess=1
      end if

      ! just give up if we have overrun the maximum no. of tries specified above
      if (itries > ntries) then
             write(*,*) 'EstCovmat.f90: Couldn''t find covmat with all positive evals after ', &
                         ntries,' tries - check starting ranges'
         stop
      end if

   end do

   ! uncomment this bit to update the propose_matrix..
    if (sucess>0) then
      U=Hess
      call Matrix_Inverse(U) ! gives the covariance matrix
      EstCovmat(1:num_params_used,1:num_params_used)=U
    end if

   call CenterParams%Clear(keep=Params)

   deallocate(step_positions,is_weakly_constrained)

 end function EstCovmat


 function PlaceGrid(CenterParams,parsteps,bestdll) result(CenterLike)
   ! Takes the input CenterParams and parsteps values and improves on them
   ! While doing so, fills the likelhood grid
   ! CenterParams and parsteps values
   Type(ParamSet) P, CenterParams
   real(mcp) parsteps(num_params_used),bestdll, CenterLike
   integer maxstepsizetries

   if (Feedback>1) write(*,*) 'Improving the grid spacings'
   P = CenterParams
   CenterLike = GetLogLike(CenterParams)
   call P%Clear(keep=CenterParams)
   if (CenterLike == LogZero) then
         write(*,*) 'ERROR: Trial parameters excluded by prior or bad likelihood'
         write(*,*) 'Try starting further away from problem regions?'
         call DoAbort
   end if

   ! Decide how big to make the grid
   maxstepsizetries=20
   call GetStepsForDChisq1(CenterParams,CenterLike,maxstepsizetries,parsteps,bestdll)
   ! See if this is a good place to center it
   !if (.any.istepsizetries>maxstepsizetries) then it failed...

 end function PlaceGrid


 subroutine GetStepsForDChisq1(CenterParams,LCenter, maxstepsizetries,parsteps,bestdll)
   ! Stays centered on CenterParams
   ! But for each paramter, iterates for a maximum of maxstepsizetries to try to find
   ! parsteps such that Delta chisq is of order 1 (when stepping only in that parameter)
   ! Inputs: CenterParams, LCenter, maxstepsizetries
   ! Outputs: parsteps
   Type(ParamSet) CenterParams
   integer, intent(in) :: maxstepsizetries
   real(mcp) parsteps(num_params_used), bestdll
   integer i, ii, tries
   Type(ParamSet) StepParams
   real(mcp) LCenter, LStep,dloglike
   real(mcp) step_min, step_max
   real(mcp) step_limit
   integer asign

   ! The center of the whole grid
   is_weakly_constrained = .false.
   !! Estimate some good step sizes to use by interating until dchisq~1
   do i=1, num_params_used
      ii= params_used(i)

      step_max = 0.
      step_min = 0.

      if (CenterParams%P(ii)- BaseParams%PMin(ii) < BaseParams%PMax(ii) - CenterParams%P(ii)) then
          asign=1
          step_limit = (BaseParams%PMax(ii) - CenterParams%P(ii))/2
      else
          asign = -1
          step_limit = (CenterParams%P(ii)- BaseParams%PMin(ii))/2
      end if
      parsteps(i) = min(parsteps(i), step_limit * 0.8)

      ! Find step sizes that give a half sensible Delta chisq eg. 0.5 to 10
      do tries = 1, maxstepsizetries

         if (parsteps(i) > step_limit) then
           is_weakly_constrained(i) = .true.
           parsteps(i) = step_limit
           if (Feedback >1 ) write(*,*) &
            ' Parameter '//trim(BaseParams%UsedParamNameOrNumber(i))//' is weakly constrained, neglect correlations'
           exit
         end if
         StepParams=CenterParams ! reset back to original starting point
         StepParams%P(ii)=CenterParams%P(ii) + parsteps(i)*asign
         LStep=GetLogLike(StepParams)
         call StepParams%Clear(keep=CenterParams)

         dloglike=abs(LStep-LCenter)
         if (dloglike > 3*bestdll) then
           step_max = parsteps(i)
           parsteps(i) = (parsteps(i) + step_min) / 2
         elseif (dloglike < bestdll/3) then
           step_min = parsteps(i)
           if (step_max==0.) then
            parsteps(i) = parsteps(i) * 2
           else
            parsteps(i) = (parsteps(i) + step_max) / 2
           end if
         else
          exit
         end if

         ! just give up if we have overrun the maximum no. of tries specified above
         if (tries == maxstepsizetries) then
            write(*,*) 'Couldn''t find a good stepsize for parameter ',ii,' with Delta chisq ~1 but continuing anyway'
         end if
      end do
      if (Feedback>1) write(*,*) ' Decided on stepsize ',parsteps(i),' for parameter '//trim(BaseParams%UsedParamNameOrNumber(i))

      if (CenterParams%P(ii) - parsteps(i) < BaseParams%PMin(ii)) then
        !Is hard prior bound below, one sided steps above only
         step_positions(:,i) = (/0,1,2/)
      else if (CenterParams%P(ii) + parsteps(i) > BaseParams%PMax(ii)) then
        !Is hard prior bound below, one sided steps below only
         step_positions(:,i) = (/-2,-1,0/)
      else
         step_positions(:,i)  = (/-1,0,1/)
      end if
   end do

 end subroutine GetStepsForDChisq1



 subroutine GetHess(CenterParams,CenterLike,parsteps,Hess)
   ! Get the curvature matrix (Hess) for the grid defined by CenterParams
   ! and parsteps.
   ! Return the answer in Hess.
   ! Input only: CenterParams, parsteps
   ! Ouput only: Hess
   Type(ParamSet) CenterParams,StepParams
   real(mcp) parsteps(num_params_used)
   real(mcp) Hess(num_params_used,num_params_used)
   real(mcp) wii,wjj
   real(mcp) LGrid(3), LGrid2(3,3), CenterLike
   integer i,j,ii,jj,istep,jstep

   ! Find likelihood values on a grid to estimate curvature matrix (Hess)
   if (Feedback>1) write(*,*) ' Finding curvature matrix ...'
   Hess = 0
   do i=1, num_params_used
      ii= params_used(i)
      if (Feedback>1) write(*,*)
      if (Feedback>1) write(*,*) ' Parameter '//trim(BaseParams%UsedParamNameOrNumber(i))

      if (is_weakly_constrained(i)) then
         Hess(i,i) = 1/((BaseParams%PMax(ii) - BaseParams%PMin(ii))/2)
      else
          wii=parsteps(i)
          do istep=1,3
            if (step_positions(istep,i)==0) then
             Lgrid(istep) = CenterLike
             cycle
            end if
            StepParams=CenterParams ! rest back to original starting point
            StepParams%P(ii)=CenterParams%P(ii) + step_positions(istep,i)*wii
            Lgrid(istep) = GetLogLike(StepParams)
            call StepParams%Clear(keep=CenterParams)
            if (Lgrid(istep) == LogZero) then
              write(*,*) 'ERROR: Trial parameters hit prior or error in function evaluation'
              write(*,*) 'Try starting further away from problem regions?'
              call DoAbort
            end if
          end do
          Hess(i,i)=1/(wii*wii) *(Lgrid(1) + Lgrid(3) - 2*LGrid(2))
      end if
       do j=(i+1), num_params_used
         jj=params_used(j)
         if (is_weakly_constrained(j)) cycle !just put zero in off diagonal for unconstrained
         if (Feedback>1) write(*,*) ' Parameter pair '//trim(UsedParamNameOrNumber(i))// &
                                    ' , '//trim(UsedParamNameOrNumber(j))
         wjj=parsteps(j)

         do istep=1,3,2
             do jstep= 1,3,2  ! ie just do jstep=1 and jstep=3
               StepParams=CenterParams ! rest back to original starting point
               StepParams%P(jj)=CenterParams%P(jj) + step_positions(jstep,j)*wjj
               StepParams%P(ii)=CenterParams%P(ii) + step_positions(istep,i)*wii
               LGrid2(istep,jstep) = GetLogLike(StepParams)
               call StepParams%Clear(keep=CenterParams)
               if (LGrid2(istep,jstep) == LogZero) then
                write(*,*) 'ERROR: Trial parameters hit prior or error in function evaluation'
                     write(*,*) 'Try starting further away from problem regions?'
                     write(*,*) '  - parameter pair is '//trim(UsedParamNameOrNumber(i))// &
                                    ' , '//trim(UsedParamNameOrNumber(j))
                call DoAbort
               end if

            end do
         end do

        Hess(i,j)=1/(wii*wjj*4)*(Lgrid2(1,1) + Lgrid2(3,3) - Lgrid2(1,3) - Lgrid2(3,1))
        Hess(j,i)=Hess(i,j)

      end do
   end do

 end subroutine GetHess

end module EstCovmatModule
