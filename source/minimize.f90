
    module minimize
    use ParamDef
    use CalcLike
    use Powell_ConstrainedOptimize
    implicit none
    private 

    Type(ParamSet) MinParams

    public FindBestFit, WriteBestFitParams

    contains

    subroutine VectToParams(vect, P)
    real(Powell_CO_prec) vect(num_params_used)
    Type(ParamSet) P
    integer i,ii

    do i=1, num_params_used
        ii=params_used(i)
        P%P(ii)=vect(i)*Scales%PWidth(ii) + Scales%center(ii)
    end do

    end subroutine VectToparams

    subroutine ParamsToVect(P,vect)
    real(Powell_CO_prec) vect(num_params_used)
    Type(ParamSet) P
    integer i,ii

    do i=1, num_params_used
        ii=params_used(i)
        vect(i)=(P%P(ii)-Scales%center(ii))/Scales%PWidth(ii)
    end do

    end subroutine ParamsToVect


    function ffn(n,vect) result(like)
    implicit none
    integer, intent(in) :: n
    real(Powell_CO_prec) :: like
    real(Powell_CO_prec) vect(n)
    Type(ParamSet) P

    P = MinParams
    call VectToParams(vect, P)

    like = GetLogLike(P)

    call AcceptReject(.true.,MinParams%Info, P%Info)
    MinParams = P !want to keep e.g. the Age calculation

    end function ffn

    function FindBestFit(Params,sigma_frac_err, max_iterations) result(best_like)
    use ParamDef
    Type(ParamSet) Params
    integer, intent(in) :: max_iterations
    real, intent(in) :: sigma_frac_err
    real best_like
    real(Powell_CO_prec) :: vect(num_params_used), XL(num_params_used), XU(num_params_used)
    real(Powell_CO_prec) rhobeg, rhoend
    integer npt

    MinParams = Params
    ! scale the params so they are all roughly the same order of magnitude
    !call ParamsToVect(MinParams, vect)
    vect = 0 !start at zero by definition (from input center values)
    !normalized parameter bounds
    XL = (Scales%PMin(params_used)-Scales%center(params_used)) / Scales%PWidth(params_used)
    XU = (Scales%PMax(params_used)-Scales%center(params_used)) / Scales%PWidth(params_used)  

    !Initial and final radius of region required (in normalized units)
    rhobeg = 0.1*sqrt(real(num_params_used))
    rhoend = sigma_frac_err*sqrt(real(num_params_used))
    npt = 2*num_params_used +1 !Hessian approximation is singular in general
    if (.not. BOBYQA (ffn, num_params_used ,npt, vect,XL,XU,rhobeg,rhoend,FeedBack+1, max_iterations)) &
       stop 'minimize: FindBestFit failed'

    best_like = ffn(num_params_used,vect)

    !back to real units
    call AcceptReject(.true.,Params%Info, MinParams%Info)
    Params = MinParams

    end function FindBestFit


    subroutine WriteParamsHumanText(aunit, P, like)
    use settings
    use cmbtypes
    use ParamDef
    implicit none
    Type(ParamSet) P
    real, intent(in), optional :: like
    integer, intent(in) :: aunit
    Type(real_pointer) :: derived
    integer numderived 
    integer CalcDerivedParams
    external CalcDerivedParams
    integer isused,i

    if (present(like)) then
        write (aunit,*) '-log(Like) = ',like
        write (aunit,*) ''
    end if

    do isused = 0,1
        do i=1, num_real_params
            if (isused==0 .and. Scales%PWidth(i)/=0 .or. isused==1 .and. Scales%PWidth(i)==0) then
                write(aunit,'(1I5,1E15.7,"   '//trim(ParamNames_name(NameMapping,i))//'")') i, P%P(i)
            end if 
        end do
        write (aunit,*) ''
    end do

    if (generic_mcmc) return

    numderived = CalcDerivedParams(P, derived)
    do i=1, numderived
        write(aunit,'(1I5,1E15.7,"   '//trim(ParamNames_name(NameMapping,num_real_params + i ))//'")') &
        num_real_params+i, derived%P(i)
    end do
    deallocate(derived%P)

    if (nuisance_params_used>0) then
        write (aunit,*) ''
        do i=1, nuisance_params_used
            write(aunit,'(1I5,1E15.7)') num_real_params+i, P%P(num_real_params+i)
        end do
    end if

    end  subroutine WriteParamsHumanText

#ifdef f2003
use, intrinsic :: iso_fortran_env, only : input_unit=>stdin, &
                                          output_unit=>stdout, &
                                          error_unit=>stderr
#else
#define stdin  5
#define stdout 6
#define stderr 0
#endif

    subroutine WriteBestFitParams(like, Params, fname)
    real like
    Type(ParamSet) Params
    character(LEN=*), intent(in) :: fname

     
    call CreateTxtFile(fname,tmp_file_unit)
    call WriteParamsHumanText(tmp_file_unit,Params, like)
    close(tmp_file_unit)
    if (Feedback>0) call WriteParamsHumanText(stdout,Params, like)
    
   
    end subroutine WriteBestFitParams
    
    !function BestFitCovmatEstimate(tol) result (M)
    !!Cannot find any way to get usable Hessian matrix from BOBYQA, even with high NPT
    !!This is not used..
    ! use MatrixUtils
    ! real M(num_params_used,num_params_used), diag(num_params_used)
    ! real tol !tol is the smallest normalized eigenvalue to allow (e..g fraction of input width)
    ! integer i
    ! 
    !  M = BOBYQA_Hessian
    !  !this may not be invertible, let's make sure all eigenvalues are positive
    !  
    !  call  Matrix_Diagonalize(M, diag, num_params_used)
    !       !Does m = U diag U^T, returning U in M
    !  if (Feedback > 0) print *,  'Un-regularized Hessian evalues: ', diag
    !  do i=1,num_params_used
    !    diag(i) = max(diag(i), tol)
    !    M(:,i) = M(:,i) / sqrt(diag(i))
    !  end do
    !  M = matmul(M, transpose(M))
    !
    !  !now go back to un-normalized units
    !  do i=1, num_params_used
    !       M(:,i) = M(:,i) * Scales%PWidth(params_used(i))
    !       M(i,:) = M(i,:) * Scales%PWidth(params_used(i))
    !  end do
    !
    !end function BestFitCovmatEstimate


    end module minimize