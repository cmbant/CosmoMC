
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
        P%P(ii)=vect(i)*Scales%PWidth(ii)
    end do

    end subroutine VectToparams

    subroutine ParamsToVect(P,vect)
    real(Powell_CO_prec) vect(num_params_used)
    Type(ParamSet) P
    integer i,ii

    do i=1, num_params_used
        ii=params_used(i)
        vect(i)=P%P(ii)/Scales%PWidth(ii)
    end do

    end subroutine ParamsToVect


    function ffn(n,vect) result(like)
    use ParamDef
    use CalcLike
    implicit none
    integer, intent(in) :: n
    real(Powell_CO_prec) :: like
    real(Powell_CO_prec) vect(n)
    Type(ParamSet) P

    P = MinParams
    call VectToParams(vect, P)

    like = GetLogLike(P)

    call AcceptReject(.true., P%Info,MinParams%Info)
    
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
    call ParamsToVect(MinParams, vect)
    !normalized parameter bounds
    XL = Scales%PMin(params_used) / Scales%PWidth(params_used)
    XU = Scales%PMax(params_used) / Scales%PWidth(params_used)  

    !Initial and final radius of region required (in normalized units)
    rhobeg = 0.1
    rhoend = sigma_frac_err  

    npt = 2*num_params_used +1
    if (.not. BOBYQA (ffn, num_params_used ,npt, vect,XL,XU,rhobeg,rhoend,FeedBack+1, max_iterations)) &
      stop 'FindBestFit failed'

    best_like = ffn(num_params_used,vect)

    !back to real units
    call VectToParams(vect, Params)
    
    end function FindBestFit


    subroutine WriteBestFitParams(like, Params, fname)
    integer i
    real like
    Type(ParamSet) Params
    character(LEN=*), intent(in) :: fname

    if (Feedback>0) write (*,*) 'Best fit parameters, -log(like) = ', like
    call CreateTxtFile(fname,tmp_file_unit)
    write (tmp_file_unit,*) 'Best fit sample -log(Like) = ',like
    write (tmp_file_unit,*) ''
    do i=1, num_params_used
        write(tmp_file_unit,'(1I5,1E15.7,"   '//trim(ParamNames_name(NameMapping,params_used(i)))//'")') &
           params_used(i), Params%P(params_used(i))
        if (Feedback>0) write (*,'(1I5,1E15.7,"   '//trim(ParamNames_name(NameMapping,params_used(i)))//'")') &
           params_used(i), Params%P(params_used(i))
    end do
    close(tmp_file_unit)

    end subroutine WriteBestFitParams

    end module minimize