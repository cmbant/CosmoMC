    module CosmoTheory
    use settings
    use cmbtypes
    use GeneralTypes
    use likelihood
    use Interpolation
    implicit none
    private

    Type, extends(TInterpGrid2D) :: TCosmoTheoryPK
        !whether power is stored as log(P_k)
        logical :: islog = .true.
    contains
    procedure :: PowerAt
    end Type TCosmoTheoryPK

    Type, extends(TTheoryPredictions) :: TCosmoTheoryPredictions
        real(mcp) cl(lmax,num_cls_tot)
        !TT, TE, EE (BB) + other C_l (e.g. lensing)  in that order
        real(mcp) sigma_8
        real(mcp) tensor_ratio_r10, tensor_ratio_02
        integer numderived
        real(mcp) derived_parameters(max_derived_parameters)
        !MPK's are interpolator objects now
        !MPK%x = logkh, MPK%y = z (redshift), MPK%z =>actual data array
        !MPK%nx = num_k, MPK%ny = num_z
        type(TCosmoTheoryPK), allocatable :: MPK
        type(TCosmoTheoryPK), allocatable :: NL_MPK
    contains
    procedure :: FreePK
    procedure :: ClsFromTheoryData
    procedure :: WriteTextCls
    !Inherited overrides
    procedure :: WriteTheory => TCosmoTheoryPredictions_WriteTheory
    procedure :: ReadTheory => TCosmoTheoryPredictions_ReadTheory
    procedure :: WriteBestFitData => TCosmoTheoryPredictions_WriteBestFitData
    end Type TCosmoTheoryPredictions

    public TCosmoTheoryPredictions, TCosmoTheoryPK
    contains
    
    function PowerAt(PK,k,z) result(outpower)
    class(TCosmoTheoryPK) PK
    real(mcp), intent(in) :: k,z
    real(mcp) :: logk
    real(mcp) :: outpower
    integer :: error
    
    logk=log(k)
    if(.not. allocated(PK%x)) then
        write(*,*) 'ERROR:  PowerAt least one of your PK arrays is not initialized:'
        write(*,*) '        Make sure you are calling a SetPk and filling your power spectra.'
        write(*,*) '        This error could also mean you are doing importance sampling'
        write(*,*) '        and need to turn on redo_pk.'
        call MPIstop()
    end if
    
    if(PK%islog) then 
        outpower = exp(PK%Value(logk,z))
    else
        outpower = PK%Value(logk,z)
    end if
    
    end function PowerAt    

    subroutine FreePK(T)
    class(TCosmoTheoryPredictions) T
    
    if(allocated(T%MPK))deallocate(T%MPK)
    if(allocated(T%NL_MPK)) deallocate(T%NL_MPK)
    
    end subroutine FreePK   

    subroutine ClsFromTheoryData(T, Cls)
    class(TCosmoTheoryPredictions) T
    real(mcp) Cls(lmax,num_cls_tot)

    Cls(2:lmax,1:num_clsS) =T%cl(2:lmax,1:num_clsS)
    if (num_cls>3 .and. num_ClsS==3) Cls(2:lmax,num_cls)=0
    if (num_cls_ext>0) then
        Cls(2:lmax,num_cls+1:num_cls_tot) =T%cl(2:lmax,num_clsS+1:num_clsS+num_cls_ext)
    end if

    end subroutine ClsFromTheoryData

    subroutine WriteTextCls(T,aname)
    class(TCosmoTheoryPredictions) T
    character (LEN=*), intent(in) :: aname
    integer l
    real(mcp) Cls(lmax,num_cls_tot), nm
    character(LEN=*), parameter :: fmt = '(1I6,*(E15.5))'
    Type(TTextFile) :: F

    call F%CreateFile(aname)
    call T%ClsFromTheoryData(Cls)
    do l = 2, lmax
        nm = 2*pi/(l*(l+1))
        if (num_cls_ext > 0) then
            write (F%unit,fmt) l, cls(l,1:num_cls)/nm, cls(l,num_cls+1:num_cls_tot)
        else
            write (F%unit,fmt) l, cls(l,:)/nm
        end if
    end do
    call F%Close()

    end subroutine WriteTextCls

    subroutine TCosmoTheoryPredictions_WriteTheory(T, F)
    Class(TCosmoTheoryPredictions) T
    class(TFileStream) :: F
    integer, parameter :: varcount = 0
    integer tmp(varcount)
    logical, save :: first = .true.

    if (first .and. new_chains) then
        first = .false.
        write(F%unit) use_LSS, compute_tensors, get_sigma8
        write(F%unit) lmax, lmax_tensor, num_cls, num_cls_ext
        write(F%unit) varcount
        write(F%unit) tmp(1:varcount)
    end if

    write(F%unit) T%numderived
    write(F%unit) T%derived_parameters(1:T%numderived)
    write(F%unit) T%cl(2:lmax,1:num_cls)
    if (num_cls_ext>0) write(F%unit) T%cl(2:lmax,num_cls+1:num_cls_tot)

    if (compute_tensors) then
        write(F%unit) T%tensor_ratio_02, T%tensor_ratio_r10
    end if

    if (get_sigma8 .or. use_LSS) write(F%unit) T%sigma_8

    if (use_LSS) then
        write(F%unit) T%MPK%nx, T%MPK%ny
        write(F%unit) T%MPK%x
        write(F%unit) T%MPK%y
        write(F%unit) T%MPK%z
        write(F%unit) use_nonlinear
        if(use_nonlinear) write(F%unit) T%NL_MPK%z
    end if

    end subroutine TCosmoTheoryPredictions_WriteTheory

    subroutine TCosmoTheoryPredictions_ReadTheory(T, F)
    Class(TCosmoTheoryPredictions) T
    class(TFileStream) :: F
    integer unused
    logical, save :: first = .true.
    logical, save :: has_sigma8, has_LSS, has_tensors
    integer, save :: almax, almaxtensor, anumcls, anumclsext, tmp(1)
    !JD 02/14 new variables for handling new pk arrays
    integer :: num_k, num_z, stat
    logical  has_nonlinear
    real(mcp), pointer :: temp(:,:) =>null()
    real(mcp), allocatable :: k(:), z(:)

    if (first) then
        first = .false.
        read(F%unit) has_LSS, has_tensors, has_sigma8
        read(F%unit) almax, almaxtensor, anumcls, anumclsext
        if (almax > lmax) call MpiStop('ReadTheory: reading file with larger lmax')
        if (anumcls /= num_cls) call MpiStop('ReadTheory: reading file with different Cls')
        if (anumclsext /= num_cls_ext) call MpiStop('ReadTheory: reading file with different ext Cls')
        read(F%unit) unused
        read(F%unit) tmp(1:unused)
    end if

    T%cl = 0
    T%derived_parameters=0
    read(F%unit) T%numderived
    read(F%unit) T%derived_parameters(1:T%numderived)
    read(F%unit) T%cl(2:almax,1:anumcls)
    if (anumclsext >0) read(F%unit) T%cl(2:almax,num_cls+1:num_cls+anumclsext)

    if (has_tensors) then
        read(F%unit) T%tensor_ratio_02, T%tensor_ratio_r10
    end if

    if (has_sigma8 .or. has_LSS) read(F%unit) T%sigma_8
    if (has_LSS) then
        call T%FreePK()
        allocate(T%MPK)
        read(F%unit) num_k, num_z
        allocate(temp(num_k,num_z))
        allocate(k(num_k))
        allocate(z(num_z)) 
        read(F%unit) k
        read(F%unit) z
        read(F%unit) temp
        call T%MPK%Init(k,z,temp)
        read(F%unit)has_nonlinear
        if(has_nonlinear) then
            allocate(T%NL_MPK)
            read(F%unit)temp
            call T%NL_MPK%Init(k,z,temp)
            if(.not. use_nonlinear) then
                write(*,*)"WARNING:  ReadTheory - Your data files have nonlinear power spectra,"
                write(*,*)"          but you are not using them. Be careful that this"
                write(*,*)"          is what you intended."
            end if
        end if
    end if

    end subroutine TCosmoTheoryPredictions_ReadTheory

    subroutine TCosmoTheoryPredictions_WriteBestFitData(Theory,fnameroot)
    class(TCosmoTheoryPredictions) Theory
    character(LEN=*), intent(in) :: fnameroot

    if (use_CMB) call Theory%WriteTextCls(fnameroot //'.bestfit_cl')

    end subroutine TCosmoTheoryPredictions_WriteBestFitData

    end module CosmoTheory
