    !Gets predicted BBN mass fraction YHe by interpolating in pre-computed grid

    !Module originally from Jan Hamann, 4/2010
    !Rewrittten with updated table from AlterBBN by AL Dec 2013
    module bbn
    use settings
    use Interpolation
    implicit none

    type, extends(TInterpGrid2D) :: TBBNPredictions
        integer :: data_col = 4
    contains
    procedure :: Error => BBNPredictions_error
    procedure :: FirstUse =>  BBN_Init
    end type TBBNPredictions

    type(TBBNPredictions), target, save :: BBN_YHe
    !Helium mass fraction (not Y_P^BBN nucleon fraction, which is column 5)

    type(TBBNPredictions), target, save :: BBN_YpBBN = TBBNPredictions(data_col=5)
    type(TBBNPredictions), target, save :: BBN_YpBBN_err = TBBNPredictions(data_col=6)
    type(TBBNPredictions), target, save :: BBN_DH = TBBNPredictions(data_col=7)
    type(TBBNPredictions), target, save :: BBN_DH_err = TBBNPredictions(data_col=8)

    character(LEN=*), parameter :: BBN_data_file = 'BBN_full_Parthenelope_880.3.dat'

    contains

    subroutine BBN_Init(this)
    class(TBBNPredictions):: this

    if (feedback > 1) print*,'Initialising BBN data '//BBN_data_file//', col:', this%data_col

    call this%InitFromFile(trim(DataDir)//BBN_data_file, xcol=1,ycol=3, zcol=this%data_col)

    if (feedback > 1) print*,'Done. Interpolation table is ', this%nx,' by ',this%ny

    end subroutine BBN_Init

    subroutine BBNPredictions_error(this,S,v1,v2)
    class(TBBNPredictions):: this
    character(LEN=*), intent(in) :: S
    class(*), intent(in), optional :: v1, v2
    
    call this%TInterpGrid2D%Error('BBN Error: '//S,v1,v2)

    end subroutine BBNPredictions_error


    end module bbn
