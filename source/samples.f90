    module SampleList
    use ObjectLists
    use settings

    integer, parameter :: sample_prec = mcp

    Type, extends(TObjectList):: TRealArrayList
    contains
    procedure :: Value
    procedure :: RealItem
    procedure :: Compare => CompareReal
    procedure :: ConfidVal
    generic :: Item => Value, RealItem
    end Type TRealArrayList

    contains

    function RealItem(L, i) result(P)
    Class(TRealArrayList) :: L
    integer, intent(in) :: i
    real(sample_prec), pointer :: P(:)
    class(*), pointer :: Item(:)

    Item => L%ArrayItem(i)
    select type (pt=>Item)
    type is (real(kind=sample_prec))
        P=> pt
        class default
        stop 'TRealArrayList: object of wrong type'
    end select

    end function RealItem

    function Value(L, i, j) result(P)
    Class(TRealArrayList) :: L
    integer, intent(in) :: i, j
    real(sample_prec) :: P
    class(*), pointer :: C

    C => L%ArrayItemIndex(i,j)
    select type (Arr=> C)
    Type is (real(sample_prec))
        P = Arr
    end select

    end function Value

    integer function CompareReal(this, R1, R2) result(comp)
    Class(TRealArrayList) :: this
    class(*) R1,R2 
    real(sample_prec) R

    select type (RR1 => R1)
    type is (real(sample_prec))
        select type (RR2 => R2)
        type is (real(sample_prec))
            R = RR1-RR2
            if (R< 0) then
                comp =-1
            elseif (R>0) then
                comp = 1
            else 
                comp = 0
            end if
            return
        end select
    end select


    end function CompareReal

    subroutine ConfidVal(L, ix, limfrac, ix1, ix2, Lower, Upper)
    !Taking the ix'th entry in each array to be a sample, value for which
    !limfrac of the items between ix1 and ix2 (inc) are above or below
    !e.g. if limfrac = 0.05 get two tail 90% confidence limits
    Class(TRealArrayList) :: L
    Type(TRealArrayList) :: SortItems
    integer, intent(IN) :: ix
    real(sample_prec), intent(IN) :: limfrac
    real(sample_prec), intent(OUT), optional :: Lower, Upper
    integer, intent(IN), optional :: ix1,ix2
    integer b,t,samps
    real(sample_prec) pos, d

    b=1
    t=L%Count
    if (present(ix1)) b = ix1
    if (present(ix2)) t = ix2
    samps = t - b + 1
    call SortItems%AssignPointers(L, b, t)
    call SortItems%SortArr(ix)

    if (present(Lower)) then
        pos = (samps-1)*limfrac + 1 
        b = max(int(pos),1)
        Lower = SortItems%Value(b, ix)
        if (b < samps .and. pos>b) then
            d = pos - b
            Lower = Lower*(1 - d) + d * SortItems%Value(b+1,ix)
        end if
    end if
    if (present(Upper)) then
        pos = (samps-1)*(1.-limfrac) + 1
        b = max(int(pos),1)
        Upper = SortItems%Value(b,ix)
        if (b < samps .and. pos>b) then
            d = pos - b
            Upper = Upper*(1 - d) + d * SortItems%Value(b+1,ix) 
        end if
    end if

    end subroutine ConfidVal

    end module SampleList

