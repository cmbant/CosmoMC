
    module Samples
    use ObjectLists
    use settings

    integer, parameter :: sample_prec = mcp

    Type, extends(TRealArrayList):: TSampleList
    contains
    procedure :: ConfidVal => TSampleList_ConfidVal
    end Type TSampleList

    contains
    
    subroutine GelmanRubinEvalues(cov, meanscov, evals, num)
     use MatrixUtils
     integer, intent(in) :: num
     real(sample_prec) :: cov(num,num), meanscov(num,num), evals(num)
     integer jj
     real(sample_prec) rot(num,num), rotmeans(num,num)
     real(sample_prec) :: sc
     
     rot = cov
     rotmeans = meanscov
     do jj=1,num
                sc = sqrt(cov(jj,jj))
                rot(jj,:) = rot(jj,:) / sc
                rot(:,jj) = rot(:,jj) / sc
                rotmeans(jj,:) = rotmeans(jj,:) /sc
                rotmeans(:,jj) = rotmeans(:,jj) /sc
     end do

     call Matrix_CholeskyRootInverse(rot)
     rotmeans =  matmul(matmul(rot, rotmeans), transpose(rot))
     call Matrix_Diagonalize(rotmeans, evals, num)
    
    end subroutine GelmanRubinEvalues
         

    subroutine TSampleList_ConfidVal(L, ix, limfrac, ix1, ix2, Lower, Upper)
    !Taking the ix'th entry in each array to be a sample, value for which
    !limfrac of the items between ix1 and ix2 (inc) are above or below
    !e.g. if limfrac = 0.05 get two tail 90% confidence limits
    Class(TSampleList) :: L
    Type(TSampleList) :: SortItems
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
    call SortItems%Clear(itemsOnly=.true.)

    end subroutine TSampleList_ConfidVal

    end module Samples

