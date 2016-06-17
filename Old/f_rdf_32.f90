module f_rdf
use iso_fortran_env
contains

subroutine pair_dist_arr(xyz, n, a)
! generate pairs of distances from an xyz matrix xyz of size (n, 3)
! 09/01/16
!    integer, parameter :: k = selected_int_kind(16)
!    integer, parameter :: dp = selected_int_kind(8)
    integer, intent(in) :: n
    real, intent(in) :: xyz(n, 3)
    real, intent(out) :: a(n*(n-1)/2)
    integer :: i, j, cnt

    cnt = 1
    do i = 1, n
        do j = i+1, n
            a(cnt) = sqrt(sum((xyz(i, :) - xyz(j, :))**2))
            cnt = cnt + 1
        enddo
    enddo
end subroutine

subroutine pair_dist_arr2(xyz, n, xyz2, n2, a)
! generate pairs of mutual distances for two xyz matrices
! 09/01/16
!    integer, parameter :: k = selected_int_kind(16)
    integer, parameter :: dp = selected_int_kind(8)
    integer, intent(in) :: n
    integer, intent(in) :: n2
    real(dp), intent(in) :: xyz(n, 3)
    real(dp), intent(in) :: xyz2(n2, 3)
    real(dp), intent(out) :: a(n*n2)
    integer(kind=int64) :: i, j, cnt

    cnt = 1
    do i = 1, n
        do j = 1, n2
            a(cnt) = sqrt(sum((xyz(i, :) - xyz2(j, :))**2))
            cnt = cnt + 1
        enddo
    enddo  
end subroutine 

subroutine histogram(a, n, n_b, hist, bins)
! Create a histogram with n_b+1 bins for vector a
! 10/01/16
    integer, intent(in) :: n, n_b
    real, intent(in) :: a(n)
    real :: dr
    integer, intent(out) :: hist(n_b)
    real, intent(out) :: bins(n_b+1)
    integer :: i, j

    bins(1) = minval(a)
    dr = (maxval(a) - minval(a))/n_b
    do i = 1, n_b+1
        bins(i+1) = bins(i) + dr
    enddo
    
    do j = 1, size(a)
        do i = 1, n_b
            if (a(j) > bins(i) .and. a(j) < bins(i+1)) then
                hist(i) = hist(i) + 1
                exit
            endif
        enddo
    enddo
end subroutine

subroutine pair_dist_hist(xyz, n, n_b, hist, bins)
! Putting functions pair_dist_mat and histogram together
! 10/01/16
    integer, intent(in) :: n, n_b
    real, intent(in) :: xyz(n, 3)
    real :: a(n*(n-1)/2)
    integer, intent(out) :: hist(n_b)
    real, intent(out) :: bins(n_b+1)

    call pair_dist_arr(xyz, n, a)
    call histogram(a, n*(n-1)/2, n_b, hist, bins)
end subroutine

end module
