module f_rdf
contains

subroutine pair_dist_mat(xyz, n, arr)
! generate pairs of distances from an xyz matrix xyz of size (n, 3)
! 09/01/16
    integer, parameter :: k = selected_int_kind(16)
    integer, parameter :: dp = selected_int_kind(8)
    integer(k), intent(in) :: n
    real(dp), intent(in) :: xyz(n, 3)
    real(dp), intent(out) :: arr(n*(n-1)/2)
    integer :: i, j, cnt

    cnt = 1
    do i = 1, n
        do j = i+1, n
            arr(cnt) = sqrt(sum((xyz(i, :) - xyz(j, :))**2))
            cnt = cnt + 1
        enddo
    enddo
end subroutine

subroutine pair_dist_mat2(xyz, n, xyz2, n2, arr)
! generate pairs of mutual distances for two xyz matrices
! 09/01/16
    integer, parameter :: k = selected_int_kind(16)
    integer, parameter :: dp = selected_int_kind(8)
    integer(k), intent(in) :: n
    integer(k), intent(in) :: n2
    real(dp), intent(in) :: xyz(n, 3)
    real(dp), intent(in) :: xyz2(n2, 3)
    real(dp), intent(out) :: arr(n*(n-1)/2)
    integer :: i, j, cnt

    cnt = 1
    do i = 1, n
        do j = 1, n2
            arr(cnt) = sqrt(sum((xyz(i, :) - xyz2(j, :))**2))
            cnt = cnt + 1
        enddo
    enddo  
end subroutine
end module
