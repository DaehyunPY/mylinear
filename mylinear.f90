module mylinear
    use kind_const
    implicit none
contains


! ==================================================
! DETERMINANT
! ==================================================
! DGETRF -------------------------------------------
function det_real(A)
    real   (dp), intent(in)  :: A(1:, 1:)
    integer(i8), allocatable :: ipiv(:)
    real   (dp) :: det_real, prod
    integer(i8) :: n, m, lda, info, i
    n     = size(A(:, 1))
    m     = size(A(1, :))
    lda   = size(A(:, 1))
    info  = 0
    if(allocated(ipiv)) deallocate(ipiv)
    allocate(ipiv(1:min(n, m)))
    call DGETRF(m, n, A, lda, ipiv, info)
    if(info /= 0) then
        print *, "Error #501: LAPACK DGETRF", info
    end if
    prod = 1.d0
    do i = 1, min(n, m)
        prod = prod*A(i, i)
        if(ipiv(i) /= i) prod = -prod
    end do
    det_real = prod
    if(allocated(ipiv)) deallocate(ipiv)
end function det_real
! ZGETRF -------------------------------------------
function det_cmplx(A)
    complex(dp), intent(in)  :: A(1:, 1:)
    integer(i8), allocatable :: ipiv(:)
    complex(dp) :: det_cmplx, prod
    integer(i8) :: n, m, lda, info, i
    n     = size(A(:, 1))
    m     = size(A(1, :))
    lda   = size(A(:, 1))
    info  = 0
    if(allocated(ipiv)) deallocate(ipiv)
    allocate(ipiv(1:min(n, m)))
    info  = 0
    call ZGETRF(m, n, A, lda, ipiv, info)
    if(info /= 0) then
        print *, "Error #502: LAPACK ZGETRF", info
    end if
    prod = 1.d0
    do i = 1, min(n, m)
        prod = prod*A(i, i)
        if(ipiv(i) /= i) prod = -prod
    end do
    det_cmplx = prod
    if(allocated(ipiv)) deallocate(ipiv)
end function det_cmplx
! end determinant ----------------------------------


! ==================================================
! INVERSE
! ==================================================
! DGETRF/DGETRI ------------------------------------
subroutine inverse_real(A)
    real   (dp), intent(inout) :: A(1:, 1:)
    real   (dp), allocatable   :: work(:)
    integer(i8), allocatable   :: ipiv(:)
    integer(i8) :: lwork, n, m, lda, info
    lwork = -1
    n     = size(A(:, 1))
    m     = size(A(1, :))
    lda   = size(A(:, 1))
    info  = 0
    if(allocated(work)) deallocate(work)
    if(allocated(ipiv)) deallocate(ipiv)
    allocate(work(1))
    allocate(ipiv(1:min(n, m)))
    call DGETRI(n, A, lda, ipiv, work, lwork, info)
    if(info /= 0) then
        print *, "Error #511: LAPACK DGETRI", info
    end if
    lwork = int(work(1))
    if(allocated(work)) deallocate(work)
    info  = 0
    call DGETRF(m, n, A, lda, ipiv, info)
    if(info /= 0) then
        print *, "Error #512: LAPACK DGETRF", info
    end if
    info  = 0
    if(allocated(work)) deallocate(work)
    allocate(work(1:lwork))
    call DGETRI(n, A, lda, ipiv, work, lwork, info)
    if(info /= 0) then
        print *, "Error #513: LAPACK DGETRI", info
    end if
    if(allocated(work)) deallocate(work)
    if(allocated(ipiv)) deallocate(ipiv)
end subroutine inverse_real
! ZGETRF/ZGETRI ------------------------------------
subroutine inverse_cmplx(A)
    complex(dp), intent(inout) :: A(1:, 1:)
    complex(dp), allocatable   :: work(:)
    integer(i8), allocatable   :: ipiv(:)
    integer(i8) :: lwork, n, m, lda
    integer(i8) :: info
    lwork = -1
    n     = size(A(:, 1))
    m     = size(A(1, :))
    lda   = size(A(:, 1))
    info  = 0
    if(allocated(work)) deallocate(work)
    if(allocated(ipiv)) deallocate(ipiv)
    allocate(work(1))
    allocate(ipiv(1:min(n, m)))
    call ZGETRI(n, A, lda, ipiv, work, lwork, info)
    if(info /= 0) then
        print *, "Error #514: LAPACK ZGETRI", info
    end if
    lwork = int(real(work(1)))
    if(allocated(work)) deallocate(work)
    info  = 0
    call ZGETRF(m, n, A, lda, ipiv, info)
    if(info /= 0) then
        print *, "Error #515: LAPACK ZGETRF", info
    end if
    info  = 0
    if(allocated(work)) deallocate(work)
    allocate(work(1:lwork))
    call ZGETRI(n, A, lda, ipiv, work, lwork, info)
    if(info /= 0) then
        print *, "Error #516: LAPACK ZGETRI", info
    end if
    if(allocated(work)) deallocate(work)
    if(allocated(ipiv)) deallocate(ipiv)
end subroutine inverse_cmplx
! end inverse --------------------------------------


! ==================================================
! SOLVE
! ==================================================
! DPBSV --------------------------------------------
subroutine solve_sym_band_vec(AB, B, info)
    character(1), parameter :: uplo = 'Upper'
    real(dp), intent(in) :: AB(1:, 1:)
    real(dp), intent(inout) :: B(1:)
    integer(i8), intent(out) :: info
    integer(i8), parameter :: nrhs = 1
    integer(i8) :: n, kd, ldab, ldb
    ldab = size(AB(:, 1))
    n    = size(AB(1, :))
    ldb  = size(B(:))
    kd   = ldab -1
    info = 0
    call DPBSV(uplo, n, kd, nrhs, AB, ldab, B, ldb, info)
    if(info < 0) then
        print *, "Error LAPACK DPBSV: There is an illegal value, ", -info
    else if(info > 0) then
        print *, "Warning LAPACK DPBSV: There is an order of A that is not positive definite, ", info
    end if
end subroutine solve_sym_band_vec
! DGBSV --------------------------------------------
subroutine solve_band_vec(AB, B, info)
    real(dp), intent(in) :: AB(1:, 1:)
    real(dp), intent(inout) :: B(1:)
    integer(i8), intent(out) :: info
    integer(i8), allocatable :: ipiv(:)
    integer(i8), parameter :: nrhs = 1
    integer(i8) :: n, kl, ku, ldab, ldb
    ldab = size(AB(:, 1))
    n    = size(AB(1, :))
    ldb  = size(B(:))
    kl   = (ldab -1)/3
    ku   = (ldab -1)/3
    info = 0
    if(allocated(ipiv)) deallocate(ipiv)
    allocate(ipiv(1:n))
    call DGBSV(n, kl, ku, nrhs, AB, ldab, ipiv, B, ldb, info)
    if(info < 0) then
        print *, "Error LAPACK DGBSV: There is an illegal value, ", -info
    else if(info > 0) then
        print *, "Warning LAPACK DGBSV: There is an element of A that is exactly zero, ", info
    end if
    if(allocated(ipiv)) deallocate(ipiv)
end subroutine solve_band_vec
! DGESV --------------------------------------------
subroutine solve_real_vec(A, B, info)
    real(dp), intent(in) :: A(1:, 1:)
    real(dp), intent(inout) :: B(1:)
    integer(i8), intent(out) :: info
    integer(i8), allocatable :: ipiv(:)
    integer(i8), parameter :: nrhs = 1
    integer(i8) :: n, lda, ldb
    lda  = size(A(:, 1))
    n    = size(A(1, :))
    ldb  = size(B(:))
    info = 0
    if(allocated(ipiv)) deallocate(ipiv)
    allocate(ipiv(1:n))
    call DGESV(n, nrhs, A, lda, ipiv, B, ldb, info)
    if(info < 0) then
        print *, "Error LAPACK DGESV: There is an illegal value, ", -info
    else if(info > 0) then
        print *, "Warning LAPACK DGESV: There is an element of A that is exactly zero, ", info
    end if
    if(allocated(ipiv)) deallocate(ipiv)
end subroutine solve_real_vec
! DGESV --------------------------------------------
subroutine solve_real_mat(A, B, info)
    real(dp), intent(in) :: A(1:, 1:)
    real(dp), intent(inout) :: B(1:, 1:)
    integer(i8), intent(out) :: info
    integer(i8), allocatable :: ipiv(:)
    integer(i8) :: n, nrhs, lda, ldb
    lda  = size(A(:, 1))
    n    = size(A(1, :))
    ldb  = size(B(:, 1))
    nrhs = size(B(1, :))
    info = 0
    if(allocated(ipiv)) deallocate(ipiv)
    allocate(ipiv(1:n))
    call DGESV(n, nrhs, A, lda, ipiv, B, ldb, info)
    if(info < 0) then
        print *, "Error LAPACK DGESV: There is an illegal value, ", -info
    else if(info > 0) then
        print *, "Warning LAPACK DGESV: There is an element of A that is exactly zero, ", info
    end if
    if(allocated(ipiv)) deallocate(ipiv)
end subroutine solve_real_mat
! ZGESV --------------------------------------------
subroutine solve_cmplx_vec(A, B, info)
    complex(dp), intent(in) :: A(1:, 1:)
    complex(dp), intent(inout) :: B(1:)
    integer(i8), intent(out) :: info
    integer(i8), allocatable :: ipiv(:)
    integer(i8), parameter :: nrhs = 1
    integer(i8) :: n, lda, ldb
    lda  = size(A(:, 1))
    n    = size(A(1, :))
    ldb  = size(B(:))
    info = 0
    if(allocated(ipiv)) deallocate(ipiv)
    allocate(ipiv(1:n))
    call ZGESV(n, nrhs, A, lda, ipiv, B, ldb, info)
    if(info < 0) then
        print *, "Error LAPACK ZGESV: There is an illegal value, ", -info
    else if(info > 0) then
        print *, "Warning LAPACK ZGESV: There is an element of A that is exactly zero, ", info
    end if
    if(allocated(ipiv)) deallocate(ipiv)
end subroutine solve_cmplx_vec
! ZGESV --------------------------------------------
subroutine solve_cmplx_mat(A, B, info)
    complex(dp), intent(in) :: A(1:, 1:)
    complex(dp), intent(inout) :: B(1:, 1:)
    integer(i8), intent(out) :: info
    integer(i8), allocatable :: ipiv(:)
    integer(i8) :: n, nrhs, lda, ldb
    lda  = size(A(:, 1))
    n    = size(A(1, :))
    ldb  = size(B(:, 1))
    nrhs = size(B(1, :))
    info = 0
    if(allocated(ipiv)) deallocate(ipiv)
    allocate(ipiv(1:n))
    call ZGESV(n, nrhs, A, lda, ipiv, B, ldb, info)
    if(info < 0) then
        print *, "Error LAPACK ZGESV: There is an illegal value, ", -info
    else if(info > 0) then
        print *, "Warning LAPACK ZGESV: There is an element of A that is exactly zero, ", info
    end if
    if(allocated(ipiv)) deallocate(ipiv)
end subroutine solve_cmplx_mat
! end solve ----------------------------------------


! ==================================================
! DIAGONALIZATION
! ==================================================
! DSBEV --------------------------------------------
subroutine diag_sym_band(AB)
    character(1), parameter :: jobz = 'No vectors', uplo = 'Upper'
    real    (dp), intent(inout) :: AB(:, :)
    real    (dp), allocatable   :: w(:), z(:, :), work(:)
    integer (i8) :: n, kd, ldab, ldz, info
    ldab = size(AB(:, 1))
    n    = size(AB(1, :))
    kd   = ldab -1
    if(allocated(w))    deallocate(w)
    if(allocated(z))    deallocate(z)
    if(allocated(work)) deallocate(work)
    allocate(w(1:n))
    allocate(z(1, 1))
    allocate(work(1:3*n -2))
    ldz  = size(z(:, 1))
    info = 0
    call DSBEV(jobz, uplo, n, kd, AB, ldab, w, z, ldz, work, info)
    if(info /= 0) then
        print *, "Error #531: LAPACK DSBEV", info
    end if
    AB(:, :) = 0.d0
    AB(1, :) = w(:)
    if(allocated(w))    deallocate(w)
    if(allocated(z))    deallocate(z)
    if(allocated(work)) deallocate(work)
end subroutine diag_sym_band
! DSYEV --------------------------------------------
subroutine diag_sym(A, W)
    character(1), parameter :: jobz = 'Vectors', uplo = 'Upper'
    real(dp), intent(inout) :: A(1:, 1:), W(1:)
    real(dp), allocatable   :: work(:)
    integer(i8) :: lwork, n, lda, info
    lwork = -1
    lda   = size(A(:, 1))
    n     = size(A(1, :))
    info  = 0
    if(allocated(work)) deallocate(work)
    allocate(work(1))
    call DSYEV(jobz, uplo, n, A, lda, W, work, lwork, info)
    if(info /= 0) then
        print *, "Error #532: LAPACK DSYEV", info
    end if
    lwork = int(work(1))
    if(allocated(work)) deallocate(work)
    info  = 0
    if(allocated(work)) deallocate(work)
    allocate(work(1:lwork))
    call DSYEV(jobz, uplo, n, A, lda, W, work, lwork, info)
    if(info /= 0) then
        print *, "Error #533: LAPACK DSYEV", info
    end if
    if(allocated(work)) deallocate(work)
end subroutine diag_sym
! ZHEEV --------------------------------------------
subroutine diag_her(A, W)
    character(1), parameter     :: jobz = 'Vectors', uplo = 'Upper'
    complex (dp), intent(inout) :: A(1:, 1:)
    real    (dp), intent(inout) :: W(1:)
    complex (dp), allocatable   :: work(:)
    real    (dp), allocatable   :: rwork(:)
    integer (i8) :: lwork, n, lda, info
    lwork = -1
    lda   = size(A(:, 1))
    n     = size(A(1, :))
    info  = 0
    if(allocated(work))  deallocate(work)
    if(allocated(rwork)) deallocate(rwork)
    allocate(work(1))
    allocate(rwork(1:max(1, 3*n -2)))
    call ZHEEV(jobz, uplo, n, A, lda, W, work, lwork, rwork, info)
    if(info /= 0) then
        print *, "Error #534: LAPACK ZHEEV", info
    end if
    lwork = int(real(work(1)))
    if(allocated(work)) deallocate(work)
    info  = 0
    if(allocated(work))  deallocate(work)
    allocate(work(1:lwork))
    call ZHEEV(jobz, uplo, n, A, lda, W, work, lwork, rwork, info)
    if(info /= 0) then
        print *, "Error #535: LAPACK ZHEEV", info
    end if
    if(allocated(work))  deallocate(work)
    if(allocated(rwork)) deallocate(rwork)
end subroutine diag_her
! end diagonalization ------------------------------
end module mylinear
