subroutine calc_gofr(pos, lbox, ngr, natom, ndim, gofr)
  integer, intent(in) :: natom, ndim, ngr
  real*8, intent(in) :: pos(natom, ndim), lbox
  real*8, intent(out) :: gofr(ngr)
  integer i, j, idx
  real*8 drij(ndim), rij, dr
  dr = lbox/2./ngr
  gofr(:) = 0
  do i=1,natom
    do j=i+1,natom
      drij = pos(i, :) - pos(j, :)
      drij = drij - lbox*nint(drij/lbox)
      rij = sqrt(dot_product(drij, drij))
      if (rij .gt. lbox/2.) cycle
      idx = nint(rij/dr)
      if ((idx .gt. ngr) .or. (idx .le. 0)) cycle
      gofr(idx) = gofr(idx) + 1
    end do
  end do
end subroutine

subroutine calc_sofk(kvecs, pos, lbox, natom, ndim, nk, sofk)
  integer, intent(in) :: natom, ndim, nk
  real*8, intent(in) :: kvecs(nk, ndim), pos(natom, ndim), lbox
  real*8, intent(out) :: sofk(nk)
  double complex :: rhok(nk)
  integer ik, j
  real*8 kdotr
  rhok(:) = 0
  do ik=1,nk
    do j=1,natom
      kdotr = dot_product(kvecs(ik, :), pos(j, :))
      rhok(ik) = rhok(ik) + exp(CMPLX(0, 1)*kdotr)
    end do
  end do
  sofk = CONJG(rhok)*rhok/natom
end subroutine
