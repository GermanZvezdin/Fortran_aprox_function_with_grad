program main
  implicit none
  !---------------------------var----------------------------------------------------------
  real(8), allocatable :: a(:), res(:)
  real(8) :: x
  real(8) :: ans, eps, c, d, p, norm, F
  integer :: n, m, n0, i
  !---------------------------end-var----------------------------------------------------------
  read(*, *) n, c, d, n0 
  m = 2
  n0 = n0 + 1
  allocate(a(1:n0))
  allocate(res(1:n0))
  res(:) = dble(2.123)
  print *, res
  !x =F(a, c, d,n0, n)
  a(:) = 0.0
  call aprox(a, c, d, n, n0, dble(0.0001), res)
  print *, res
  stop
end program main


subroutine grad(a, c, d, n, k , gradient)
	implicit none 
	real(8), intent(in) :: a(1:n), c, d
	integer, intent(in) :: n, k
	real(8), intent(inout) :: gradient(1:n)
	real(8) :: y, p, x
	integer :: i, j
	gradient(:) = 0.0
	do j=1,n,1
		x = c
		do i=0, k,1
			gradient(j) = gradient(j) - 2 * (y(x) - p(a, x, n)) * (x ** (j - 1))
			x = x + (d - c) / k
		end do
	end do
end subroutine grad

function y(x)
	real(8), intent(in) :: x
	real(8) :: y
	y = cos(x / 10) * atan(x)
	return 
end function y
function norm(x, n)
	implicit none
	real(8), intent(in) :: x(1:n)
	integer, intent(in) :: n
	real(8) :: norm 
	norm = sqrt(DOT_PRODUCT(x, x))
	return
end function norm

function p(a, x, n)
  implicit none
  real(8), intent(in) :: a(1:n)
  real(8), intent(in) :: x
  integer, intent(in) :: n
  real(8) :: p
  integer :: i
  p = 0.0
  do i = 1, n, 1
	p = p + a(i) * x**(i - 1)
  end do
  return
end function p

function F(a, c, d, n, k)
	implicit none 
	real(8), intent(in) :: a(1:n), c, d
	integer, intent(in) :: n, k
	real(8) :: y, p, x, F
	integer :: i
	x = c
	print *, x
	F = 0.0
	do i=1, k + 1, 1
		F = F + (y(x) - p(a, x, n)) ** 2
		x = x + (d - c) / k
	end do
	return
end function F
function del(cur_a, old_a, n)
	real(8), intent(in) :: cur_a(1:n), old_a(1:n)
	integer, intent(in) :: n
	integer :: i
	real(8) :: del, m
	del = 0.0
	do i=1, n, 1
		if (cur_a(i) /= 0) then 
			if (abs((cur_a(i) - old_a(i)) / cur_a(i)) > del) then
				del = abs((cur_a(i) - old_a(i)) / cur_a(i))
			end if
		end if
	end do
end function del
function lambda(a, c, d, n, k, g)
	implicit none 
	real(8), intent(in) :: a(1:n), g(1:n), c, d
	integer, intent(in) :: n, k
	real(8) :: y, p, x, F, lambda
	lambda = (F(a - g, c, d, n, k) - F(a + g, c, d, n, k)) / (2 * (F(a - g, c, d, n, k) - 2 * F(a, c, d, n, k) + F(a + g, c, d, n, k)))
	return
end function lambda 
subroutine aprox(cur_a, c, d, k, n, eps, res)
	implicit none 
	real(8), intent(in) :: c, d, eps
	integer, intent(in) :: n, k
	real(8), intent(inout) :: res(1:n), cur_a(1:n)
	real(8), allocatable :: old_a(:), old_grad(:),cur_grad(:), v(:)
	real(8) :: y, p, lambda, betta, delta, norm, del
	integer :: i, j
	allocate(old_a(1:n))
	allocate(old_grad(1:n))
	allocate(cur_grad(1:n))
	allocate(v(1:n))
	delta = 1
	j = 0
	call grad(cur_a, c, d, n, k, v)
	v = -v
	do while (delta > eps)
		call grad(cur_a, c, d, n, k, old_grad)
		old_a = cur_a
		cur_a = old_a + lambda(old_a, c, d, n, k, v) * v
		call grad(cur_a, c, d, n, k, cur_grad)
		betta = (norm(cur_grad, n) ** 2) / (norm(old_grad, n) ** 2)
		v = -cur_grad + betta * v
		delta = del(cur_a, old_a, n)
		j = j + 1
		print *, norm(cur_grad, n)
	end do
	print *, 'Iter', j 
	res = cur_a
end subroutine aprox


















