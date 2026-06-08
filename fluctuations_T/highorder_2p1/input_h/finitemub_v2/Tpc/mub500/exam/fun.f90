subroutine Fnb(x,T,nb)          
!boson distribution function
  implicit none

  real(16) x(18), T
  real(16) nb(18)
  real(16) over(18)
  integer i

  over=x/T
  do i=1,18
  if(over(i) > 5.Q+03)then
    nb(i)=0.Q+00
  else
    nb(i)=1.Q+00/(exp(over(i))-1.Q+00)
  end if
  end do
end

subroutine Fnf0(x,T,l,lb,nf0,nf1,nf2)
!fermion distribution function
  implicit none

  real(16) x(3),T,l,lb
  real(16) nf0(3),nf1(3),nf2(3)
  real(16) over(3)
  integer i

  over=x/T
  do i=1,3
  if(over(i) > 5.Q+03)then
    nf0(i)=0.Q+00
    nf1(i)=0.Q+00
    nf2(i)=0.Q+00
  else
    nf0(i)=1.Q+00/(1.Q+00+3.Q+00*lb*exp(over(i))+3.Q+00*l*exp(over(i))**2+exp(over(i))**3)
    nf1(i)=(2.Q+00*exp(over(i)))/(1.Q+00+3.Q+00*lb*exp(over(i))+3.Q+00*l*exp(over(i))**2+exp(over(i))**3)
    nf2(i)=(exp(over(i))**2)/(1.Q+00+3.Q+00*lb*exp(over(i))+3.Q+00*l*exp(over(i))**2+exp(over(i))**3)
  end if
  end do
end
