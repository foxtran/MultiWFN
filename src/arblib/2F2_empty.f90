real(8) function hyp2F2(a1, a2, b1, b2, z) bind(C, name="hyp2F2")
  real(8), intent(in), value :: a1, a2, b1, b2, z
  error stop "Arb Library is not linked!"
end function hyp2F2
