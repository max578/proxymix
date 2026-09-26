parameters { vector[2] x; }
model {
  real z2 = x[2] - 0.5 * (x[1]^2 - 1);
  target += -0.5 * (x[1]^2 + z2^2);
}
