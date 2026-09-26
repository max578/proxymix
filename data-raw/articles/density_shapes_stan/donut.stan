parameters { vector[2] x; }
model { target += -square(sqrt(dot_self(x)) - 2.5) / (2 * 0.25); }
