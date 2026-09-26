parameters { vector<lower=-1, upper=1>[1] x; }
model { target += log1m(square(x[1])); }
