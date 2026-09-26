parameters { vector[2] x; }
model {
  target += log_sum_exp([
    log(0.3) + multi_normal_lpdf(x | [-2, -2]', [[0.6, 0], [0, 0.6]]),
    log(0.4) + multi_normal_lpdf(x | [0, 0]', [[0.5, 0.2], [0.2, 0.5]]),
    log(0.3) + multi_normal_lpdf(x | [2, 2]', [[0.4, -0.1], [-0.1, 0.4]])]');
}
