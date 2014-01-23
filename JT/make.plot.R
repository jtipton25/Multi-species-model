make.plot <- function(out){
  layout(matrix(1:9, nrow = 3))
  hist(out$N.save[floor(n.mcmc/10) : n.mcmc], breaks = 20)
  abline(v = N, col = 'red')
  #
  plot(out$alpha.p.save[floor(n.mcmc/10) : n.mcmc], type = 'l', main = paste('accept rate', round(out$p.accept, 2)), ylab = 'alpha_p')
  abline(h = alpha.p, col = 'red')
  #
  plot(out$beta.p.save[floor(n.mcmc/10) : n.mcmc], type = 'l', main = paste('accept rate', round(out$p.accept, 2)), ylab = 'beta_p')
  abline(h = beta.p, col = 'red')
  #
  plot(out$alpha.psi.save[floor(n.mcmc/10) : n.mcmc], type = 'l', main = paste('accept rate', round(out$psi.accept, 2)), ylab = 'alpha_psi')
  abline(h = alpha.psi, col = 'red')
  #
  plot(out$beta.psi.save[floor(n.mcmc/10) : n.mcmc], type = 'l', main = paste('accept rate', round(out$psi.accept, 2)), ylab = 'beta_psi')
  abline(h = beta.psi, col = 'red')
}