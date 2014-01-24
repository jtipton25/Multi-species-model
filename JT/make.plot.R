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
  #
  curve(dbeta(x, alpha.p, beta.p))
  curve(dbeta(x, mean(out$alpha.p.save[(n.burn + 1):n.mcmc]), mean(out$beta.p.save[(n.burn + 1):n.mcmc])), add = TRUE, col = 'blue')
  #
  curve(dbeta(x, alpha.psi, beta.psi))
  curve(dbeta(x, mean(out$alpha.psi.save[(n.burn + 1):n.mcmc]), mean(out$beta.psi.save[(n.burn + 1):n.mcmc])), add = TRUE, col = 'blue')
  
}