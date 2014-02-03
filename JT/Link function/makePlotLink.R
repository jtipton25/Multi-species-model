make.plot <- function(out){
  layout(matrix(1:9, nrow = 3))
  hist(out$N.save[n.burn:n.mcmc], breaks = 20, main = "Species Richness", xlab = "N")
  abline(v = N, col = 'red')
  #
  plot(out$alpha.p.save[n.burn:n.mcmc], type = 'l', main = paste('accept rate', round(out$p.accept, 2)), ylab = 'alpha_p', xlab = 'MCMC iteration')
  abline(h = alpha.p, col = 'red')
  #
  plot(out$beta.p.save[n.burn:n.mcmc], type = 'l', main = paste('accept rate', round(out$p.accept, 2)), ylab = 'beta_p', xlab = 'MCMC iteration')
  abline(h = beta.p, col = 'red')
  #
  plot(out$alpha.psi.save[n.burn:n.mcmc], type = 'l', main = paste('accept rate', round(out$psi.accept, 2)), ylab = 'alpha_psi', xlab = 'MCMC iteration')
  abline(h = alpha.psi, col = 'red')
  #
  plot(out$beta.psi.save[n.burn:n.mcmc], type = 'l', main = paste('accept rate', round(out$psi.accept, 2)), ylab = 'beta_psi',, xlab = 'MCMC iteration')
  abline(h = beta.psi, col = 'red')
  #
  curve(dbeta(x, alpha.p, beta.p))
  curve(dbeta(x, mean(out$alpha.p.save[n.burn:n.mcmc]), mean(out$beta.p.save[(n.burn + 1):n.mcmc])), add = TRUE, col = 'blue')
  #
  curve(dbeta(x, alpha.psi, beta.psi))
  curve(dbeta(x, mean(out$alpha.psi.save[(n.burn + 1):n.mcmc]), mean(out$beta.psi.save[(n.burn + 1):n.mcmc])), add = TRUE, col = 'blue')
  
}