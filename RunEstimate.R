source("GenerateData.R")  # simulate_etas_modified を入れてある想定

theta_true <- list(lambda0=2, alpha1=-0.01,
                   A=0.30, alpha=1.0,
                   p=1.2, c=0.01, D_trig=3, beta_m=1.2)

mu_true <- rep(1/(4*2^8), 4*2^8)

dat <- simulate_etas_modified(
        T_end=100,
        lambda0=theta_true$lambda0, alpha1=theta_true$alpha1,
        mu_vec=mu_true, b_loc=8,
        beta_m=theta_true$beta_m,
        A=theta_true$A, alpha=theta_true$alpha,
        p=theta_true$p, c=theta_true$c, D_trig=theta_true$D_trig,
        max_events=50000, seed=1
)

res <- fit_etas_mu_Aapc(
        t=dat$time_t, x=dat$moduleID_x, loc_norm=dat$LOCnorm_y, m=dat$difficulty_m,
        lambda0=theta_true$lambda0, alpha1=theta_true$alpha1,
        b_loc=8, D_trig=theta_true$D_trig, D_mu=3,
        A_init=theta_true$A, alpha_init=theta_true$alpha,
        p_init=theta_true$p, c_init=theta_true$c,
        max_iter=30, verbose=TRUE
)

saveRDS(list(theta_true=theta_true, dat=dat, res=res), "run1.rds")
write.csv(res$history, "fit_history.csv", row.names=FALSE)
write.csv(cbind(res$state_space, mu=res$mu), "mu_hat.csv", row.names=FALSE)
write.csv(dat, "sim_data.csv", row.names=FALSE)

#背景確率可視化
b_loc <- 8
n_loc <- 2^b_loc
mu_mat <- matrix(res$mu, nrow = 4, ncol = n_loc, byrow = TRUE)
# そのまま描く（画像）
image(t(mu_mat),
      xlab = "moduleID (0..3)", ylab = "LOC bin (0..255)")
# 軸ラベルをちゃんと出したい場合
axis(1, at=seq(0,1,length.out=4), labels=0:3)
axis(2, at=c(0,0.5,1), labels=c(0,128,255))
#ビンを０１に
loc_grid <- seq(0, 1, length.out = n_loc)
image(x = 0:3, y = loc_grid, z = mu_mat,
      xlab = "moduleID", ylab = "LOCnorm (0..1)")
axis(1, at = 0:3, labels = 0:3)
#各モデュールの折れ線
cols <- 1:4  # matplotのデフォルト色と同じ並び（黒/赤/緑/青系）
matplot(loc_grid, t(mu_mat), type = "l", lty = 1, col = cols,
        xlab = "LOCnorm (0..1)", ylab = "mu")
legend("topright",
       legend = paste0("module ", 0:3),
       col = cols, lty = 1, bty = "n")
#４パネル
plot_mu_by_module(res, b_loc = 8, file = "mu_by_module_4panel.png")
# 収束サマリ出力
# summarize_history(res0$res$history)
# 確率最大セル
# mu <- res0$res$mu           # 例：長さ 1024 のベクトル
# ss <- res0$res$state_space  # 例：moduleID と loc_bin 等の表（同じ行数）
# idx <- which.max(mu)
# mu_max <- mu[idx]
# ss_max <- ss[idx, , drop=FALSE]
# list(idx=idx, mu_max=mu_max, cell=ss_max)
