# --- 便利関数 ---
U_cum <- function(t, lambda0, alpha1){
        if(abs(alpha1) < 1e-12) return(lambda0 * t)
        (lambda0/alpha1) * (exp(alpha1 * t) - 1)
}

G_cum <- function(tau, p, c){
        out <- rep(0, length(tau))
        ok <- tau > 0
        out[ok] <- 1 - (1 + tau[ok]/c)^(1 - p)
        out
}

k_childrenNum <- function(A, alpha, m, Mc=3.0){
        ifelse(m >= Mc, A * exp(alpha * (m - Mc)), 0)
}

# --- Λ_total(t) を「データ(t_i,m_i)とパラメータ」から計算 ---
Lambda_total <- function(t_grid, dat, theta, Mc=3.0){
        t_i <- dat$time_t
        m_i <- dat$difficulty_m
        
        # 背景
        bg <- U_cum(t_grid, theta$lambda0, theta$alpha1)
        
        # 自己励起：各時刻tで sum_j k(m_j) G(t - t_j)
        km <- k_childrenNum(theta$A, theta$alpha, m_i, Mc=Mc)
        
        trig <- vapply(t_grid, function(tt){
                idx <- which(t_i < tt)
                if(length(idx)==0) return(0)
                sum(km[idx] * G_cum(tt - t_i[idx], theta$p, theta$c))
        }, numeric(1))
        
        bg + trig
}

# --- メイン：保存結果を読んでΛ(t)とCIを作る ---
expected_count_from_saved <- function(
                out_dir="output/boot_with_lambda",
                t_grid=seq(0, 100, by=1),
                Mc=3.0
){
        rds_path <- file.path(out_dir, "bootstrap_result_with_lambda.rds")
        dat_path <- file.path(out_dir, "data_original.csv")
        boot_path <- file.path(out_dir, "boot_theta.csv")
        
        obj <- readRDS(rds_path)
        dat <- read.csv(dat_path)
        
        # 点推定（theta_hat）
        theta_hat <- obj$theta_hat
        Lam_hat <- Lambda_total(t_grid, dat, theta_hat, Mc=Mc)
        
        # ブートパラメータ標本（CI用）
        boot <- read.csv(boot_path)
        ok <- !boot$failed & is.finite(boot$A) & is.finite(boot$alpha) & is.finite(boot$p) & is.finite(boot$c)
        boot_ok <- boot[ok, c("A","alpha","p","c")]
        
        # Λ(t) を各ブート標本で計算（t_grid×B）
        Lam_mat <- sapply(seq_len(nrow(boot_ok)), function(b){
                th <- theta_hat
                th$A <- boot_ok$A[b]
                th$alpha <- boot_ok$alpha[b]
                th$p <- boot_ok$p[b]
                th$c <- boot_ok$c[b]
                Lambda_total(t_grid, dat, th, Mc=Mc)
        })
        
        Lam_lo <- apply(Lam_mat, 1, quantile, probs=0.025, na.rm=TRUE)
        Lam_hi <- apply(Lam_mat, 1, quantile, probs=0.975, na.rm=TRUE)
        
        out <- data.frame(t=t_grid, Lambda_hat=Lam_hat, Lambda_lo=Lam_lo, Lambda_hi=Lam_hi)
        
        # CSV保存
        write.csv(out, file.path(out_dir, "Lambda_total_CI.csv"), row.names=FALSE)
        
        # 簡易プロット（base）
        png(file.path(out_dir, "Lambda_total_CI.png"), width=1200, height=800, res=150)
        ymax <- max(out$Lambda_hi, na.rm = TRUE)
        ymin <- min(out$Lambda_lo, na.rm = TRUE)
        pad  <- 0.03 * (ymax - ymin)  # 3% 余白（好みで 0.05 など）
        
        plot(out$t, out$Lambda_hat, type="l",
             xlab="t", ylab="Expected count (Lambda)",
             main="Cumulative expected number of events",
             ylim = c(ymin, ymax + pad))
        lines(out$t, out$Lambda_lo, lty=2)
        lines(out$t, out$Lambda_hi, lty=2)
        legend("topleft", legend=c("point estimate","95% CI"),
               lty=c(1,2), bty="n")
        dev.off()
        
        # 簡易プロット（base）: PDF保存版
        pdf(file.path(out_dir, "Lambda_total_CI.pdf"),
            width = 8, height = 5, family = "Helvetica")  # サイズは好みで調整
        
        op <- par(mar = c(4.2, 4.2, 2.5, 1.0))  # 余白調整（任意）
        
        plot(out$t, out$Lambda_hat, type="l",
             xlab="t", ylab="Expected count (Lambda)",
             main="Cumulative expected number of events")
        lines(out$t, out$Lambda_lo, lty=2)
        lines(out$t, out$Lambda_hi, lty=2)
        legend("topleft", legend=c("point estimate","95% CI"),
               lty=c(1,2), bty="n")
        
        par(op)
        dev.off()
        
        out
}

# ---- 実行例 ----
Lam_res <- expected_count_from_saved(
        out_dir="output/boot_with_lambda",
        t_grid=seq(0, 100, by=1),
        Mc=3.0
)
head(Lam_res)
