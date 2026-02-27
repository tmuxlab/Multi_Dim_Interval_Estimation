#s(m)
b_beta_from_MLE <- function(M, Mc = 3.0, dM = 0.1){
        M <- M[M >= Mc]
        b    <- log10(exp(1)) / (mean(M) - (Mc - dM/2))
        beta <- 1 / (mean(M) - (Mc - dM/2))   # = b*log(10)
        list(b = b, beta = beta, n = length(M))
}

s_gr_trunc_pmf <- function(m, beta, Mc = 3.0, Mmax = 6.0, dM = 0.1){
        # m は 0.1刻みのベクトル想定（例: 3.0, 3.1, ...）
        Z <- 1 - exp(-beta * (Mmax - Mc))
        
        # 各丸め値 m に対応するビン境界（端ははみ出さないようクリップ）
        lo <- pmax(Mc,   m - dM/2)
        hi <- pmin(Mmax, m + dM/2)
        
        # ビンが無効なら 0
        out <- rep(0, length(m))
        ok <- lo < hi
        
        # truncated exponential の CDF: F(x)=(1-exp(-beta*(x-Mc)))/Z
        F <- function(x) (1 - exp(-beta * (x - Mc))) / Z
        
        out[ok] <- F(hi[ok]) - F(lo[ok])
        out
}

# -------------------------
# 背景時間レート u(t)
# -------------------------
u_time <- function(t, lambda0, alpha1){
        lambda0 * exp(alpha1 * t)
}

#beta_hat <- b_beta_from_MLE(difficulty, Mc = 3.0, dM = 0.1)$beta
beta_hat <- 1.2
#s_i <- s_gr_trunc_pmf(difficulty, beta = beta_hat, Mc=3.0, Mmax=6.0, dM=0.1)

#k(m)
k_childrenNum <- function(A, alpha, m, m_c = 3.0){
        ifelse(m >= m_c, A * exp(alpha * (m - m_c)), 0)
}

#g(t)
g_time <- function(t, p, c){
        stopifnot(is.finite(p), is.finite(c), p > 1, c > 0)
        out <- rep(0, length(t))
        ok  <- (t > 0)
        out[ok] <- (p-1)/c*(1+t[ok]/c)^(-p)
        out
}

#f(x,y|m)
encode_loc <- function(y, b_loc = 8){
        # y: 0..1 を想定（外れたらクリップ）
        L <- 2^b_loc - 1
        y <- pmin(pmax(y, 0), 1)
        as.integer(round(y * L))
}
encode_xy <- function(x, y_loc_int, b_loc = 8){
        stopifnot(all(x %in% 0:3))
        bitwShiftL(as.integer(x), b_loc) + as.integer(y_loc_int)
}
vec1024_to_mat256x4 <- function(v){
        stopifnot(length(v) == 1024)
        matrix(v, nrow = 256, ncol = 4)  # rows: loc_int 0..255, cols: module 0..3
}
make_z_code <- function(x, loc_norm, b_loc = 8){
        encode_xy(x, encode_loc(loc_norm, b_loc), b_loc)
}
# module固定・loc_gridから z_vec(code) を作る
make_zvec_module_loc <- function(module_id, loc_grid, b_loc = 8){
        L <- 2^b_loc - 1
        loc_int <- as.integer(round(pmin(pmax(loc_grid,0),1) * L))
        module_id * 2^b_loc + loc_int
}
# 0..255 の popcount ルックアップ表
.pop8 <- as.integer(colSums(matrix(as.integer(rawToBits(as.raw(0:255))),nrow = 8)))

# 32bit 整数の popcount（xorした値）を byte ごとに数える
popcount32 <- function(z){
        z <- as.integer(z)
        b1 <- bitwAnd(z, 255L)
        b2 <- bitwAnd(bitwShiftR(z,  8L), 255L)
        b3 <- bitwAnd(bitwShiftR(z, 16L), 255L)
        b4 <- bitwAnd(bitwShiftR(z, 24L), 255L)
        .pop8[b1 + 1L] + .pop8[b2 + 1L] + .pop8[b3 + 1L] + .pop8[b4 + 1L] #1Lはinteger(0)を防ぐ
}

# ハミング距離（ベクトル対応）
hamming_code <- function(code1, code2){
        popcount32(bitwXor(as.integer(code1), as.integer(code2)))
}

# 正規化済みハミング距離カーネル f
f_space_hamming_norm <- function(x1, loc1_norm, x2, loc2_norm,
                                 b_loc = 8, D = 3){
        # エンコード
        c1 <- encode_xy(x1, encode_loc(loc1_norm, b_loc), b_loc)
        c2 <- encode_xy(x2, encode_loc(loc2_norm, b_loc), b_loc)
        
        # ハミング距離
        dH <- hamming_code(c1, c2)
        
        # 正規化定数（B bit の完全空間を仮定）
        B <- 2 + b_loc
        Z <- (1 + exp(-1 / D))^B
        
        # 正規化済み kernel
        exp(-dH / D) / Z
}
## 別方式 hamming_kernel_vec_normと同じだった．
f_hamming_norm_code <- function(code_vec, code_j, D = 3, b_loc = 8){
        Bbits <- 2 + b_loc
        Z <- (1 + exp(-1 / D))^Bbits
        dH <- hamming_code(code_vec, code_j)
        exp(-dH / D) / Z
}
## lambda_z at single time t_eval for ALL z cells (0..1023)
lambda_z_at_time_allz <- function(t_eval, dat, theta, mu_vec,
                                  b_loc = 8, D_trig = 3, Mc = 3.0){
        # 条件付きH_t dat: time_t, moduleID_x, LOCnorm_y, difficulty_m
        t <- dat$time_t
        x <- dat$moduleID_x
        y <- dat$LOCnorm_y
        m <- dat$difficulty_m
        
        z_j <- make_z_code(x, y, b_loc)
        z_all <- 0:(length(mu_vec)-1)
        
        # background part
        lam_bg <- u_time(t_eval, theta$lambda0, theta$alpha1) * mu_vec[z_all + 1L]
        
        # trigger part (only parents with t_j < t_eval)
        idx <- which(t < t_eval)
        if(length(idx) == 0) return(lam_bg)
        
        tau <- t_eval - t[idx]
        gj  <- g_time(tau, theta$p, theta$c)
        kj  <- k_childrenNum(theta$A, theta$alpha, m[idx], Mc)
        
        lam_trig <- rep(0, length(z_all))
        for(a in seq_along(idx)){
                j <- idx[a]
                fj <- f_hamming_norm_code(z_all, z_j[j], D = D_trig, b_loc = b_loc)
                lam_trig <- lam_trig + (kj[a] * gj[a]) * fj
        }
        
        lam_bg + lam_trig
}
# 任意の z_vec だけに対する λz(t,z)（3D用に高速）
lambda_z_at_time_zvec <- function(t_eval, dat, theta, mu_vec, z_vec,
                                  b_loc = 8, D_trig = 3, Mc = 3.0){
        t <- dat$time_t
        x <- dat$moduleID_x
        y <- dat$LOCnorm_y
        m <- dat$difficulty_m
        z_j <- make_z_code(x, y, b_loc)
        
        lam_bg <- u_time(t_eval, theta$lambda0, theta$alpha1) * mu_vec[z_vec + 1L]
        idx <- which(t < t_eval)
        if(length(idx) == 0) return(lam_bg)
        
        tau <- t_eval - t[idx]
        gj  <- g_time(tau, theta$p, theta$c)
        kj  <- k_childrenNum(theta$A, theta$alpha, m[idx], Mc)
        
        lam_trig <- rep(0, length(z_vec))
        for(a in seq_along(idx)){
                j <- idx[a]
                fj <- f_hamming_norm_code(z_vec, z_j[j], D = D_trig, b_loc = b_loc)
                lam_trig <- lam_trig + (kj[a] * gj[a]) * fj
        }
        lam_bg + lam_trig
}
# あるデータ(dat)・パラメータ(theta,mu)で module_id の (t×loc) 格子 λz を計算
calc_lambdaZ_grid <- function(dat, theta, mu_vec, module_id,
                              t_grid, loc_grid,
                              b_loc=8, D_trig=3, Mc=3.0){
        z_vec <- make_zvec_module_loc(module_id, loc_grid, b_loc)
        out <- matrix(NA_real_, nrow=length(t_grid), ncol=length(loc_grid))
        for(i in seq_along(t_grid)){
                out[i,] <- lambda_z_at_time_zvec(
                        t_eval=t_grid[i], dat=dat, theta=theta, mu_vec=mu_vec, z_vec=z_vec,
                        b_loc=b_loc, D_trig=D_trig, Mc=Mc
                )
        }
        out
}
#推定後４パネル
plot_mu_by_module <- function(res, b_loc = 8,
                              file = NULL, width = 900, height = 600){
        n_loc <- 2^b_loc
        loc_grid <- seq(0, 1, length.out = n_loc)
        
        mu_mat <- matrix(res$mu, nrow = 4, ncol = n_loc, byrow = TRUE)
        
        if(!is.null(file)){
                png(file, width = width, height = height)
                on.exit(dev.off(), add = TRUE)
        }
        
        op <- par(mfrow = c(2,2), mar = c(4,4,3,1))
        on.exit(par(op), add = TRUE)
        
        for(mod in 0:3){
                plot(loc_grid, mu_mat[mod+1, ],
                     type = "l", lty = 1,
                     xlab = "LOCnorm (0..1)",
                     ylab = expression(mu),
                     main = paste0("Background mu: module ", mod))
                grid()
        }
        
        invisible(mu_mat)
}
#推定後収束サマリ確認
summarize_history <- function(history, tail_k = 5){
        diff_end  <- tail(history$diff, 1)
        diff_start<- head(history$diff, 1)
        ll_end    <- tail(history$loglik, 1)
        ll_start  <- head(history$loglik, 1)
        ll_gain   <- ll_end - ll_start
        ll_inc    <- diff(tail(history$loglik, tail_k + 1))
        list(
                iter_end = tail(history$iter, 1),
                diff_start = diff_start,
                diff_end   = diff_end,
                diff_ratio = diff_end / diff_start,
                loglik_start = ll_start,
                loglik_end   = ll_end,
                loglik_gain  = ll_gain,
                loglik_inc_last = tail(ll_inc, 1),
                loglik_inc_mean_tail = mean(ll_inc)
        )
}