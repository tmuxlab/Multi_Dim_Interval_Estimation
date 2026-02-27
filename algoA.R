source("funclist.R")
# -------------------------
# 状態空間（全セル）とイベントのセルindex
# -------------------------
make_state_space <- function(b_loc = 8){
        loc_int <- 0:(2^b_loc - 1)
        x <- 0:3
        grid <- expand.grid(x = x, loc_int = loc_int)
        grid$code <- bitwShiftL(as.integer(grid$x), b_loc) + as.integer(grid$loc_int)
        grid
}

event_state_index <- function(x, loc_norm, b_loc = 8){
        loc_i <- encode_loc(loc_norm, b_loc)
        as.integer(x) * (2^b_loc) + as.integer(loc_i) + 1L
}

# -------------------------
# ハミング距離カーネル K（全状態へのベクトル）: 正規化済み
# -------------------------
hamming_kernel_vec_norm <- function(state_codes, center_code, b_loc = 8, D = 3){
        B <- 2 + b_loc
        Z <- (1 + exp(-1 / D))^B
        dH <- hamming_code(state_codes, center_code)
        exp(-dH / D) / Z
}

# -------------------------
# g の累積 G(t)=∫_0^t g(s) ds  (p>1, c>0)
# -------------------------
G_time <- function(t, p, c){
        out <- rep(0, length(t))
        ok <- (t > 0)
        out[ok] <- 1 - (1 + t[ok]/c)^(1 - p)
        out
}
# -------------------------
# ∫_0^T u(t) dt  （u(t)=lambda0*exp(alpha1*t)）
# -------------------------
U_time_integral <- function(T_end, lambda0, alpha1){
        if(abs(alpha1) < 1e-12){
                lambda0 * T_end
        } else {
                lambda0 / alpha1 * (exp(alpha1 * T_end) - 1)
        }
}

# -------------------------
# 観測対数尤度（時空間部分）
# l = Σ log λ_i - U(T) - Σ k(m_j) G(T-t_j)
# -------------------------
loglik_obs_ts <- function(t, m, E, lambda0, alpha1, A, alpha, p, c, m_c = 3.0){
        T_end <- max(t)
        k_j   <- k_childrenNum(A, alpha, m, m_c = m_c)
        integ_bg   <- U_time_integral(T_end, lambda0, alpha1)
        integ_trig <- sum(k_j * G_time(T_end - t, p, c))
        sum(log(pmax(E$lambda, 1e-300))) - integ_bg - integ_trig
}

# -------------------------
# Q（期待完全対数尤度：時空間部分）
# Q = Σ φ_i( log u(t_i) + log μ(s_i) )
#   + Σ_{i>j} φ_ij( log k(m_j) + log g(dt) + log f(s_i|s_j) )
#   - U(T) - Σ k(m_j)G(T-t_j)
# -------------------------
Q_ts <- function(t, x, loc_norm, m, mu_vec, E,
                 lambda0, alpha1, A, alpha, p, c,
                 b_loc = 8, D_trig = 3, m_c = 3.0){
        
        n     <- length(t)
        T_end <- max(t)
        
        # 背景項
        u_i  <- u_time(t, lambda0, alpha1)
        mu_i <- mu_vec[E$s_idx]
        term_bg <- sum(E$phi_bg * (log(pmax(u_i, 1e-300)) + log(pmax(mu_i, 1e-300))))
        
        # 親ごとの期待子数 w_j = Σ_{i>j} φ_ij
        wj  <- expected_children_per_parent(E$phi_list)
        dm  <- pmax(0, m - m_c)
        k_j <- k_childrenNum(A, alpha, m, m_c = m_c)
        
        term_k <- sum(wj * log(pmax(k_j, 1e-300)))
        
        # g と f の項（iごとに合算）
        logZ <- (2 + b_loc) * log1p(exp(-1 / D_trig))  # log( (1+exp(-1/D))^B )
        term_gf <- 0
        
        for(i in 2:n){
                w <- E$phi_list[[i]]
                if(length(w) == 0) next
                if(sum(w) <= 0) next
                
                dt <- t[i] - t[1:(i-1)]
                # log g(dt) = log(p-1) - log(c) - p*log(1+dt/c)
                lg <- (log(p - 1) - log(c)) - p * log1p(dt / c)
                
                # log f = -dH/D - logZ（dHはハミング距離）
                c1 <- encode_xy(x[i], encode_loc(loc_norm[i], b_loc), b_loc)
                c2 <- encode_xy(x[1:(i-1)], encode_loc(loc_norm[1:(i-1)], b_loc), b_loc)
                dH <- hamming_code(c1, c2)
                lf <- -(dH / D_trig) - logZ
                
                term_gf <- term_gf + sum(w * (lg + lf))
        }
        
        # 積分項
        integ_bg   <- U_time_integral(T_end, lambda0, alpha1)
        integ_trig <- sum(k_j * G_time(T_end - t, p, c))
        
        term_bg + term_k + term_gf - integ_bg - integ_trig
}
# -------------------------
# E-step: phi_i（背景）と phi_{ij}（親ごと）を計算して保存
# O(n^2) です
# -------------------------
etas_E_step_storePhi <- function(t, x, loc_norm, m,
                                 mu_vec,
                                 lambda0, alpha1,
                                 A, alpha,
                                 p, c,
                                 b_loc = 8,
                                 D_trig = 3,
                                 m_c = 3.0){
        n <- length(t)
        stopifnot(length(x)==n, length(loc_norm)==n, length(m)==n)
        
        s_idx <- event_state_index(x, loc_norm, b_loc)
        u_i   <- u_time(t, lambda0, alpha1)
        k_j   <- k_childrenNum(A, alpha, m, m_c = m_c)
        
        lambda_all <- numeric(n)
        phi_bg     <- numeric(n)
        phi_list   <- vector("list", n)  # phi_list[[i]]: length i-1 のベクトル（親ごと）
        
        for(i in seq_len(n)){
                bg <- u_i[i] * mu_vec[s_idx[i]]
                
                trig_vec <- numeric(0)
                trig_sum <- 0
                
                if(i > 1){
                        dt <- t[i] - t[1:(i-1)]
                        g  <- g_time(dt, p, c)
                        
                        f  <- f_space_hamming_norm(
                                x1 = x[i], loc1_norm = loc_norm[i],
                                x2 = x[1:(i-1)], loc2_norm = loc_norm[1:(i-1)],
                                b_loc = b_loc, D = D_trig
                        )
                        
                        trig_vec <- k_j[1:(i-1)] * g * f
                        trig_sum <- sum(trig_vec)
                }
                
                lam <- bg + trig_sum
                lambda_all[i] <- lam
                
                if(lam > 0){
                        phi_bg[i] <- bg / lam
                        if(i > 1) phi_list[[i]] <- trig_vec / lam
                        else      phi_list[[i]] <- numeric(0)
                } else {
                        phi_bg[i]   <- 0
                        phi_list[[i]] <- if(i>1) rep(0, i-1) else numeric(0)
                }
        }
        
        list(lambda=lambda_all, phi_bg=phi_bg, phi_list=phi_list, s_idx=s_idx)
}

# -------------------------
# M-step: mu を平滑更新（Kは正規化済み）
# -------------------------
update_mu_smooth <- function(phi_bg, s_idx, b_loc = 8, D_mu = 3,
                             state_space = NULL){
        if(is.null(state_space)) state_space <- make_state_space(b_loc)
        S <- nrow(state_space)
        state_codes <- state_space$code
        
        denom <- sum(phi_bg)
        if(denom <= 0) return(rep(1/S, S))  # 退避：一様
        
        mu_new <- rep(0, S)
        for(i in seq_along(phi_bg)){
                center_code <- state_codes[s_idx[i]]
                Kvec <- hamming_kernel_vec_norm(state_codes, center_code, b_loc = b_loc, D = D_mu)
                mu_new <- mu_new + phi_bg[i] * Kvec
        }
        mu_new / denom
}

# -------------------------
# 親ごとの期待子数 w_j = Σ_{i>j} phi_{ij} を作る
# -------------------------
expected_children_per_parent <- function(phi_list){
        n <- length(phi_list)
        w <- rep(0, n)
        if(n <= 1) return(w)
        for(i in 2:n){
                w[1:(i-1)] <- w[1:(i-1)] + phi_list[[i]]
        }
        w
}

# -------------------------
# M-step: A と alpha を更新（phi固定）
#  alpha は1次元方程式を uniroot で解く（ダメなら最小化）
# -------------------------
update_A_alpha <- function(t, m, phi_list, p, c, A_old, alpha_old,
                           T_end = max(t), m_c = 3.0){
        n <- length(t)
        wj <- expected_children_per_parent(phi_list)      # 期待子数（親jごと）
        S1 <- sum(wj)                                     # 期待トリガ総数
        dm <- pmax(0, m - m_c)                             # ここでは m>=m_c 想定
        GJ <- G_time(T_end - t, p, c)                      # 窓内での時間積分
        
        if(S1 <= 0){
                # トリガがほぼゼロなら、Aを小さく保ち、alphaは据え置き
                return(list(A = 0, alpha = alpha_old, wj = wj))
        }
        
        # alpha を解くための補助
        S1m <- sum(wj * dm)
        
        f_alpha <- function(a){
                w  <- exp(a * dm) * GJ
                S2 <- sum(w)
                S2m <- sum(w * dm)
                # A = S1/S2 を代入した条件
                S1m - (S1 / S2) * S2m
        }
        
        # uniroot の範囲を広めに試す
        lo <- -10; hi <- 10
        val_lo <- f_alpha(lo); val_hi <- f_alpha(hi)
        
        if(is.finite(val_lo) && is.finite(val_hi) && val_lo * val_hi < 0){
                alpha_new <- uniroot(f_alpha, lower = lo, upper = hi)$root
        } else {
                # 符号が変わらない等：二乗誤差を最小化
                alpha_new <- optimize(function(a) f_alpha(a)^2, interval = c(lo, hi))$minimum
        }
        
        S2 <- sum(exp(alpha_new * dm) * GJ)
        A_new <- S1 / S2
        
        list(A = A_new, alpha = alpha_new, wj = wj)
}

# -------------------------
# M-step: p,c を数値最適化で更新（phi固定）
#  制約: p>1, c>0 を変換で担保
# -------------------------
update_p_c <- function(t, x, loc_norm, m,
                       phi_list,
                       A, alpha,
                       b_loc = 8,
                       D_trig = 3,
                       T_end = max(t),
                       m_c = 3.0,
                       p_init = 1.2,
                       c_init = 0.01){
        
        n <- length(t)
        k_j <- k_childrenNum(A, alpha, m, m_c = m_c)
        
        # 目的関数：-(Qのトリガ部分)
        # Q = Σ_{i>j} phi_ij log g(dt;p,c) - Σ_j k_j * G(T-t_j;p,c)
        obj <- function(eta){
                p <- 1 + exp(eta[1])
                c <- exp(eta[2])
                
                # 1) Σ phi_ij log g
                term1 <- 0
                for(i in 2:n){
                        w <- phi_list[[i]]
                        if(length(w) == 0) next
                        
                        dt <- t[i] - t[1:(i-1)]
                        # log g(dt) = log(p-1) - log(c) - p*log(1+dt/c)
                        lg <- (log(p - 1) - log(c)) - p * log1p(dt / c)
                        
                        # dt<=0 は理論上無いが、念のため
                        lg[dt <= 0] <- -Inf
                        
                        # wが0のところは無視してOK
                        term1 <- term1 + sum(w * lg, na.rm = TRUE)
                }
                
                # 2) Σ k_j * G(T-t_j)
                term2 <- sum(k_j * G_time(T_end - t, p, c))
                
                # 最大化したいのは (term1 - term2)
                -(term1 - term2)
        }
        
        eta0 <- c(log(p_init - 1), log(c_init))
        fit <- optim(eta0, obj, method = "Nelder-Mead",
                     control = list(maxit = 2000))
        
        p_new <- 1 + exp(fit$par[1])
        c_new <- exp(fit$par[2])
        
        list(p = p_new, c = c_new, optim = fit)
}

# -------------------------
# 全体反復：mu（平滑） + A,alpha,p,c を更新
# -------------------------
fit_etas_mu_Aapc <- function(t, x, loc_norm, m,
                             lambda0 = 2, alpha1 = -0.01,
                             b_loc = 8,
                             D_trig = 3,
                             D_mu   = 3,
                             m_c = 3.0,
                             A_init = 0.3, alpha_init = 1.0,
                             p_init = 1.2, c_init = 0.01,
                             max_iter = 30,
                             tol = 1e-6,
                             verbose = TRUE){
        
        state_space <- make_state_space(b_loc)
        S <- nrow(state_space)
        
        mu <- rep(1/S, S)
        A <- A_init; alpha <- alpha_init
        p <- p_init; c <- c_init
        
        hist <- data.frame(
                iter = integer(0),
                loglik = numeric(0),
                Q = numeric(0),
                A = numeric(0),
                alpha = numeric(0),
                p = numeric(0),
                c = numeric(0),
                diff = numeric(0)
        )
        
        for(it in seq_len(max_iter)){
                
                # --- E-step（現パラメータで責任度など）---
                E <- etas_E_step_storePhi(
                        t, x, loc_norm, m,
                        mu_vec = mu,
                        lambda0 = lambda0, alpha1 = alpha1,
                        A = A, alpha = alpha,
                        p = p, c = c,
                        b_loc = b_loc,
                        D_trig = D_trig,
                        m_c = m_c
                )
                
                # --- 収束モニタ用：観測対数尤度とQ（いずれも時空間部分）---
                ll_now <- loglik_obs_ts(t, m, E, lambda0, alpha1, A, alpha, p, c, m_c)
                Q_now  <- Q_ts(t, x, loc_norm, m, mu, E, lambda0, alpha1, A, alpha, p, c,
                               b_loc = b_loc, D_trig = D_trig, m_c = m_c)
                
                # --- M-step(1): mu 平滑更新 ---
                mu_new <- update_mu_smooth(E$phi_bg, E$s_idx, b_loc = b_loc, D_mu = D_mu,
                                           state_space = state_space)
                
                # --- M-step(2): p,c 更新（数値最適化） ---
                pc <- update_p_c(t, x, loc_norm, m,
                                 phi_list = E$phi_list,
                                 A = A, alpha = alpha,
                                 b_loc = b_loc,
                                 D_trig = D_trig,
                                 T_end = max(t),
                                 m_c = m_c,
                                 p_init = p, c_init = c)
                p_new <- pc$p
                c_new <- pc$c
                
                # --- M-step(3): A,alpha 更新（閉形式＋1次元） ---
                Aa <- update_A_alpha(t, m,
                                     phi_list = E$phi_list,
                                     p = p_new, c = c_new,
                                     A_old = A, alpha_old = alpha,
                                     T_end = max(t),
                                     m_c = m_c)
                A_new <- Aa$A
                alpha_new <- Aa$alpha
                
                # パラメータ変化量（収束用）
                diff <- sum(abs(mu_new - mu)) +
                        abs(A_new - A) + abs(alpha_new - alpha) +
                        abs(p_new - p) + abs(c_new - c)
                
                # 履歴保存（※この行の A,p... は “E-step時点” の値＝ ll_now/Q_now と対応）
                hist <- rbind(hist, data.frame(
                        iter = it, loglik = ll_now, Q = Q_now,
                        A = A, alpha = alpha, p = p, c = c,
                        diff = diff
                ))
                
                if(verbose){
                        cat(sprintf("iter=%d  loglik=%.3f  Q=%.3f  diff=%.3e  (A=%.3f alpha=%.3f p=%.3f c=%.4f)\n",
                                    it, ll_now, Q_now, diff, A, alpha, p, c))
                }
                
                # 更新
                mu <- mu_new
                A <- A_new; alpha <- alpha_new
                p <- p_new; c <- c_new
                
                if(diff < tol) break
        }
        
        list(mu = mu, state_space = state_space,
             A = A, alpha = alpha, p = p, c = c,
             last_E = E,
             history = hist)
}
