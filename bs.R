#RunEstimate後
parametric_bootstrap_CI <- function(
                dat0,                       # 元データ（data.frame: time_t, moduleID_x, LOCnorm_y, difficulty_m）
                theta_true = NULL,          # 真値（list）※NULLなら包含判定しない
                T_end = 100,
                b_loc = 8,
                D_mu = 3,
                B = 200,
                seed = 1,
                lambda0 = 2,
                alpha1  = -0.01,
                # ブート生成に使う固定値（必要なら引数で変える）
                beta_m  = 1.2,
                Mc = 3.0,
                Mmax = 6.0,
                dM = 0.1,
                m_set = c(3.0, 4.5, 6.0), # 最終時刻で断面出すm
                D_trig  = 3,
                # 初期値（NULLなら最初は適当、ブートは前回推定値を初期値に使う）
                init = list(A=0.3, alpha=1.0, p=1.2, c=0.01),
                max_iter = 30,
                tol = 1e-6,
                verbose = TRUE,
                # 保存
                out_dir = NULL
){
        
        set.seed(seed)
        
        # ---- 0) 入力整形 ----
        dat0 <- dat0[order(dat0$time_t), ]
        t0   <- dat0$time_t
        x0   <- dat0$moduleID_x
        loc0 <- dat0$LOCnorm_y
        m0   <- dat0$difficulty_m
        
        # ---- 1) 元データのフィット（res0）----
        if(verbose) cat("Fitting on original data...\n")
        res0 <- fit_etas_mu_Aapc(
                t = t0, x = x0, loc_norm = loc0, m = m0,
                lambda0 = lambda0, alpha1 = alpha1,
                b_loc = b_loc, D_trig = D_trig, D_mu = D_mu,
                A_init = init$A, alpha_init = init$alpha,
                p_init = init$p, c_init = init$c,
                max_iter = max_iter, tol = tol,
                verbose = verbose
        )
        
        # 返り値構造の差を吸収（res0$mu or res0$res$mu どちらでもOK）
        get_mu <- function(obj){
                if(!is.null(obj$mu)) return(obj$mu)
                if(!is.null(obj$res) && !is.null(obj$res$mu)) return(obj$res$mu)
                stop("mu not found in fit result.")
        }
        get_hist <- function(obj){
                if(!is.null(obj$history)) return(obj$history)
                if(!is.null(obj$res) && !is.null(obj$res$history)) return(obj$res$history)
                return(NULL)
        }
        
        mu_hat <- get_mu(res0)
        
        # ブート生成に使う“推定真値”
        theta_hat <- list(
                lambda0 = lambda0, alpha1 = alpha1,
                A = res0$A, alpha = res0$alpha,
                p = res0$p, c = res0$c,
                D_trig = D_trig, Mc=Mc
        )
        mu_hat <- res0$mu
        
        if(verbose){
                cat(sprintf("theta_hat: A=%.4f alpha=%.4f p=%.4f c=%.6f\n",
                            theta_hat$A, theta_hat$alpha, theta_hat$p, theta_hat$c))
        }
        
        ## ---- λz(T,z) for original data (point estimate) ----
        dat0_code <- dat0
        lambdaZ_hat_T <- lambda_z_at_time_allz(
                t_eval = T_end, dat = dat0_code, theta = theta_hat, mu_vec = mu_hat,
                b_loc = b_loc, D_trig = D_trig, Mc = Mc
        )
        
        ## ---- s(m) fixed at beta_true ----
        s_set <- s_gr_trunc_pmf(m_set, beta = beta_m, Mc = Mc, Mmax = Mmax, dM = dM)
        
        
        # ---- 2) ブートループ ----
        boot_mat <- matrix(NA_real_, nrow = B, ncol = 4)
        colnames(boot_mat) <- c("A","alpha","p","c")
        
        n_events <- integer(B)
        conv_it  <- integer(B)  # 収束した反復回数（historyの最終iter）
        failed   <- logical(B)
        
        lambdaZ_boot <- matrix(NA_real_, nrow = B, ncol = length(mu_hat))
        
        # ブートでの初期値は「直前の推定値」を使うと安定
        A_init_b     <- theta_hat$A
        alpha_init_b <- theta_hat$alpha
        p_init_b     <- theta_hat$p
        c_init_b     <- theta_hat$c
        
        for(b in 1:B){
                if(verbose && (b %% max(1, floor(B/10)) == 0)) cat("bootstrap", b, "/", B, "\n")
                
                sim <- simulate_etas_modified(
                        T_end = T_end,
                        lambda0 = theta_hat$lambda0, alpha1 = theta_hat$alpha1,
                        mu_vec = mu_hat, b_loc = b_loc,
                        beta_m = beta_m,
                        A = theta_hat$A, alpha = theta_hat$alpha,
                        p = theta_hat$p, c = theta_hat$c,
                        D_trig = theta_hat$D_trig,
                        max_events = 50000
                )
                
                n_events[b] <- nrow(sim)
                
                # フィット
                fitb <- try(
                        fit_etas_mu_Aapc(
                                t = sim$time_t, x = sim$moduleID_x, loc_norm = sim$LOCnorm_y, m = sim$difficulty_m,
                                lambda0 = theta_hat$lambda0, alpha1 = theta_hat$alpha1,
                                b_loc = b_loc, D_trig = theta_hat$D_trig, D_mu = D_mu,
                                A_init = A_init_b, alpha_init = alpha_init_b,
                                p_init = p_init_b, c_init = c_init_b,
                                max_iter = max_iter, tol = tol,
                                verbose = FALSE
                        ),
                        silent = TRUE
                )
                
                if(inherits(fitb, "try-error")){
                        failed[b] <- TRUE
                        next
                }
                
                boot_mat[b, ] <- c(fitb$A, fitb$alpha, fitb$p, fitb$c)
                
                # 次回初期値を更新（連鎖的に安定させる）
                A_init_b     <- fitb$A
                alpha_init_b <- fitb$alpha
                p_init_b     <- fitb$p
                c_init_b     <- fitb$c
                
                hb <- get_hist(fitb)
                if(!is.null(hb)) conv_it[b] <- tail(hb$iter, 1)
                
                # λz(T,z) を保存
                mu_b <- get_mu(fitb)
                theta_b <- list(
                        lambda0 = lambda0, alpha1 = alpha1,
                        A = fitb$A, alpha = fitb$alpha,
                        p = fitb$p, c = fitb$c,
                        D_trig = D_trig,
                        Mc = Mc
                )
                lambdaZ_boot[b, ] <- lambda_z_at_time_allz(
                        t_eval = T_end, dat = sim, theta = theta_b, mu_vec = mu_b,
                        b_loc = b_loc, D_trig = D_trig, Mc = Mc
                )
                
                # 収束回数（historyがある前提）
                if(!is.null(fitb$history)){
                        conv_it[b] <- tail(fitb$history$iter, 1)
                }
        }
        
        boot_df <- as.data.frame(boot_mat)
        boot_df$N <- n_events
        boot_df$failed <- failed
        boot_df$conv_iter <- conv_it
        
        ok <- !failed & is.finite(boot_df$A) & is.finite(boot_df$alpha) & is.finite(boot_df$p) & is.finite(boot_df$c) & apply(lambdaZ_boot, 1, function(r) all(is.finite(r)))
        # ---- 3) CI（percentile）----
        CI <- apply(boot_df[ok, c("A","alpha","p","c")], 2,
                    quantile, probs = c(0.025, 0.975), na.rm = TRUE)
        
        # ---- 4) 真値包含判定（任意）----
        cover <- NULL
        if(!is.null(theta_true)){
                cover <- c(
                        A     = (CI[1,"A"]     <= theta_true$A     && theta_true$A     <= CI[2,"A"]),
                        alpha = (CI[1,"alpha"] <= theta_true$alpha && theta_true$alpha <= CI[2,"alpha"]),
                        p     = (CI[1,"p"]     <= theta_true$p     && theta_true$p     <= CI[2,"p"]),
                        c     = (CI[1,"c"]     <= theta_true$c     && theta_true$c     <= CI[2,"c"])
                )
        }
        
        ## ---- 5) λz(T,z) CI ----
        Lz_lo <- apply(lambdaZ_boot[ok, , drop=FALSE], 2, quantile, probs = 0.025, na.rm = TRUE)
        Lz_hi <- apply(lambdaZ_boot[ok, , drop=FALSE], 2, quantile, probs = 0.975, na.rm = TRUE)
        
        lambdaZ_CI_T <- list(lo = Lz_lo, hi = Lz_hi)
        
        ## ---- 6) mark-slice at final time: λ(T,z,m)=s(m)*λz(T,z)
        # point estimate
        lambdaZM_hat_T <- lapply(seq_along(m_set), function(k) s_set[k] * lambdaZ_hat_T)
        names(lambdaZM_hat_T) <- paste0("m=", m_set)
        
        # CI (scale)
        lambdaZM_CI_T <- lapply(seq_along(m_set), function(k){
                list(lo = s_set[k] * Lz_lo,
                     hi = s_set[k] * Lz_hi)
        })
        names(lambdaZM_CI_T) <- paste0("m=", m_set)
        
        # ---- 7) 保存（任意）----
        if(!is.null(out_dir)){
                dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
                saveRDS(list(res0=res0, theta_hat=theta_hat, CI=CI, cover=cover, boot=boot_df,ok_rate = mean(ok),
                             lambdaZ_hat_T = lambdaZ_hat_T, lambdaZ_CI_T  = lambdaZ_CI_T, m_set = m_set, s_set = s_set,
                             lambdaZM_hat_T = lambdaZM_hat_T,lambdaZM_CI_T  = lambdaZM_CI_T),
                        file = file.path(out_dir, "bootstrap_result_with_lambda.rds"))
                write.csv(res0$history, file.path(out_dir, "fit_history.csv"), row.names = FALSE)
                write.csv(cbind(res0$state_space, mu=res0$mu), file.path(out_dir, "mu_hat.csv"), row.names = FALSE)
                write.csv(boot_df, file.path(out_dir, "boot_theta.csv"), row.names = FALSE)
                # 元データも残す
                write.csv(dat0, file.path(out_dir, "data_original.csv"), row.names = FALSE)
                
                # optional: λz CI export
                write.csv(data.frame(code=0:(length(mu_hat)-1),
                                     lambdaZ_hat_T = lambdaZ_hat_T,
                                     lambdaZ_lo_T = Lz_lo,
                                     lambdaZ_hi_T = Lz_hi),
                          file.path(out_dir, "lambdaZ_T_CI.csv"), row.names = FALSE)
                
                # optional: mark slices export (m_set only)
                for(k in seq_along(m_set)){
                        nm <- paste0("m_", gsub("\\.","_", as.character(m_set[k])))
                        write.csv(data.frame(code=0:(length(mu_hat)-1),
                                             lambdaZM_hat_T = lambdaZM_hat_T[[k]],
                                             lambdaZM_lo_T  = lambdaZM_CI_T[[k]]$lo,
                                             lambdaZM_hi_T  = lambdaZM_CI_T[[k]]$hi),
                                  file.path(out_dir, paste0("lambdaZM_T_CI_", nm, ".csv")),
                                  row.names = FALSE)
                }
        }
        
        list(
                res0 = res0,
                theta_hat = theta_hat,
                CI = CI,
                cover = cover,
                boot = boot_df,
                ok_rate = mean(ok),
                
                lambdaZ_hat_T = lambdaZ_hat_T,
                lambdaZ_CI_T  = lambdaZ_CI_T,
                m_set = m_set,
                s_set = s_set,
                lambdaZM_hat_T = lambdaZM_hat_T,
                lambdaZM_CI_T  = lambdaZM_CI_T
        )
}

# ブートストラップCI
bs <- parametric_bootstrap_CI(
        dat0 = dat,
        theta_true = theta_true,
        T_end = 100, b_loc = 8,
        lambda0 = theta_true$lambda0,
        alpha1  = theta_true$alpha1,
        beta_m  = theta_true$beta_m,
        m_set = c(3.0, 4.5, 6.0),
        D_trig  = theta_true$D_trig,
        B = 1000,
        out_dir = "output/boot_with_lambda", #旧BS_out
        verbose = TRUE
)

bs$CI
bs$cover
bs$ok_rate
