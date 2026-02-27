bs <- readRDS("output/boot_with_lambda/bootstrap_result_with_lambda.rds")
#3D格子
parametric_bootstrap_CI <- function(
                dat0,
                theta_true = NULL,
                T_end = 100,
                b_loc = 8,
                D_mu = 3,
                B = 200,
                seed = 1,
                lambda0 = 2,
                alpha1  = -0.01,
                
                beta_m  = 1.2,   # s(m)は固定
                Mc = 3.0,
                Mmax = 6.0,
                dM = 0.1,
                m_set = c(3.0, 4.5, 6.0),
                
                D_trig  = 3,
                
                init = list(A=0.3, alpha=1.0, p=1.2, c=0.01),
                max_iter = 30,
                tol = 1e-6,
                verbose = TRUE,
                
                # ---- 3D (t,m,LOCnorm) 用の格子出力 ----
                do_grid = TRUE,
                grid_modules = 0:3,
                t_grid  = seq(0, 100, length.out = 40),
                loc_grid = seq(0, 1, length.out = 64),
                m_grid  = seq(3.0, 6.0, by = 0.1),
                
                out_dir = NULL
){
        
        set.seed(seed)
        
        ## ---- 0) 入力整形 ----
        dat0 <- dat0[order(dat0$time_t), ]
        t0   <- dat0$time_t
        x0   <- dat0$moduleID_x
        loc0 <- dat0$LOCnorm_y
        m0   <- dat0$difficulty_m
        
        ## ---- 1) 元データのフィット（res0）----
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
        
        # 返り値構造の差を吸収
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
        get_state_space <- function(obj){
                if(!is.null(obj$state_space)) return(obj$state_space)
                if(!is.null(obj$res) && !is.null(obj$res$state_space)) return(obj$res$state_space)
                return(NULL)
        }
        
        mu_hat <- get_mu(res0)
        
        theta_hat <- list(
                lambda0 = lambda0, alpha1 = alpha1,
                A = res0$A, alpha = res0$alpha,
                p = res0$p, c = res0$c,
                D_trig = D_trig, Mc = Mc
        )
        
        if(verbose){
                cat(sprintf("theta_hat: A=%.4f alpha=%.4f p=%.4f c=%.6f\n",
                            theta_hat$A, theta_hat$alpha, theta_hat$p, theta_hat$c))
        }
        
        ## ---- s(m) fixed ----
        s_set  <- s_gr_trunc_pmf(m_set,  beta = beta_m, Mc = Mc, Mmax = Mmax, dM = dM)
        s_grid <- s_gr_trunc_pmf(m_grid, beta = beta_m, Mc = Mc, Mmax = Mmax, dM = dM)
        
        ## ---- λz(T,z) for original data ----
        lambdaZ_hat_T <- lambda_z_at_time_allz(
                t_eval = T_end, dat = dat0, theta = theta_hat, mu_vec = mu_hat,
                b_loc = b_loc, D_trig = D_trig, Mc = Mc
        )
        
        ## ---- (option) grid point estimate on original data: λz(t,loc) by module ----
        lambdaZ_hat_grid <- NULL
        if(do_grid){
                lambdaZ_hat_grid <- lapply(grid_modules, function(mod){
                        calc_lambdaZ_grid(dat0, theta_hat, mu_hat, module_id = mod,
                                          t_grid = t_grid, loc_grid = loc_grid,
                                          b_loc = b_loc, D_trig = D_trig, Mc = Mc)
                })
                names(lambdaZ_hat_grid) <- paste0("module=", grid_modules)
        }
        
        ## ---- 2) bootstrap loop ----
        boot_mat <- matrix(NA_real_, nrow = B, ncol = 4)
        colnames(boot_mat) <- c("A","alpha","p","c")
        
        n_events <- integer(B)
        conv_it  <- integer(B)
        failed   <- logical(B)
        
        lambdaZ_boot_T <- matrix(NA_real_, nrow = B, ncol = length(mu_hat))
        
        # grid保存：moduleごとに (B × nt × nloc)
        lambdaZ_boot_grid <- NULL
        if(do_grid){
                lambdaZ_boot_grid <- lapply(grid_modules, function(mod){
                        array(NA_real_, dim = c(B, length(t_grid), length(loc_grid)))
                })
                names(lambdaZ_boot_grid) <- paste0("module=", grid_modules)
        }
        
        # 初期値（連鎖更新）
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
                
                # 次回初期値を更新（安定化）
                A_init_b     <- fitb$A
                alpha_init_b <- fitb$alpha
                p_init_b     <- fitb$p
                c_init_b     <- fitb$c
                
                hb <- get_hist(fitb)
                if(!is.null(hb) && !is.null(hb$iter)) conv_it[b] <- tail(hb$iter, 1)
                
                # λz(T,z) 保存
                mu_b <- get_mu(fitb)
                theta_b <- list(
                        lambda0 = lambda0, alpha1 = alpha1,
                        A = fitb$A, alpha = fitb$alpha,
                        p = fitb$p, c = fitb$c,
                        D_trig = D_trig, Mc = Mc
                )
                
                lambdaZ_boot_T[b, ] <- lambda_z_at_time_allz(
                        t_eval = T_end, dat = sim, theta = theta_b, mu_vec = mu_b,
                        b_loc = b_loc, D_trig = D_trig, Mc = Mc
                )
                
                # (option) gridも保存
                if(do_grid){
                        for(ii in seq_along(grid_modules)){
                                mod <- grid_modules[ii]
                                lambdaZ_boot_grid[[ii]][b,,] <- calc_lambdaZ_grid(
                                        dat = sim, theta = theta_b, mu_vec = mu_b, module_id = mod,
                                        t_grid = t_grid, loc_grid = loc_grid,
                                        b_loc = b_loc, D_trig = D_trig, Mc = Mc
                                )
                        }
                }
        }
        
        boot_df <- as.data.frame(boot_mat)
        boot_df$N <- n_events
        boot_df$failed <- failed
        boot_df$conv_iter <- conv_it
        
        ok <- !failed &
                is.finite(boot_df$A) & is.finite(boot_df$alpha) & is.finite(boot_df$p) & is.finite(boot_df$c) &
                apply(lambdaZ_boot_T, 1, function(r) all(is.finite(r)))
        
        ## ---- 3) parameter CI ----
        CI <- apply(boot_df[ok, c("A","alpha","p","c")], 2,
                    quantile, probs = c(0.025, 0.975), na.rm = TRUE)
        
        ## ---- 4) cover check ----
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
        Lz_lo_T <- apply(lambdaZ_boot_T[ok, , drop=FALSE], 2, quantile, probs = 0.025, na.rm = TRUE)
        Lz_hi_T <- apply(lambdaZ_boot_T[ok, , drop=FALSE], 2, quantile, probs = 0.975, na.rm = TRUE)
        lambdaZ_CI_T <- list(lo = Lz_lo_T, hi = Lz_hi_T)
        
        ## ---- 6) mark-slice at final time: λ(T,z,m)=s(m)*λz(T,z) ----
        lambdaZM_hat_T <- lapply(seq_along(m_set), function(k) s_set[k] * lambdaZ_hat_T)
        names(lambdaZM_hat_T) <- paste0("m=", m_set)
        
        lambdaZM_CI_T <- lapply(seq_along(m_set), function(k){
                list(lo = s_set[k] * Lz_lo_T,
                     hi = s_set[k] * Lz_hi_T)
        })
        names(lambdaZM_CI_T) <- paste0("m=", m_set)
        
        ## ---- 7) (option) grid CI: moduleごとに λz(t,loc) のCI ----
        lambdaZ_CI_grid <- NULL
        if(do_grid){
                ok_grid <- ok
                # NA/Inf混入をさらに除外
                for(ii in seq_along(grid_modules)){
                        arr <- lambdaZ_boot_grid[[ii]]
                        ok_grid <- ok_grid & apply(arr, 1, function(a) all(is.finite(a)))
                }
                
                lambdaZ_CI_grid <- lapply(seq_along(grid_modules), function(ii){
                        arr <- lambdaZ_boot_grid[[ii]][ok_grid,,,drop=FALSE]  # reps × nt × nloc
                        lo <- apply(arr, c(2,3), quantile, probs=0.025, na.rm=TRUE)
                        hi <- apply(arr, c(2,3), quantile, probs=0.975, na.rm=TRUE)
                        list(lo=lo, hi=hi, ok_rate_grid = mean(ok_grid))
                })
                names(lambdaZ_CI_grid) <- paste0("module=", grid_modules)
        }
        
        ## ---- 8) (option) grid → λ(t,m,LOCnorm) のCIを作るための“材料”を返す ----
        # ここでは巨大な 4次元( t×loc×m )配列は作らず、
        # 使う時に s_grid を掛けて作れるようにして返す（軽量）
        grid_bundle <- NULL
        if(do_grid){
                grid_bundle <- list(
                        t_grid = t_grid,
                        loc_grid = loc_grid,
                        m_grid = m_grid,
                        s_grid = s_grid,
                        modules = grid_modules,
                        lambdaZ_hat_grid = lambdaZ_hat_grid,   # module -> (nt×nloc)
                        lambdaZ_CI_grid  = lambdaZ_CI_grid     # module -> list(lo(nt×nloc), hi(nt×nloc))
                )
        }
        
        ## ---- 9) 保存 ----
        if(!is.null(out_dir)){
                dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
                
                saveRDS(
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
                                lambdaZM_CI_T  = lambdaZM_CI_T,
                                
                                grid = grid_bundle
                        ),
                        file = file.path(out_dir, "bootstrap_result_with_lambda.rds")
                )
                
                h0 <- get_hist(res0)
                if(!is.null(h0)) write.csv(h0, file.path(out_dir, "fit_history.csv"), row.names = FALSE)
                
                ss <- get_state_space(res0)
                if(!is.null(ss)){
                        write.csv(cbind(ss, mu = mu_hat), file.path(out_dir, "mu_hat.csv"), row.names = FALSE)
                } else {
                        write.csv(data.frame(code = 0:(length(mu_hat)-1), mu = mu_hat),
                                  file.path(out_dir, "mu_hat.csv"), row.names = FALSE)
                }
                
                write.csv(boot_df, file.path(out_dir, "boot_theta.csv"), row.names = FALSE)
                write.csv(dat0, file.path(out_dir, "data_original.csv"), row.names = FALSE)
                
                # λz(T,z) export
                write.csv(data.frame(code=0:(length(mu_hat)-1),
                                     lambdaZ_hat_T = lambdaZ_hat_T,
                                     lambdaZ_lo_T  = Lz_lo_T,
                                     lambdaZ_hi_T  = Lz_hi_T),
                          file.path(out_dir, "lambdaZ_T_CI.csv"), row.names = FALSE)
                
                # λ(T,z,m_set) export
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
                lambdaZM_CI_T  = lambdaZM_CI_T,
                
                grid = grid_bundle
        )
}

## =========================================================
## 3Dプロット：moduleごとに x=t, y=m, z=LOCnorm（色=λ, hover=CI）
## =========================================================
plot_lambdaZM_3d <- function(bs, module_id = 0, which = c("hat","lo","hi"),
                             use_CI = TRUE, log_scale = TRUE){
        which <- match.arg(which)
        if(is.null(bs$grid)) stop("bs$grid が NULL です。parametric_bootstrap_CI(do_grid=TRUE) で実行してください。")
        if(!requireNamespace("plotly", quietly = TRUE)) stop("plotly を install.packages('plotly') してください。")
        
        g <- bs$grid
        key <- paste0("module=", module_id)
        if(is.null(g$lambdaZ_hat_grid[[key]])) stop("その module_id の格子がありません。grid_modules を確認してください。")
        
        t_grid  <- g$t_grid
        loc_grid <- g$loc_grid
        m_grid  <- g$m_grid
        s_grid  <- g$s_grid
        
        Z_hat <- g$lambdaZ_hat_grid[[key]]
        
        if(use_CI){
                CIg <- g$lambdaZ_CI_grid[[key]]
                Z_lo <- CIg$lo
                Z_hi <- CIg$hi
        } else {
                Z_lo <- Z_hat
                Z_hi <- Z_hat
        }
        
        # λ(t,loc,m) = s(m) * λz(t,loc)
        # 4Dを作らず df で展開（軽量）
        df <- expand.grid(t = t_grid, loc = loc_grid, m = m_grid)
        nt <- length(t_grid); nl <- length(loc_grid); nm <- length(m_grid)
        
        idx_t  <- rep(seq_len(nt), times = nl * nm)
        idx_loc <- rep(rep(seq_len(nl), each = nt), times = nm)
        idx_m  <- rep(seq_len(nm), each = nt * nl)
        
        lam_hat <- s_grid[idx_m] * Z_hat[cbind(idx_t, idx_loc)]
        lam_lo  <- s_grid[idx_m] * Z_lo [cbind(idx_t, idx_loc)]
        lam_hi  <- s_grid[idx_m] * Z_hi [cbind(idx_t, idx_loc)]
        
        df$hat <- lam_hat
        df$lo  <- lam_lo
        df$hi  <- lam_hi
        
        val <- switch(which, hat=df$hat, lo=df$lo, hi=df$hi)
        colval <- if(log_scale) log10(df$hat + 1e-12) else df$hat
        
        plotly::plot_ly(
                df,
                x = ~t, y = ~m, z = ~loc,
                type = "scatter3d", mode = "markers",
                marker = list(size = 2),
                color = colval,
                text = ~sprintf("module=%d<br>t=%.3g loc=%.3g m=%.1f<br>hat=%.3g<br>CI=[%.3g, %.3g]",
                                module_id, t, loc, m, hat, lo, hi),
                hoverinfo = "text"
        ) |>
                plotly::layout(
                        title = paste0("module ", module_id, " : λ(t,z,m)  (color = ",
                                       if(log_scale) "log10(hat)" else "hat", ")"),
                        scene = list(
                                xaxis = list(title = "t"),
                                yaxis = list(title = "m"),
                                zaxis = list(title = "LOCnorm")
                        )
                )
}

plot_lambdaZM_band_surface <- function(
                bs,
                module_id = 0,
                m_slices = c(3.0, 4.5, 6.0),
                use_log10 = TRUE,
                z_as = c("loc", "m"),     # z軸に LOCnorm を置くか、m を置くか
                show_hat = TRUE,
                opacity_band = 0.35,
                opacity_hat  = 0.85
){
        z_as <- match.arg(z_as)
        
        if(is.null(bs$grid)) stop("bs$grid が NULL。parametric_bootstrap_CI(do_grid=TRUE) で実行してください。")
        if(!requireNamespace("plotly", quietly = TRUE)) stop("plotly を install.packages('plotly') してください。")
        
        g <- bs$grid
        key <- paste0("module=", module_id)
        if(is.null(g$lambdaZ_hat_grid[[key]])) stop("その module_id の格子がありません。grid_modules を確認してください。")
        
        # 格子
        t_grid   <- g$t_grid
        loc_grid <- g$loc_grid
        m_grid   <- g$m_grid
        s_grid   <- g$s_grid
        
        # m_slices を m_grid に合わせて最近傍に丸める
        nearest_m <- function(mv){
                m_grid[which.min(abs(m_grid - mv))]
        }
        m_slices2 <- vapply(m_slices, nearest_m, numeric(1))
        m_slices2 <- unique(m_slices2)
        
        # base: λz(t,loc)
        Z_hat <- g$lambdaZ_hat_grid[[key]]
        CIg   <- g$lambdaZ_CI_grid[[key]]
        Z_lo  <- CIg$lo
        Z_hi  <- CIg$hi
        
        # plotly の surface は z が matrix(nrow=length(x), ncol=length(y)) を期待
        # 今 Z_* は nt×nloc なので、x=t, y=loc にそのまま渡せる
        x <- t_grid
        y <- loc_grid
        
        # 値変換
        transf <- function(M){
                if(!use_log10) return(M)
                log10(pmax(M, 1e-12))
        }
        
        p <- plotly::plot_ly()
        
        for(ms in m_slices2){
                # s(m)
                s_m <- s_grid[which(m_grid == ms)]
                if(length(s_m) != 1) next
                
                Lhat <- s_m * Z_hat
                Llo  <- s_m * Z_lo
                Lhi  <- s_m * Z_hi
                
                # 表示を log10 にする場合
                Lhat2 <- transf(Lhat)
                Llo2  <- transf(Llo)
                Lhi2  <- transf(Lhi)
                
                # 3枚のサーフェスをどう配置するか：
                # - z_as="m": z軸に m を置く（面が m 方向に並ぶ）、color で λ
                # - z_as="loc": z軸に LOCnorm を置く（通常の t×loc 面）、color で λ、m はタイトル表示
                #
                # 注意：surface は z が“高さ”で、color は surfacecolor
                # ここでは「帯」を見せたいので、高さ=z軸は (t,loc) 平面を保ちつつ、
                #   surfacecolor に λ、そして lo/hi/hat を重ねる
                # とすると重なりが見やすい。
                #
                # ただし user の希望「CIの上下を2枚の3Dサーフェスで挟む」に忠実にするなら、
                # 高さ=z に λ（lo/hat/hi）を置く（⇒本当に上下に挟める）。
                # その場合、軸は x=t, y=loc, z=λ になる。
                #
                # ここでは “挟む” を優先して z=λ（変換後）で描く。
                
                z_hat <- Lhat2
                z_lo  <- Llo2
                z_hi  <- Lhi2
                
                # --- lo surface ---
                p <- p |>
                        plotly::add_surface(
                                x = x, y = y, z = z_lo,
                                opacity = opacity_band,
                                showscale = FALSE,
                                name = paste0("lo m=", ms),
                                hovertemplate = paste0(
                                        "module=", module_id,
                                        "<br>m=", ms,
                                        "<br>t=%{x:.3g}",
                                        "<br>LOC=%{y:.3g}",
                                        "<br>lo(z)=%{z:.3g}<extra></extra>"
                                )
                        )
                
                # --- hi surface ---
                p <- p |>
                        plotly::add_surface(
                                x = x, y = y, z = z_hi,
                                opacity = opacity_band,
                                showscale = FALSE,
                                name = paste0("hi m=", ms),
                                hovertemplate = paste0(
                                        "module=", module_id,
                                        "<br>m=", ms,
                                        "<br>t=%{x:.3g}",
                                        "<br>LOC=%{y:.3g}",
                                        "<br>hi(z)=%{z:.3g}<extra></extra>"
                                )
                        )
                
                # --- hat surface（任意） ---
                if(show_hat){
                        p <- p |>
                                plotly::add_surface(
                                        x = x, y = y, z = z_hat,
                                        opacity = opacity_hat,
                                        showscale = FALSE,
                                        name = paste0("hat m=", ms),
                                        hovertemplate = paste0(
                                                "module=", module_id,
                                                "<br>m=", ms,
                                                "<br>t=%{x:.3g}",
                                                "<br>LOC=%{y:.3g}",
                                                "<br>hat(z)=%{z:.3g}<extra></extra>"
                                        )
                                )
                }
        }
        
        ztitle <- if(use_log10) "log10(λ(t,z,m))" else "λ(t,z,m)"
        
        p |>
                plotly::layout(
                        title = paste0(
                                "CI band surfaces (module ", module_id, "), m slices: ",
                                paste(m_slices2, collapse=", "),
                                if(use_log10) "  [z=log10 λ]" else "  [z=λ]"
                        ),
                        scene = list(
                                xaxis = list(title = "t"),
                                yaxis = list(title = "LOCnorm"),
                                zaxis = list(title = ztitle)
                        )
                )
}

## =========================================================
## 共通：格子データ取得（LOC切り取りはここで実行）
## =========================================================
.get_lambdaZ_grid_for_plot <- function(bs, module_id, m_slice,
                                       loc_range = c(0,1),
                                       use_log10 = TRUE){
        if(is.null(bs$grid)) stop("bs$grid が NULL。parametric_bootstrap_CI(do_grid=TRUE) で実行してください。")
        g <- bs$grid
        key <- paste0("module=", module_id)
        if(is.null(g$lambdaZ_hat_grid[[key]])) stop("その module_id の格子がありません。grid_modules を確認してください。")
        CIg <- g$lambdaZ_CI_grid[[key]]
        if(is.null(CIg$lo) || is.null(CIg$hi)) stop("lambdaZ_CI_grid が見つかりません。do_grid=TRUE を確認。")
        
        # m を grid に合わせて最近傍
        m_grid <- g$m_grid
        s_grid <- g$s_grid
        ms <- m_grid[which.min(abs(m_grid - m_slice))]
        s_m <- s_grid[which(m_grid == ms)]
        if(length(s_m) != 1) stop("m_slice に対応する s(m) が取れません。m_grid を確認してください。")
        
        # LOCnorm 範囲でフィルタ（列を切り取り）
        loc_grid <- g$loc_grid
        loc_min <- loc_range[1]; loc_max <- loc_range[2]
        if(loc_min > loc_max) { tmp <- loc_min; loc_min <- loc_max; loc_max <- tmp }
        keep <- which(loc_grid >= loc_min & loc_grid <= loc_max)
        if(length(keep) < 2) stop("loc_range が狭すぎます（点が2つ未満）。")
        
        t_grid <- g$t_grid
        Zhat <- g$lambdaZ_hat_grid[[key]][, keep, drop=FALSE]   # nt×nloc_keep
        Zlo  <- CIg$lo[, keep, drop=FALSE]
        Zhi  <- CIg$hi[, keep, drop=FALSE]
        
        # λ(t,z,m)=s(m)*λz(t,z)
        Lhat <- s_m * Zhat
        Llo  <- s_m * Zlo
        Lhi  <- s_m * Zhi
        
        transf <- function(M){
                if(!use_log10) return(M)
                log10(pmax(M, 1e-12))
        }
        
        list(
                module_id = module_id,
                m_used = ms,
                s_m = s_m,
                t = t_grid,
                loc_full = loc_grid,        # 0..1 の全格子（軸固定に使う）
                loc = loc_grid[keep],       # 切り取った格子
                loc_range = c(loc_min, loc_max),
                hat = transf(Lhat),
                lo  = transf(Llo),
                hi  = transf(Lhi),
                use_log10 = use_log10
        )
}

## =========================================================
## (1) 1パネル：lo/hat/hi を 3D surface で表示（m固定）
## =========================================================
plot_band_surface_panel <- function(bs, module_id=0, m_slice=3.0,
                                    loc_range=c(0,1),
                                    use_log10=TRUE,
                                    show_hat=TRUE,
                                    opacity_band=0.35,
                                    opacity_hat=0.85,
                                    showscale=FALSE){
        if(!requireNamespace("plotly", quietly=TRUE)) stop("plotly を install.packages('plotly') してください。")
        
        d <- .get_lambdaZ_grid_for_plot(bs, module_id, m_slice, loc_range, use_log10)
        ztitle <- if(d$use_log10) "log10(λ(t,z,m))" else "λ(t,z,m)"
        
        p <- plotly::plot_ly() |>
                plotly::add_surface(
                        x=d$t, y=d$loc, z=d$lo,
                        opacity=opacity_band, showscale=FALSE, name=paste0("lo m=", d$m_used),
                        hovertemplate=paste0("module=",module_id,"<br>m=",d$m_used,
                                             "<br>t=%{x:.3g}<br>LOC=%{y:.3g}<br>lo=%{z:.3g}<extra></extra>")
                ) |>
                plotly::add_surface(
                        x=d$t, y=d$loc, z=d$hi,
                        opacity=opacity_band, showscale=FALSE, name=paste0("hi m=", d$m_used),
                        hovertemplate=paste0("module=",module_id,"<br>m=",d$m_used,
                                             "<br>t=%{x:.3g}<br>LOC=%{y:.3g}<br>hi=%{z:.3g}<extra></extra>")
                )
        
        if(show_hat){
                p <- p |>
                        plotly::add_surface(
                                x=d$t, y=d$loc, z=d$hat,
                                opacity=opacity_hat, showscale=showscale, name=paste0("hat m=", d$m_used),
                                hovertemplate=paste0("module=",module_id,"<br>m=",d$m_used,
                                                     "<br>t=%{x:.3g}<br>LOC=%{y:.3g}<br>hat=%{z:.3g}<extra></extra>")
                        )
        }
        
        p |> plotly::layout(
                title = paste0("Surface band (module ", module_id, ", m=", d$m_used,
                               ", LOC∈[", min(d$loc), ",", max(d$loc), "])"),
                scene = list(
                        xaxis=list(title="t"),
                        yaxis=list(title="LOCnorm"),
                        zaxis=list(title=ztitle)
                ),
                showlegend = FALSE
        )
}

## =========================================================
## (1) サブプロット：m_slices ごとに別パネルで表示（surface band）
## =========================================================
plot_band_surface_subplots <- function(bs, module_id=0,
                                       m_slices=c(3.0,4.5,6.0),
                                       loc_range=c(0,1),
                                       use_log10=TRUE,
                                       show_hat=TRUE,
                                       ncols=2,
                                       opacity_band=0.35,
                                       opacity_hat=0.85){
        if(!requireNamespace("plotly", quietly=TRUE)) stop("plotly を install.packages('plotly') してください。")
        
        plist <- lapply(m_slices, function(ms){
                plot_band_surface_panel(
                        bs, module_id=module_id, m_slice=ms,
                        loc_range=loc_range, use_log10=use_log10,
                        show_hat=show_hat,
                        opacity_band=opacity_band, opacity_hat=opacity_hat,
                        showscale=FALSE
                )
        })
        
        plotly::subplot(plist, nrows=ceiling(length(plist)/ncols), shareX=FALSE, shareY=FALSE) |>
                plotly::layout(
                        title = paste0("Surface bands by m (module ", module_id, ")"),
                        showlegend = FALSE
                )
}

## =========================================================
## (3) “リボン” mesh3d：lo–hi の間を側面付きで埋めた 3D
##  - x=t, y=LOCnorm, z=λ（log10可）
##  - 1つの m_slice 用（mごとは subplot で並べるのがおすすめ）
## =========================================================
.build_band_mesh3d <- function(x, y, zlo, zhi){
        nt <- length(x); nl <- length(y)
        stopifnot(is.matrix(zlo), is.matrix(zhi), nrow(zlo)==nt, ncol(zlo)==nl)
        
        idx <- function(i,j) (j-1L)*nt + i
        nV <- nt*nl
        
        X <- rep(x, times=nl)
        Y <- rep(y, each=nt)
        Zlo <- as.vector(zlo)
        Zhi <- as.vector(zhi)
        
        vx <- c(X, X)
        vy <- c(Y, Y)
        vz <- c(Zlo, Zhi)
        
        I <- integer(0); J <- integer(0); K <- integer(0)
        add_tri <- function(a,b,c){
                I <<- c(I, a-1L); J <<- c(J, b-1L); K <<- c(K, c-1L)
        }
        
        # top/bottom faces
        for(j in 1:(nl-1)){
                for(i in 1:(nt-1)){
                        v11 <- idx(i,   j)
                        v21 <- idx(i+1, j)
                        v12 <- idx(i,   j+1)
                        v22 <- idx(i+1, j+1)
                        
                        # hi face
                        add_tri(v11+nV, v21+nV, v22+nV)
                        add_tri(v11+nV, v22+nV, v12+nV)
                        
                        # lo face (reverse)
                        add_tri(v11, v22, v21)
                        add_tri(v11, v12, v22)
                }
        }
        
        # side faces
        j <- 1
        for(i in 1:(nt-1)){
                v1 <- idx(i, j); v2 <- idx(i+1, j)
                add_tri(v1, v2, v2+nV)
                add_tri(v1, v2+nV, v1+nV)
        }
        j <- nl
        for(i in 1:(nt-1)){
                v1 <- idx(i, j); v2 <- idx(i+1, j)
                add_tri(v1, v2+nV, v2)
                add_tri(v1, v1+nV, v2+nV)
        }
        i <- 1
        for(j in 1:(nl-1)){
                v1 <- idx(i, j); v2 <- idx(i, j+1)
                add_tri(v1, v2, v2+nV)
                add_tri(v1, v2+nV, v1+nV)
        }
        i <- nt
        for(j in 1:(nl-1)){
                v1 <- idx(i, j); v2 <- idx(i, j+1)
                add_tri(v1, v2+nV, v2)
                add_tri(v1, v1+nV, v2+nV)
        }
        
        list(vx=vx, vy=vy, vz=vz, i=I, j=J, k=K)
}

plot_band_ribbon_panel <- function(bs, module_id=0, m_slice=3.0,
                                   loc_range=c(0,1),
                                   use_log10=TRUE,
                                   show_hat_surface=TRUE,
                                   opacity_ribbon=0.35,
                                   opacity_hat=0.9){
        if(!requireNamespace("plotly", quietly=TRUE)) stop("plotly を install.packages('plotly') してください。")
        
        d <- .get_lambdaZ_grid_for_plot(bs, module_id, m_slice, loc_range, use_log10)
        ztitle <- if(d$use_log10) "log10(λ(t,z,m))" else "λ(t,z,m)"
        
        mesh <- .build_band_mesh3d(d$t, d$loc, d$lo, d$hi)
        
        p <- plotly::plot_ly() |>
                plotly::add_mesh3d(
                        x = mesh$vx, y = mesh$vy, z = mesh$vz,
                        i = mesh$i, j = mesh$j, k = mesh$k,
                        opacity = opacity_ribbon,
                        name = paste0("band m=", d$m_used),
                        hoverinfo = "skip"
                )
        
        if(show_hat_surface){
                p <- p |>
                        plotly::add_surface(
                                x=d$t, y=d$loc, z=d$hat,
                                opacity=opacity_hat, showscale=FALSE,
                                name=paste0("hat m=", d$m_used),
                                hovertemplate=paste0("module=",module_id,"<br>m=",d$m_used,
                                                     "<br>t=%{x:.3g}<br>LOC=%{y:.3g}<br>hat=%{z:.3g}<extra></extra>")
                        )
        }
        
        p |> plotly::layout(
                title = paste0("Ribbon band (mesh3d) (module ", module_id, ", m=", d$m_used,
                               ", LOC∈[", min(d$loc), ",", max(d$loc), "])"),
                scene = list(
                        xaxis=list(title="t"),
                        yaxis=list(title="LOCnorm"),
                        zaxis=list(title=ztitle)
                ),
                showlegend = FALSE
        )
}

## =========================================================
## (3) mごとに別パネル（subplot）：リボン band を並べる
## =========================================================
plot_band_ribbon_subplots <- function(bs, module_id=0,
                                      m_slices=c(3.0,4.5,6.0),
                                      loc_range=c(0,1),
                                      use_log10=TRUE,
                                      show_hat_surface=TRUE,
                                      ncols=2,
                                      opacity_ribbon=0.35,
                                      opacity_hat=0.9){
        if(!requireNamespace("plotly", quietly=TRUE)) stop("plotly を install.packages('plotly') してください。")
        
        plist <- lapply(m_slices, function(ms){
                plot_band_ribbon_panel(
                        bs, module_id=module_id, m_slice=ms,
                        loc_range=loc_range, use_log10=use_log10,
                        show_hat_surface=show_hat_surface,
                        opacity_ribbon=opacity_ribbon, opacity_hat=opacity_hat
                )
        })
        
        plotly::subplot(plist, nrows=ceiling(length(plist)/ncols), shareX=FALSE, shareY=FALSE) |>
                plotly::layout(
                        title = paste0("Ribbon bands (mesh3d) by m (module ", module_id, ")"),
                        showlegend = FALSE
                )
}

#実行
bs <- parametric_bootstrap_CI(
  dat0 = dat,
  theta_true = theta_true,
  T_end = 100, b_loc = 8,
  lambda0 = theta_true$lambda0,
  alpha1  = theta_true$alpha1,
  beta_m  = theta_true$beta_m,
  m_set = c(3.0, 4.5, 6.0),
  D_trig  = theta_true$D_trig,
  B = 200,
  do_grid = TRUE,
  grid_modules = 0:3,
  t_grid = seq(0, 100, length.out = 40),
  loc_grid = seq(0, 1, length.out = 64),
  m_grid = seq(3.0, 6.0, by = 0.1),
  out_dir = "output/boot_with_lambda",
  verbose = TRUE
)
bs$CI
bs$cover
bs$ok_rate
# 3D（module 0）
plot_lambdaZM_3d(bs, module_id = 0)
plot_lambdaZM_3d(bs, module_id = 1)
plot_lambdaZM_3d(bs, module_id = 2)
plot_lambdaZM_3d(bs, module_id = 3)

# 事前に bs <- parametric_bootstrap_CI(..., do_grid=TRUE, ...) を実行しておく
plot_lambdaZM_band_surface(
        bs,
        module_id = 0,
        m_slices = c(3.0, 4.5, 6.0),
        use_log10 = TRUE,
        show_hat = TRUE
)
# surface band（(1)）: mごとに別パネル
plot_band_surface_subplots(
        bs,
        module_id = 0,
        m_slices = c(3.0, 4.5, 6.0),
        loc_range = c(0, 1)     # ←ここを c(0.2,0.8) などに
)

# ribbon band（(3)）: mごとに別パネル
plot_band_ribbon_subplots(
        bs,
        module_id = 0,
        m_slices = c(3.0, 4.5, 6.0),
        loc_range = c(0, 1),
        use_log10 = TRUE
)

########################
plot_band_surface_panel2 <- function(bs, module_id=0, m_slice=3.0,
                                     loc_range=c(0,1),
                                     use_log10=TRUE,
                                     show_hat=TRUE,
                                     opacity_band=0.35,
                                     opacity_hat=0.85){
        if(!requireNamespace("plotly", quietly=TRUE)) stop("install.packages('plotly') が必要です。")
        
        d <- .get_lambdaZ_grid_for_plot(bs, module_id, m_slice, loc_range, use_log10)
        ztitle <- if(d$use_log10) "log10(λ(t,z,m))" else "λ(t,z,m)"
        
        p <- plotly::plot_ly()
        
        # lo
        p <- plotly::add_surface(
                p, x=d$t, y=d$loc, z=d$lo,
                opacity=opacity_band, showscale=FALSE, name=paste0("lo m=", d$m_used)
        )
        # hi
        p <- plotly::add_surface(
                p, x=d$t, y=d$loc, z=d$hi,
                opacity=opacity_band, showscale=FALSE, name=paste0("hi m=", d$m_used)
        )
        # hat
        if(show_hat){
                p <- plotly::add_surface(
                        p, x=d$t, y=d$loc, z=d$hat,
                        opacity=opacity_hat, showscale=FALSE, name=paste0("hat m=", d$m_used)
                )
        }
        
        plotly::layout(
                p,
                title = list(text = paste0("Surface band (module ", module_id, ", m=", d$m_used,
                                           ", LOC∈[", min(d$loc), ", ", max(d$loc), "])")),
                showlegend = FALSE,
                scene = list(
                        xaxis=list(title="t"),
                        yaxis=list(title="LOCnorm"),
                        zaxis=list(title=ztitle)
                )
        )
}

plot_band_surface_subplots2 <- function(bs, module_id=0,
                                        m_slices=c(3.0,4.5,6.0),
                                        loc_range=c(0,1),
                                        use_log10=TRUE,
                                        show_hat=TRUE,
                                        ncols=2,
                                        opacity_band=0.35,
                                        opacity_hat=0.85){
        if(!requireNamespace("plotly", quietly=TRUE)) stop("install.packages('plotly') が必要です。")
        
        plist <- lapply(m_slices, function(ms){
                plot_band_surface_panel2(
                        bs, module_id=module_id, m_slice=ms,
                        loc_range=loc_range, use_log10=use_log10,
                        show_hat=show_hat,
                        opacity_band=opacity_band, opacity_hat=opacity_hat
                )
        })
        
        nrows <- ceiling(length(plist)/ncols)
        
        p <- plotly::subplot(plist, nrows=nrows, shareX=FALSE, shareY=FALSE, margin=0.03)
        
        plotly::layout(
                p,
                title = list(text = paste0("Surface bands by m (module ", module_id, ")")),
                showlegend = FALSE
        )
}

plot_band_surface_subplots2(
        bs,
        module_id = 0,
        m_slices = c(3.0, 4.5, 6.0),
        loc_range = c(0, 1)   # ←ここを変える（例 c(0.2,0.8)）
)
###############
plot_band_ribbon_panel2 <- function(bs, module_id=0, m_slice=3.0,
                                    loc_range=c(0,1),
                                    use_log10=TRUE,
                                    show_hat_surface=TRUE,
                                    opacity_ribbon=0.35,
                                    opacity_hat=0.9){
        if(!requireNamespace("plotly", quietly=TRUE)) stop("install.packages('plotly') が必要です。")
        
        d <- .get_lambdaZ_grid_for_plot(bs, module_id, m_slice, loc_range, use_log10)
        ztitle <- if(d$use_log10) "log10(λ(t,z,m))" else "λ(t,z,m)"
        
        mesh <- .build_band_mesh3d(d$t, d$loc, d$lo, d$hi)
        
        p <- plotly::plot_ly()
        
        # ribbon (mesh3d)
        p <- plotly::add_trace(
                p,
                type = "mesh3d",
                x = mesh$vx, y = mesh$vy, z = mesh$vz,
                i = mesh$i, j = mesh$j, k = mesh$k,
                opacity = opacity_ribbon,
                hoverinfo = "skip",
                name = paste0("band m=", d$m_used)
        )
        
        # hat surface
        if(show_hat_surface){
                p <- plotly::add_surface(
                        p, x=d$t, y=d$loc, z=d$hat,
                        opacity=opacity_hat, showscale=FALSE,
                        name=paste0("hat m=", d$m_used)
                )
        }
        
        plotly::layout(
                p,
                title = list(text = paste0("Ribbon band (mesh3d) (module ", module_id, ", m=", d$m_used,
                                           ", LOC∈[", min(d$loc), ", ", max(d$loc), "])")),
                showlegend = FALSE,
                scene = list(
                        xaxis=list(title="t"),
                        yaxis=list(title="LOCnorm"),
                        zaxis=list(title=ztitle)
                )
        )
}

plot_band_ribbon_subplots2 <- function(bs, module_id=0,
                                       m_slices=c(3.0,4.5,6.0),
                                       loc_range=c(0,1),
                                       use_log10=TRUE,
                                       show_hat_surface=TRUE,
                                       ncols=2,
                                       opacity_ribbon=0.35,
                                       opacity_hat=0.9){
        if(!requireNamespace("plotly", quietly=TRUE)) stop("install.packages('plotly') が必要です。")
        
        plist <- lapply(m_slices, function(ms){
                plot_band_ribbon_panel2(
                        bs, module_id=module_id, m_slice=ms,
                        loc_range=loc_range, use_log10=use_log10,
                        show_hat_surface=show_hat_surface,
                        opacity_ribbon=opacity_ribbon, opacity_hat=opacity_hat
                )
        })
        
        nrows <- ceiling(length(plist)/ncols)
        p <- plotly::subplot(plist, nrows=nrows, shareX=FALSE, shareY=FALSE, margin=0.03)
        
        plotly::layout(
                p,
                title = list(text = paste0("Ribbon bands (mesh3d) by m (module ", module_id, ")")),
                showlegend = FALSE
        )
}

plot_band_ribbon_subplots2(
        bs,
        module_id = 0,
        m_slices = c(3.0, 4.5, 6.0),
        loc_range = c(0, 1)
)
## =========================================================
## (1) m別subplot：surfaceで lo/hat/hi を重ねる（軸固定オプション付き）
## =========================================================
plot_band_surface_panel3 <- function(bs, module_id=0, m_slice=3.0,
                                     loc_range=c(0,1),
                                     fix_loc_axis=TRUE,
                                     use_log10=TRUE,
                                     show_hat=TRUE,
                                     opacity_band=0.35,
                                     opacity_hat=0.85){
        if(!requireNamespace("plotly", quietly=TRUE)) stop("install.packages('plotly') が必要です。")
        
        d <- .get_lambdaZ_grid_for_plot(bs, module_id, m_slice, loc_range, use_log10)
        ztitle <- if(d$use_log10) "log10(λ(t,z,m))" else "λ(t,z,m)"
        
        yaxis <- list(title="LOCnorm")
        if(isTRUE(fix_loc_axis)) yaxis$range <- c(0, 1)
        
        p <- plotly::plot_ly()
        
        p <- plotly::add_surface(p, x=d$t, y=d$loc, z=d$lo,
                                 opacity=opacity_band, showscale=FALSE, name="lo")
        p <- plotly::add_surface(p, x=d$t, y=d$loc, z=d$hi,
                                 opacity=opacity_band, showscale=FALSE, name="hi")
        if(show_hat){
                p <- plotly::add_surface(p, x=d$t, y=d$loc, z=d$hat,
                                         opacity=opacity_hat, showscale=FALSE, name="hat")
        }
        
        plotly::layout(
                p,
                title = list(text = paste0(
                        "Surface band (module ", module_id,
                        ", m=", d$m_used,
                        ", slice LOC∈[", d$loc_range[1], ", ", d$loc_range[2], "]",
                        if(fix_loc_axis) "  (axis fixed to [0,1])" else ""
                )),
                showlegend = FALSE,
                scene = list(
                        xaxis=list(title="t"),
                        yaxis=yaxis,
                        zaxis=list(title=ztitle)
                )
        )
}

plot_band_surface_subplots3 <- function(bs, module_id=0,
                                        m_slices=c(3.0,4.5,6.0),
                                        loc_range=c(0,1),
                                        fix_loc_axis=TRUE,
                                        use_log10=TRUE,
                                        show_hat=TRUE,
                                        ncols=2,
                                        opacity_band=0.35,
                                        opacity_hat=0.85){
        if(!requireNamespace("plotly", quietly=TRUE)) stop("install.packages('plotly') が必要です。")
        
        plist <- lapply(m_slices, function(ms){
                plot_band_surface_panel3(
                        bs, module_id=module_id, m_slice=ms,
                        loc_range=loc_range, fix_loc_axis=fix_loc_axis,
                        use_log10=use_log10, show_hat=show_hat,
                        opacity_band=opacity_band, opacity_hat=opacity_hat
                )
        })
        
        nrows <- ceiling(length(plist)/ncols)
        p <- plotly::subplot(plist, nrows=nrows, shareX=FALSE, shareY=FALSE, margin=0.03)
        
        plotly::layout(
                p,
                title = list(text = paste0("Surface bands by m (module ", module_id, ")")),
                showlegend = FALSE
        )
}
plot_band_ribbon_panel3 <- function(bs, module_id=0, m_slice=3.0,
                                    loc_range=c(0,1),
                                    fix_loc_axis=TRUE,
                                    use_log10=TRUE,
                                    show_hat_surface=TRUE,
                                    opacity_ribbon=0.35,
                                    opacity_hat=0.9){
        if(!requireNamespace("plotly", quietly=TRUE)) stop("install.packages('plotly') が必要です。")
        
        d <- .get_lambdaZ_grid_for_plot(bs, module_id, m_slice, loc_range, use_log10)
        ztitle <- if(d$use_log10) "log10(λ(t,z,m))" else "λ(t,z,m)"
        
        yaxis <- list(title="LOCnorm")
        if(isTRUE(fix_loc_axis)) yaxis$range <- c(0, 1)
        
        mesh <- .build_band_mesh3d(d$t, d$loc, d$lo, d$hi)
        
        p <- plotly::plot_ly()
        
        p <- plotly::add_trace(
                p,
                type="mesh3d",
                x=mesh$vx, y=mesh$vy, z=mesh$vz,
                i=mesh$i, j=mesh$j, k=mesh$k,
                opacity=opacity_ribbon,
                hoverinfo="skip",
                name="band"
        )
        
        if(show_hat_surface){
                p <- plotly::add_surface(
                        p, x=d$t, y=d$loc, z=d$hat,
                        opacity=opacity_hat, showscale=FALSE, name="hat"
                )
        }
        
        plotly::layout(
                p,
                title = list(text = paste0(
                        "Ribbon band (mesh3d) (module ", module_id,
                        ", m=", d$m_used,
                        ", slice LOC∈[", d$loc_range[1], ", ", d$loc_range[2], "]",
                        if(fix_loc_axis) "  (axis fixed to [0,1])" else ""
                )),
                showlegend = FALSE,
                scene = list(
                        xaxis=list(title="t"),
                        yaxis=yaxis,
                        zaxis=list(title=ztitle)
                )
        )
}

plot_band_ribbon_subplots3 <- function(bs, module_id=0,
                                       m_slices=c(3.0,4.5,6.0),
                                       loc_range=c(0,1),
                                       fix_loc_axis=TRUE,
                                       use_log10=TRUE,
                                       show_hat_surface=TRUE,
                                       ncols=2,
                                       opacity_ribbon=0.35,
                                       opacity_hat=0.9){
        if(!requireNamespace("plotly", quietly=TRUE)) stop("install.packages('plotly') が必要です。")
        
        plist <- lapply(m_slices, function(ms){
                plot_band_ribbon_panel3(
                        bs, module_id=module_id, m_slice=ms,
                        loc_range=loc_range, fix_loc_axis=fix_loc_axis,
                        use_log10=use_log10,
                        show_hat_surface=show_hat_surface,
                        opacity_ribbon=opacity_ribbon,
                        opacity_hat=opacity_hat
                )
        })
        
        nrows <- ceiling(length(plist)/ncols)
        p <- plotly::subplot(plist, nrows=nrows, shareX=FALSE, shareY=FALSE, margin=0.03)
        
        plotly::layout(
                p,
                title = list(text = paste0("Ribbon bands (mesh3d) by m (module ", module_id, ")")),
                showlegend = FALSE
        )
}

plot_band_ribbon_subplots3(
        bs,
        module_id = 0,
        m_slices = c(3.0, 4.5, 6.0),
        loc_range = c(0.0, 1.0),
        fix_loc_axis = TRUE
)

plot_band_ribbon_subplots3(
        bs,
        module_id = 0,
        m_slices = c(3.0, 4.5, 6.0),
        loc_range = c(0.2, 0.8),
        fix_loc_axis = TRUE
)
