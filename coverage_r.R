#bs後
binom_ci_wilson <- function(x, n, conf=0.95){
        # x successes out of n
        z <- stats::qnorm(1 - (1-conf)/2)
        phat <- x/n
        denom <- 1 + z^2/n
        center <- (phat + z^2/(2*n))/denom
        half <- z*sqrt(phat*(1-phat)/n + z^2/(4*n^2))/denom
        c(lo = center-half, hi = center+half)
}
parametric_bootstrap_CI_lambda_slices <- function(
                dat0,
                theta_true = NULL,
                T_end = 100,
                t_eval_vec = c(0.5*T_end, T_end),   # ★ 50%時点と最終時刻
                b_loc = 8,
                D_mu = 3,
                B = 200,
                seed = 1,
                lambda0 = 2,
                alpha1  = -0.01,
                beta_m  = 1.2,
                Mc = 3.0,
                Mmax = 6.0,
                dM = 0.1,
                m_set = c(3.0, 4.5, 6.0),
                D_trig  = 3,
                init = list(A=0.3, alpha=1.0, p=1.2, c=0.01),
                max_iter = 30,
                tol = 1e-6,
                verbose = TRUE
){
        set.seed(seed)
        
        dat0 <- dat0[order(dat0$time_t), ]
        t0   <- dat0$time_t
        x0   <- dat0$moduleID_x
        loc0 <- dat0$LOCnorm_y
        m0   <- dat0$difficulty_m
        
        # ---- 1) fit on original ----
        if(verbose) cat("Fitting on original data...\n")
        res0 <- fit_etas_mu_Aapc(
                t=t0, x=x0, loc_norm=loc0, m=m0,
                lambda0=lambda0, alpha1=alpha1,
                b_loc=b_loc, D_trig=D_trig, D_mu=D_mu,
                A_init=init$A, alpha_init=init$alpha, p_init=init$p, c_init=init$c,
                max_iter=max_iter, tol=tol, verbose=verbose
        )
        
        get_mu <- function(obj){
                if(!is.null(obj$mu)) return(obj$mu)
                if(!is.null(obj$res) && !is.null(obj$res$mu)) return(obj$res$mu)
                stop("mu not found in fit result.")
        }
        
        mu_hat <- get_mu(res0)
        
        theta_hat <- list(
                lambda0=lambda0, alpha1=alpha1,
                A=res0$A, alpha=res0$alpha, p=res0$p, c=res0$c,
                D_trig=D_trig, Mc=Mc
        )
        
        # ---- s(m) fixed ----
        s_set <- s_gr_trunc_pmf(m_set, beta=beta_m, Mc=Mc, Mmax=Mmax, dM=dM)
        names(s_set) <- paste0("m=", m_set)
        
        # ---- point estimate λz at multiple times ----
        lambdaZ_hat <- lapply(t_eval_vec, function(tt){
                lambda_z_at_time_allz(
                        t_eval=tt, dat=dat0, theta=theta_hat, mu_vec=mu_hat,
                        b_loc=b_loc, D_trig=D_trig, Mc=Mc
                )
        })
        names(lambdaZ_hat) <- paste0("t=", t_eval_vec)
        
        # ---- bootstrap storage ----
        nZ <- length(mu_hat)
        nT <- length(t_eval_vec)
        
        boot_mat <- matrix(NA_real_, nrow=B, ncol=4)
        colnames(boot_mat) <- c("A","alpha","p","c")
        failed <- logical(B)
        
        # ★ λz(t,z) を時点ごとに保存：配列 B×nT×nZ
        lambdaZ_boot <- array(NA_real_, dim=c(B, nT, nZ))
        
        # init chain
        A_init_b <- theta_hat$A
        alpha_init_b <- theta_hat$alpha
        p_init_b <- theta_hat$p
        c_init_b <- theta_hat$c
        
        for(b in 1:B){
                if(verbose && (b %% max(1, floor(B/10)) == 0)) cat("bootstrap", b, "/", B, "\n")
                
                sim <- simulate_etas_modified(
                        T_end=T_end,
                        lambda0=theta_hat$lambda0, alpha1=theta_hat$alpha1,
                        mu_vec=mu_hat, b_loc=b_loc,
                        beta_m=beta_m,
                        A=theta_hat$A, alpha=theta_hat$alpha, p=theta_hat$p, c=theta_hat$c,
                        D_trig=theta_hat$D_trig,
                        max_events=50000
                )
                
                fitb <- try(
                        fit_etas_mu_Aapc(
                                t=sim$time_t, x=sim$moduleID_x, loc_norm=sim$LOCnorm_y, m=sim$difficulty_m,
                                lambda0=theta_hat$lambda0, alpha1=theta_hat$alpha1,
                                b_loc=b_loc, D_trig=theta_hat$D_trig, D_mu=D_mu,
                                A_init=A_init_b, alpha_init=alpha_init_b, p_init=p_init_b, c_init=c_init_b,
                                max_iter=max_iter, tol=tol,
                                verbose=FALSE
                        ),
                        silent=TRUE
                )
                
                if(inherits(fitb, "try-error")){
                        failed[b] <- TRUE
                        next
                }
                
                boot_mat[b,] <- c(fitb$A, fitb$alpha, fitb$p, fitb$c)
                
                # update chain
                A_init_b <- fitb$A
                alpha_init_b <- fitb$alpha
                p_init_b <- fitb$p
                c_init_b <- fitb$c
                
                mu_b <- get_mu(fitb)
                theta_b <- list(
                        lambda0=lambda0, alpha1=alpha1,
                        A=fitb$A, alpha=fitb$alpha, p=fitb$p, c=fitb$c,
                        D_trig=D_trig, Mc=Mc
                )
                
                # ★ 複数時点の λz を保存
                for(k in seq_along(t_eval_vec)){
                        tt <- t_eval_vec[k]
                        lambdaZ_boot[b, k, ] <- lambda_z_at_time_allz(
                                t_eval=tt, dat=sim, theta=theta_b, mu_vec=mu_b,
                                b_loc=b_loc, D_trig=D_trig, Mc=Mc
                        )
                }
        }
        
        boot_df <- as.data.frame(boot_mat)
        boot_df$failed <- failed
        
        ok <- !failed &
                is.finite(boot_df$A) & is.finite(boot_df$alpha) & is.finite(boot_df$p) & is.finite(boot_df$c) &
                apply(lambdaZ_boot, 1, function(v) all(is.finite(v)))
        
        # ---- parameter CI ----
        CI_param <- apply(boot_df[ok, c("A","alpha","p","c")], 2, quantile,
                          probs=c(0.025,0.975), na.rm=TRUE)
        
        # ---- λz CI at each time ----
        lambdaZ_CI <- lapply(seq_along(t_eval_vec), function(k){
                mat <- lambdaZ_boot[ok, k, , drop=FALSE]  # n_ok × 1 × nZ
                # drop to matrix n_ok × nZ
                mat <- matrix(mat, ncol=nZ)
                lo <- apply(mat, 2, quantile, probs=0.025, na.rm=TRUE)
                hi <- apply(mat, 2, quantile, probs=0.975, na.rm=TRUE)
                list(lo=lo, hi=hi)
        })
        names(lambdaZ_CI) <- paste0("t=", t_eval_vec)
        
        # ---- λ(t,z,m)=s(m)*λz(t,z) : point & CI ----
        lambdaZM_hat <- list()
        lambdaZM_CI  <- list()
        for(k in seq_along(t_eval_vec)){
                nm_t <- paste0("t=", t_eval_vec[k])
                lambdaZM_hat[[nm_t]] <- lapply(seq_along(m_set), function(j){
                        s_set[j] * lambdaZ_hat[[nm_t]]
                })
                names(lambdaZM_hat[[nm_t]]) <- paste0("m=", m_set)
                
                lambdaZM_CI[[nm_t]] <- lapply(seq_along(m_set), function(j){
                        list(
                                lo = s_set[j] * lambdaZ_CI[[nm_t]]$lo,
                                hi = s_set[j] * lambdaZ_CI[[nm_t]]$hi
                        )
                })
                names(lambdaZM_CI[[nm_t]]) <- paste0("m=", m_set)
        }
        
        list(
                res0=res0,
                theta_hat=theta_hat,
                CI_param=CI_param,
                boot=boot_df,
                ok_rate=mean(ok),
                t_eval_vec=t_eval_vec,
                m_set=m_set,
                s_set=s_set,
                lambdaZ_hat=lambdaZ_hat,
                lambdaZ_CI=lambdaZ_CI,
                lambdaZM_hat=lambdaZM_hat,
                lambdaZM_CI=lambdaZM_CI
        )
}
coverage_experiment_lambda_slices <- function(
                theta_true,
                mu_true,
                R=30,
                B=200,
                T_end=100,
                t_eval_vec = c(0.5*T_end, T_end),  # ★ 50%時点と最終
                b_loc=8,
                D_mu=3,
                seed=1,
                max_events=50000,
                init=list(A=0.3, alpha=1.0, p=1.2, c=0.01),
                max_iter=30,
                tol=1e-6,
                Mc=3.0,
                Mmax=6.0,
                dM=0.1,
                m_set=c(3.0, 4.5, 6.0),
                verbose=TRUE,
                out_dir=NULL
){
        stopifnot(is.list(theta_true))
        set.seed(seed)
        
        if(!is.null(out_dir)) dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)
        
        nZ <- length(mu_true)
        nT <- length(t_eval_vec)
        nM <- length(m_set)
        
        # カウントで持つ（巨大配列を避ける）
        hit_ZM <- array(0L, dim=c(nT, nM, nZ))  # 入った回数
        tot_rep <- 0L
        
        failed <- logical(R)
        Nvec   <- integer(R)
        ok_rate_inner <- numeric(R)
        
        for(r in 1:R){
                if(verbose) cat(sprintf("\n=== Lambda-coverage replicate %d / %d ===\n", r, R))
                
                dat <- simulate_etas_modified(
                        T_end=T_end,
                        lambda0=theta_true$lambda0,
                        alpha1 =theta_true$alpha1,
                        mu_vec =mu_true,
                        b_loc  =b_loc,
                        beta_m =theta_true$beta_m,
                        A      =theta_true$A,
                        alpha  =theta_true$alpha,
                        p      =theta_true$p,
                        c      =theta_true$c,
                        D_trig =theta_true$D_trig,
                        max_events=max_events,
                        seed=sample.int(1e9,1)
                )
                Nvec[r] <- nrow(dat)
                
                bs <- try(
                        parametric_bootstrap_CI_lambda_slices(
                                dat0=dat,
                                T_end=T_end,
                                t_eval_vec=t_eval_vec,
                                b_loc=b_loc, D_mu=D_mu, B=B,
                                seed=sample.int(1e9,1),
                                lambda0=theta_true$lambda0,
                                alpha1 =theta_true$alpha1,
                                beta_m =theta_true$beta_m,
                                Mc=Mc, Mmax=Mmax, dM=dM,
                                m_set=m_set,
                                D_trig=theta_true$D_trig,
                                init=init,
                                max_iter=max_iter, tol=tol,
                                verbose=FALSE
                        ),
                        silent=TRUE
                )
                if(inherits(bs, "try-error")){
                        failed[r] <- TRUE
                        if(verbose) cat("  -> failed\n")
                        next
                }
                ok_rate_inner[r] <- bs$ok_rate
                
                # このrepの「真の」λ(t,z,m|H) を同じ履歴(dat)で計算
                thetaT <- list(
                        lambda0=theta_true$lambda0, alpha1=theta_true$alpha1,
                        A=theta_true$A, alpha=theta_true$alpha, p=theta_true$p, c=theta_true$c,
                        D_trig=theta_true$D_trig, Mc=Mc
                )
                
                # s(m) は真値beta固定（あなたの前提）
                s_set_true <- s_gr_trunc_pmf(m_set, beta=theta_true$beta_m, Mc=Mc, Mmax=Mmax, dM=dM)
                
                # 判定：各 (t,m,z) で真値がCIに入ったか
                for(k in seq_along(t_eval_vec)){
                        nm_t <- paste0("t=", t_eval_vec[k])
                        
                        Lz_true <- lambda_z_at_time_allz(
                                t_eval=t_eval_vec[k], dat=dat, theta=thetaT, mu_vec=mu_true,
                                b_loc=b_loc, D_trig=theta_true$D_trig, Mc=Mc
                        )
                        
                        for(j in seq_along(m_set)){
                                L_true <- s_set_true[j] * Lz_true
                                lo <- bs$lambdaZM_CI[[nm_t]][[paste0("m=", m_set[j])]]$lo
                                hi <- bs$lambdaZM_CI[[nm_t]][[paste0("m=", m_set[j])]]$hi
                                hit_ZM[k, j, ] <- hit_ZM[k, j, ] + as.integer(lo <= L_true & L_true <= hi)
                        }
                }
                
                tot_rep <- tot_rep + 1L
                if(verbose) cat(sprintf("  N=%d  ok_rate=%.3f\n", Nvec[r], ok_rate_inner[r]))
        }
        
        # 被覆率（zごと）
        cov_ZM <- hit_ZM / max(1L, tot_rep)
        
        # 全z平均（論文で1行にしやすい）
        cov_mean <- apply(cov_ZM, c(1,2), mean, na.rm=TRUE)
        dimnames(cov_mean) <- list(paste0("t=", t_eval_vec), paste0("m=", m_set))
        
        out <- list(
                theta_true=theta_true,
                mu_true=mu_true,
                R=R, B=B,
                success_reps=tot_rep,
                success_rate=tot_rep/R,
                n_events=Nvec,
                failed=failed,
                ok_rate_inner=ok_rate_inner,
                t_eval_vec=t_eval_vec,
                m_set=m_set,
                coverage_by_z = cov_ZM,        # [t,m,z]
                coverage_mean = cov_mean       # [t,m]
        )
        
        if(!is.null(out_dir)){
                saveRDS(out, file=file.path(out_dir, "coverage_lambda_slices.rds"))
                write.csv(cov_mean, file.path(out_dir, "coverage_mean_tm.csv"), row.names=TRUE)
                write.csv(data.frame(rep=1:R, N=Nvec, failed=failed, ok_rate_inner=ok_rate_inner),
                          file.path(out_dir, "rep_stats.csv"), row.names=FALSE)
        }
        
        out
}
# ---- cov_res から (t,m) の mat4xL (4×256) を作る ----
coverage_heatmap_module_loc <- function(cov_res, t_eval, m_val, b_loc=8){
        cov_ZM <- cov_res$coverage_by_z
        t_eval_vec <- cov_res$t_eval_vec
        m_set <- cov_res$m_set
        
        kt <- which.min(abs(t_eval_vec - t_eval))
        km <- which.min(abs(m_set - m_val))
        
        v <- cov_ZM[kt, km, ]  # length 1024
        
        L <- 2^b_loc
        mat <- matrix(v, nrow=L, ncol=4)   # locbin × module（z=module*L+locbin の並びを仮定）
        t(mat)  # 4 × L (module × locbin)
}

plot_heatmap_module_loc_one_with_legend <- function(mat4xL, main="coverage",
                                                    zlim=c(0,1), ncol_legend=100){
        loc_norm <- seq(0, 1, length.out = ncol(mat4xL))
        moduleID <- 0:3
        Z <- t(mat4xL)  # 256×4
        
        old <- par(c("mar","mfrow","oma","mgp","las","xaxs","yaxs"))
        on.exit(par(old), add=TRUE)
        
        layout(matrix(c(1,2), nrow=1), widths=c(1,0.18))
        par(mar=c(4,4,3,1))
        
        pal <- grDevices::colorRampPalette(c("white","yellow","orange","red","darkred"))(ncol_legend)
        
        image(x=loc_norm, y=moduleID, z=Z,
              col=pal, zlim=zlim,
              xlab="LOCnorm", ylab="moduleID",
              main=main)
        axis(2, at=moduleID, labels=moduleID, las=1)
        box()
        
        # legend
        par(mar=c(4,2,3,4))
        zseq <- seq(zlim[1], zlim[2], length.out=ncol_legend)
        image(x=1, y=zseq, z=matrix(zseq, nrow=1),
              col=pal, xlab="", ylab="coverage", xaxt="n",
              main="legend")
        axis(4, las=1)
        box()
}


plot_heatmap_module_loc_with_legend <- function(mat4xL, main="coverage",
                                                zlim=c(0,1), ncol_legend=100){
        loc_norm <- seq(0, 1, length.out = ncol(mat4xL))
        y <- c(0, 1)
        
        old <- par(c("mar","mfrow","oma","mgp","las","xaxs","yaxs"))
        on.exit(par(old), add=TRUE)
        
        layout(matrix(c(1,2,5,
                        3,4,5), nrow=2, byrow=TRUE),
               widths=c(1,1,0.18), heights=c(1,1))
        
        par(mar=c(4,4,3,1))
        pal <- grDevices::colorRampPalette(c("white","yellow","orange","red","darkred"))(ncol_legend)
        
        for(mod in 1:4){
                v <- mat4xL[mod, ]
                z <- cbind(v, v)  # 256×2 の帯
                image(x=loc_norm, y=y, z=z,
                      col=pal, zlim=zlim,
                      xlab="LOCnorm", ylab="", yaxt="n",
                      main=paste0(main, "  module=", mod-1))
                box()
        }
        
        # legend
        par(mar=c(4,2,3,4))
        zseq <- seq(zlim[1], zlim[2], length.out=ncol_legend)
        image(x=1, y=zseq, z=matrix(zseq, nrow=1),
              col=pal, xlab="", ylab="coverage", xaxt="n",
              main="legend")
        axis(4, las=1)
        box()
}

# ---- PNG保存ラッパー：tやmを複数指定して一括出力 ----
save_coverage_heatmaps_png <- function(
                cov_res,
                out_dir = "COV_png",
                t_list = c(50, 100),
                m_list = c(3.0, 4.5, 6.0),
                b_loc = 8,
                zlim = c(0,1),
                width = 1600,   # ★少し大きめ
                height = 1200,  # ★少し大きめ
                res = 150,
                save_one = TRUE,
                save_4panel = TRUE
){
        dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
        
        for(tt in t_list){
                for(mm in m_list){
                        mat <- coverage_heatmap_module_loc(cov_res, t_eval=tt, m_val=mm, b_loc=b_loc)
                        tag <- paste0("t_", gsub("\\.","_", as.character(tt)),
                                      "__m_", gsub("\\.","_", as.character(mm)))
                        
                        if(save_one){
                                f1 <- file.path(out_dir, paste0("coverage_one_", tag, ".png"))
                                grDevices::png(filename=f1, width=width, height=height, res=res)
                                plot_heatmap_one_bottom_legend(
                                        mat,
                                        main=paste0("pointwise coverage (t=", tt, ", m=", mm, ")"),
                                        zlim=zlim
                                )
                                grDevices::dev.off()
                        }
                        
                        if(save_4panel){
                                f2 <- file.path(out_dir, paste0("coverage_4panel_", tag, ".png"))
                                grDevices::png(filename=f2, width=width, height=height, res=res)
                                plot_heatmap_4panel_bottom_legend(
                                        mat,
                                        main=paste0("pointwise coverage (t=", tt, ", m=", mm, ")"),
                                        zlim=zlim
                                )
                                grDevices::dev.off()
                        }
                }
        }
        
        invisible(TRUE)
}

plot_heatmap_module_loc <- function(mat4xL, main="coverage", b_loc=8){
        loc_norm <- seq(0, 1, length.out = ncol(mat4xL))  # 256点
        y <- c(0, 1)                                      # ★長さ2にする（帯の高さ）
        
        old <- par(no.readonly=TRUE); on.exit(par(old))
        par(mfrow=c(2,2), mar=c(4,4,3,1))
        
        for(mod in 1:4){
                v <- mat4xL[mod, ]                 # length 256
                z <- cbind(v, v)                   # ★256×2 にして「帯」にする
                
                image(x=loc_norm, y=y, z=z,
                      xlab="LOCnorm", ylab="", yaxt="n",
                      main=paste0(main, "  module=", mod-1))
                box()
        }
}
plot_heatmap_module_loc_one <- function(mat4xL, main="coverage"){
        loc_norm <- seq(0, 1, length.out = ncol(mat4xL))  # 256
        moduleID <- 0:3                                   # 4
        
        Z <- t(mat4xL)  # ★ 256×4 にする（x=loc, y=module）
        
        image(x=loc_norm, y=moduleID, z=Z,
              xlab="LOCnorm", ylab="moduleID",
              main=main)
        box()
}
plot_heatmap_one_bottom_legend <- function(mat4xL, main="coverage",
                                           zlim=c(0,1), ncol_legend=100){
        loc_norm <- seq(0, 1, length.out = ncol(mat4xL))
        moduleID <- 0:3
        Z <- t(mat4xL)  # 256×4
        
        pal <- grDevices::colorRampPalette(c("white","yellow","orange","red","darkred"))(ncol_legend)
        
        # 上：本体、下：凡例
        layout(matrix(c(1,2), nrow=2), heights=c(0.85, 0.15))
        
        # 本体（余白を小さめに）
        par(mar=c(4,4,3,1))
        image(x=loc_norm, y=moduleID, z=Z,
              col=pal, zlim=zlim,
              xlab="LOCnorm", ylab="moduleID",
              main=main)
        axis(2, at=moduleID, labels=moduleID, las=1)
        box()
        
        # 凡例（横）
        par(mar=c(3,4,1,1))
        zseq <- seq(zlim[1], zlim[2], length.out=ncol_legend)
        image(x=zseq, y=1, z=matrix(zseq, nrow=length(zseq), ncol=1),
              col=pal, xlab="coverage", ylab="", yaxt="n")
        ticks <- seq(zlim[1], zlim[2], by=0.2)
        axis(1, at=ticks, labels=sprintf("%.1f", ticks))
        box()
}

plot_heatmap_4panel_bottom_legend <- function(mat4xL, main="coverage",
                                              zlim=c(0,1), ncol_legend=100){
        loc_norm <- seq(0, 1, length.out = ncol(mat4xL))
        y <- c(0, 1)
        
        pal <- grDevices::colorRampPalette(c("white","yellow","orange","red","darkred"))(ncol_legend)
        
        # 2×2 + 下凡例（計5パネル）
        layout(matrix(c(1,2,
                        3,4,
                        5,5), nrow=3, byrow=TRUE),
               heights=c(0.45,0.45,0.10))
        
        par(mar=c(4,4,3,1))
        for(mod in 1:4){
                v <- mat4xL[mod, ]
                z <- cbind(v, v)  # 256×2
                image(x=loc_norm, y=y, z=z,
                      col=pal, zlim=zlim,
                      xlab="LOCnorm", ylab="", yaxt="n",
                      main=paste0(main, "  module=", mod-1))
                box()
        }
        
        # 凡例（横）
        par(mar=c(3,4,1,1))
        zseq <- seq(zlim[1], zlim[2], length.out=ncol_legend)
        image(x=zseq, y=1, z=matrix(zseq, nrow=length(zseq), ncol=1),
              col=pal, xlab="coverage", ylab="", yaxt="n")
        ticks <- seq(zlim[1], zlim[2], by=0.2)
        axis(1, at=ticks, labels=sprintf("%.1f", ticks))
        box()
}
