source("algoA.R")     # fit_etas_mu_Aapc を使う
##疑似データ生成
r_bg_times_custom <- function(T_end, lambda, alpha1){
        if(lambda <= 0) stop("lambda must be > 0")
        if(alpha1 >= 0) stop("This generator requires alpha1 < 0 (decaying).")
        sita <- -lambda / alpha1
        beta <- -alpha1
        
        get_nhpp_realization <- function(sita, beta, tn){
                i <- 1; t <- 0
                m <- rpois(1, sita)
                if(m <= 0) return(numeric(0))
                A <- numeric(0)
                while(i <= m){
                        t <- rexp(1, rate = 1) / beta / (m - i + 1) + t
                        if(t > tn) break
                        A <- c(A, t)
                        i <- i + 1
                }
                A
        }
        
        sort(get_nhpp_realization(sita, beta, T_end))
}

r_mark_trunc_gr_discrete <- function(n, beta, Mc=3.0, Mmax=6.0, dM=0.1){
        grid <- seq(Mc, Mmax, by=dM)
        prob <- s_gr_trunc_pmf(grid, beta=beta, Mc=Mc, Mmax=Mmax, dM=dM)
        prob <- prob / sum(prob)
        sample(grid, size=n, replace=TRUE, prob=prob)
}

r_trunc_omori <- function(n, tau_max, p, c){
        # G(tau)=1-(1+tau/c)^(1-p)
        Gmax <- 1 - (1 + tau_max/c)^(1 - p)
        u <- runif(n)
        # 1-(1+tau/c)^(1-p) = u*Gmax から tau を解く
        tau <- c * ((1 - u*Gmax)^(1/(1 - p)) - 1)
        tau
}

decode_xy <- function(code, b_loc=8){
        mask <- as.integer(2^b_loc - 1)
        x <- bitwAnd(bitwShiftR(as.integer(code), b_loc), 3L)       # 上位2bit
        loc_int <- bitwAnd(as.integer(code), mask)                  # 下位b_loc bit
        list(x = x, loc_int = loc_int, loc_norm = loc_int / mask)
}

r_child_state_hamming <- function(parent_code, b_loc=8, D=3){
        B <- 2 + b_loc
        q <- exp(-1/D) / (1 + exp(-1/D))    # flip確率
        code <- as.integer(parent_code)
        for(k in 0:(B-1)){
                if(runif(1) < q){
                        code <- bitwXor(code, bitwShiftL(1L, k))
                }
        }
        decode_xy(code, b_loc=b_loc)
}

simulate_etas_modified <- function(
                T_end = 100,
                # background u(t)=lambda0*exp(alpha1*t)
                lambda0 = 2, alpha1 = -0.01,
                # mu over state space (length S = 4*2^b_loc)
                mu_vec = NULL,
                b_loc = 8,
                # marks
                beta_m = 1.2, Mc = 3.0, Mmax = 6.0, dM = 0.1,
                # triggering
                A = 0.3, alpha = 1.0, p = 1.2, c = 0.01, D_trig = 3,
                # safety
                max_events = 50000,
                seed = NULL
){
        if(!is.null(seed)) set.seed(seed)
        
        state_space <- make_state_space(b_loc)
        S <- nrow(state_space)
        if(is.null(mu_vec)){
                mu_vec <- rep(1/S, S)
        } else {
                mu_vec <- pmax(mu_vec, 0); mu_vec <- mu_vec / sum(mu_vec)
        }
        
        # --- 1) background events ---
        t_bg <- r_bg_times_custom(T_end, lambda=lambda0, alpha1=alpha1)
        n_bg <- length(t_bg)
        if(n_bg == 0){
                return(data.frame(index=integer(0), time_t=numeric(0), moduleID_x=integer(0),
                                  LOCnorm_y=numeric(0), difficulty_m=numeric(0)))
        }
        
        # state: sample code by mu
        s_idx_bg <- sample.int(S, size=n_bg, replace=TRUE, prob=mu_vec)
        code_bg  <- state_space$code[s_idx_bg]
        
        xy_bg <- lapply(code_bg, decode_xy, b_loc=b_loc)
        x_bg   <- vapply(xy_bg, `[[`, integer(1), "x")
        loc_bg <- vapply(xy_bg, `[[`, numeric(1), "loc_norm")
        
        # mark
        m_bg <- r_mark_trunc_gr_discrete(n_bg, beta=beta_m, Mc=Mc, Mmax=Mmax, dM=dM)
        
        # event table (will grow)
        t_all   <- t_bg
        code_all<- code_bg
        x_all   <- x_bg
        loc_all <- loc_bg
        m_all   <- m_bg
        
        # --- 2) branching (queue) ---
        q_ptr <- 1
        while(q_ptr <= length(t_all)){
                if(length(t_all) >= max_events) break
                
                t0 <- t_all[q_ptr]
                if(t0 >= T_end){ q_ptr <- q_ptr + 1; next }
                
                # expected children inside window: k(m)*G(T-t0)
                k0 <- k_childrenNum(A, alpha, m_all[q_ptr], m_c = Mc)  # funclistの定義を利用:contentReference[oaicite:7]{index=7}
                Gmax <- 1 - (1 + (T_end - t0)/c)^(1 - p)               # = G_time(T_end-t0) と同じ:contentReference[oaicite:8]{index=8}
                mean_off <- k0 * Gmax
                n_off <- rpois(1, mean_off)
                
                if(n_off > 0){
                        # child times
                        tau <- r_trunc_omori(n_off, tau_max = (T_end - t0), p=p, c=c)
                        t_child <- t0 + tau
                        
                        # child states via hamming kernel sampling
                        parent_code <- code_all[q_ptr]
                        child_xy <- replicate(n_off, r_child_state_hamming(parent_code, b_loc=b_loc, D=D_trig), simplify=FALSE)
                        x_child   <- vapply(child_xy, `[[`, integer(1), "x")
                        loc_child <- vapply(child_xy, `[[`, numeric(1), "loc_norm")
                        code_child<- vapply(child_xy, function(z) encode_xy(z$x, encode_loc(z$loc_norm, b_loc), b_loc),
                                            integer(1))
                        
                        # child marks (independent)
                        m_child <- r_mark_trunc_gr_discrete(n_off, beta=beta_m, Mc=Mc, Mmax=Mmax, dM=dM)
                        
                        # append
                        t_all    <- c(t_all, t_child)
                        code_all <- c(code_all, code_child)
                        x_all    <- c(x_all, x_child)
                        loc_all  <- c(loc_all, loc_child)
                        m_all    <- c(m_all, m_child)
                }
                
                q_ptr <- q_ptr + 1
        }
        
        # --- 3) sort by time ---
        o <- order(t_all)
        out <- data.frame(
                index = seq_along(o),
                time_t = t_all[o],
                moduleID_x = x_all[o],
                LOCnorm_y = pmin(pmax(loc_all[o], 0), 1),
                difficulty_m = m_all[o]
        )
        out
}
