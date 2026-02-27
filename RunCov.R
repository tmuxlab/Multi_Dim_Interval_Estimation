source("GenerateData.R")
source("coverage_r.R")

theta_true <- list(lambda0=2, alpha1=-0.01,
                   A=0.30, alpha=1.0,
                   p=1.2, c=0.01, D_trig=3, beta_m=1.2)

mu_true <- rep(1/(4*2^8), 4*2^8)

covL <- coverage_experiment_lambda_slices(
        theta_true=theta_true, mu_true=mu_true,
        R=50, B=300, T_end=100,
        t_eval_vec=c(50,100),          # 50%時点=50, 最終=100（T_end=100なら）
        m_set=c(3.0,4.5,6.0),
        out_dir="COV_lambda"
)

# 平均被覆率（論文の表向け）
covL$coverage_mean

# ヒートマップ（t=50, m=4.5）
#mat <- coverage_heatmap_module_loc(covL, t_eval=50, m_val=4.5, b_loc=8)
#plot_heatmap_module_loc(mat, main="pointwise coverage (t=50, m=4.5)", b_loc=8)

# 例：t=50 と t=100、m=4.5 だけ欲しいなら
# save_coverage_heatmaps_png(
#         cov_res = covL,              # ← coverage_experiment_lambda_slices の結果
#         out_dir = "COV_png",
#         t_list = c(50, 100),
#         m_list = c(4.5),
#         b_loc = 8,
#         save_one = TRUE,
#         save_4panel = FALSE          # 1枚版だけ保存
# )

# 例：t=50,100 × m=3.0,4.5,6.0 を全部保存（1枚版+4パネル版）
save_coverage_heatmaps_png(
        cov_res = covL,
        out_dir = "COV_png",
        t_list = c(50, 100),
        m_list = c(3.0, 4.5, 6.0),
        b_loc = 8,
        save_one = TRUE,
        save_4panel = TRUE
)
