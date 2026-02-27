##
##funclist.R
##
c("b_beta_from_MLE     :mに基づくβ推定")
c("s_gr_trunc_pmf      :倍率s(m)累積分布関数")
c("u_time              :GS-法,u(t)の計算式")
c("k_childrenNum       :子数,k(m)の計算式")
c("g_time              :時間間隔分布,g(t)の計算式")
c("encode_loc          :y=Locnormの整数化0:2^b_loc-1")
c("encode_xy           :例x=moduleID[0:3]を0,256,512,768,+encode_loc(Locnorm,8)")
c("vec1024_to_mat256x4 :ベクトル1024からマトリクス256*4へ")
c("make_z_code         :encode_xy(x,encode_loc(y,b),b)")
c("make_zvec_module_loc:自作関数要らずでmake_z_code,0~1023")
c(".pop8               :0..255の表8bitにある1の数")
c("bitwAnd             :引数2進数化，共に1になる数")
c("bitwShiftR          :数値bitを右シフト=2^bで割る,4byteまで")
c("popcount32          :2進数10bitにしたとき1の数")
c("hamming_code        :二つの引数0~1023よりハミング距離取得")
c("f_space_hamming_norm:正規化ハミング距離カーネル関数")
c("f_hamming_norm_code :0~1023コード２つにカーネル関数")
#条件付きλ計算：パラメータと履歴必要
c("lambda_z_at_time_allz : 全空間z0~1023にて条件付き強度計算")
c("lambda_z_at_time_zvec : 軽量化，特定のzベクトルに")
c("calc_lambdaZ_grid     : 履歴dat,パラメータθ,μで特定時刻の特定局所瞬間発見率λ計算")
#プロット関数
c("plot_mu_by_module     : 反復結果からμの４パネルプロット")
##
##AlgoA.R 反復アルゴリズム
##
c("make_state_space            : 内部情報にx,y,z=0:3,0:255,0:1023を持つ全セルgridベクトルを生成")
C("event_state_index           : x,loc_normをｚに変換,encode_xyと同じ？")
c("hamming_kernel_vec_norm     : f_hamming_norm_codeと全部同じ")
c("G_time                      : gの累積関数")
c("U_time_integral             : GS-の累積関数")
c("logllik_obs_ts              : 観測対数尤度(時空間部分の計算関数)")
c("Q_ts                        : 期待完全対数尤度の計算")
c("etas_E_step_storePhi        : 背景確率Φと親確率Φ_ijを計算，保存")
c("update_mu_smooth            : μの平滑更新，全セルで背景確率Φ*カーネルの総和/背景Φ和")
c("expected_children_per_parent: 親確率Φ_ijより親毎の期待個数計算")
c("update_A_alpha              : パラメータθ=(A,alpha,p,c)のA,alpha更新")
c("update_p_c                  : パラメータθ=(A,alpha,p,c)のp,c更新")
c("fit_etas_mu_Aapc            : 上記更新関数の反復実行，背景空間成分μとパラメータの更新")
##
##GenerateData.R
##
c("r_bg_times_custom       :GS-によるサンプル生成")
c("r_mark_trunc_gr_discrete:n個のリサンプルmを生成,βに従う")
c("r_trunc_omori           :時間間隔関数G")
c("decode_xy               :z0~1023をx0~3,y_loc_int0~255,loc_norm0~1に戻す")
c("r_child_state_hamming   :親zを元に,反転確率に基づき子のデコード(x,y,Locnorm)を返す")
c("simulate_etas_modified  :提案モデルに基づくメトリクスデータ生成，再現性アリ")
