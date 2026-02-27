# Multi Dim Interval Estimationについて

製作者: Kento Oiwa

**Multi Dim Interval Estimation**は，ソフトウェアフォールト発見メトリクスデータから，パラメトリック・ブートストラップ方法を利用して区間推定を実行するプログラムです．

# Multi Dim Interval Estimationの使い方(反復推定部分)
## funclist.R
funclist.R は、提案モデル（ETAS 風自己励起型点過程）を構成する **基本関数群（時間・マーク・擬似空間カーネル、強度計算の補助、可視化）**をまとめたファイルです。
推定スクリプト（例：algoA.R）から source("funclist.R") して利用します。

使い方
 ・このリポジトリを clone した後、同じディレクトリで以下を実行して関数を読み込みます。

注意点
 ・moduleID は 0〜3 を想定しています（encode_xy() 内でチェックあり）。
 ・b_loc を変えると擬似空間セル数（2^{b_loc}×4）が変わるため、mu_vec の長さも整合している必要があります。

## AlgoA.R
algoA.R は、提案している 自己励起型点過程（ETAS 風）モデルのパラメータ推定（EM 風反復）を実行するスクリプトです。
内部で funclist.R を source("funclist.R") しているため、同じディレクトリに funclist.R が必要です。
依存関係
 funclist.R（必須）
 ・encode_loc(), hamming_code() など、擬似空間の符号化・距離計算系を利用します。

## GenerateData.R
GenerateData.R は、提案している 自己励起型点過程（ETAS 風）モデルに基づいて、疑似データ（イベント列）を生成するためのスクリプトです。
出力は、推定側（例：algoA.R の fit_etas_mu_Aapc()）がそのまま受け取れる形式のデータフレームになります。
出力列（イベントデータ）
time_t：イベント時刻
moduleID_x：モジュールID（0–3想定）
LOCnorm_y：0–1正規化された位置指標（擬似空間用）
difficulty_m：難易度（マグニチュード相当）
依存関係
 algoA.R（GenerateData.R の先頭で source("algoA.R") しています）
 ・さらに algoA.R 側で funclist.R を利用する想定です（同じディレクトリ推奨）。
使い方
 ・同じディレクトリで以下を実行し、疑似データ生成関数を読み込みます。

## bs.R
bs.R は、提案モデル（ETAS 風自己励起型点過程）に対して パラメトリック・ブートストラップを行い、
・パラメータ(A,α,p,c) の パーセンタイル信頼区間
・最終時刻Tにおける擬似空間強度λz(T,z) の 信頼区間
・さらに difficulty（マーク）で切ったλ(T,z,m)=s(m)λz(T,z) の 断面CI
をまとめて出力・保存するスクリプトです。
依存関係
 bs.R は以下の関数を使うため、事前に読み込める状態が必要です。
 ・fit_etas_mu_Aapc()（推定）：algoA.R
 ・simulate_etas_modified()（疑似データ生成）：GenerateData.R
 ・lambda_z_at_time_allz(), s_gr_trunc_pmf() など：funclist.R
実運用では、algoA.R / GenerateData.R 側で source() されている想定でもOKですが、同一ディレクトリに揃えておくのが安全。

## bs_after_plot.R
bs_after_plot.R は、bs.R（パラメトリック・ブートストラップ）で保存した結果を読み込み、推定結果・信頼区間・被覆状況などを図として確認するための後処理スクリプトです。
「推定→ブート→CI作成」まで終えた後に、確認用のプロットをまとめて出す用途で使います。
前提（先に実行しておくもの）
bs.R を実行して、結果（例：bootstrap_result_with_lambda.rds や各種 .csv）が out_dir に保存されていること
描画に使うため、推定・カーネル・強度計算の関数が必要な場合があります（環境により以下を source()）
 ・funclist.R
 ・algoA.R
 ・GenerateData.R
使い方
・保存済みのブート結果（RDS や CSV）を読み込んでプロットします。
（スクリプト内で readRDS() / read.csv() している箇所の out_dir を、あなたの保存先に合わせて変更してください。）
このスクリプトで確認できること（例）
パラメータ (A,α,p,c) の
・ブートストラップ分布（ヒストグラム等）
・percentile CI の位置
・真値（シミュレーション真値がある場合）との比較
最終時刻 T における
・λz(T,z) の点推定と CI（擬似空間セル z ごと）
・m で切ったλ(T,z,m)=s(m)λz(T,z) の断面 CI
（シミュレーション実験の場合）
・被覆（coverage） の確認（CI が真値を含むかどうか）
入出力の対応
入力：bs.R が保存した
・bootstrap_result_with_lambda.rds
・lambdaZ_T_CI.csv
・lambdaZM_T_CI_m_*.csv
・boot_theta.csv など
出力：図（画面表示、または png() などでファイル保存）

## bs_after_Number.R
bs_after_Number.R は、bs.R（パラメトリック・ブートストラップ）で作成した推定結果を使って、イベント個数に関する量（例：累積期待個数 Λ(t)、最終時刻までの期待個数、またはそれに準ずる集計量）を計算・整理するための後処理スクリプトです。
強度 λ の CI だけでなく、**「時刻 t までにどれくらい発生するか」**の不確実性（CI）を見たいときに使います。
前提（先に実行しておくもの）
bs.R を実行して、out_dir に結果が保存されていること（RDS / CSV）
例：bootstrap_result_with_lambda.rds、boot_theta.csv、mu_hat.csv など
必要に応じて以下の関数群を読み込める状態（環境により source()）
・funclist.R（強度計算・カーネル等）
・algoA.R（推定器の返り値構造など）
・GenerateData.R（シミュレーション真値がある場合など）
使い方
保存済み結果を読み込み、個数系の指標を計算します。
（スクリプト内の out_dir をあなたの保存先に合わせて変更してください。）
入出力の対応
入力：bs.R が保存した
・bootstrap_result_with_lambda.rds
・boot_theta.csv
・mu_hat.csv
（必要なら）元データ data_original.csv
出力：
・代表時点でのΛ(t) 点推定＋CI（CSVやRオブジェクト）
・coverage や幅（interval width）の集計（必要なら）

## bs3D.R
bs3D.R は、bs.R のパラメトリック・ブートストラップを拡張して、最終時刻だけでなく
・時刻 t
・擬似空間の位置 LOCnorm
・difficulty（マーク）m の 
3D格子上でλ(T,z,m)=s(m)λz(t,z) などの値を計算し、**分位点CI（percentile）**を作ったり、3Dサーフェスや帯（band）として可視化するためのスクリプトです。
前提（依存関係）
中で以下の関数を使う想定です（同一ディレクトリ推奨）。
・fit_etas_mu_Aapc()（推定）：algoA.R
・simulate_etas_modified()（疑似データ生成）：GenerateData.R
・calc_lambdaZ_grid(), lambda_z_at_time_allz(), s_gr_trunc_pmf() など：funclist.R

## RunEstimate.R
RunEstimate.R は、提案モデル（ETAS 風自己励起型点過程）の 推定処理をまとめて実行するための実行用スクリプトです。
「データ準備 → 推定（EM風反復） → 結果保存（mu・history・推定パラメータ）」の流れを、1本のスクリプトとして回す用途で使います。
前提（依存関係）
RunEstimate.R は、推定本体やカーネル計算を利用するため、同一ディレクトリに以下がある前提で使うのが安全です。
・algoA.R（推定：fit_etas_mu_Aapc()）
・funclist.R（擬似空間符号化・カーネル・強度計算など）
（データを生成してから推定する場合）GenerateData.R
使い方
1) データを用意する
実データを読み込む、または GenerateData.R で疑似データを作成します。
データ形式は以下を想定します：
・time_t
・moduleID_x
・LOCnorm_y
・difficulty_m
2) 推定を実行する
RunEstimate.R を実行すると、内部で fit_etas_mu_Aapc() を呼び出して推定します。

## coverage_r.R
coverage_r.R は、シミュレーション実験で作った信頼区間（CI）について、真値（真のパラメータ／真の強度など）をどれだけ含めているかを定量評価するためのスクリプトです。
特に、パラメトリック・ブートストラップ（例：bs.R, bs3D.R）で構成した CI の coverageを計算・集計します。
前提（先に用意しておくもの）
シミュレーションで 真値（true） が分かっていること
・例：theta_true（A, α, p, c 等）、真の μ、真の λ/Λ など
ブートストラップ結果（CI）が保存されていること
・例：bootstrap_result_with_lambda.rds、boot_theta.csv、lambdaZ_T_CI.csv、lambdaZM_T_CI_m_*.csv など
（どれを使うかはスクリプト内の読み込み設定に依存）
使い方
・保存済みのCI結果を読み込んで、真値がCIに入っているかを判定し、カバレッジを計算します。
（スクリプト内の out_dir やファイル名を、あなたの保存先に合わせてください。）

## RunCov.R
RunCov.R は、シミュレーション実験における coverage（被覆率）評価をまとめて実行するための実行用スクリプトです。
主に coverage_r.R の関数や処理を呼び出して、
・どの設定（例：パラメータ設定、乱数種、B回数、格子設定など）で
・どの出力フォルダ（out_dir）の結果に対して
・どの被覆（パラメータ／強度／断面／3D格子など）を
評価するかを指定し、一括で計算・保存します。
前提（依存関係）
coverage_r.R（計算の本体）
・RunCov.R から source("coverage_r.R") して使う想定です。
（評価対象の結果が必要）
・bs.R / bs3D.R 等で作ったブートストラップ結果（RDS/CSV）が out_dir に保存されていること
（真値が必要）
・theta_true や真の μ、真の λ/Λ 等が RunCov.R 内で用意されている、または読み込めること
使い方
・coverage を取りたい結果フォルダや設定を RunCov.R 内で指定し、実行します。
