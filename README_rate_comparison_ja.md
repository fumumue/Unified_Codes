# BBGS+Goppa と RS+RM の同率比較

## 今回の比較で実際に使う復号法

- **BBGS+Goppa**: 内側は binary Goppa 復号（Patterson 系）。外側は Newton syndrome + **FR-VOTE**。この連接経路は `decode_fr_vote()` を直接呼ぶため、通常の `-vote`（BMS majority voting）ではありません。現状の連接経路に通常 VOTE/FR/BMS を選択する CLI はありません。
- **RS+RM**: 内側 RM 復号。内側で復号できなかったブロックは消失として RS/GRS 外復号に渡します。
- **インターリーブ**: 両方とも `D=128`。連接 runner の同じ D-way matrix interleaver を使用します。バーストは送信列の先頭から連続して注入します。
- **比較パラメータ**: BBGS outer `[1008,k_AG]` over GF(64) × Goppa `[60,6]`。対して RS outer `[945,k_RS]` over GF(4096) × RM(2,6) `[64,22]` を12情報ビットに短縮。`k_RS=k_AG/2` として N=60480 と K=6*k_AG=12*k_RS を一致させます。

## 実行

```sh
cc -O2 -std=gnu11 -fopenmp -Wall -Wextra \
  unified_codes_matched_errors.c -lm -o unified_codes_matched_errors

# 基準点 M=260、符号化率 8.63095% の比較
TRIALS=100 ./compare_bbgs_goppa_vs_rmrs.sh

# Mを振り、各MでRS+RMのKを合わせる。上限側のM=246..260も含む
TRIALS=100 ./sweep_bbgs_goppa_vs_rsrm_rate.sh
```

短い動作確認:

```sh
TRIALS=2 M_LIST="246 260" BER_LIST="5 8" BURST_LIST="1000" \
  ./sweep_bbgs_goppa_vs_rsrm_rate.sh
```

結果は `summary.csv` と各ログに保存されます。

## 符号化率の範囲と注意

この固定構成（BBGS q=4,m=6,a=1,b=2、内符号 Goppa[60,6]）でMを動かせる範囲は、おおむね 2.44%〜8.77%です。M=246で outer k=884、全体 K=5304、R=8.77%。M=260で K=5220、R=8.63%。したがって、この sweep は現構成の上限付近までを含みますが、20%や50%のような一般的な高レート符号の比較ではありません。

また、これは異なる外符号方式の比較です。BBGS+Goppa は内側失敗をゼロシンボルに置換し、外側では通常のシンボル誤りとしてFR-VOTEに渡します。一方、RS+RMは内側失敗位置を消失としてRSに渡します。この違いは実装上の設計差として結果に影響します。誤り率だけでなく、inner failures、outer symbol errors、detected failures、miscorrectionsも併せて確認してください。

このスクリプトは通常VOTEを比較していません。FR-VOTEと通常VOTEを同一構成で比較するには、連接側に復号法選択を実装する追加改修が必要です。
