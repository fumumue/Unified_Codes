# unified_codes: Srivastava / Generalized Srivastava 追加版

この版は、既存の `unified_codes` に **Srivastava符号** と **Generalized Srivastava (GS) 符号**を追加した実験用Cプログラムです。既存のBCH、Goppa、Alternant、AG符号などの実装はそのまま残しています。

## 追加した符号

### Original Srivastava

`--code srivastava` で選択します。Generalized Srivastava の multiplicity 1 に相当し、重みを

\[
z_i=\alpha_i^{\mu}
\]

として、パリティ検査行列のブロックを

\[
H_{j,i}=\frac{z_i}{\alpha_i-w_j}
\]

で構成します。設計距離の下界は \(d\ge s+1\)、保証訂正誤り数は \(\lfloor s/2\rfloor\) です。

### Generalized Srivastava (GS)

`--code gs` で選択します。相異なる極 \(w_1,\ldots,w_s\)、評価点 \(\alpha_i\)、非零重み \(z_i\) を使い、行列要素は

\[
H_{(j,\ell),i}=\frac{z_i}{(\alpha_i-w_j)^\ell},
\qquad 1\le j\le s,\quad 1\le\ell\le t.
\]

ここで `--s` が極の数、`--t` が各極の multiplicity です。設計距離の下界と保証訂正誤り数は

\[
d\ge st+1,\qquad t_{\rm err}=\left\lfloor\frac{st}{2}\right\rfloor.
\]

実装では、この検査行列を一般Alternant形式に変換してから、既存のAlternant用二元核・符号化処理を利用します。具体的には

\[
y_i=\frac{z_i}{\prod_{j=1}^{s}(\alpha_i-w_j)^t}
\]

とし、次数 \(<st\) の多項式基底を使う表現です。したがって、有限体演算、GF(2)への展開、系統生成行列の構築は既存コードを再利用します。

## 復号器

- デフォルト：Berlekamp–Massey (BM) 法
- `-eea`：一般Alternant鍵方程式を拡張ユークリッド互除法 (EEA) で解く復号器
- `-bm`：BM法を明示選択

一般Alternant/Srivastava/GSのEEAは、シンドローム多項式に対して

\[
\Lambda(x)S(x)\equiv\Omega(x)\pmod{x^r}
\]

を解き、Chien探索とForney公式で誤り位置・誤り値を求めます。EEAモードは現在、消失訂正をサポートしません。消失を含む場合はデフォルトのBMモードを使用してください。

Goppa符号の既存EEA復号は従来どおり別経路で、\(g(x)^2\) を法とする鍵方程式を使います。

## ビルド

```sh
gcc -O2 -std=gnu11 -fopenmp -Wall -Wextra \\
  unified_codes_srivastava.c -lm -o unified_codes
```

## 使用例

Original Srivastava、BM法：

```sh
./unified_codes --family classical --code srivastava \\
  --m 6 --s 8 --mu 0 --errors 4 -n 200 123
```

Original Srivastava、EEA：

```sh
./unified_codes --family classical --code srivastava \\
  --m 6 --s 8 --mu 0 --errors 4 -eea -n 200 123
```

Generalized Srivastava、BM法：

```sh
./unified_codes --family classical --code gs \\
  --m 6 --s 3 --t 2 --errors 3 -n 200 123
```

Generalized Srivastava、EEA：

```sh
./unified_codes --family classical --code gs \\
  --m 8 --s 4 --t 3 --errors 6 -eea -n 100 777
```

## オプション

| オプション | 意味 |
|---|---|
| `--code srivastava` | Original Srivastava |
| `--code gs` | Generalized Srivastava |
| `--m m` | 拡大体 GF(2^m)、既存実装の範囲は 2 ≤ m ≤ 12 |
| `--s s` | 極の数 |
| `--t t` | GSの極ごとのmultiplicity |
| `--mu mu` | Original Srivastavaの重み指数 \(z_i=\alpha_i^\mu\) |
| `--len n` | 符号長 |
| `--errors e` | 注入するビット誤り数 |
| `--erasures u` | 消失数。BMモードのみ |
| `-eea` | EEA復号 |
| `-bm` | BM復号 |
| `-n N` | 試行回数 |
| `-v` | 各試行を表示 |
| `--list` | パラメータだけ表示 |

`--family classical --help` で古典符号側のヘルプを表示します。

## 実装上の制約・注意

1. 現在の実装は **二元符号**（GF(2^m)上で検査行列を構成し、GF(2)に制限した符号）です。
2. 既存Forney/Chien処理が \(\alpha_i^{-1}\) を使うため、評価点 \(\alpha_i\) は非零に限定しています。また極 \(w_j\) と評価点は相異なるように生成します。このため自動生成時の最大長は \(2^m-1-s\) です。
3. GSの \(z_i\) は非零のランダム値です。Original Srivastavaは指定された \(\mu\) に対して \(z_i=\alpha_i^\mu\) を使います。
4. `--len` を指定しない場合、極を除いた利用可能な非零評価点をすべて使います。
5. `k` はGF(2)に展開したパリティ検査行列の厳密な階数から計算します。理論下界 \(k\ge n-mst\) と実際の次元が一致するとは限りません。
6. `-eea` は現在エラーのみを扱い、消失訂正はBMモードに限ります。
7. これは符号構成・復号の実験実装であり、暗号用途の安全性を保証するものではありません。公開鍵暗号としての採用には、鍵構造の漏えい・既知攻撃・パラメータ選定を別途評価する必要があります。

## 実行検証

この追加版を `gcc -O2 -std=gnu11 -fopenmp -Wall -Wextra` でコンパイルし、次を確認しました。

| テスト | 試行 | 結果 |
|---|---:|---:|
| Srivastava GF(2^6), s=8, BM, 4 errors | 200 | 200/200 PASS |
| Srivastava GF(2^6), s=8, EEA, 4 errors | 200 | 200/200 PASS |
| GS GF(2^6), s=3, t=2, BM, 3 errors | 200 | 200/200 PASS |
| GS GF(2^6), s=3, t=2, EEA, 3 errors | 200 | 200/200 PASS |
| GS GF(2^8), s=4, t=3, EEA, 6 errors | 100 | 100/100 PASS |
| GS GF(2^6), s=3, t=2, BM, 2 errors + 2 erasures | 100 | 100/100 PASS |
| Original Srivastava GF(2^7), s=9, mu=3, EEA, 4 errors | 100 | 100/100 PASS |
| 既存Goppa GF(2^8), deg g=10, EEA, 10 errors | 30 | 30/30 PASS |
| 既存Alternant GF(2^6), r=6, BM, 3 errors | 100 | 100/100 PASS |

コンパイルは成功しました。既存AGコード側の未使用関数に関する `-Wunused-function` 警告が4件出ます。今回追加したコードに由来するコンパイルエラーはありません。
