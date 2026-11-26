# Unified Codes

**Unified Codes** は、代数幾何（AG）符号と古典的な誤り訂正符号を、1つの
C ソースファイル・1つの実行ファイルから実験するためのプログラムです。

AG 符号側と古典符号側の実装を統合し、コマンドラインの `--family`
で切り替えられるようにしています。内部のアルゴリズムやデータ構造は各実装に残してあり、すべてを共通の符号
API に書き換えたものではありません。

## 主な機能

### 代数幾何符号（`ag`）

-   BBGS 曲線
-   Hermitian 曲線
-   Suzuki 曲線（専用エンジン）
-   `--terms` による平面曲線モデルの指定
-   1 点 AG 符号と多点 AG 符号
-   Feng--Rao（FR）、BMS、VOTE、BASIC などの復号経路
-   有限体上のシンドローム計算・誤り位置の探索
-   `-encode` による符号語生成を含めた復号確認

### 古典符号（`classical`）

-   BCH 符号
-   Golay 符号（通常形・拡張形）
-   Goppa 符号
-   Alternant 符号
-   Reed--Muller 符号
-   Hadamard 符号
-   Reed–Solomon（RS）符号
-   Generalized Reed–Solomon（GRS）符号
-   Gabidulin（rank-metric）符号
-   積符号（Product code）
-   連接符号（Concatenated code）

Reed--Muller 符号では、`RM(r,m)` を扱います。`RM(1,m)` は `-fht`
を指定すると、高速 Walsh--Hadamard
変換（FWHT）による最尤復号経路を使用します。`r > 1` の場合は通常の
Reed--Muller 復号経路を使用します。

## 必要環境

-   C コンパイラ（GCC / Clang など）
-   C11 互換のコンパイル環境
-   `libm`（`-lm`）
-   OpenMP は任意（利用できる環境では `-fopenmp` を付けます）

## ビルド

`unified_codes.c` と `Makefile` をリポジトリのルートに置き、次を実行します。

```sh
# OpenMP を有効にしてビルド（既定）
make

# OpenMP を無効にしてビルド
make OPENMP=

# 生成した実行ファイルを削除
make clean
```

Makefile を使わず直接ビルドする場合：

```sh
cc -O2 -std=gnu11 -Wall -Wextra -fopenmp \
  -o unified_codes unified_codes.c -lm
```

OpenMP を使わない場合：

```sh
cc -O2 -std=gnu11 -Wall -Wextra \
  -o unified_codes unified_codes.c -lm
```

`make` の既定設定は `CC=cc`、`CFLAGS=-O2 -std=gnu11 -Wall -Wextra`、
`OPENMP=-fopenmp`、`LDLIBS=-lm` です。必要に応じて `CC`、`CFLAGS`、
`OPENMP`、`LDFLAGS`、`LDLIBS` を make のコマンドラインで上書きできます。

未使用関数に関する警告が表示される場合があります。OpenMP を無効にした場合は、
OpenMP プラグマを無視する警告も表示されます。いずれも警告であり、
コンパイルエラーとは区別してください。

## ヘルプ

`--help` はコマンドラインの概要を表示するヘルプです。対話式メニューではありません。`--help-all` で AG 側と古典符号側の詳細ヘルプを続けて表示できます。

``` sh
# 統合ヘルプ
./unified_codes --help

# AG 符号のオプション
./unified_codes --family ag --help

# 古典符号のオプション
./unified_codes --family classical --help

# 両ファミリーの詳細ヘルプ
./unified_codes --help-all
```

## ファミリーの選択

AG 符号は既定のファミリーです。

``` sh
./unified_codes --family ag [AG オプション]
```

古典符号は次のように指定します。

``` sh
./unified_codes --family classical --code [符号名] [オプション]
```

`--code` を先頭に置く形式では、古典符号ファミリーが自動選択されます。

``` sh
./unified_codes --code rm --r 1 --m 8 -fht
```

## 実行例

### Reed--Muller 符号

`RM(1,8)` に対して Walsh--Hadamard 変換を使う例：

``` sh
./unified_codes --family classical \
  --code rm --r 1 --m 8 --errors 20 -fht -n 100
```

`RM(2,8)` の例：

``` sh
./unified_codes --family classical \
  --code rm --r 2 --m 8 --errors 7 -n 100
```

### BCH 符号

``` sh
./unified_codes --family classical \
  --code bch --m 6 --t 3 --errors 3 -n 100
```

### Golay 符号

``` sh
./unified_codes --family classical \
  --code golay --ext --errors 3 -n 100
```

### Goppa 符号

``` sh
./unified_codes --family classical \
  --code goppa --m 8 --t 10 --errors 10 -n 20
```

### Reed–Solomon（RS）符号とGeneralized Reed–Solomon（GRS）符号

RS 符号は GRS 符号の特別な場合です。符号語は
\((v_1f(a_1),\ldots,v_nf(a_n))\) として構成されます。

- `--code rs`：列乗数をすべて \(v_i=1\) に固定する通常の RS 符号
- `--code grs`：各列乗数に非零のランダムな有限体要素を使う GRS 符号
- `-shuffle`：評価点の順序をランダム化

RS 符号の例：

```sh
./unified_codes --family classical \
  --code rs --m 8 --len 255 --k 223 \
  --errors 10 --erasures 12
```

GRS 符号の例：

```sh
./unified_codes --family classical \
  --code grs --m 8 --len 255 --k 223 \
  --errors 10 --erasures 12
```

どちらも誤り数 \(e\) と消失数 \(s\) が \(2e+s\le n-k\) を満たす範囲を
一意復号の保証範囲とします。同じ `m`、`len`、`k` を指定すれば、RS と GRS の
符号長・次元・最小距離の下界を揃えて比較できます。
### Alternant 符号

```sh
./unified_codes --code alt --m 8 --len 200 --t 10 --errors 10 -n 100
```

### Gabidulin（rank-metric）符号

``` sh
./unified_codes --family classical \
  --code rank --m 8 --len 8 --k 4 --errors 2
```

### 積符号

``` sh
./unified_codes --family classical \
  --code product --c1 bch:5:2 --c2 golay --errors 8
```

### 連接符号

``` sh
./unified_codes --family classical \
  --code concat --m 12 --len 100 --k 60 \
  --inner golay --errors 100
```

### BBGS 多点 AG 符号

``` sh
./unified_codes --family ag \
  --m 3 --q 2 --errors 2 \
  --pole '#0:6,#1:9' -basic -n 20
```

AG
符号では曲線・有限体・極因子・シンドローム境界などの組み合わせに制約があります。最初に
`-limits`
を実行し、そのビルドと曲線で扱えるパラメータ範囲を確認してください。

``` sh
./unified_codes --family ag -limits
```

## 主なオプション

オプションの全一覧と各符号固有の指定は、必ず実行ファイルのヘルプを参照してください。

### 共通的な実験オプション

  ----------------------------------------------------------------------------------------------
  オプション                          意味
  ----------------------------------- ----------------------------------------------------------
  `--errors N`                        注入する誤り数

  `-n N`                              ランダム試行回数

  `-v`                                各試行の詳細表示（各ファミリーの実装に依存）

  `-show-values`                      誤り値やベクトルなどを表示

  `--poly P`                          有限体の既約多項式を指定（対応する符号・ファミリーのみ）
  ----------------------------------------------------------------------------------------------

### 古典符号側

  オプション       意味
  ---------------- ----------------------------------------------
  `--code NAME`    符号の種類
  `--m M`          有限体の次数パラメータ
  `--len N`        符号長
  `--k K`          符号次元（対応する符号のみ）
  `--t T`          訂正能力・設計パラメータ（対応する符号のみ）
  `--r R`          Reed--Muller の次数
  `--erasures S`   消失数（対応する符号のみ）
  `-fht`           `RM(1,m)` の FWHT 復号
  `--list`         符号パラメータを表示して終了

### AG 符号側

  オプション         意味
  ------------------ -------------------------------------------------
  `--curve NAME`     `bbgs` / `hermitian` / `suzuki`
  `--q Q`            曲線パラメータ
  `--a A --b B`      BBGS パラメータ
  `--terms SPEC`     平面曲線の項を指定
  `--pole SPEC`      多点極因子を指定
  `--M M`            シンドローム境界・極次数
  `-fr -bms -vote`   1 点復号器の選択
  `-frvote`          Feng--Rao 多数決復号経路
  `-basic`           多点 Niederreiter 型シンドローム復号テスト
  `-encode`          `-basic` と併用し、符号語生成も含めて復号を検証
  `-limits`          ビルド・曲線のパラメータ上限を表示

## 実験結果の読み方

このプログラムは、指定された誤りをランダムに注入し、復号結果を集計する実験プログラムです。

-   試行回数を増やすと、ランダム実験の成功率をより多くの試行で確認できます。
-   `--errors`
    を増やした結果が成功しても、それだけで一般的な訂正保証が証明されたことにはなりません。
-   理論上の最小距離・訂正能力と、実験で観測した成功率は分けて評価してください。
-   実行時間を比較するときは、同じコンパイラ設定・同じパラメータ・同じ試行回数で測定してください。
-   AG
    符号では、曲線・有限体・極因子・基底の構成条件によって利用可能なパラメータが変わります。

## 実装上の注意

-   本プログラムは研究・実験用です。暗号用途や、本番システムの信頼性を保証するものではありません。
-   AG
    側と古典符号側は、単一の実行ファイルに統合されていますが、内部の実装は独立しています。
-   `-fht` は `RM(1,m)` 用です。一般の `RM(r,m)`（`r > 1`）に FWHT
    最尤復号を適用する指定ではありません。
-   AG 側の `-encode` は `-basic` と組み合わせて使用します。
-   メモリ上限・有限体次数・曲線の条件などはソース内の実装制限に従います。`-limits`
    と `--help-all` を確認してください。

## ライセンス

このリポジトリにライセンスを設定する場合は、利用条件に合ったライセンスファイル（例：`LICENSE`）を別途追加してください。現時点でこの
README は特定のオープンソースライセンスを宣言しません。

## 開発・検証状況

この README は `unified_codes.c`
のコマンドラインヘルプと実装内容に基づいています。コンパイル可否や各符号の復号性能は、利用する環境・パラメータごとに個別に検証してください。未実施のベンチマークや未確認の訂正性能を保証するものではありません。

## ファイル入出力（古典二元符号）

ファイル入出力は、現在 `bch` / `golay` / `goppa` / `alt` / `rm` /
`hadamard` に対応しています。`grs` / `rank`（Gabidulin）/ `product` /
`concat` および AG 符号側は、このファイル入出力モードの対象外です。

入力ファイルはメッセージビット列として読み込まれ、符号化、必要に応じた誤り付加、復号を行います。

| オプション | 意味 |
|---|---|
| `--input FILE` | 入力ファイル |
| `--input-format binary` | 任意のバイナリデータとして読む（既定値） |
| `--input-format text` | `0` と `1` からなるビット列テキストとして読む。空白・改行は無視 |
| `--no-errors` | 誤りを付加しない |
| `--errors N` | 各符号ブロックに N 個のランダムなビット誤りを付加 |
| `--output FILE` | 復号結果の出力先 |
| `--output-format binary` | 復号結果をバイナリで保存 |
| `--output-format text` | 復号結果を `0`/`1` のテキストで保存 |
| `--output-format screen` | 復号結果を画面に表示（既定値） |

バイナリ入力は各バイトの最上位ビットから順にビット列へ展開します。入力ビット列は符号の情報長 `k` ごとに分割されます。最後のブロックが `k` ビットに満たない場合は、符号化時に末尾を 0 で埋め、復号後は元の入力ビット長に切り詰めます。そのため、バイナリ入力では誤りなしの往復時に元のファイル長が保たれます。テキスト入力は `0` と `1` および空白文字のみを受け付けます。テキスト入力を `binary` 出力する場合、最終バイトの未使用ビットは 0 で埋められます。

`--errors` を省略したファイルモードは誤りなし（0 個）です。`--errors N` はファイル全体ではなく、符号化した各ブロックに N 個ずつ誤りを付加します。現段階ではファイルモードの消失訂正（`--erasures`）は未対応です。`-n` もファイルモードでは使用しません。各ブロックの復号成功・検出失敗・誤訂正の件数と処理時間は標準エラー出力に表示されます。

### 例1：バイナリファイルを誤りなしで往復

```sh
./unified_codes --family classical --code golay \
  --input input.bin --input-format binary \
  --no-errors \
  --output decoded.bin --output-format binary
cmp input.bin decoded.bin
```

`cmp` が何も表示せず終了コード 0 なら、入力と復号結果はバイト単位で一致しています。

### 例2：バイナリファイルに各ブロック2個の誤りを付加

```sh
./unified_codes --family classical --code golay \
  --input input.bin --input-format binary \
  --errors 2 \
  --output decoded.bin --output-format binary
cmp input.bin decoded.bin
```

### 例3：ビット列テキストを読み込んで復号結果をテキスト保存

入力 `message.txt` の例：

```text
01010110
11001010
```

実行：

```sh
./unified_codes --family classical --code bch --m 5 --t 2 \
  --input message.txt --input-format text \
  --errors 2 \
  --output decoded.txt --output-format text
```

### 例4：復号結果を画面に表示

```sh
./unified_codes --family classical --code rm --r 1 --m 5 \
  --input message.txt --input-format text \
  --no-errors --output-format screen
```

画面表示の場合は `--output FILE` を指定しません。

### 例5：ほかの符号のパラメータを確認

```sh
./unified_codes --family classical --code golay --list
./unified_codes --family classical --code bch --m 6 --t 3 --list
./unified_codes --family classical --code rm --r 2 --m 6 --list
```

## 積符号・Gabidulin符号の復号

### 積符号 (`--code product`)

積符号は、二つの二元符号 `C1` と `C2` の直積として構成します。配列の各行に `C1`、各列に `C2` の制約を課し、行方向・列方向の有界距離復号を交互に繰り返します。

```sh
# Golay x Golay。各試行で24個のランダムなビット誤りを加える
./unified_codes --family classical --code product \
  --c1 golay --c2 golay --errors 24 --iters 10 -n 500 123

# BCH x Golay
./unified_codes --family classical --code product \
  --c1 bch:5:2 --c2 golay --errors 8 --iters 10 -n 100 123
```

`--iters` は行・列の復号を交互に行う最大スイープ回数です。最終段階で全行が `C1`、全列が `C2` の符号語になっているかを再検査します。出力に表示される `code-theoretic unique-decoding radius` は、積符号の最小距離 `d` から計算した理論上の半径 `floor((d-1)/2)` です。これは、この行・列の反復復号器がその半径内で常に成功するという保証ではありません。反復復号後の制約検査を通過しても、別の積符号語への誤訂正は理論上あり得ます。

### Gabidulin符号 (`--code rank`)

Gabidulin符号は `GF(2^m)` 上のランク距離符号です。長さ `n` は `n <= m` を満たす必要があります。復号器はシンドロームから線形化誤り位置多項式の係数を求め、その根空間とパリティ検査行列を用いて誤りを復元します。

```sh
# GF(2^8)、長さ6、次元2。ランク重み2の誤りを訂正
./unified_codes --family classical --code rank \
  --m 8 --len 6 --k 2 --errors 2 -n 1000 123

# GF(2^12)、長さ12、次元6
./unified_codes --family classical --code rank \
  --m 12 --len 12 --k 6 --errors 3 -n 100 123
```

理論上の一意復号半径はランク重み `t = floor((n-k)/2)` です。`--errors` はハミング重みではなく、誤りベクトルのランク重みを指定します。訂正半径を超える誤りについては、復号失敗だけでなく誤訂正も起こり得ます。

### 実装上の確認

- 積符号：反復復号後に全行・全列の符号制約を検査します。
- Gabidulin：シンドローム復号後に全パリティ検査式を再評価します。
- これら二つの方式は現在、乱数による符号語生成・誤り付加・復号実験に対応します。汎用ファイル入出力にはまだ接続していません。
