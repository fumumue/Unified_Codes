# unified_codes — 統合誤り訂正符号実験ツール

## 1. 概要

`unified_codes` は、代数幾何符号（AG符号）と古典的な誤り訂正符号を、単一のCプログラムで生成・符号化・復号・評価するための実験用ツールです。

主な目的は次のとおりです。

- 複数の符号構成法を同じ環境で比較する
- 符号長・次元・符号化率・訂正能力を調べる
- 復号アルゴリズムの正確性と処理速度を測定する
- ランダム誤りを加え、復号結果を検証する
- ファイル単位で符号化・誤り付加・復号を行う
- 符号の理論的な訂正能力と実際の復号性能を比較する

本プログラムは研究・実験用です。すべての符号や復号方式が、実用的な暗号方式として安全性を保証するものではありません。

---

## 2. コンパイル

基本的なコンパイル方法：

```bash
gcc -O3 -march=native -mtune=native \
    -std=gnu11 -fopenmp -Wall -Wextra \
    unified_codes.c -lm -o unified_codes
```

OpenMPが利用できない環境では `-fopenmp` を外してください。

```bash
gcc -O3 -march=native -mtune=native \
    -std=gnu11 -Wall -Wextra \
    unified_codes.c -lm -o unified_codes
```

実行ファイルのヘルプ：

```bash
./unified_codes --help
```

全オプション：

```bash
./unified_codes --help-all
```

---

## 3. モジュール構成

本プログラムは大きく2つの符号群に分かれています。

### 3.1 AG符号系

代数幾何符号と、そのシンドローム復号に関する実験機能です。

対象：

- BBGS曲線
- Hermitian曲線
- Suzuki曲線
- Klein曲線
- 独自平面曲線
- 一点AG符号
- 多点AG符号

主な復号方式：

| 方式 | 概要 |
|---|---|
| FR | Feng-Rao系の復号 |
| BMS | BMS型のシンドローム復号 |
| VOTE | 多数決・投票を利用する復号 |
| FR-VOTE | Feng-Raoカーネルを用いた投票型復号 |
| ECP | Error-Correcting Pairによる復号 |
| BASIC | 多点AG符号の基本的なシンドローム復号 |

### 3.2 古典符号系

次の符号を扱います。

| 符号 | 概要 |
|---|---|
| BCH | BCH符号 |
| Golay | Golay符号 |
| Goppa | 二元Goppa符号 |
| Alternant | Alternant符号 |
| Srivastava | Srivastava符号 |
| GS | Generalized Srivastava系 |
| RM | Reed-Muller符号 |
| dRM | RM符号の反復構成 |
| Hadamard | Hadamard符号 |
| LDPC | Gallager型LDPC符号 |
| RS | Reed-Solomon符号 |
| GRS | Generalized Reed-Solomon符号 |
| Gabidulin | Rank metric符号 |
| Product | 積符号 |
| Concatenated | 連接符号 |

---

## 4. 実行モード

### 4.1 AG符号系

AG符号系は標準では次の形式で実行します。

```bash
./unified_codes --family ag [options]
```

例：

```bash
./unified_codes --family ag --curve bbgs
```

### 4.2 古典符号系

```bash
./unified_codes --family classical [options]
```

また、`--code` を指定すると古典符号系が自動選択されます。

```bash
./unified_codes --code bch --m 6 --t 3
```

---

# 5. AG符号系の使い方

## 5.1 曲線の選択

```bash
--curve bbgs
--curve hermitian
--curve suzuki
--curve klein
```

例：

```bash
./unified_codes --family ag --curve bbgs
```

## 5.2 有限体の指定

```bash
--m 6
--poly 0x43
```

`--m` は基礎有限体 GF(2^m) の次数です。

`--poly` は有限体の既約多項式・原始多項式の指定に使います。指定可能な値は実装上の検証条件に従います。

## 5.3 BBGSパラメータ

```bash
--q auto
--q 8
--a 1
--b 2
```

BBGS曲線では、互いに素な正整数 a,b を指定します。

\[
c=a+b
\]

有限体は GF(q^c) となります。

例：

```bash
./unified_codes --family ag \
    --curve bbgs --q 8 --a 1 --b 2
```

## 5.4 独自平面曲線

```bash
--terms "16,0;0,1;15,4;12,16"
```

各項は \(X^iY^j\) を表します。

係数も指定する場合：

```bash
--terms "2,0,1;0,1,1;1,2,3"
```

具体的な曲線の次数や係数は、有限体と曲線の条件に合わせて指定してください。

## 5.5 シンドローム復号パラメータ

```bash
--errors 10
--M 100
--place 0
--len 100
```

| オプション | 意味 |
|---|---|
| `--errors t` | 注入する誤りの個数 |
| `--M M` | シンドロームの上限・極次数 |
| `--place k` | 使用する場所のインデックス |
| `--len L` | 使用する評価点を先頭L個に制限 |
| `-half` | 評価点を半分に制限 |
| `--fixed-error-value v` | 誤り値を同じ非零有限体要素に固定 |

`--M` は復号可能性や基底サイズなどに関係する重要なパラメータです。

## 5.6 多点AG符号

```bash
--pole "E2:140,E1:85"
```

または、

```bash
--pole "#3:100"
```

`E<k>` は辺、`#k` は軌道を表します。

基底の列順序：

```bash
--basis legacy
--basis newton
```

デフォルトは `legacy` です