# 使用方法
- 標準入力に以下の形式でデータ入力

  widthの種類数

  width[um] β_even β_odd β_1 β_2

  ...


- 伝搬方向zに対する複素振幅を標準エラーに出力

- コンパイルに当たって，ライブラリ"Eigen"が必要

  事前にダウンロード & パスを通しておくこと

- "make all"でmake clean & make

- サンプル例

  ./cmt <input_1.550.pre 2>output
  
  で、実行

- 計算結果の妥当性チェックはしてある(落合のMDλテーパADCで確認済み)

- LI.classで伝搬定数の線形補間

# 今後
- CMT未完成

- 今後は、「Structure」と「複素振幅の初期状態」と「伝搬定数」の三つを入力して、resultに、「zに対する複素振幅」と「最終の複素振幅」を出力するようにする

- Solverディレクトリからアクセスできるようにする