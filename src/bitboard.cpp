/*
  Stockfish, a UCI chess playing engine derived from Glaurung 2.1
  Copyright (C) 2004-2008 Tord Romstad (Glaurung author)
  Copyright (C) 2008-2014 Marco Costalba, Joona Kiiski, Tord Romstad

  Stockfish is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  Stockfish is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <algorithm>
#include <cstring> // For memset

#include "bitboard.h"
#include "bitcount.h"
#include "rkiss.h"

CACHE_LINE_ALIGNMENT

Bitboard RMasks[SQUARE_NB];
Bitboard RMagics[SQUARE_NB];
Bitboard* RAttacks[SQUARE_NB];
unsigned RShifts[SQUARE_NB];

Bitboard BMasks[SQUARE_NB];
Bitboard BMagics[SQUARE_NB];
Bitboard* BAttacks[SQUARE_NB];
unsigned BShifts[SQUARE_NB];

Bitboard SquareBB[SQUARE_NB];
Bitboard FileBB[FILE_NB];
Bitboard RankBB[RANK_NB];
Bitboard AdjacentFilesBB[FILE_NB];
Bitboard InFrontBB[COLOR_NB][RANK_NB];
/*
非飛び駒（bishop,rook,queenを除く駒）が駒をとれる位置に
ビットが立ったbitboardの配列
*/
Bitboard StepAttacksBB[PIECE_NB][SQUARE_NB];
Bitboard BetweenBB[SQUARE_NB][SQUARE_NB];
Bitboard LineBB[SQUARE_NB][SQUARE_NB];
Bitboard DistanceRingsBB[SQUARE_NB][8];
Bitboard ForwardBB[COLOR_NB][SQUARE_NB];
Bitboard PassedPawnMask[COLOR_NB][SQUARE_NB];
Bitboard PawnAttackSpan[COLOR_NB][SQUARE_NB];
/*
ROOK,BISHOP,QUEENの座標ごとの利きを入れておく
他に駒がない状態の時の利き
*/
Bitboard PseudoAttacks[PIECE_TYPE_NB][SQUARE_NB];

int SquareDistance[SQUARE_NB][SQUARE_NB];

namespace {

  // De Bruijn sequences. See chessprogramming.wikispaces.com/BitScan

	const uint64_t DeBruijn_64 = 0x3F79D71B4CB0A89ULL;
	const uint32_t DeBruijn_32 = 0x783A9B23;

	CACHE_LINE_ALIGNMENT

	int MS1BTable[256];
	Square BSFTable[SQUARE_NB];
	Bitboard RTable[0x19000]; // Storage space for rook attacks
	Bitboard BTable[0x1480];  // Storage space for bishop attacks

	typedef unsigned (Fn)(Square, Bitboard);

	void init_magics(Bitboard table[], Bitboard* attacks[], Bitboard magics[],
	Bitboard masks[], unsigned shifts[], Square deltas[], Fn index);

	/*
	bit scanしているが仕組みがよくわからん
	ここに解説らしきものが書かれているが読めない
	https://chessprogramming.wikispaces.com/De+Bruijn+Sequence+Generator
	*/
	FORCE_INLINE unsigned bsf_index(Bitboard b) {

		// Matt Taylor's folding for 32 bit systems, extended to 64 bits by Kim Walisch
		b ^= (b - 1);
		return Is64Bit ? (b * DeBruijn_64) >> 58
			: ((unsigned(b) ^ unsigned(b >> 32)) * DeBruijn_32) >> 26;
	}
}

/// lsb()/msb() finds the least/most significant bit in a non-zero bitboard.
/// pop_lsb() finds and clears the least significant bit in a non-zero bitboard.

#ifndef USE_BSFQ

Square lsb(Bitboard b) { return BSFTable[bsf_index(b)]; }

Square pop_lsb(Bitboard* b) {

  Bitboard bb = *b;
  *b = bb & (bb - 1);
  return BSFTable[bsf_index(bb)];
}

Square msb(Bitboard b) {

  unsigned b32;
  int result = 0;

  if (b > 0xFFFFFFFF)
  {
      b >>= 32;
      result = 32;
  }

  b32 = unsigned(b);

  if (b32 > 0xFFFF)
  {
      b32 >>= 16;
      result += 16;
  }

  if (b32 > 0xFF)
  {
      b32 >>= 8;
      result += 8;
  }

  return Square(result + MS1BTable[b32]);
}

#endif // ifndef USE_BSFQ


/// Bitboards::pretty() returns an ASCII representation of a bitboard to be
/// printed to standard output. This is sometimes useful for debugging.
/*
引数のbiboardを表示する
ただしstackfishの内部盤データは行が通常のchess盤と上下逆さまになっている
のでここで元に戻している
＜表示例＞
+---+---+---+---+---+---+---+---+
|   |   |   |   |   |   |   |   |
+---+---+---+---+---+---+---+---+
|   |   |   |   |   |   |   |   |
+---+---+---+---+---+---+---+---+
|   |   |   |   |   |   |   |   |
+---+---+---+---+---+---+---+---+
|   |   |   |   |   |   |   |   |
+---+---+---+---+---+---+---+---+
|   |   |   |   |   |   |   |   |
+---+---+---+---+---+---+---+---+
| X | X | X |   |   |   |   |   |
+---+---+---+---+---+---+---+---+
|   |   | X |   |   |   |   |   |
+---+---+---+---+---+---+---+---+
|   |   | X |   |   |   |   |   |
+---+---+---+---+---+---+---+---+
*/
const std::string Bitboards::pretty(Bitboard b) {

  std::string s = "+---+---+---+---+---+---+---+---+\n";

  for (Rank r = RANK_8; r >= RANK_1; --r)
  {
      for (File f = FILE_A; f <= FILE_H; ++f)
          s.append(b & make_square(f, r) ? "| X " : "|   ");

      s.append("|\n+---+---+---+---+---+---+---+---+\n");
  }

  return s;
}


/// Bitboards::init() initializes various bitboard tables. It is called at
/// startup and relies on global objects to be already zero-initialized.

void Bitboards::init() {

  for (Square s = SQ_A1; s <= SQ_H8; ++s)
      BSFTable[bsf_index(SquareBB[s] = 1ULL << s)] = s;

  for (Bitboard b = 1; b < 256; ++b)
      MS1BTable[b] = more_than_one(b) ? MS1BTable[b - 1] : lsb(b);
  /*
  FILE_A
  +---+---+---+---+---+---+---+---+
  | X |   |   |   |   |   |   |   |
  +---+---+---+---+---+---+---+---+
  | X |   |   |   |   |   |   |   |
  +---+---+---+---+---+---+---+---+
  | X |   |   |   |   |   |   |   |
  +---+---+---+---+---+---+---+---+
  | X |   |   |   |   |   |   |   |
  +---+---+---+---+---+---+---+---+
  | X |   |   |   |   |   |   |   |
  +---+---+---+---+---+---+---+---+
  | X |   |   |   |   |   |   |   |
  +---+---+---+---+---+---+---+---+
  | X |   |   |   |   |   |   |   |
  +---+---+---+---+---+---+---+---+
  | X |   |   |   |   |   |   |   |
  +---+---+---+---+---+---+---+---+  */
  for (File f = FILE_A; f <= FILE_H; ++f)
      FileBB[f] = f > FILE_A ? FileBB[f - 1] << 1 : FileABB;
  /*
  RANK1
  +---+---+---+---+---+---+---+---+
  | X | X | X | X | X | X | X | X |
  +---+---+---+---+---+---+---+---+
  |   |   |   |   |   |   |   |   |
  +---+---+---+---+---+---+---+---+
  |   |   |   |   |   |   |   |   |
  +---+---+---+---+---+---+---+---+
  |   |   |   |   |   |   |   |   |
  +---+---+---+---+---+---+---+---+
  |   |   |   |   |   |   |   |   |
  +---+---+---+---+---+---+---+---+
  |   |   |   |   |   |   |   |   |
  +---+---+---+---+---+---+---+---+
  |   |   |   |   |   |   |   |   |
  +---+---+---+---+---+---+---+---+
  |   |   |   |   |   |   |   |   |
  +---+---+---+---+---+---+---+---+
  */
  for (Rank r = RANK_1; r <= RANK_8; ++r)
      RankBB[r] = r > RANK_1 ? RankBB[r - 1] << 8 : Rank1BB;
  /*
  FILE_B
  隣の列のbitboardを返す
  +---+---+---+---+---+---+---+---+
  | X |   | X |   |   |   |   |   |
  +---+---+---+---+---+---+---+---+
  | X |   | X |   |   |   |   |   |
  +---+---+---+---+---+---+---+---+
  | X |   | X |   |   |   |   |   |
  +---+---+---+---+---+---+---+---+
  | X |   | X |   |   |   |   |   |
  +---+---+---+---+---+---+---+---+
  | X |   | X |   |   |   |   |   |
  +---+---+---+---+---+---+---+---+
  | X |   | X |   |   |   |   |   |
  +---+---+---+---+---+---+---+---+
  | X |   | X |   |   |   |   |   |
  +---+---+---+---+---+---+---+---+
  | X |   | X |   |   |   |   |   |
  +---+---+---+---+---+---+---+---+  */
  for (File f = FILE_A; f <= FILE_H; ++f)
      AdjacentFilesBB[f] = (f > FILE_A ? FileBB[f - 1] : 0) | (f < FILE_H ? FileBB[f + 1] : 0);
  /*
  InFrontBB[WHITE][RANK_2]
  この絵は下側がWHITEでRANK_2は下から２番目の行となっている
  つまりInFrontBBは指定した行から敵陣側のbitboardを返す
  +---+---+---+---+---+---+---+---+
  | X | X | X | X | X | X | X | X |
  +---+---+---+---+---+---+---+---+
  | X | X | X | X | X | X | X | X |
  +---+---+---+---+---+---+---+---+
  | X | X | X | X | X | X | X | X |
  +---+---+---+---+---+---+---+---+
  | X | X | X | X | X | X | X | X |
  +---+---+---+---+---+---+---+---+
  | X | X | X | X | X | X | X | X |
  +---+---+---+---+---+---+---+---+
  | X | X | X | X | X | X | X | X |
  +---+---+---+---+---+---+---+---+
  |   |   |   |   |   |   |   |   |
  +---+---+---+---+---+---+---+---+
  |   |   |   |   |   |   |   |   |
  +---+---+---+---+---+---+---+---+
  */
  for (Rank r = RANK_1; r < RANK_8; ++r)
      InFrontBB[WHITE][r] = ~(InFrontBB[BLACK][r + 1] = InFrontBB[BLACK][r] | RankBB[r]);
  /*
  ForwardBB
	自分の位置から前方（Forward）のbitboardを返す
	＜例＞
	color=0 sq=0
	 +---+---+---+---+---+---+---+---+
	| X |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	| X |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	| X |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	| X |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	| X |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	| X |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	| X |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	| s |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+  
	color=0 sq=10
	+---+---+---+---+---+---+---+---+
	|   |   | X |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   | X |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   | X |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   | X |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   | X |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   | X |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   | s |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	color=1(BLACK) sq=41
	カラーが変わるとFORWORDの方向が変わる
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   | S |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   | X |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   | X |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   | X |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   | X |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   | X |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+	
	PawnAttackSpanの確認
	＜例＞
	  color=0 sq=1
	  PAWNが駒を取れる位置のbitboardを返す
	 +---+---+---+---+---+---+---+---+
	| X |   | X |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	| X |   | X |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	| X |   | X |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	| X |   | X |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	| X |   | X |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	| X |   | X |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	| X |   | X |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   | s |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	カラーが変わると利きの方向も変わる
	color=1 sq=17
	 +---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   | s |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	| X |   | X |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	| X |   | X |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+	
	PassedPawnMaskの確認
	PassedPawnMaskはForwardBBビットボードとPawnAttackSpanビットボードのOR演算なので
	PAWN駒の３列にbitboardが立つことになる。用途不明
	＜例＞
	color=0 sq=2
	 +---+---+---+---+---+---+---+---+
	|   | X | X | X |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   | X | X | X |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   | X | X | X |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   | X | X | X |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   | X | X | X |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   | X | X | X |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   | X | X | X |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   | s |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+	
	*/
  for (Color c = WHITE; c <= BLACK; ++c)
      for (Square s = SQ_A1; s <= SQ_H8; ++s)
      {
          ForwardBB[c][s]      = InFrontBB[c][rank_of(s)] & FileBB[file_of(s)];
          PawnAttackSpan[c][s] = InFrontBB[c][rank_of(s)] & AdjacentFilesBB[file_of(s)];
          PassedPawnMask[c][s] = ForwardBB[c][s] | PawnAttackSpan[c][s];
      }
	/*
	SquareDistance変数の初期化
		SquareDistance変数には座標s1とs2の距離を入れている
		但し列と行の距離の大きい方をいれる。
	DistanceRingsBB変数の初期化
		DistanceRingsBBはよくわからない

	*/
  for (Square s1 = SQ_A1; s1 <= SQ_H8; ++s1)
      for (Square s2 = SQ_A1; s2 <= SQ_H8; ++s2)
          if (s1 != s2)
          {
              SquareDistance[s1][s2] = std::max(file_distance(s1, s2), rank_distance(s1, s2));
              DistanceRingsBB[s1][SquareDistance[s1][s2] - 1] |= s2;
          }

  int steps[][9] = { {}, { 7, 9 }, { 17, 15, 10, 6, -6, -10, -15, -17 },
                     {}, {}, {}, { 9, 7, -7, -9, 8, 1, -1, -8 } };
	/*
	StepAttacksBB変数の初期設定
	宣言 Bitboard StepAttacksBB[PIECE_NB][SQUARE_NB];
	非飛び駒（bishop,rook,queenを除く駒）が駒をとれる位置に
	ビットが立ったbitboardの配列
	<例>
	piece=W_KING,sq=SQ_B2(内部の盤データ構造と上下が
	逆になっていることに注意
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	| X | X | X |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	| X |   | X |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	| X | X | X |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	*/
	for (Color c = WHITE; c <= BLACK; ++c)
		for (PieceType pt = PAWN; pt <= KING; ++pt)
			for (Square s = SQ_A1; s <= SQ_H8; ++s)
				for (int i = 0; steps[pt][i]; ++i)
				{
					Square to = s + Square(c == WHITE ? steps[pt][i] : -steps[pt][i]);

					if (is_ok(to) && square_distance(s, to) < 3)
						StepAttacksBB[make_piece(c, pt)][s] |= to;
				}

	Square RDeltas[] = { DELTA_N,  DELTA_E,  DELTA_S,  DELTA_W  };
	Square BDeltas[] = { DELTA_NE, DELTA_SE, DELTA_SW, DELTA_NW };

	init_magics(RTable, RAttacks, RMagics, RMasks, RShifts, RDeltas, magic_index<ROOK>);
	init_magics(BTable, BAttacks, BMagics, BMasks, BShifts, BDeltas, magic_index<BISHOP>);

	/*
	PseudoAttacks[ROOK],PseudoAttacks[BISHOP],PseudoAttacks[QUEEN]に他の駒が
	ない状態の利きbitboardをいれる。

	PseudoAttacks[BISHOP][s1] & s2　->　bitboardとsquareのAND演算は演算子オーバーロードに
	定義してある。（bitboard.hに定義してある　PseudoAttacks[BISHOP][s1] & s2は
	PseudoAttacks[BISHOP][s1] & SquareBB[s2]となる）

	LineBB[s1][s2]はs1,s2のあいだにRook,bishopの利きがあるのbitboardを返す
	分かりにくいので例を示す
	<例>
	s1=SQ_A1,s2=SQ_B2のLineBBは
	駒の配置
	+---+---+---+---+---+---+---+---+
	|s1 |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |s2 |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	返すbitboard
	+---+---+---+---+---+---+---+---+
	| X |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   | X |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   | X |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   | X |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   | X |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   | X |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   | X |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   | X |
	+---+---+---+---+---+---+---+---+

	駒の配置
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |s1 |   |   |   |s2 |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	返すbitboard
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	| X | X | X | X | X | X | X | X |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+

	駒の配置
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |s1 |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |S2 |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	返すbitboard,s1とs2の間に利きが（Rook,bishop)が
	存在しない場合は空のbitboardを返す
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+


	BetweenBB[s1][s2]はs1,s2のあいだにRook,bishopの利きがあるのbitboardを返す
	分かりにくいので例を示す
	<例>
	s1=SQ_A1,s2=SQ_B2のBetweenBBは
	駒の配置
	+---+---+---+---+---+---+---+---+
	|s1 |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |s2 |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	返すbitboards1とs2の間には遠方利きがないでの
	駒同士には利きが利いているが空きの升には
	利きがないのでこのようなbitboardをかえす
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+

	駒の配置
	+---+---+---+---+---+---+---+---+
	|s1 |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |s2 |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	s1とs2の間にある利きを返す
	Betweenは何何のあいだという意味がある
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   | x |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+

	駒の配置
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |s1 |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |s2 |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	このような駒配置の場合rook,bishopの利きが
	存在しないのでbitboardもbitがたたない（表示省略）

	*/
	for (Square s1 = SQ_A1; s1 <= SQ_H8; ++s1)
	{
		PseudoAttacks[QUEEN][s1]  = PseudoAttacks[BISHOP][s1] = attacks_bb<BISHOP>(s1, 0);
		PseudoAttacks[QUEEN][s1] |= PseudoAttacks[  ROOK][s1] = attacks_bb<  ROOK>(s1, 0);

		for (Square s2 = SQ_A1; s2 <= SQ_H8; ++s2)
		{
			Piece pc = (PseudoAttacks[BISHOP][s1] & s2) ? W_BISHOP :
			(PseudoAttacks[ROOK][s1]   & s2) ? W_ROOK   : NO_PIECE;

			if (pc == NO_PIECE)
				continue;

			LineBB[s1][s2] = (attacks_bb(pc, s1, 0) & attacks_bb(pc, s2, 0)) | s1 | s2;
			BetweenBB[s1][s2] = attacks_bb(pc, s1, SquareBB[s2]) & attacks_bb(pc, s2, SquareBB[s1]);
		}
	}
}


namespace {

  Bitboard sliding_attack(Square deltas[], Square sq, Bitboard occupied) {

    Bitboard attack = 0;

    for (int i = 0; i < 4; ++i)
        for (Square s = sq + deltas[i];
             is_ok(s) && square_distance(s, s - deltas[i]) == 1;
             s += deltas[i])
        {
            attack |= s;

            if (occupied & s)
                break;
        }

    return attack;
  }


  // init_magics() computes all rook and bishop attacks at startup. Magic
  // bitboards are used to look up attacks of sliding pieces. As a reference see
  // chessprogramming.wikispaces.com/Magic+Bitboards. In particular, here we
  // use the so called "fancy" approach.
	/*
	table=>RTable[0x19000]配列はなぜこんなに大きな配列なのか
	magic tableなので sにRookがいる場合 X に駒がある/ないのパターンは
	6bit＝64とおり（最後の升目は不要）なので
	横×縦＝64*64=4096とおり、これが４隅あるので
	4096*4=16,384
	+---+---+---+---+---+---+---+---+
	| S | X | X | X | X | X | X | ○ |
	+---+---+---+---+---+---+---+---+
	| X |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	| X |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	| X |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	| X |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	| X |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	| X |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	| ○ |   |   |   |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	辺に駒がいた場合	
	横方向は5bitなので32とおり縦は6bitのなので
	64とおり,32*64=2,048
	辺に駒をおける場所（隅を除く）２４か所なので
	2048*24=49,152とおり
	+---+---+---+---+---+---+---+---+
	| ○ | X | X | S | X | X | X | ○ |
	+---+---+---+---+---+---+---+---+
	|   |   |   | X |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   | X |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   | X |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   | X |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   | X |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   | X |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   | ○ |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	隅、辺を除く升目は36
	縦、横とも5bitなので32とおり
	32*32*36=36,864
	+---+---+---+---+---+---+---+---+
	|   |   |   | ○ |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   | X |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   | X |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	| ○ | X | X | s | X | X | X | ○ |
	+---+---+---+---+---+---+---+---+
	|   |   |   | X |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   | X |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   | X |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	|   |   |   | ○ |   |   |   |   |
	+---+---+---+---+---+---+---+---+
	16,384+49,152+36,864=102,400
	これを16進数に変換すると0x19000となる

	bishopの配列はBitboard BTable[0x1480]となって
	Rook程ではないが大きい配列
	bishopは座標によって利き数が大きく異なる
	升目に直接入れた数字がその升目での利き数
	+---+---+---+---+---+---+---+---+
	|128|///| 32|///| 32|///| 32|///|
	+---+---+---+---+---+---+---+---+
	|///| 32|///| 32|///| 32|///| 32|
	+---+---+---+---+---+---+---+---+
	| 32|///|128|///|128|///| 32|///|
	+---+---+---+---+---+---+---+---+
	|///| 32|///|512|///|128|///| 32|
	+---+---+---+---+---+---+---+---+
	| 32|///|128|///|512|///| 32|///|
	+---+---+---+---+---+---+---+---+
	|///| 32|///|128|///|128|///| 32|
	+---+---+---+---+---+---+---+---+
	| 32|///| 32|///| 32|///| 32|///|
	+---+---+---+---+---+---+---+---+
	|///| 32|///| 32|///| 32|///|128|
	+---+---+---+---+---+---+---+---+
	32が22＝704
	64が2＝128
	128が6＝768
	512が2＝1024
	2624となる
	色違いがあるので
	2624*2=5248=>１６進法に変換すると0x1480
	*/
  void init_magics(Bitboard table[], Bitboard* attacks[], Bitboard magics[],
                   Bitboard masks[], unsigned shifts[], Square deltas[], Fn index) {

    int MagicBoosters[][8] = { {  969, 1976, 2850,  542, 2069, 2852, 1708,  164 },
                               { 3101,  552, 3555,  926,  834,   26, 2131, 1117 } };

    RKISS rk;
    Bitboard occupancy[4096], reference[4096], edges, b;
    int i, size, booster;

    // attacks[s] is a pointer to the beginning of the attacks table for square 's'
    attacks[SQ_A1] = table;

    for (Square s = SQ_A1; s <= SQ_H8; ++s)
    {
        // Board edges are not considered in the relevant occupancies
        edges = ((Rank1BB | Rank8BB) & ~rank_bb(s)) | ((FileABB | FileHBB) & ~file_bb(s));

        // Given a square 's', the mask is the bitboard of sliding attacks from
        // 's' computed on an empty board. The index must be big enough to contain
        // all the attacks for each possible subset of the mask and so is 2 power
        // the number of 1s of the mask. Hence we deduce the size of the shift to
        // apply to the 64 or 32 bits word to get the index.
        masks[s]  = sliding_attack(deltas, s, 0) & ~edges;
        shifts[s] = (Is64Bit ? 64 : 32) - popcount<Max15>(masks[s]);

        // Use Carry-Rippler trick to enumerate all subsets of masks[s] and
        // store the corresponding sliding attack bitboard in reference[].
        b = size = 0;
        do {
            occupancy[size] = b;
            reference[size] = sliding_attack(deltas, s, b);

            if (HasPext)
                attacks[s][_pext_u64(b, masks[s])] = reference[size];

            size++;
            b = (b - masks[s]) & masks[s];
        } while (b);

        // Set the offset for the table of the next square. We have individual
        // table sizes for each square with "Fancy Magic Bitboards".
        if (s < SQ_H8)
            attacks[s + 1] = attacks[s] + size;

        if (HasPext)
            continue;

        booster = MagicBoosters[Is64Bit][rank_of(s)];

        // Find a magic for square 's' picking up an (almost) random number
        // until we find the one that passes the verification test.
        do {
            do magics[s] = rk.magic_rand<Bitboard>(booster);
            while (popcount<Max15>((magics[s] * masks[s]) >> 56) < 6);

            std::memset(attacks[s], 0, size * sizeof(Bitboard));

            // A good magic must map every possible occupancy to an index that
            // looks up the correct sliding attack in the attacks[s] database.
            // Note that we build up the database for square 's' as a side
            // effect of verifying the magic.
            for (i = 0; i < size; ++i)
            {
                Bitboard& attack = attacks[s][index(s, occupancy[i])];

                if (attack && attack != reference[i])
                    break;

                assert(reference[i]);

                attack = reference[i];
            }
        } while (i < size);
    }
  }
}
