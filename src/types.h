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

#ifndef TYPES_H_INCLUDED
#define TYPES_H_INCLUDED

/// For Linux and OSX configuration is done automatically using Makefile. To get
/// started type 'make help'.
///
/// For Windows, part of the configuration is detected automatically, but some
/// switches need to be set manually:
///
/// -DNDEBUG      | Disable debugging mode. Always use this.
///
/// -DNO_PREFETCH | Disable use of prefetch asm-instruction. A must if you want
///               | the executable to run on some very old machines.
///
/// -DUSE_POPCNT  | Add runtime support for use of popcnt asm-instruction. Works
///               | only in 64-bit mode. For compiling requires hardware with
///               | popcnt support.

#include <cassert>
#include <cctype>
#include <climits>
#include <cstdlib>

#include "platform.h"

/*
用途不明
*/
#define unlikely(x) (x) // For code annotation purposes

/*
*/
#if defined(_WIN64) && !defined(IS_64BIT)
#  include <intrin.h> // MSVC popcnt and bsfq instrinsics
#  define IS_64BIT
#  define USE_BSFQ
#endif

/*
POPCNT命令を使いたいがINTEL COMPILERではない
*/
#if defined(USE_POPCNT) && defined(_MSC_VER) && defined(__INTEL_COMPILER)
#  include <nmmintrin.h> // Intel header for _mm_popcnt_u64() intrinsic
#endif

/*
PEXT命令を使いたいがCPUがサポートしていない
PEXT命令をとは下記URL参照のこと
BMI2 (Bit Manipulation Instruction Set 2)
http://en.wikipedia.org/wiki/Bit_Manipulation_Instruction_Sets#BMI2_.28Bit_Manipulation_Instruction_Set_2.29
*/
#if defined(USE_PEXT)
#  include <immintrin.h> // Header for _pext_u64() intrinsic
#else
#  define _pext_u64(b, m) (0)
#endif

/*
先行読みをしたいので<xmmintrin.h>を読み込んでおく
*/
#  if !defined(NO_PREFETCH) && (defined(__INTEL_COMPILER) || defined(_MSC_VER))
#   include <xmmintrin.h> // Intel and Microsoft header for _mm_prefetch()
#  endif


#define CACHE_LINE_SIZE 64
#if defined(_MSC_VER) || defined(__INTEL_COMPILER)
#  define CACHE_LINE_ALIGNMENT __declspec(align(CACHE_LINE_SIZE))
#else
#  define CACHE_LINE_ALIGNMENT  __attribute__ ((aligned(CACHE_LINE_SIZE)))
#endif

/*
inline関数の展開の仕方？
*/
#ifdef _MSC_VER
#  define FORCE_INLINE  __forceinline
#elif defined(__GNUC__)
#  define FORCE_INLINE  inline __attribute__((always_inline))
#else
#  define FORCE_INLINE  inline
#endif

/*
POPCNT命令を使いたかったらHasPopCntをtrueに
*/
#ifdef USE_POPCNT
const bool HasPopCnt = true;
#else
const bool HasPopCnt = false;
#endif

/*
PEXT命令を使いたいならHasPextをtrueに
*/
#ifdef USE_PEXT
const bool HasPext = true;
#else
const bool HasPext = false;
#endif

/*
このIS_64BITはマシンのことか。コンパイラのことか？
（追記）
OSが６４bitならtrueになる
付属のMakeFileに記述してある
*/
#ifdef IS_64BIT
const bool Is64Bit = true;
#else
const bool Is64Bit = false;
#endif

typedef uint64_t Key;
typedef uint64_t Bitboard;

/*
着手リストの最大長さ
*/
const int MAX_MOVES      = 256;
/*
最大探索深さ
*/
const int MAX_PLY        = 120;
const int MAX_PLY_PLUS_6 = MAX_PLY + 6;

/// A move needs 16 bits to be stored
///
/// bit  0- 5: destination square (from 0 to 63)
/// bit  6-11: origin square (from 0 to 63)
/// bit 12-13: promotion piece type - 2 (from KNIGHT-2 to QUEEN-2)
/// bit 14-15: special move flag: promotion (1), en passant (2), castling (3)
/// NOTE: EN-PASSANT bit is set only when a pawn can be captured
///
/// Special cases are MOVE_NONE and MOVE_NULL. We can sneak these in because in
/// any normal move destination square is always different from origin square
/// while MOVE_NONE and MOVE_NULL have the same origin and destination square.
/*
着手のデータ構造
MOVE_NONEは手がないことを表す
MOVE_NULLは用途不明
*/
enum Move {
  MOVE_NONE,
  MOVE_NULL = 65
};

/*
着手種別
  PROMOTION = 成り
  ENPASSANT = アンパッサン
  CASTLING  = キャスリング

*/
enum MoveType {
  NORMAL,
  PROMOTION = 1 << 14,
  ENPASSANT = 2 << 14,
  CASTLING  = 3 << 14
};

/*
陣営、NO_COLORは使っていない
COLOR_NBはカラーの数,配列の数を決める
*/
enum Color {
  WHITE, BLACK, NO_COLOR, COLOR_NB = 2
};

/*
キャスリングの種別（KING側にキャスリングする、Queen側にキャスリングする２種類ある）
CASTLING_SIDE_NBは種別の数
*/
enum CastlingSide {
  KING_SIDE, QUEEN_SIDE, CASTLING_SIDE_NB = 2
};

/*
用途不明
*/
enum CastlingRight {  // Defined as in PolyGlot book hash key
  NO_CASTLING,
  WHITE_OO,
  WHITE_OOO   = WHITE_OO << 1,
  BLACK_OO    = WHITE_OO << 2,
  BLACK_OOO   = WHITE_OO << 3,
  ANY_CASTLING = WHITE_OO | WHITE_OOO | BLACK_OO | BLACK_OOO,
  CASTLING_RIGHT_NB = 16
};

/*
用途不明
*/
template<Color C, CastlingSide S> struct MakeCastling {
  static const CastlingRight
  right = C == WHITE ? S == QUEEN_SIDE ? WHITE_OOO : WHITE_OO
                     : S == QUEEN_SIDE ? BLACK_OOO : BLACK_OO;
};

/*
用途不明
評価値に関するなにか？
*/
enum Phase {
  PHASE_ENDGAME,
  PHASE_MIDGAME = 128,
  MG = 0, EG = 1, PHASE_NB = 2
};

/*
用途不明
評価値に関するなにか？
*/
enum ScaleFactor {
  SCALE_FACTOR_DRAW    = 0,
  SCALE_FACTOR_ONEPAWN = 48,
  SCALE_FACTOR_NORMAL  = 64,
  SCALE_FACTOR_MAX     = 128,
  SCALE_FACTOR_NONE    = 255
};

/*
用途不明
評価値に関するなにか？
評価値の上限、下限
*/
enum Bound {
  BOUND_NONE,
  BOUND_UPPER,
  BOUND_LOWER,
  BOUND_EXACT = BOUND_UPPER | BOUND_LOWER
};

/*
用途不明
評価値に関するなにか？
*/
enum Value {
  VALUE_ZERO      = 0,
  VALUE_DRAW      = 0,
  VALUE_KNOWN_WIN = 10000,
  VALUE_MATE      = 32000,
  VALUE_INFINITE  = 32001,
  VALUE_NONE      = 32002,

  VALUE_MATE_IN_MAX_PLY  =  VALUE_MATE - MAX_PLY,
  VALUE_MATED_IN_MAX_PLY = -VALUE_MATE + MAX_PLY,

  VALUE_ENSURE_INTEGER_SIZE_P = INT_MAX,
  VALUE_ENSURE_INTEGER_SIZE_N = INT_MIN,
  /*
  駒の評価値？
  */
  PawnValueMg   = 198,   PawnValueEg   = 258,
  KnightValueMg = 817,   KnightValueEg = 846,
  BishopValueMg = 836,   BishopValueEg = 857,
  RookValueMg   = 1270,  RookValueEg   = 1278,
  QueenValueMg  = 2521,  QueenValueEg  = 2558,
  /*
  用途不明
  */
  MidgameLimit  = 15581, EndgameLimit  = 3998
};

/*
駒種
*/
enum PieceType {
  NO_PIECE_TYPE, PAWN, KNIGHT, BISHOP, ROOK, QUEEN, KING,
  ALL_PIECES = 0,
  PIECE_TYPE_NB = 8
};

/*
カラーを含めた駒コード
*/
enum Piece {
  NO_PIECE,
  W_PAWN = 1, W_KNIGHT, W_BISHOP, W_ROOK, W_QUEEN, W_KING,
  B_PAWN = 9, B_KNIGHT, B_BISHOP, B_ROOK, B_QUEEN, B_KING,
  PIECE_NB = 16
};

/*
探索深さ
この探索では１手はONE_PLY＝２
*/
enum Depth {

  ONE_PLY = 2,

  DEPTH_ZERO          =  0 * ONE_PLY,
  DEPTH_QS_CHECKS     =  0 * ONE_PLY,
  DEPTH_QS_NO_CHECKS  = -1 * ONE_PLY,
  DEPTH_QS_RECAPTURES = -5 * ONE_PLY,

  DEPTH_NONE = -127 * ONE_PLY
};

/*
座標番号
0-63
*/
enum Square {
  SQ_A1, SQ_B1, SQ_C1, SQ_D1, SQ_E1, SQ_F1, SQ_G1, SQ_H1,
  SQ_A2, SQ_B2, SQ_C2, SQ_D2, SQ_E2, SQ_F2, SQ_G2, SQ_H2,
  SQ_A3, SQ_B3, SQ_C3, SQ_D3, SQ_E3, SQ_F3, SQ_G3, SQ_H3,
  SQ_A4, SQ_B4, SQ_C4, SQ_D4, SQ_E4, SQ_F4, SQ_G4, SQ_H4,
  SQ_A5, SQ_B5, SQ_C5, SQ_D5, SQ_E5, SQ_F5, SQ_G5, SQ_H5,
  SQ_A6, SQ_B6, SQ_C6, SQ_D6, SQ_E6, SQ_F6, SQ_G6, SQ_H6,
  SQ_A7, SQ_B7, SQ_C7, SQ_D7, SQ_E7, SQ_F7, SQ_G7, SQ_H7,
  SQ_A8, SQ_B8, SQ_C8, SQ_D8, SQ_E8, SQ_F8, SQ_G8, SQ_H8,
  SQ_NONE,

  SQUARE_NB = 64,
  /*
  方向子、周りに番兵は行い配列になっている
  DELTA_NNDELTA_SS,は使用していない
  */
  DELTA_N =  8,
  DELTA_E =  1,
  DELTA_S = -8,
  DELTA_W = -1,

  DELTA_NN = DELTA_N + DELTA_N,
  DELTA_NE = DELTA_N + DELTA_E,
  DELTA_SE = DELTA_S + DELTA_E,
  DELTA_SS = DELTA_S + DELTA_S,
  DELTA_SW = DELTA_S + DELTA_W,
  DELTA_NW = DELTA_N + DELTA_W
};

/*
列
*/
enum File {
  FILE_A, FILE_B, FILE_C, FILE_D, FILE_E, FILE_F, FILE_G, FILE_H, FILE_NB
};

/*
行
*/
enum Rank {
  RANK_1, RANK_2, RANK_3, RANK_4, RANK_5, RANK_6, RANK_7, RANK_8, RANK_NB
};


/// The Score enum stores a middlegame and an endgame value in a single integer
/// (enum). The least significant 16 bits are used to store the endgame value
/// and the upper 16 bits are used to store the middlegame value. The compiler
/// is free to choose the enum type as long as it can store the data, so we
/// ensure that Score is an integer type by assigning some big int values.
/*
stackfishでは評価値はミドルゲーム（mg 中盤という意味だと思う）での駒評価値と
エンドゲーム（eg 終盤だと思う）での駒評価値と別々にある。
それを一つのint変数に格納している
32bit変数の上位16bitにミドルゲーム評価値を下位16bitにエンドゲーム評価値を入れる
そのためのデータ構造体がScoreView,そしてそのデータ構造体を作るのがmake_score関数
 (int16_t)(mg - (uint16_t(eg) >> 15));の部分が今一つ理解できない

 SCORE_ENSURE_INTEGER_SIZE_P、SCORE_ENSURE_INTEGER_SIZE_Nは使っていない？
 （追記）
 上のコードはScoreView v変数がunit型なのでわざわざこのようなことをしている
 のだと思うmg変数をv.half.mgに代入する時点で32bitの上位に代入され下位の16bitに
 干渉するので(uint16_t(eg) >> 15)を引いているものと思われる。
*/
enum Score {
  SCORE_ZERO,
  SCORE_ENSURE_INTEGER_SIZE_P = INT_MAX,
  SCORE_ENSURE_INTEGER_SIZE_N = INT_MIN
};

typedef union {
  uint32_t full;
  struct { int16_t eg, mg; } half;
} ScoreView;

inline Score make_score(int mg, int eg) {
  ScoreView v;
  v.half.mg = (int16_t)(mg - (uint16_t(eg) >> 15));
  v.half.eg = (int16_t)eg;
  return Score(v.full);
}

/*
駒評価値の内ミドルゲーム評価値を取り出す
上位ビットに収められている

（追記）
駒評価値ではないでのは
進行度のようなもの？
引数はScore型で戻り値はValue型
*/
inline Value mg_value(Score s) {
  ScoreView v;
  v.full = s;
  return Value(v.half.mg + (uint16_t(v.half.eg) >> 15));
}

/*
駒評価値の内エンドゲーム評価値を取り出す
下位ビットに収められている

（追記）
駒評価値ではないでのは
進行度のようなもの？
引数はScore型で戻り値はValue型
*/
inline Value eg_value(Score s) {
  ScoreView v;
  v.full = s;
  return Value(v.half.eg);
}
/*
stackfishではいろいろな型(Value,PieceType,Piece,Color
Depth,Square,Filem,Rank)が定義
してあるのでその型ごとの演算（四則演算、ビット演算など）をここで
再定義している
*/
#define ENABLE_BASE_OPERATORS_ON(T)                                         \
inline T operator+(const T d1, const T d2) { return T(int(d1) + int(d2)); } \
inline T operator-(const T d1, const T d2) { return T(int(d1) - int(d2)); } \
inline T operator*(int i, const T d) { return T(i * int(d)); }              \
inline T operator*(const T d, int i) { return T(int(d) * i); }              \
inline T operator-(const T d) { return T(-int(d)); }                        \
inline T& operator+=(T& d1, const T d2) { return d1 = d1 + d2; }            \
inline T& operator-=(T& d1, const T d2) { return d1 = d1 - d2; }            \
inline T& operator*=(T& d, int i) { return d = T(int(d) * i); }

ENABLE_BASE_OPERATORS_ON(Score)

#define ENABLE_FULL_OPERATORS_ON(T)                                         \
ENABLE_BASE_OPERATORS_ON(T)                                                 \
inline T& operator++(T& d) { return d = T(int(d) + 1); }                    \
inline T& operator--(T& d) { return d = T(int(d) - 1); }                    \
inline T operator/(const T d, int i) { return T(int(d) / i); }              \
inline T& operator/=(T& d, int i) { return d = T(int(d) / i); }

ENABLE_FULL_OPERATORS_ON(Value)
ENABLE_FULL_OPERATORS_ON(PieceType)
ENABLE_FULL_OPERATORS_ON(Piece)
ENABLE_FULL_OPERATORS_ON(Color)
ENABLE_FULL_OPERATORS_ON(Depth)
ENABLE_FULL_OPERATORS_ON(Square)
ENABLE_FULL_OPERATORS_ON(File)
ENABLE_FULL_OPERATORS_ON(Rank)

#undef ENABLE_FULL_OPERATORS_ON
#undef ENABLE_BASE_OPERATORS_ON

/// Additional operators to add integers to a Value
inline Value operator+(Value v, int i) { return Value(int(v) + i); }
inline Value operator-(Value v, int i) { return Value(int(v) - i); }
inline Value& operator+=(Value& v, int i) { return v = v + i; }
inline Value& operator-=(Value& v, int i) { return v = v - i; }

/// Only declared but not defined. We don't want to multiply two scores due to
/// a very high risk of overflow. So user should explicitly convert to integer.
inline Score operator*(Score s1, Score s2);

/// Division of a Score must be handled separately for each term
inline Score operator/(Score s, int i) {
  return make_score(mg_value(s) / i, eg_value(s) / i);
}

/*
駒の評価値、
PieceValue[0][] ミドルゲーム用駒評価値
PieceValue[1][] エンドゲーム用駒評価値
position.cppで初期化
*/
extern Value PieceValue[PHASE_NB][PIECE_NB];

/*
着手リストのデータ構造体
着手リストはmovepick.hに定義してあるMovePickerクラスの
中に入っている。
*/
struct ExtMove {
  Move move;
  Value value;
};

/*
ExtMove同士の比較演算子のオーバライド
*/
inline bool operator<(const ExtMove& f, const ExtMove& s) {
  return f.value < s.value;
}

/*
COLOR=WHITE=0,BLACK=1なので
0 ^ 1 = 1
1 ^ 1 = 0
となってカラーを切り替えることができる
*/
inline Color operator~(Color c) {
  return Color(c ^ BLACK);
}

/*
A1->A8
A2->A7
...
A8->A1

B1->B8
B2->B7
...
B8->B1
*/
inline Square operator~(Square s) {
  return Square(s ^ SQ_A8); // Vertical flip SQ_A1 -> SQ_A8
}

/*
キャステイングのための演算子
のようだが用途不明
*/
inline CastlingRight operator|(Color c, CastlingSide s) {
  return CastlingRight(WHITE_OO << ((s == QUEEN_SIDE) + 2 * c));
}

/*
seach.cppで使用されているが用途不明
*/
inline Value mate_in(int ply) {
  return VALUE_MATE - ply;
}

/*
seach.cppで使用されているが用途不明
*/
inline Value mated_in(int ply) {
  return -VALUE_MATE + ply;
}

/*
列と行から盤座標を計算する
*/
inline Square make_square(File f, Rank r) {
  return Square((r << 3) | f);
}

/*
駒種とカラーで駒コードを算出する
*/
inline Piece make_piece(Color c, PieceType pt) {
  return Piece((c << 3) | pt);
}

/*
駒コードから駒種を取り出す
*/
inline PieceType type_of(Piece pc)  {
  return PieceType(pc & 7);
}

/*
駒コードからカラーを取り出す
*/
inline Color color_of(Piece pc) {
  assert(pc != NO_PIECE);
  return Color(pc >> 3);
}

/*
座標値が盤内に入っているかチエック
*/
inline bool is_ok(Square s) {
  return s >= SQ_A1 && s <= SQ_H8;
}

/*
座標値の下位３bitを抜き出すことで列番号を取り出せる
*/
inline File file_of(Square s) {
  return File(s & 7);
}

/*
座標値から行番号を抜き出す
*/
inline Rank rank_of(Square s) {
  return Rank(s >> 3);
}

/*
行を入れ替える、カラーがWHITEのときはそのまま
列は入れ替えない
relative_square(WHITE,sq))の実行例
  0  1  2  3  4  5  6  7
  8  9 10 11 12 13 14 15
 16 17 18 19 20 21 22 23
 24 25 26 27 28 29 30 31
 32 33 34 35 36 37 38 39
 40 41 42 43 44 45 46 47
 48 49 50 51 52 53 54 55
 56 57 58 59 60 61 62 63
relative_square(BLACK,sq)の実行例
 56 57 58 59 60 61 62 63
 48 49 50 51 52 53 54 55
 40 41 42 43 44 45 46 47
 32 33 34 35 36 37 38 39
 24 25 26 27 28 29 30 31
 16 17 18 19 20 21 22 23
  8  9 10 11 12 13 14 15
  0  1  2  3  4  5  6  7
*/
inline Square relative_square(Color c, Square s) {
  return Square(s ^ (c * 56));
}

/*
行番号を入れ替える、カラーがWHITEのときはそのまま
relative_rank(WHITE,sq))の実行例
  0  0  0  0  0  0  0  0
  1  1  1  1  1  1  1  1
  2  2  2  2  2  2  2  2
  3  3  3  3  3  3  3  3
  4  4  4  4  4  4  4  4
  5  5  5  5  5  5  5  5
  6  6  6  6  6  6  6  6
  7  7  7  7  7  7  7  7
relative_rank(BLACK,sq)の実行例
  7  7  7  7  7  7  7  7
  6  6  6  6  6  6  6  6
  5  5  5  5  5  5  5  5
  4  4  4  4  4  4  4  4
  3  3  3  3  3  3  3  3
  2  2  2  2  2  2  2  2
  1  1  1  1  1  1  1  1
  0  0  0  0  0  0  0  0
*/
inline Rank relative_rank(Color c, Rank r) {
  return Rank(r ^ (c * 7));
}

/*
上と一緒だが引数に座標値を使う
行番号を入れ替える、カラーがWHITEのときはそのまま
*/
inline Rank relative_rank(Color c, Square s) {
  return relative_rank(c, rank_of(s));
}

/*
s2の位置のカラーを返す
（自駒の位置と同じカラーなら０
反対なら１を返す）
chess boardの市松模様をイメージすると
良い。
opposite_colors(sq,SQ_A1)
 0  1  0  1  0  1  0  1
 1  0  1  0  1  0  1  0
 0  1  0  1  0  1  0  1
 1  0  1  0  1  0  1  0
 0  1  0  1  0  1  0  1
 1  0  1  0  1  0  1  0
 0  1  0  1  0  1  0  1
 1  0  1  0  1  0  1  0
 opposite_colors(sq,SQ_A2)
  1  0  1  0  1  0  1  0
  0  1  0  1  0  1  0  1
  1  0  1  0  1  0  1  0
  0  1  0  1  0  1  0  1
  1  0  1  0  1  0  1  0
  0  1  0  1  0  1  0  1
  1  0  1  0  1  0  1  0
  0  1  0  1  0  1  0  1

*/
inline bool opposite_colors(Square s1, Square s2) {
  int s = int(s1) ^ int(s2);
  return ((s >> 3) ^ s) & 1;
}

/*
列番号の文字を返す
File_A -> a or A
File_B -> b or B
...
*/
inline char to_char(File f, bool tolower = true) {
  return char(f - FILE_A + (tolower ? 'a' : 'A'));
}

/*
行番号の文字を返す
RANK_1 -> 1
RANK_2 -> 2
...
*/
inline char to_char(Rank r) {
  return char(r - RANK_1 + '1');
}

/*
PAWNのカラーから移動方向をDELTA_NかDELTA_Sに決定している
*/
inline Square pawn_push(Color c) {
  return c == WHITE ? DELTA_N : DELTA_S;
}

/*
着手データ構造から移動元座標値を取り出す
*/
inline Square from_sq(Move m) {
  return Square((m >> 6) & 0x3F);
}

/*
着手データ構造から移動先座標値を取り出す
*/
inline Square to_sq(Move m) {
  return Square(m & 0x3F);
}

/*
着手の種別を取り出している
promotion (1), en passant (2), castling (3)
*/
inline MoveType type_of(Move m) {
  return MoveType(m & (3 << 14));
}

/*
成った駒の駒種
promotion piece type - 2 (from KNIGHT-2 to QUEEN-2)
*/
inline PieceType promotion_type(Move m) {
  return PieceType(((m >> 12) & 3) + 2);
}

/*
移動元(6bit)、移動先の座標(6bit)を
Move型にパックしている
*/
inline Move make_move(Square from, Square to) {
  return Move(to | (from << 6));
}

/*
テンプレート関数
Moveデータ構造にするための関数
*/
template<MoveType T>
inline Move make(Square from, Square to, PieceType pt = KNIGHT) {
  return Move(to | (from << 6) | T | ((pt - KNIGHT) << 12));
}

/*
移動のチエック関数
*/
inline bool is_ok(Move m) {
  return from_sq(m) != to_sq(m); // Catches also MOVE_NULL and MOVE_NONE
}

#include <string>

/*
座標値を座標文字列にして返す
*/
inline const std::string to_string(Square s) {
  char ch[] = { to_char(file_of(s)), to_char(rank_of(s)), 0 };
  return ch;
}

#endif // #ifndef TYPES_H_INCLUDED
