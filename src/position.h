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

#ifndef POSITION_H_INCLUDED
#define POSITION_H_INCLUDED

#include <cassert>
#include <cstddef>

#include "bitboard.h"
#include "types.h"


/// The checkInfo struct is initialized at c'tor time and keeps info used
/// to detect if a move gives check.
class Position;
struct Thread;

/*
用途不明
*/
struct CheckInfo {

  explicit CheckInfo(const Position&);

  Bitboard dcCandidates;
  Bitboard pinned;
  Bitboard checkSq[PIECE_TYPE_NB];
  Square ksq;
};


/// The StateInfo struct stores information needed to restore a Position
/// object to its previous state when we retract a move. Whenever a move
/// is made on the board (by calling Position::do_move), a StateInfo
/// object must be passed as a parameter.

/*
用途不明
（追記）
positionは局面を保持するクラスで、do_moveによって
変更を加えられるが（局面更新）もとの局面に戻す時に
使用される情報を入れておく構造体かな
*/
struct StateInfo {
  Key pawnKey, materialKey;
  Value npMaterial[COLOR_NB];
  int castlingRights, rule50, pliesFromNull;
  Score psq;
  Square epSquare;

  Key key;
  Bitboard checkersBB;
  PieceType capturedType;
  StateInfo* previous;
};


/// When making a move the current StateInfo up to 'key' excluded is copied to
/// the new one. Here we calculate the quad words (64bits) needed to be copied.
/*
offsetofの機能
    構造体のメンバのバイト位置を示す整数を返す
用途不明
*/
const size_t StateCopySize64 = offsetof(StateInfo, key) / sizeof(uint64_t) + 1;


/// The Position class stores the information regarding the board representation
/// like pieces, side to move, hash keys, castling info, etc. The most important
/// methods are do_move() and undo_move(), used by the search to update node info
/// when traversing the search tree.
/*
局面を表現するクラス
*/
class Position {
public:
  Position() {}
  Position(const Position& pos, Thread* t) { *this = pos; thisThread = t; }
  Position(const std::string& f, bool c960, Thread* t) { set(f, c960, t); }
  Position& operator=(const Position&);
  static void init();

  // Text input/output
  /*
  fenStrは局面を文字列で表現したもの
	＜例＞
	"rnbqkbnr/pppppppp/8/8/8/8/PPPPPPPP/RNBQKBNR w KQkq - 0 1",
	R--rook,N--night,B--bishop,Q--queen,K--king,P--pawn,大文字はwhite（先手）
	r--rook,n--night,b--bishop,q--queen,k--king,p--pawn,小文字はblack（後手）
	数字は空白の数,/は行の終わり 局面を表現する文字列のあと空白を入れてこの局面で次に
	指すカラーをw/bで表現、その次のKQkqは不明 -も不明　0 1も不明
	（追記）
	KQkqはキャスリングに関係するなにか
  */
	/*
	set関数はPositionクラスのコンストラクタから呼ばれ
	fenStrを解析して内部データを更新している。
	更新される内部データは
		board[s]
		byTypeBB[ALL_PIECES]
	  byColorBB[c]
		index[s]
		pieceList[c][pt][index[s]]
	他にもいろいろ局面保持、局面更新に必要なものを初期化しているが
	詳細不明
	*/
  void set(const std::string& fenStr, bool isChess960, Thread* th);
	/*
	現在の局面のfenStr文字列に変換する
	*/
  const std::string fen() const;
	/*
	指した手と現在の局面を文字列にして返す
	*/
  const std::string pretty(Move m = MOVE_NONE) const;

  // Position representation
	/*
	用途不明
	（追記）
	全ての駒のビットがたったbitboardを返す
	関数の宣言のあとについているconstはこの
	positionクラスのメンバー変数は変更できない
	ことを示している。（bitboardを読みだして返す
	関数なので変更することがない）
	*/
  Bitboard pieces() const;
	/*
	指定した駒種が立ったbitboardを返す
	*/
  Bitboard pieces(PieceType pt) const;
	/*
	pt1,pt2で指定した駒種が立ったbitboardを返す
	*/
  Bitboard pieces(PieceType pt1, PieceType pt2) const;
	/*
	指定したカラーの駒が立っているbitboardを返す
	*/
  Bitboard pieces(Color c) const;
	/*
	指定した駒種、指定したカラーの駒が立っているbitboardを返す
	*/
  Bitboard pieces(Color c, PieceType pt) const;
	/*
	指定した駒種pt1,pt2、指定したカラーの駒が立っているbitboardを返す
	*/
  Bitboard pieces(Color c, PieceType pt1, PieceType pt2) const;
	/*
	メンバー変数board[64]の座標を指定して駒コードを返す
	*/
  Piece piece_on(Square s) const;
	/*
	カラーを指定してpieceList[COLOR][PIECE_TYPE][16]変数からkingの座標を返す
	*/
  Square king_square(Color c) const;
	/*
	用途不明
	*/
  Square ep_square() const;
	/*
	指定した座標に駒がなかったらtrueを返す
	*/
  bool empty(Square s) const;
	/*
	カラーと駒種を指定すればpieceCount[COLOR][PIECE_TYPE]を
	返す。pieceCount配列がなんなのか分かっていないが
	多分駒種、カラーごとの駒数だと思われる
	ところでtemplate<PieceType Pt>とはなんだろう
	おそらく駒種ごとにcount関数をテンプレート化しているのでは
	それだけ高速に動作する？
	*/
  template<PieceType Pt> int count(Color c) const;
	/*
	カラーと駒種を指定すればpieceList[COLOR][PIECE_TYPE]の
	座標リスト（配列）を返してくる
	この関数も駒種ごとにテンプレート化されている？
	*/
  template<PieceType Pt> const Square* list(Color c) const;

  // Castling
	/*
	用途不明
	多分、キャステイングのことだと思う
	*/
  int can_castle(Color c) const;
	/*
	用途不明
	多分、キャステイングのことだと思う
	*/
  int can_castle(CastlingRight cr) const;
	/*
	用途不明
	多分、キャステイングのことだと思う
	*/
  bool castling_impeded(CastlingRight cr) const;
	/*
	用途不明
	多分、キャステイングのことだと思う
	*/
  Square castling_rook_square(CastlingRight cr) const;

  // Checking
	/*
	用途不明
	StateInfoのメンバーcheckersBBを返してくる
	*/
  Bitboard checkers() const;
	/*
	用途不明
	多分、kingへのチエックに関するなにか
	*/
  Bitboard discovered_check_candidates() const;
	/*
	用途不明
	多分、kingへのチエックに関するなにか
	*/
  Bitboard pinned_pieces(Color c) const;

  // Attacks to/from a given square
	/*
	用途不明
	利き計算だとおもう
	*/
  Bitboard attackers_to(Square s) const;
  Bitboard attackers_to(Square s, Bitboard occ) const;
  Bitboard attacks_from(Piece pc, Square s) const;
  template<PieceType> Bitboard attacks_from(Square s) const;
  template<PieceType> Bitboard attacks_from(Square s, Color c) const;

  // Properties of moves
	/*
	渡された手が合法手かチエックする
	*/
  bool legal(Move m, Bitboard pinned) const;
	/*
	用途不明
	*/
  bool pseudo_legal(const Move m) const;
	/*
	着手が駒を取る手ならtrueを返す
	*/
  bool capture(Move m) const;
	/*
	用途不明
	着手の種別をチエックしているが
	関数名から判断するに駒をとる手
	もしくは成る手を判定しているようだが
	？
	*/
  bool capture_or_promotion(Move m) const;
	/*
	用途不明
	合法手の判定？
	*/
  bool gives_check(Move m, const CheckInfo& ci) const;
	/*
	用途不明
	PAWNに関することのようであるが不明
	*/
  bool advanced_pawn_push(Move m) const;
	/*
	着手データから駒コードを取得する
	*/
  Piece moved_piece(Move m) const;
	/*
	用途不明
	取る駒種？
	*/
  PieceType captured_piece_type() const;

  // Piece specific
	/*
	用途不明
	*/
  bool pawn_passed(Color c, Square s) const;
	/*
	PAWNが相手陣地の最下段にいるかチエック、いたら
	trueを返す、PAWNが成るためのチエックに使用しているのかも？
	*/
  bool pawn_on_7th(Color c) const;
	/*
	用途不明
	同じカラーのbishopが２個以上あって（つまりとられていない）
	boardのカラー（市松模様の色）が異なっていればtrue
	ただ味方同士のbishopは互いに異なるboardカラーに配置されるので
	trueが返ってくるのは当たり前のような気がするが、chess960では
	違うのかも？
	*/
  bool bishop_pair(Color c) const;
	/*
	WHITE、BLACK側のbishopが互い違いのboardカラーに位置すれば
	trueを返す
	つまり違うboardカラーの場合お互いに取り合うことはないと
	判断できる
	*/
  bool opposite_bishops() const;

  // Doing and undoing moves
  /*
	局面更新、詳細不明
	ここからいろいろ引数を追加して下の
	do_moveを呼んでいるラッパー
	*/
	void do_move(Move m, StateInfo& st);
	/*
	局面更新、詳細不明
	*/
  void do_move(Move m, StateInfo& st, const CheckInfo& ci, bool moveIsCheck);
	/*
	局面復元、詳細不明
	*/
  void undo_move(Move m);
	/*
	null moveを動かす,詳細不明
	*/
  void do_null_move(StateInfo& st);
	/*
	do_null_moveの復元
	*/
  void undo_null_move();

  // Static exchange evaluation
	/*
	用途不明
	SEEの名前からして静止探索？
	*/
  Value see(Move m) const;
	/*
	用途不明
	*/
  Value see_sign(Move m) const;

  // Accessing hash keys
	/*
	用途不明
	StateInfoのメンバーkeyを返すだけ
	*/
  Key key() const;
	/*
	用途不明
	*/
  Key exclusion_key() const;
	/*
	用途不明
	StateInfoのメンバーpawnKeyを返すだけ
	*/
  Key pawn_key() const;
	/*
	用途不明
	StateInfoのメンバーmaterialKeyを返すだけ
	*/
  Key material_key() const;

  // Incremental piece-square evaluation
	/*
	用途不明
	StateInfoのメンバーpsqを返す
	*/
  Score psq_score() const;
	/*
	用途不明
	StateInfoのメンバーnpMaterial
	*/
  Value non_pawn_material(Color c) const;

  // Other properties of the position
	/*
	メンバー変数sideToMoveを返す
	*/
  Color side_to_move() const;
	/*
	ゲームの何手目かを返す
	*/
  int game_ply() const;
	/*
	チェス960（Chess 960）は、変則チェスの一種
	chess960かどうかを返す
	*/
  bool is_chess960() const;
	/*
	用途不明
	メンバー変数thisThreadを返す
	*/
  Thread* this_thread() const;
	/*
	用途不明
	メンバー変数nodesを返す
	*/
  uint64_t nodes_searched() const;
	/*
	用途不明
	*/
  void set_nodes_searched(uint64_t n);
	/*
	用途不明
	引き分けの判定をしている？
	*/
  bool is_draw() const;

  // Position consistency check, for debugging
	/*
	局面の不整合などをチエックしている
	詳細不明なところもある
	*/
  bool pos_is_ok(int* step = NULL) const;
	/*
	ここまで、position.cppにコメント入れ中
	*/
  void flip();

private:
  // Initialization helpers (used while setting up a position)
  void clear();
  void set_castling_right(Color c, Square rfrom);
  void set_state(StateInfo* si) const;

  // Helper functions
  Bitboard check_blockers(Color c, Color kingColor) const;
  void put_piece(Square s, Color c, PieceType pt);
  void remove_piece(Square s, Color c, PieceType pt);
  void move_piece(Square from, Square to, Color c, PieceType pt);
  template<bool Do>
  void do_castling(Square from, Square& to, Square& rfrom, Square& rto);

  // Board and pieces
  Piece board[SQUARE_NB];
  Bitboard byTypeBB[PIECE_TYPE_NB];
  Bitboard byColorBB[COLOR_NB];
  int pieceCount[COLOR_NB][PIECE_TYPE_NB];
  Square pieceList[COLOR_NB][PIECE_TYPE_NB][16];
  int index[SQUARE_NB];

  // Other info
  int castlingRightsMask[SQUARE_NB];
  Square castlingRookSquare[CASTLING_RIGHT_NB];
  Bitboard castlingPath[CASTLING_RIGHT_NB];
  StateInfo startState;
  uint64_t nodes;
  int gamePly;
  Color sideToMove;
  Thread* thisThread;
  StateInfo* st;
  bool chess960;
};

inline uint64_t Position::nodes_searched() const {
  return nodes;
}

inline void Position::set_nodes_searched(uint64_t n) {
  nodes = n;
}

inline Piece Position::piece_on(Square s) const {
  return board[s];
}

/*
着手データから駒データを取得する
*/
inline Piece Position::moved_piece(Move m) const {
  return board[from_sq(m)];
}

inline bool Position::empty(Square s) const {
  return board[s] == NO_PIECE;
}

inline Color Position::side_to_move() const {
  return sideToMove;
}

inline Bitboard Position::pieces() const {
  return byTypeBB[ALL_PIECES];
}

inline Bitboard Position::pieces(PieceType pt) const {
  return byTypeBB[pt];
}

inline Bitboard Position::pieces(PieceType pt1, PieceType pt2) const {
  return byTypeBB[pt1] | byTypeBB[pt2];
}

inline Bitboard Position::pieces(Color c) const {
  return byColorBB[c];
}

inline Bitboard Position::pieces(Color c, PieceType pt) const {
  return byColorBB[c] & byTypeBB[pt];
}

inline Bitboard Position::pieces(Color c, PieceType pt1, PieceType pt2) const {
  return byColorBB[c] & (byTypeBB[pt1] | byTypeBB[pt2]);
}

template<PieceType Pt> inline int Position::count(Color c) const {
  return pieceCount[c][Pt];
}

template<PieceType Pt> inline const Square* Position::list(Color c) const {
  return pieceList[c][Pt];
}

inline Square Position::ep_square() const {
  return st->epSquare;
}

inline Square Position::king_square(Color c) const {
  return pieceList[c][KING][0];
}

inline int Position::can_castle(CastlingRight cr) const {
  return st->castlingRights & cr;
}

inline int Position::can_castle(Color c) const {
  return st->castlingRights & ((WHITE_OO | WHITE_OOO) << (2 * c));
}

inline bool Position::castling_impeded(CastlingRight cr) const {
  return byTypeBB[ALL_PIECES] & castlingPath[cr];
}

inline Square Position::castling_rook_square(CastlingRight cr) const {
  return castlingRookSquare[cr];
}

template<PieceType Pt>
inline Bitboard Position::attacks_from(Square s) const {

  return  Pt == BISHOP || Pt == ROOK ? attacks_bb<Pt>(s, byTypeBB[ALL_PIECES])
        : Pt == QUEEN  ? attacks_from<ROOK>(s) | attacks_from<BISHOP>(s)
        : StepAttacksBB[Pt][s];
}

template<>
inline Bitboard Position::attacks_from<PAWN>(Square s, Color c) const {
  return StepAttacksBB[make_piece(c, PAWN)][s];
}

inline Bitboard Position::attacks_from(Piece pc, Square s) const {
  return attacks_bb(pc, s, byTypeBB[ALL_PIECES]);
}

inline Bitboard Position::attackers_to(Square s) const {
  return attackers_to(s, byTypeBB[ALL_PIECES]);
}

inline Bitboard Position::checkers() const {
  return st->checkersBB;
}

inline Bitboard Position::discovered_check_candidates() const {
  return check_blockers(sideToMove, ~sideToMove);
}

inline Bitboard Position::pinned_pieces(Color c) const {
  return check_blockers(c, c);
}

inline bool Position::pawn_passed(Color c, Square s) const {
  return !(pieces(~c, PAWN) & passed_pawn_mask(c, s));
}

inline bool Position::advanced_pawn_push(Move m) const {
  return   type_of(moved_piece(m)) == PAWN
        && relative_rank(sideToMove, from_sq(m)) > RANK_4;
}

inline Key Position::key() const {
  return st->key;
}

inline Key Position::pawn_key() const {
  return st->pawnKey;
}

inline Key Position::material_key() const {
  return st->materialKey;
}

inline Score Position::psq_score() const {
  return st->psq;
}

inline Value Position::non_pawn_material(Color c) const {
  return st->npMaterial[c];
}

/*
ゲームの何手目かを返す
*/
inline int Position::game_ply() const {
  return gamePly;
}

/*
WHITE、BLACK側のbishopが互い違いのboardカラーに位置すれば
trueを返す
つまり違うboardカラーの場合お互いに取り合うことはないと
判断できる
*/
inline bool Position::opposite_bishops() const {

  return   pieceCount[WHITE][BISHOP] == 1
        && pieceCount[BLACK][BISHOP] == 1
        && opposite_colors(pieceList[WHITE][BISHOP][0], pieceList[BLACK][BISHOP][0]);
}

/*
用途不明
同じカラーのbishopが２個以上あって（つまりとられていない）
boardのカラー（市松模様の色）が異なっていればtrue
ただ味方同士のbishopは互いに異なるboardカラーに配置されるので
trueが返ってくるのは当たり前のような気がするが、chess960では
違うのかも？
*/
inline bool Position::bishop_pair(Color c) const {

  return   pieceCount[c][BISHOP] >= 2
        && opposite_colors(pieceList[c][BISHOP][0], pieceList[c][BISHOP][1]);
}

inline bool Position::pawn_on_7th(Color c) const {
  return pieces(c, PAWN) & rank_bb(relative_rank(c, RANK_7));
}

inline bool Position::is_chess960() const {
  return chess960;
}

inline bool Position::capture_or_promotion(Move m) const {

  assert(is_ok(m));
  return type_of(m) != NORMAL ? type_of(m) != CASTLING : !empty(to_sq(m));
}

inline bool Position::capture(Move m) const {

  // Note that castling is encoded as "king captures the rook"
  assert(is_ok(m));
  return (!empty(to_sq(m)) && type_of(m) != CASTLING) || type_of(m) == ENPASSANT;
}

inline PieceType Position::captured_piece_type() const {
  return st->capturedType;
}

inline Thread* Position::this_thread() const {
  return thisThread;
}

inline void Position::put_piece(Square s, Color c, PieceType pt) {

  board[s] = make_piece(c, pt);
  byTypeBB[ALL_PIECES] |= s;
  byTypeBB[pt] |= s;
  byColorBB[c] |= s;
  index[s] = pieceCount[c][pt]++;
  pieceList[c][pt][index[s]] = s;
}

inline void Position::move_piece(Square from, Square to, Color c, PieceType pt) {

  // index[from] is not updated and becomes stale. This works as long
  // as index[] is accessed just by known occupied squares.
  Bitboard from_to_bb = SquareBB[from] ^ SquareBB[to];
  byTypeBB[ALL_PIECES] ^= from_to_bb;
  byTypeBB[pt] ^= from_to_bb;
  byColorBB[c] ^= from_to_bb;
  board[from] = NO_PIECE;
  board[to] = make_piece(c, pt);
  index[to] = index[from];
  pieceList[c][pt][index[to]] = to;
}

inline void Position::remove_piece(Square s, Color c, PieceType pt) {

  // WARNING: This is not a reversible operation. If we remove a piece in
  // do_move() and then replace it in undo_move() we will put it at the end of
  // the list and not in its original place, it means index[] and pieceList[]
  // are not guaranteed to be invariant to a do_move() + undo_move() sequence.
  byTypeBB[ALL_PIECES] ^= s;
  byTypeBB[pt] ^= s;
  byColorBB[c] ^= s;
  /* board[s] = NO_PIECE; */ // Not needed, will be overwritten by capturing
  Square lastSquare = pieceList[c][pt][--pieceCount[c][pt]];
  index[lastSquare] = index[s];
  pieceList[c][pt][index[lastSquare]] = lastSquare;
  pieceList[c][pt][pieceCount[c][pt]] = SQ_NONE;
}

#endif // #ifndef POSITION_H_INCLUDED
