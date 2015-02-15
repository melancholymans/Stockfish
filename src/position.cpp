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
#include <cassert>
#include <cstring>
#include <iomanip>
#include <sstream>

#include "bitcount.h"
#include "movegen.h"
#include "notation.h"
#include "position.h"
#include "psqtab.h"
#include "rkiss.h"
#include "thread.h"
#include "tt.h"

using std::string;

static const string PieceToChar(" PNBRQK  pnbrqk");

CACHE_LINE_ALIGNMENT

Value PieceValue[PHASE_NB][PIECE_NB] = {
{ VALUE_ZERO, PawnValueMg, KnightValueMg, BishopValueMg, RookValueMg, QueenValueMg },
{ VALUE_ZERO, PawnValueEg, KnightValueEg, BishopValueEg, RookValueEg, QueenValueEg } };

static Score psq[COLOR_NB][PIECE_TYPE_NB][SQUARE_NB];

/*
chessなどの局面の状態を１つのハッシュ値で代表させる方法
参考HP:http://hackemdown.blogspot.jp/2014/06/zobrist-hashing.html
初期化しているのはPosition::init()関数内
*/
namespace Zobrist {

  Key psq[COLOR_NB][PIECE_TYPE_NB][SQUARE_NB];
  Key enpassant[FILE_NB];
  Key castling[CASTLING_RIGHT_NB];
  Key side;
  Key exclusion;
}

/*
用途不明
*/
Key Position::exclusion_key() const
{ 
	return st->key ^ Zobrist::exclusion;
}

namespace {

// min_attacker() is a helper function used by see() to locate the least
// valuable attacker for the side to move, remove the attacker we just found
// from the bitboards and scan for new X-ray attacks behind it.
/*
min_attacker関数はsee関数のヘルパー関数
でsee関数はおそらく静止探索だと思う
詳細不明
*/
template<int Pt> FORCE_INLINE
PieceType min_attacker(const Bitboard* bb, const Square& to, const Bitboard& stmAttackers,
                       Bitboard& occupied, Bitboard& attackers) {

  Bitboard b = stmAttackers & bb[Pt];
  if (!b)
      return min_attacker<Pt+1>(bb, to, stmAttackers, occupied, attackers);

  occupied ^= b & ~(b - 1);
  /*
  ここの処理不明
  */
  if (Pt == PAWN || Pt == BISHOP || Pt == QUEEN)
      attackers |= attacks_bb<BISHOP>(to, occupied) & (bb[BISHOP] | bb[QUEEN]);

  if (Pt == ROOK || Pt == QUEEN)
      attackers |= attacks_bb<ROOK>(to, occupied) & (bb[ROOK] | bb[QUEEN]);

  attackers &= occupied; // After X-ray that may add already processed pieces
  return (PieceType)Pt;
}

/*
テンプレート関数（KINGの明示化）
*/
template<> FORCE_INLINE
PieceType min_attacker<KING>(const Bitboard*, const Square&, const Bitboard&, Bitboard&, Bitboard&) {
  return KING; // No need to update bitboards: it is the last cycle
}

} // namespace


/// CheckInfo c'tor
/*
このCheckInfoクラスはコンストラクタしかない
局面クラスpositionを受け取って現局面で王手をかけている駒種ごとのbitboard（もちろんチエックがかかっていない場合は0）
敵KINGに対してpin付けされている駒のbitboardを返す
*/
CheckInfo::CheckInfo(const Position& pos) {

  Color them = ~pos.side_to_move();
  ksq = pos.king_square(them);

  pinned = pos.pinned_pieces(pos.side_to_move());	
  dcCandidates = pos.discovered_check_candidates();	//dcCandidatesは相手のKINGへの利きの邪魔になっている自陣側の駒bitboardを返す

  checkSq[PAWN]   = pos.attacks_from<PAWN>(ksq, them);
  checkSq[KNIGHT] = pos.attacks_from<KNIGHT>(ksq);
  checkSq[BISHOP] = pos.attacks_from<BISHOP>(ksq);
  checkSq[ROOK]   = pos.attacks_from<ROOK>(ksq);
  checkSq[QUEEN]  = checkSq[BISHOP] | checkSq[ROOK];
  checkSq[KING]   = 0;
}


/// Position::init() initializes at startup the various arrays used to compute
/// hash keys and the piece square tables. The latter is a two-step operation:
/// Firstly, the white halves of the tables are copied from PSQT[] tables.
/// Secondly, the black halves of the tables are initialized by flipping and
/// changing the sign of the white scores.
/*
Zobristを初期化している
そのあと評価値（駒評価値と位置評価値）を初期している
*/
void Position::init() {

  RKISS rk;

  /*
  升目、駒種ごとに乱数をあらかじめ設定しておき
  局面の状態に応じて１意のハッシュ値（異局面で同一のハッシュ値がでる可能性はある）
  */
  for (Color c = WHITE; c <= BLACK; ++c)
      for (PieceType pt = PAWN; pt <= KING; ++pt)
          for (Square s = SQ_A1; s <= SQ_H8; ++s)
              Zobrist::psq[c][pt][s] = rk.rand<Key>();
  /*
  アンパッサンのハッシュ値？
  用途不明
  */
  for (File f = FILE_A; f <= FILE_H; ++f)
      Zobrist::enpassant[f] = rk.rand<Key>();
  /*
  キャスリングのハッシュ値？
  用途不明
  */
  for (int cf = NO_CASTLING; cf <= ANY_CASTLING; ++cf)
  {
      Bitboard b = cf;
      while (b)
      {
          Key k = Zobrist::castling[1ULL << pop_lsb(&b)];
          Zobrist::castling[cf] ^= k ? k : rk.rand<Key>();
      }
  }
  /*
  用途不明
  */
  Zobrist::side = rk.rand<Key>();
  Zobrist::exclusion  = rk.rand<Key>();

  /*
  WHITE側の駒評価値は直接設定している（このposition.cppの冒頭部分）
  ここではWHITE側の評価値をBLACK側にコピーしている
  */
  for (PieceType pt = PAWN; pt <= KING; ++pt)
  {
      PieceValue[MG][make_piece(BLACK, pt)] = PieceValue[MG][pt];
      PieceValue[EG][make_piece(BLACK, pt)] = PieceValue[EG][pt];

      Score v = make_score(PieceValue[MG][pt], PieceValue[EG][pt]);
	  /*
	  PSQTはpsqtab.hに定義してある配列でScore変数が（32bitの上位16bitにミドルゲーム駒評価値を、下位16bitにエンドゲーム駒評価値を設定してある）
	  盤位置に応じて格納されている。位置評価値の基本位置評価値と言える
	  psq[BLACK][pt][~s]の~sは演算子のオーバーロードで座標変換している (例A1->A8,B2->B7）
	  white側（先手）がプラス、BLACK側が（後手）マイナスをもつ
	  基本位置評価値に駒評価値を加算してpsq配列を初期化している
	  ｐｓｑ配列はstatic Score psq[COLOR_NB][PIECE_TYPE_NB][SQUARE_NB];と宣言されている
	  おそらく初期化のされかた、ネーミングから駒評価値と位置評価値を組み合わせたもの
	  */
      for (Square s = SQ_A1; s <= SQ_H8; ++s)
      {
         psq[WHITE][pt][ s] =  (v + PSQT[pt][s]);
         psq[BLACK][pt][~s] = -(v + PSQT[pt][s]);
      }
  }
}


/// Position::operator=() creates a copy of 'pos'. We want the new born Position
/// object to not depend on any external data so we detach state pointer from
/// the source one.
/*
局面クラスpositionをコピーする演算子のオーバーロード
startStateは用途不明
*/
Position& Position::operator=(const Position& pos) {

  std::memcpy(this, &pos, sizeof(Position));
  startState = *st;
  st = &startState;
  nodes = 0;

  assert(pos_is_ok());

  return *this;
}


/// Position::clear() erases the position object to a pristine state, with an
/// empty board, white to move, and no castling rights.
/*
positionクラスをクリアにする、
startStateは用途不明
*/
void Position::clear() {

  std::memset(this, 0, sizeof(Position));
  startState.epSquare = SQ_NONE;
  st = &startState;

  for (int i = 0; i < PIECE_TYPE_NB; ++i)
      for (int j = 0; j < 16; ++j)
          pieceList[WHITE][i][j] = pieceList[BLACK][i][j] = SQ_NONE;
}


/// Position::set() initializes the position object with the given FEN string.
/// This function is not very robust - make sure that input FENs are correct,
/// this is assumed to be the responsibility of the GUI.
/*
FEN stringを読み込んで局面を設定している
*/
void Position::set(const string& fenStr, bool isChess960, Thread* th) {
/*
   A FEN string defines a particular position using only the ASCII character set.

   A FEN string contains six fields separated by a space. The fields are:

   1) Piece placement (from white's perspective). Each rank is described, starting
      with rank 8 and ending with rank 1. Within each rank, the contents of each
      square are described from file A through file H. Following the Standard
      Algebraic Notation (SAN), each piece is identified by a single letter taken
      from the standard English names. White pieces are designated using upper-case
      letters ("PNBRQK") whilst Black uses lowercase ("pnbrqk"). Blank squares are
      noted using digits 1 through 8 (the number of blank squares), and "/"
      separates ranks.

   2) Active color. "w" means white moves next, "b" means black.

   3) Castling availability. If neither side can castle, this is "-". Otherwise,
      this has one or more letters: "K" (White can castle kingside), "Q" (White
      can castle queenside), "k" (Black can castle kingside), and/or "q" (Black
      can castle queenside).

   4) En passant target square (in algebraic notation). If there's no en passant
      target square, this is "-". If a pawn has just made a 2-square move, this
      is the position "behind" the pawn. This is recorded regardless of whether
      there is a pawn in position to make an en passant capture.

   5) Halfmove clock. This is the number of halfmoves since the last pawn advance
      or capture. This is used to determine if a draw can be claimed under the
      fifty-move rule.

   6) Fullmove number. The number of the full move. It starts at 1, and is
      incremented after Black's move.
*/

  char col, row, token;
  size_t idx;
  Square sq = SQ_A8;
  std::istringstream ss(fenStr);

  clear();
  //空白文字をスキップさせない設定
  ss >> std::noskipws;

  // 1. Piece placement
  //FEN stringのスキャンはA8->B8->...->H8
  //A7->B7..->H7
  //A1->B1->..H1と読み取ってい行く
  while ((ss >> token) && !isspace(token))
  {
	  //数字は空白を表すので数値だけ座標を加算する
      if (isdigit(token))
          sq += Square(token - '0'); // Advance the given number of files
	  //FEN stringの見本
	  //"rnbqkbnr/pppppppp/8/8/8/8/PPPPPPPP/RNBQKBNR w KQkq - 0 1",
	  //
      else if (token == '/')
          sq -= Square(16);
	  //ここはFEN stringを駒コードに変換してput_pieceを呼んで内部データを更新している
      else if ((idx = PieceToChar.find(token)) != string::npos)
      {
          put_piece(sq, color_of(Piece(idx)), type_of(Piece(idx)));
          ++sq;
      }
  }

  // 2. Active color
  //手番を設定している
  ss >> token;
  sideToMove = (token == 'w' ? WHITE : BLACK);
  ss >> token;

  // 3. Castling availability. Compatible with 3 standards: Normal FEN standard,
  // Shredder-FEN that uses the letters of the columns on which the rooks began
  // the game instead of KQkq and also X-FEN standard that, in case of Chess960,
  // if an inner rook is associated with the castling right, the castling tag is
  // replaced by the file letter of the involved rook, as for the Shredder-FEN.
  //キャスリングに関することのようだが詳細不明
  while ((ss >> token) && !isspace(token))
  {
      Square rsq;
      Color c = islower(token) ? BLACK : WHITE;

      token = char(toupper(token));

      if (token == 'K')
          for (rsq = relative_square(c, SQ_H1); type_of(piece_on(rsq)) != ROOK; --rsq) {}

      else if (token == 'Q')
          for (rsq = relative_square(c, SQ_A1); type_of(piece_on(rsq)) != ROOK; ++rsq) {}

      else if (token >= 'A' && token <= 'H')
          rsq = make_square(File(token - 'A'), relative_rank(c, RANK_1));

      else
          continue;

      set_castling_right(c, rsq);
  }

  // 4. En passant square. Ignore if no pawn capture is possible
  //用途不明
  if (   ((ss >> col) && (col >= 'a' && col <= 'h'))
      && ((ss >> row) && (row == '3' || row == '6')))
  {
      st->epSquare = make_square(File(col - 'a'), Rank(row - '1'));

      if (!(attackers_to(st->epSquare) & pieces(sideToMove, PAWN)))
          st->epSquare = SQ_NONE;
  }

  // 5-6. Halfmove clock and fullmove number
  //用途不明
  ss >> std::skipws >> st->rule50 >> gamePly;

  // Convert from fullmove starting from 1 to ply starting from 0,
  // handle also common incorrect FEN with fullmove = 0.
  //用途不明
  gamePly = std::max(2 * (gamePly - 1), 0) + (sideToMove == BLACK);
  //用途不明
  chess960 = isChess960;
  thisThread = th;
  set_state(st);

  assert(pos_is_ok());
}


/// Position::set_castling_right() is a helper function used to set castling
/// rights given the corresponding color and the rook starting square.
/*
多分キャスリングに関するなにか
用途不明
*/
void Position::set_castling_right(Color c, Square rfrom) {

  Square kfrom = king_square(c);
  CastlingSide cs = kfrom < rfrom ? KING_SIDE : QUEEN_SIDE;
  CastlingRight cr = (c | cs);

  st->castlingRights |= cr;
  castlingRightsMask[kfrom] |= cr;
  castlingRightsMask[rfrom] |= cr;
  castlingRookSquare[cr] = rfrom;

  Square kto = relative_square(c, cs == KING_SIDE ? SQ_G1 : SQ_C1);
  Square rto = relative_square(c, cs == KING_SIDE ? SQ_F1 : SQ_D1);

  for (Square s = std::min(rfrom, rto); s <= std::max(rfrom, rto); ++s)
      if (s != kfrom && s != rfrom)
          castlingPath[cr] |= s;

  for (Square s = std::min(kfrom, kto); s <= std::max(kfrom, kto); ++s)
      if (s != kfrom && s != rfrom)
          castlingPath[cr] |= s;
}


/// Position::set_state() computes the hash keys of the position, and other
/// data that once computed is updated incrementally as moves are made.
/// The function is only used when a new position is set up, and to verify
/// the correctness of the StateInfo data when running in debug mode.

/*
Position::set関数から呼ばれており、Positionクラスの主要変数の初期化を
していると思われる。
Zobristはchessなどに使用される局面の状態を１つのハッシュ値で代表させる。
参考HP:http://hackemdown.blogspot.jp/2014/06/zobrist-hashing.html

*/
void Position::set_state(StateInfo* si) const {

  si->key = si->pawnKey = si->materialKey = 0;
  si->npMaterial[WHITE] = si->npMaterial[BLACK] = VALUE_ZERO;
  si->psq = SCORE_ZERO;

  si->checkersBB = attackers_to(king_square(sideToMove)) & pieces(~sideToMove);
  /*
  すべての駒を巡回する
  */
  for (Bitboard b = pieces(); b; )
  {
      Square s = pop_lsb(&b);
      Piece pc = piece_on(s);
	  /*
	  現局面でのハッシュ値
	  */
      si->key ^= Zobrist::psq[color_of(pc)][type_of(pc)][s];
	  /*
	  psqは駒種ごとの位置評価値なので、si->psqは現局面での位置評価値の集計
	  */
      si->psq += psq[color_of(pc)][type_of(pc)][s];
  }
  /*
  詳細不明であるがアンパッサンになにか関係ある？
  */
  if (ep_square() != SQ_NONE)
      si->key ^= Zobrist::enpassant[file_of(ep_square())];
  /*
  カラーごとのハッシュ値をつける
  */
  if (sideToMove == BLACK)
      si->key ^= Zobrist::side;
  /*
  詳細不明であるがキャスリングごとのハッシュ値をつける
  */
  si->key ^= Zobrist::castling[st->castlingRights];
  /*
  PAWNだけ巡回して、PAWNだけのハッシュ値を求める
  */
  for (Bitboard b = pieces(PAWN); b; )
  {
      Square s = pop_lsb(&b);
      si->pawnKey ^= Zobrist::psq[color_of(piece_on(s))][PAWN][s];
  }
  /*
  カラー、駒種、駒数によるハッシュ値、最後のcntは升目のハッシュ値で代用している
  */
  for (Color c = WHITE; c <= BLACK; ++c)
      for (PieceType pt = PAWN; pt <= KING; ++pt)
          for (int cnt = 0; cnt < pieceCount[c][pt]; ++cnt)
              si->materialKey ^= Zobrist::psq[c][pt][cnt];

  for (Color c = WHITE; c <= BLACK; ++c)
      for (PieceType pt = KNIGHT; pt <= QUEEN; ++pt)
          si->npMaterial[c] += pieceCount[c][pt] * PieceValue[MG][pt];
}


/// Position::fen() returns a FEN representation of the position. In case of
/// Chess960 the Shredder-FEN notation is used. This is mainly a debugging function.
//内部データよりFEN stringを作り上げる
const string Position::fen() const {

  int emptyCnt;
  std::ostringstream ss;

  for (Rank r = RANK_8; r >= RANK_1; --r)
  {
      for (File f = FILE_A; f <= FILE_H; ++f)
      {
          for (emptyCnt = 0; f <= FILE_H && empty(make_square(f, r)); ++f)
              ++emptyCnt;

          if (emptyCnt)
              ss << emptyCnt;

          if (f <= FILE_H)
              ss << PieceToChar[piece_on(make_square(f, r))];
      }

      if (r > RANK_1)
          ss << '/';
  }

  ss << (sideToMove == WHITE ? " w " : " b ");

  if (can_castle(WHITE_OO))
      ss << (chess960 ? to_char(file_of(castling_rook_square(WHITE |  KING_SIDE)), false) : 'K');

  if (can_castle(WHITE_OOO))
      ss << (chess960 ? to_char(file_of(castling_rook_square(WHITE | QUEEN_SIDE)), false) : 'Q');

  if (can_castle(BLACK_OO))
      ss << (chess960 ? to_char(file_of(castling_rook_square(BLACK |  KING_SIDE)),  true) : 'k');

  if (can_castle(BLACK_OOO))
      ss << (chess960 ? to_char(file_of(castling_rook_square(BLACK | QUEEN_SIDE)),  true) : 'q');

  if (!can_castle(WHITE) && !can_castle(BLACK))
      ss << '-';

  ss << (ep_square() == SQ_NONE ? " - " : " " + to_string(ep_square()) + " ")
     << st->rule50 << " " << 1 + (gamePly - (sideToMove == BLACK)) / 2;

  return ss.str();
}


/// Position::pretty() returns an ASCII representation of the position to be
/// printed to the standard output together with the move's san notation.
/*
内部データ（board[]など）と渡された指し手情報を表示する
*/
const string Position::pretty(Move m) const {

  std::ostringstream ss;

  if (m)
      ss << "\nMove: " << (sideToMove == BLACK ? ".." : "")
         << move_to_san(*const_cast<Position*>(this), m);

  ss << "\n +---+---+---+---+---+---+---+---+\n";

  for (Rank r = RANK_8; r >= RANK_1; --r)
  {
      for (File f = FILE_A; f <= FILE_H; ++f)
          ss << " | " << PieceToChar[piece_on(make_square(f, r))];

      ss << " |\n +---+---+---+---+---+---+---+---+\n";
  }

  ss << "\nFen: " << fen() << "\nKey: " << std::hex << std::uppercase
     << std::setfill('0') << std::setw(16) << st->key << "\nCheckers: ";

  for (Bitboard b = checkers(); b; )
      ss << to_string(pop_lsb(&b)) << " ";

  ss << "\nLegal moves: ";
  for (MoveList<LEGAL> it(*this); *it; ++it)
      ss << move_to_san(*const_cast<Position*>(this), *it) << " ";

  return ss.str();
}


/// Position::check_blockers() returns a bitboard of all the pieces with color
/// 'c' that are blocking check on the king with color 'kingColor'. A piece
/// blocks a check if removing that piece from the board would result in a
/// position where the king is in check. A check blocking piece can be either a
/// pinned or a discovered check piece, according if its color 'c' is the same
/// or the opposite of 'kingColor'.

/*
pop_lsbはなにか
	LSB（最下位bit）からスキャンして1が立っているindexを返す
	indexは0から始まる
pin付けされた駒のbitboardを返す
Colorに自陣カラー kingColorに自陣カラーを指定すれば
自KINGに釘付けされた駒（自陣側の駒、敵の駒ではない）を返す
但しpinされた駒は１つだけで２つ挟まっているとpinとは判断されない
check_blockers関数がこの役割を行う
＝つまり自KINGへのcheckをblockしている駒

Colorに自陣カラー kingColorに敵陣カラーを指定すれば
敵KINGに釘付けされた駒（自陣側の駒、敵陣の駒ではない）を返す
但しpinされた駒は１つだけで２つ挟まっているとpinとは判断されない
つまり敵のKINGへの利きを邪魔している自陣側の駒を返す
discovered_check_candidates関数がこの役割を行っている
*/
Bitboard Position::check_blockers(Color c, Color kingColor) const {

  Bitboard b, pinners, result = 0;
  Square ksq = king_square(kingColor);

  // Pinners are sliders that give check when a pinned piece is removed
  /*
  KINGに影の利きをしている敵陣側ROOK,BISHOP,QUEENのbitboardをpinnersに入れている
  */
  pinners = (  (pieces(  ROOK, QUEEN) & PseudoAttacks[ROOK  ][ksq])
             | (pieces(BISHOP, QUEEN) & PseudoAttacks[BISHOP][ksq])) & pieces(~kingColor);

  while (pinners)
  {
	  //ksqと飛び駒との間の利きと全駒との＆演算してpinされた駒があるか調べている
	  //あればbに保存される
      b = between_bb(ksq, pop_lsb(&pinners)) & pieces();
	  //bが１bitしかなかったら（つまり飛び駒とKINGの間にある駒（味方、敵関係なく）が
	  //１枚しかなかったら飛び駒の利きをブロックしている駒として登録して返す
	  //どちらのカラーの駒かは引数Color cで決まる
      if (!more_than_one(b))
          result |= b & pieces(c);
  }
  return result;
}


/// Position::attackers_to() computes a bitboard of all pieces which attack a
/// given square. Slider attacks use the occ bitboard to indicate occupancy.
/*
attacks_from<>関数の概要
	指定した駒コード、盤座標にある駒の利きのbitboardを返す。ただしPAWNは進む方向と
	駒を取る利きが違う、attacks_fromはあくまで駒をとる利きのみかえす
	attacks_from関数は３つオーバーロードしている
	attacks_from(Square s)はテンプレート引数で駒種を指定して、
		指定した座標から利きのbitboardを返す。
		飛び駒だけでなく非飛び駒の利きbitboardも返せる。
	attacks_from<PAWN>(Square s, Color c)
		駒種はPAWNだけで座標とカラーが指定できて利きbitboardを返す
	attacks_from(Piece pc, Square s)
		指定した座標、指定した駒種から利きのbitboardを返す。
		非飛び駒は対応していない
	
attackers_to関数の機能
	指定した座標に利いている全ての駒（カラーに関係なく）を検出してビットを立てた
	bitboardを返す
*/
Bitboard Position::attackers_to(Square s, Bitboard occ) const {

  return  (attacks_from<PAWN>(s, BLACK) & pieces(WHITE, PAWN))
        | (attacks_from<PAWN>(s, WHITE) & pieces(BLACK, PAWN))
        | (attacks_from<KNIGHT>(s)      & pieces(KNIGHT))
        | (attacks_bb<ROOK>(s, occ)     & pieces(ROOK, QUEEN))
        | (attacks_bb<BISHOP>(s, occ)   & pieces(BISHOP, QUEEN))
        | (attacks_from<KING>(s)        & pieces(KING));
}


/// Position::legal() tests whether a pseudo-legal move is legal
/*
引数Move mが合法手か検査する合法手かどうかは
アンパッサンだったら
	用途不明
動いた駒がKINGだったら
	移動先に敵の利きが利いていたらNG、キャスリングはOK
pinがかかっていないこと,pinがかかっていてもpinがはずれない動きならOK
*/
bool Position::legal(Move m, Bitboard pinned) const {

  assert(is_ok(m));
  assert(pinned == pinned_pieces(sideToMove));

  Color us = sideToMove;
  Square from = from_sq(m);

  assert(color_of(moved_piece(m)) == us);
  assert(piece_on(king_square(us)) == make_piece(us, KING));

  // En passant captures are a tricky special case. Because they are rather
  // uncommon, we do it simply by testing whether the king is attacked after
  // the move is made.
  if (type_of(m) == ENPASSANT)
  {
      Square ksq = king_square(us);
      Square to = to_sq(m);
      Square capsq = to - pawn_push(us);
      Bitboard occ = (pieces() ^ from ^ capsq) | to;

      assert(to == ep_square());
      assert(moved_piece(m) == make_piece(us, PAWN));
      assert(piece_on(capsq) == make_piece(~us, PAWN));
      assert(piece_on(to) == NO_PIECE);

      return   !(attacks_bb<  ROOK>(ksq, occ) & pieces(~us, QUEEN, ROOK))
            && !(attacks_bb<BISHOP>(ksq, occ) & pieces(~us, QUEEN, BISHOP));
  }

  // If the moving piece is a king, check whether the destination
  // square is attacked by the opponent. Castling moves are checked
  // for legality during move generation.
  if (type_of(piece_on(from)) == KING)
      return type_of(m) == CASTLING || !(attackers_to(to_sq(m)) & pieces(~us));

  // A non-king move is legal if and only if it is not pinned or it
  // is moving along the ray towards or away from the king.
  /*
  BISHOP,ROOKが仮にPINされていてもPINを外さないうごきならOK

  */
  return   !pinned
        || !(pinned & from)
        ||  aligned(from, to_sq(m), king_square(us));
}


/// Position::pseudo_legal() takes a random move and tests whether the move is
/// pseudo legal. It is used to validate moves from TT that can be corrupted
/// due to SMP concurrent access or hash position key aliasing.
/*
用途不明
*/
bool Position::pseudo_legal(const Move m) const {

  Color us = sideToMove;
  Square from = from_sq(m);
  Square to = to_sq(m);
  Piece pc = moved_piece(m);

  // Use a slower but simpler function for uncommon cases
  if (type_of(m) != NORMAL)
      return MoveList<LEGAL>(*this).contains(m);

  // Is not a promotion, so promotion piece must be empty
  if (promotion_type(m) - 2 != NO_PIECE_TYPE)
      return false;

  // If the 'from' square is not occupied by a piece belonging to the side to
  // move, the move is obviously not legal.
  if (pc == NO_PIECE || color_of(pc) != us)
      return false;

  // The destination square cannot be occupied by a friendly piece
  if (pieces(us) & to)
      return false;

  // Handle the special case of a pawn move
  if (type_of(pc) == PAWN)
  {
      // We have already handled promotion moves, so destination
      // cannot be on the 8th/1st rank.
      if (rank_of(to) == relative_rank(us, RANK_8))
          return false;

      if (   !(attacks_from<PAWN>(from, us) & pieces(~us) & to) // Not a capture

          && !((from + pawn_push(us) == to) && empty(to))       // Not a single push

          && !(   (from + 2 * pawn_push(us) == to)              // Not a double push
               && (rank_of(from) == relative_rank(us, RANK_2))
               && empty(to)
               && empty(to - pawn_push(us))))
          return false;
  }
  else if (!(attacks_from(pc, from) & to))
      return false;

  // Evasions generator already takes care to avoid some kind of illegal moves
  // and legal() relies on this. We therefore have to take care that the same
  // kind of moves are filtered out here.
  if (checkers())
  {
      if (type_of(pc) != KING)
      {
          // Double check? In this case a king move is required
          if (more_than_one(checkers()))
              return false;

          // Our move must be a blocking evasion or a capture of the checking piece
          if (!((between_bb(lsb(checkers()), king_square(us)) | checkers()) & to))
              return false;
      }
      // In case of king moves under check we have to remove king so as to catch
      // invalid moves like b1a1 when opposite queen is on c1.
      else if (attackers_to(to, pieces() ^ from) & pieces(~us))
          return false;
  }

  return true;
}


/// Position::gives_check() tests whether a pseudo-legal move gives a check
/*
指し手が王手であればtrueを返す
*/
bool Position::gives_check(Move m, const CheckInfo& ci) const {

  assert(is_ok(m));
  assert(ci.dcCandidates == discovered_check_candidates());
  assert(color_of(moved_piece(m)) == sideToMove);

  Square from = from_sq(m);
  Square to = to_sq(m);
  PieceType pt = type_of(piece_on(from));

  // Is there a direct check?
  //駒の移動先に移動した結果王手が掛けれる位置ならtrue
  if (ci.checkSq[pt] & to)
      return true;

  // Is there a discovered check?
  /*
  移動することで王手が掛けれるならtrueを返す
  */
  if (   unlikely(ci.dcCandidates) && (ci.dcCandidates & from) && !aligned(from, to, ci.ksq))
      return true;
  /*
  指し手パターンによって判断、NORMALなら即false
  指し手パターンが成りで成った先で王手ができるようならtrue
  */
  switch (type_of(m))
  {
  case NORMAL:
      return false;

  case PROMOTION:
      return attacks_bb(Piece(promotion_type(m)), to, pieces() ^ from) & ci.ksq;

  // En passant capture with check? We have already handled the case
  // of direct checks and ordinary discovered check, so the only case we
  // need to handle is the unusual case of a discovered check through
  // the captured pawn.
  /*
  アンパッサン関係のようだが詳細不明
  */
  case ENPASSANT:
  {
      Square capsq = make_square(file_of(to), rank_of(from));
      Bitboard b = (pieces() ^ from ^ capsq) | to;

      return  (attacks_bb<  ROOK>(ci.ksq, b) & pieces(sideToMove, QUEEN, ROOK))
            | (attacks_bb<BISHOP>(ci.ksq, b) & pieces(sideToMove, QUEEN, BISHOP));
  }
  /*
  キャスリングのようだが詳細不明
  */
  case CASTLING:
  {
      Square kfrom = from;
      Square rfrom = to; // Castling is encoded as 'King captures the rook'
      Square kto = relative_square(sideToMove, rfrom > kfrom ? SQ_G1 : SQ_C1);
      Square rto = relative_square(sideToMove, rfrom > kfrom ? SQ_F1 : SQ_D1);

      return   (PseudoAttacks[ROOK][rto] & ci.ksq)
            && (attacks_bb<ROOK>(rto, (pieces() ^ kfrom ^ rfrom) | rto | kto) & ci.ksq);
  }
  default:
      assert(false);
      return false;
  }
}


/// Position::do_move() makes a move, and saves all information necessary
/// to a StateInfo object. The move is assumed to be legal. Pseudo-legal
/// moves should be filtered out before this function is called.
/*
局面を更新する関数、下の関数do_moveのオーバーライド
*/
void Position::do_move(Move m, StateInfo& newSt) {

  CheckInfo ci(*this);
  do_move(m, newSt, ci, gives_check(m, ci));
}

/*
局面を更新する唯一の関数

*/
void Position::do_move(Move m, StateInfo& newSt, const CheckInfo& ci, bool moveIsCheck) {

  assert(is_ok(m));
  assert(&newSt != st);

  /*
  展開したノード数をカウントしている
  think関数で探索が終了したあとnodes_searched関数を呼び出し
  このnodes数を表示させる
  */
  ++nodes;
  /*
  StateInfo.keyに局面のハッシュ値が記録されている
  */
  Key k = st->key;

  // Copy some fields of the old state to our new StateInfo object except the
  // ones which are going to be recalculated from scratch anyway and then switch
  // our state pointer to point to the new (ready to be updated) state.
  /*
  StateCopySize64はStateInfo構造体のなかでkeyアイテムまでのオフセット（byte単位）数を返す
  つまりStateInfo構造体の一部だけnewStにコピーする（何故全部コピーしないのかは不明）
  */
  std::memcpy(&newSt, st, StateCopySize64 * sizeof(uint64_t));
  /*
  StateInfoをつないでいる
  */
  newSt.previous = st;
  st = &newSt;

  // Update side to move
  /*
  局面のハッシュ値を更新している
  */
  k ^= Zobrist::side;

  // Increment ply counters. In particular, rule50 will be reset to zero later on
  // in case of a capture or a pawn move.
  /*
  gamePlyはゲーム手数のカウントアップ
  rule50のためのカウントアップ
  pliesFromNullは今のところ不明do_null_move関数では0に初期化する
  */
  ++gamePly;
  ++st->rule50;
  ++st->pliesFromNull;

  Color us = sideToMove;
  Color them = ~us;
  Square from = from_sq(m);
  Square to = to_sq(m);
  Piece pc = piece_on(from);
  PieceType pt = type_of(pc);
  PieceType captured = type_of(m) == ENPASSANT ? PAWN : type_of(piece_on(to));

  assert(color_of(pc) == us);
  assert(piece_on(to) == NO_PIECE || color_of(piece_on(to)) == them || type_of(m) == CASTLING);
  assert(captured != KING);

  if (type_of(m) == CASTLING)
  {
      assert(pc == make_piece(us, KING));

      Square rfrom, rto;
      do_castling<true>(from, to, rfrom, rto);

      captured = NO_PIECE_TYPE;
      st->psq += psq[us][ROOK][rto] - psq[us][ROOK][rfrom];
      k ^= Zobrist::psq[us][ROOK][rfrom] ^ Zobrist::psq[us][ROOK][rto];
  }
  /*
  capturedは取った駒の駒種、もし取る手ではなければcapturedは0
  */
  if (captured)
  {
      Square capsq = to;

      // If the captured piece is a pawn, update pawn hash key, otherwise
      // update non-pawn material.
	  /*
	  もしとった駒種がPAWNで取り方がアンパッサンなら
	  とった駒は移動先の真後ろになる（アンパッサンのルール確認）ので
	  capsqはtoではなくto+pawn_push(them)となる
	  */
      if (captured == PAWN)
      {
          if (type_of(m) == ENPASSANT)
          {
              capsq += pawn_push(them);

              assert(pt == PAWN);
              assert(to == st->epSquare);
              assert(relative_rank(us, to) == RANK_6);
              assert(piece_on(to) == NO_PIECE);
              assert(piece_on(capsq) == make_piece(them, PAWN));

              board[capsq] = NO_PIECE;
          }

          st->pawnKey ^= Zobrist::psq[them][PAWN][capsq];
      }
      else
          st->npMaterial[them] -= PieceValue[MG][captured];

      // Update board and piece lists
	  /*
	  駒を取ったことによって変更になる
	  byTypeBB,byTypeBB,byColorBBを更新
	  index[],pieceList[][][],pieceCount[][]配列を更新する関数
	  */
      remove_piece(capsq, them, captured);

      // Update material hash key and prefetch access to materialTable
	  /*
	  局面のハッシュ値から取られた駒のハッシュ値を除去している
	  同時にmaterialKeyも変更している
	  */
      k ^= Zobrist::psq[them][captured][capsq];
      st->materialKey ^= Zobrist::psq[them][captured][pieceCount[them][captured]];
      prefetch((char*)thisThread->materialTable[st->materialKey]);

      // Update incremental scores
	  /*
	  st->psqは位置評価値の集計値なので取られた駒の差分をしている
	  */
      st->psq -= psq[them][captured][capsq];

      // Reset rule 50 counter
	  /*
	  駒を取ったのでrule50は一旦キャンセルとなる
	  */
      st->rule50 = 0;
  }
  /*
  これ以降は駒を取られていない指し手の更新
  */
  // Update hash key
  /*
  移動前のハッシュ値を除去し、移動後のハッシュ値を更新している
  */
  k ^= Zobrist::psq[us][pt][from] ^ Zobrist::psq[us][pt][to];

  // Reset en passant square
  /*
  アンパッサン関係だと思うが詳細不明
  */
  if (st->epSquare != SQ_NONE)
  {
      k ^= Zobrist::enpassant[file_of(st->epSquare)];
      st->epSquare = SQ_NONE;
  }

  // Update castling rights if needed
  /*
  キャスリング関係かな、詳細不明
  */
  if (st->castlingRights && (castlingRightsMask[from] | castlingRightsMask[to]))
  {
      int cr = castlingRightsMask[from] | castlingRightsMask[to];
      k ^= Zobrist::castling[st->castlingRights & cr];
      st->castlingRights &= ~cr;
  }

  // Prefetch TT access as soon as we know the new hash key
  prefetch((char*)TT.first_entry(k));

  // Move the piece. The tricky Chess960 castling is handled earlier
  /*
  move_piece関数の機能は
  byTypeBB,byTypeBB,byColorBBの更新
  board,index[],pieceList[][][]配列の更新
  駒取りはないのでpieceCount[]配列の更新はない
  */
  if (type_of(m) != CASTLING)
      move_piece(from, to, us, pt);

  // If the moving piece is a pawn do some special extra work
  if (pt == PAWN)
  {
      // Set en-passant square if the moved pawn can be captured
	  /*
	  (int(to) ^ int(from)) == 16となるPAWNの動きは２段とびのみ
	  でかつアンパッサンできる条件が成立している場合の処理
	  */
      if (   (int(to) ^ int(from)) == 16
          && (attacks_from<PAWN>(from + pawn_push(us), us) & pieces(them, PAWN)))
      {
          st->epSquare = Square((from + to) / 2);
          k ^= Zobrist::enpassant[file_of(st->epSquare)];
      }
	  /*
	  PAWNがなる場合の処理
	  */
      else if (type_of(m) == PROMOTION)
      {
          PieceType promotion = promotion_type(m);

          assert(relative_rank(us, to) == RANK_8);
          assert(promotion >= KNIGHT && promotion <= QUEEN);
		  /*
		  一旦PAWNを除去する処理
		  */
          remove_piece(to, us, PAWN);
		  /*
		  なった駒を移動先に置く処理
		  */
          put_piece(to, us, promotion);

          // Update hash keys
		  /*
		  局面のハッシュ値をを更新、pawn専用のハッシュ値も更新
		  */
          k ^= Zobrist::psq[us][PAWN][to] ^ Zobrist::psq[us][promotion][to];
          st->pawnKey ^= Zobrist::psq[us][PAWN][to];
          st->materialKey ^=  Zobrist::psq[us][promotion][pieceCount[us][promotion]-1]
                            ^ Zobrist::psq[us][PAWN][pieceCount[us][PAWN]];

          // Update incremental score
		  //位置評価値もも更新
          st->psq += psq[us][promotion][to] - psq[us][PAWN][to];

          // Update material
		  /*
		  駒評価値も更新
		  */
          st->npMaterial[us] += PieceValue[MG][promotion];
      }

      // Update pawn hash key and prefetch access to pawnsTable
      st->pawnKey ^= Zobrist::psq[us][PAWN][from] ^ Zobrist::psq[us][PAWN][to];
      prefetch((char*)thisThread->pawnsTable[st->pawnKey]);

      // Reset rule 50 draw counter
	  /*
	  引き分け条件をクリア
	  */
      st->rule50 = 0;
  }
  /*
  これ以降は駒を取らない指し手でPAWN以外の駒種の処理、共通処理かな
  */
  // Update incremental scores
  /*
  位置評価値の更新
  */
  st->psq += psq[us][pt][to] - psq[us][pt][from];

  // Set capture piece
  /*
  とった駒種
  */
  st->capturedType = captured;

  // Update the key with the final value
  /*
  最終ハッシュ値を登録
  */
  st->key = k;

  // Update checkers bitboard: piece must be already moved due to attacks_from()
  st->checkersBB = 0;
  /*
  moveIsCheckはdo_move関数の引数の１つ
  王手の手があるならtrue
  */
  if (moveIsCheck)
  {
      if (type_of(m) != NORMAL)
          st->checkersBB = attackers_to(king_square(them)) & pieces(us);
      else
      {
          // Direct checks
		  //ci.checkSq[pt]には駒種ごとに敵KINGに王手をかけることで出来るbitboardがはいっている
		  //今回のMoveによってその場所に移動できたかをci.checkSq[pt] & toでチエックしている
		  //そしてチエックが可能であればcheckerBBに追加している
          if (ci.checkSq[pt] & to)
              st->checkersBB |= to;

          // Discovered checks
		  //ROOKまたはBISHOPではない駒が動いたことでROOK,BISHOPの利きが敵KINGに届いたのではないかチエックしている
          if (unlikely(ci.dcCandidates) && (ci.dcCandidates & from))
          {
              if (pt != ROOK)
                  st->checkersBB |= attacks_from<ROOK>(king_square(them)) & pieces(us, QUEEN, ROOK);

              if (pt != BISHOP)
                  st->checkersBB |= attacks_from<BISHOP>(king_square(them)) & pieces(us, QUEEN, BISHOP);
          }
      }
  }
  /*
  手番の変更
  */
  sideToMove = ~sideToMove;

  assert(pos_is_ok());
}


/// Position::undo_move() unmakes a move. When it returns, the position should
/// be restored to exactly the same state as before the move was made.
/*
do_move関数にくらべすごくコード量が少ない
*/
void Position::undo_move(Move m) {

  assert(is_ok(m));

  sideToMove = ~sideToMove;

  Color us = sideToMove;
  Square from = from_sq(m);
  Square to = to_sq(m);
  PieceType pt = type_of(piece_on(to));

  assert(empty(from) || type_of(m) == CASTLING);
  assert(st->capturedType != KING);

  if (type_of(m) == PROMOTION)
  {
      assert(pt == promotion_type(m));
      assert(relative_rank(us, to) == RANK_8);
      assert(promotion_type(m) >= KNIGHT && promotion_type(m) <= QUEEN);
	  /*
	  remove_piece関数は駒を取り除く時の処理
	  ひとます、なった駒をもとに戻し通常の移動の処理と共通化する
	  */
      remove_piece(to, us, promotion_type(m));
      put_piece(to, us, PAWN);
      pt = PAWN;
  }
  /*
  キャスリング関係の戻しかな
  */
  if (type_of(m) == CASTLING)
  {
      Square rfrom, rto;
      do_castling<false>(from, to, rfrom, rto);
  }
  else
  {
	  /*
	  移動の戻し（from,toをテレコにしている）
	  */
      move_piece(to, from, us, pt); // Put the piece back at the source square
      /*
	  st->capturedTypeの使い方がよくわからん
	  */
      if (st->capturedType)
      {
          Square capsq = to;

          if (type_of(m) == ENPASSANT)
          {
              capsq -= pawn_push(us);

              assert(pt == PAWN);
              assert(to == st->previous->epSquare);
              assert(relative_rank(us, to) == RANK_6);
              assert(piece_on(capsq) == NO_PIECE);
          }

          put_piece(capsq, ~us, st->capturedType); // Restore the captured piece
      }
  }

  // Finally point our state pointer back to the previous state
  st = st->previous;
  --gamePly;

  assert(pos_is_ok());
}


/// Position::do_castling() is a helper used to do/undo a castling move. This
/// is a bit tricky, especially in Chess960.
/*
キャスリング関係、詳細不明
*/
template<bool Do>
void Position::do_castling(Square from, Square& to, Square& rfrom, Square& rto) {

  bool kingSide = to > from;
  rfrom = to; // Castling is encoded as "king captures friendly rook"
  rto = relative_square(sideToMove, kingSide ? SQ_F1 : SQ_D1);
  to  = relative_square(sideToMove, kingSide ? SQ_G1 : SQ_C1);

  // Remove both pieces first since squares could overlap in Chess960
  remove_piece(Do ?  from :  to, sideToMove, KING);
  remove_piece(Do ? rfrom : rto, sideToMove, ROOK);
  board[Do ? from : to] = board[Do ? rfrom : rto] = NO_PIECE; // Since remove_piece doesn't do it for us
  put_piece(Do ?  to :  from, sideToMove, KING);
  put_piece(Do ? rto : rfrom, sideToMove, ROOK);
}


/// Position::do(undo)_null_move() is used to do(undo) a "null move": It flips
/// the side to move without executing any move on the board.
/*
ヌルムーブ用の局面更新関数
*/
void Position::do_null_move(StateInfo& newSt) {

  assert(!checkers());

  std::memcpy(&newSt, st, sizeof(StateInfo)); // Fully copy here

  newSt.previous = st;
  st = &newSt;

  if (st->epSquare != SQ_NONE)
  {
      st->key ^= Zobrist::enpassant[file_of(st->epSquare)];
      st->epSquare = SQ_NONE;
  }

  st->key ^= Zobrist::side;
  prefetch((char*)TT.first_entry(st->key));

  ++st->rule50;
  st->pliesFromNull = 0;

  sideToMove = ~sideToMove;

  assert(pos_is_ok());
}
/*
ヌルムーブ用の局面復元関数
*/
void Position::undo_null_move() {

  assert(!checkers());

  st = st->previous;
  sideToMove = ~sideToMove;
}


/// Position::see() is a static exchange evaluator: It tries to estimate the
/// material gain or loss resulting from a move.
/*
静止探索
*/
Value Position::see_sign(Move m) const {

  assert(is_ok(m));

  // Early return if SEE cannot be negative because captured piece value
  // is not less then capturing one. Note that king moves always return
  // here because king midgame value is set to 0.
  /*
  取る駒の駒価値より取られる駒の駒価値が大きければ価値10000を返す
  そうではなければ（価値の低い駒しかとれなければ）see関数を呼ぶ
  */
  if (PieceValue[MG][moved_piece(m)] <= PieceValue[MG][piece_on(to_sq(m))])
      return VALUE_KNOWN_WIN;

  return see(m);
}

/*

*/
Value Position::see(Move m) const {

  Square from, to;
  Bitboard occupied, attackers, stmAttackers;
  Value swapList[32];
  int slIndex = 1;
  PieceType captured;
  Color stm;

  assert(is_ok(m));

  from = from_sq(m);
  to = to_sq(m);
  swapList[0] = PieceValue[MG][piece_on(to)];
  stm = color_of(piece_on(from));
  occupied = pieces() ^ from;	//occupiedはformにいる駒以外の全ての駒のbitboard

  // Castling moves are implemented as king capturing the rook so cannot be
  // handled correctly. Simply return 0 that is always the correct value
  // unless in the rare case the rook ends up under attack.
  /*
  キャスリング関係なら価値0で返る
  */
  if (type_of(m) == CASTLING)
      return VALUE_ZERO;
  /*
  指し手パターンがアンパッサンならswapListにpawnを入れておく
  */
  if (type_of(m) == ENPASSANT)
  {
      occupied ^= to - pawn_push(stm); // Remove the captured pawn
      swapList[0] = PieceValue[MG][PAWN];
  }

  // Find all attackers to the destination square, with the moving piece
  // removed, but possibly an X-ray attacker added behind it.
  /*
  指した手の移動先に利いている駒のbitboardをattackersに入れる（カラーに関係なく）
  */
  attackers = attackers_to(to, occupied) & occupied;

  // If the opponent has no attackers we are finished
  /*
  stmAttackersに敵の駒だけを入れる
  もし敵の駒がなければ（つまり取り合いがなければ）
  そこで中断
  */
  stm = ~stm;
  stmAttackers = attackers & pieces(stm);
  if (!stmAttackers)
      return swapList[0];

  // The destination square is defended, which makes things rather more
  // difficult to compute. We proceed by building up a "swap list" containing
  // the material gain or loss at each stop in a sequence of captures to the
  // destination square, where the sides alternately capture, and always
  // capture with the least valuable piece. After each capture, we look for
  // new X-ray attacks from behind the capturing piece.
  captured = type_of(piece_on(from));
  /*
  slIndexは１から始まる
  swapList[]は取り合い駒リスト
  */
  do {
      assert(slIndex < 32);

      // Add the new entry to the swap list
	  /*
	  取り合いになっていくがその駒評価値をswapListに記録する
	  */
      swapList[slIndex] = -swapList[slIndex - 1] + PieceValue[MG][captured];

      // Locate and remove the next least valuable attacker
	  /*
	  まず取り合いはPAWNからおこない、to座標にきている駒を取り合い
	  その駒種を返す
	  */
      captured = min_attacker<PAWN>(byTypeBB, to, stmAttackers, occupied, attackers);

      // Stop before processing a king capture
	  /*
	  取り合いになり徐々に駒種がPAWNから上がっていくがKNIGになったらそこでやめる
	  */
      if (captured == KING)
      {
          if (stmAttackers == attackers)
              ++slIndex;

          break;
      }

      stm = ~stm;
      stmAttackers = attackers & pieces(stm);
      ++slIndex;

  } while (stmAttackers);

  // Having built the swap list, we negamax through it to find the best
  // achievable score from the point of view of the side to move.
  /*
  swapListを遡り最小の評価値を得る？
  */
  while (--slIndex)
      swapList[slIndex - 1] = std::min(-swapList[slIndex], swapList[slIndex - 1]);

  return swapList[0];
}


/// Position::is_draw() tests whether the position is drawn by material, 50 moves
/// rule or repetition. It does not detect stalemates.
/*
ドロー（引き分け）
	次の場合は、「自動的」にドローとなる。
	ステイルメイト ： 自分の手番で、自分のキングにチェックされてはいないが、合法手がない状況を指す。
	ドロー・オファー： 片方がドローを提案し、もう片方がそれを承諾した場合。
	デッド・ポジション[8]： 駒の兵力不足のため、双方が相手のキングをチェックメイトできなくなった状況を指す。次の駒の組合せの時は、たとえ敵の駒がキング一つだけであってもチェックメイトすることはできない。[9]
	キング + ビショップ1個
	キング + ナイト1個
	（キング + ナイト2個

次の場合、一方のプレーヤーの「申請（クレーム）」によりドローとなる
50手ルール ： 50手連続して両者ともポーンが動かず、またお互いに駒を取らない場合。
スリーフォールド・レピティション（同形三復）： 同一の局面が3回現れた場合。
*/

bool Position::is_draw() const {

  if ( !pieces(PAWN)
      && (non_pawn_material(WHITE) + non_pawn_material(BLACK) <= BishopValueMg))
      return true;

  if (st->rule50 > 99 && (!checkers() || MoveList<LEGAL>(*this).size()))
      return true;

  StateInfo* stp = st;
  for (int i = 2, e = std::min(st->rule50, st->pliesFromNull); i <= e; i += 2)
  {
      stp = stp->previous->previous;

      if (stp->key == st->key)
          return true; // Draw at first repetition
  }

  return false;
}


/// Position::flip() flips position with the white and black sides reversed. This
/// is only useful for debugging e.g. for finding evaluation symmetry bugs.

static char toggle_case(char c) {
  return char(islower(c) ? toupper(c) : tolower(c));
}

/*
getline(is,str,delim);
	is:文字列の抽出元となる入力ストリーム。
	str:入力ストリームから抽出した文字の読み込み先となる文字列。
	delim:行の区切り記号。
	splitの代わりみたいな関数
string.insert(index,string)
	indexの位置に文字列を挿入する
std::transform(start,end,result,func)
	startからendまでの範囲に関数funcを適用してresultに結果を返す
	関数toggle_caseは文字（文字列ではない）を受け取りそれが小文字なら
	大文字にして返す、大文字だったら小文字にして返す
*/
void Position::flip() {

  string f, token;
  std::stringstream ss(fen());

	/*
	fen文字列を逆にしている
	rnbqkbnr/pppppppp/8/8/8/8/PPPPPPPP/RNBQKBNRが入ってくる文字
	fには"RNBQKBNR/PPPPPPPP/8/8/8/8/pppppppp/rnbqkbnr "
	と反対に構築する
	*/

  for (Rank r = RANK_8; r >= RANK_1; --r) // Piece placement
  {
      std::getline(ss, token, r > RANK_1 ? '/' : ' ');
      f.insert(0, token + (f.empty() ? " " : "/"));
  }
	/*
	カラーを変えている
	*/
  ss >> token; // Active color
  f += (token == "w" ? "B " : "W "); // Will be lowercased later

  ss >> token; // Castling availability
  f += token + " ";
	/*
	大文字を小文字に、小文字を大文字に変換
	つまりWHITEをBLACKにBLACKをWHITEにする
	*/
  std::transform(f.begin(), f.end(), f.begin(), toggle_case);

  ss >> token; // En passant square
  f += (token == "-" ? token : token.replace(1, 1, token[1] == '3' ? "6" : "3"));

  std::getline(ss, token); // Half and full moves
  f += token;

  set(f, is_chess960(), this_thread());

  assert(pos_is_ok());
}


/// Position::pos_is_ok() performs some consistency checks for the position object.
/// This is meant to be helpful when debugging.
/*
デバック用のチエックするもの？
*/
bool Position::pos_is_ok(int* step) const {

  // Which parts of the position should be verified?
  const bool all = false;

  const bool testBitboards       = all || false;
  const bool testState           = all || false;
  const bool testKingCount       = all || false;
  const bool testKingCapture     = all || false;
  const bool testPieceCounts     = all || false;
  const bool testPieceList       = all || false;
  const bool testCastlingSquares = all || false;

  if (step)
      *step = 1;
	/*
	sideToMoveがWHITE,BLACKになっているか
	ちゃんとking_square配列に正しくkingの座標が入っているか
	ep_square() != SQ_NONEは用途不明
	relative_rank(sideToMove, ep_square()) != RANK_6)は用途不明
	*/
  if (   (sideToMove != WHITE && sideToMove != BLACK)
      || piece_on(king_square(WHITE)) != W_KING
      || piece_on(king_square(BLACK)) != B_KING
      || (   ep_square() != SQ_NONE
          && relative_rank(sideToMove, ep_square()) != RANK_6))
      return false;
	/*

	*/
  if (step && ++*step, testBitboards)
  {
      // The intersection of the white and black pieces must be empty
			/*
			お互いのカラーのANDをとると必ず0になるはず、成らなかったら
			おかしい
			*/
      if (pieces(WHITE) & pieces(BLACK))
          return false;

      // The union of the white and black pieces must be equal to all
      // occupied squares
			/*
			WHITEとBLACKのORは全ての駒と一緒のはず
			そうでなければおかしい
			*/
      if ((pieces(WHITE) | pieces(BLACK)) != pieces())
          return false;

      // Separate piece type bitboards must have empty intersections
			/*
			駒種が違うもの同士は重ならない
			重なったらおかしい
			*/
      for (PieceType p1 = PAWN; p1 <= KING; ++p1)
          for (PieceType p2 = PAWN; p2 <= KING; ++p2)
              if (p1 != p2 && (pieces(p1) & pieces(p2)))
                  return false;
  }
	/*
	用途不明
	*/
  if (step && ++*step, testState)
  {
      StateInfo si;
      set_state(&si);
      if (   st->key != si.key
          || st->pawnKey != si.pawnKey
          || st->materialKey != si.materialKey
          || st->npMaterial[WHITE] != si.npMaterial[WHITE]
          || st->npMaterial[BLACK] != si.npMaterial[BLACK]
          || st->psq != si.psq
          || st->checkersBB != si.checkersBB)
          return false;
  }
	/*
	std::countは、配列やコンテナの特定の要素がいくつ含まれているか返すテンプレート関数です。
	std::countの第１引数には、配列の最初を指定します。 std::countの第２引数には、配列の終わりを指定します。 std::countの第３引数には、探したい対象を指定します。
	つまりboard配列の中を調べてW_KING,B_KINGが１以外あるのはおかしいことを検出している
	*/
  if (step && ++*step, testKingCount)
      if (   std::count(board, board + SQUARE_NB, W_KING) != 1
          || std::count(board, board + SQUARE_NB, B_KING) != 1)
          return false;
	/*
	ここまで
	*/
  if (step && ++*step, testKingCapture)
      if (attackers_to(king_square(~sideToMove)) & pieces(sideToMove))
          return false;

  if (step && ++*step, testPieceCounts)
      for (Color c = WHITE; c <= BLACK; ++c)
          for (PieceType pt = PAWN; pt <= KING; ++pt)
              if (pieceCount[c][pt] != popcount<Full>(pieces(c, pt)))
                  return false;

  if (step && ++*step, testPieceList)
      for (Color c = WHITE; c <= BLACK; ++c)
          for (PieceType pt = PAWN; pt <= KING; ++pt)
              for (int i = 0; i < pieceCount[c][pt];  ++i)
                  if (   board[pieceList[c][pt][i]] != make_piece(c, pt)
                      || index[pieceList[c][pt][i]] != i)
                      return false;

  if (step && ++*step, testCastlingSquares)
      for (Color c = WHITE; c <= BLACK; ++c)
          for (CastlingSide s = KING_SIDE; s <= QUEEN_SIDE; s = CastlingSide(s + 1))
          {
              if (!can_castle(c | s))
                  continue;

              if (  (castlingRightsMask[king_square(c)] & (c | s)) != (c | s)
                  || piece_on(castlingRookSquare[c | s]) != make_piece(c, ROOK)
                  || castlingRightsMask[castlingRookSquare[c | s]] != (c | s))
                  return false;
          }

  return true;
}
