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
�p�r�s��
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
�p�r�s��
�i�ǋL�j
position�͋ǖʂ�ێ�����N���X�ŁAdo_move�ɂ����
�ύX���������邪�i�ǖʍX�V�j���Ƃ̋ǖʂɖ߂�����
�g�p�����������Ă����\���̂���
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
offsetof�̋@�\
    �\���̂̃����o�̃o�C�g�ʒu������������Ԃ�
�p�r�s��
*/
const size_t StateCopySize64 = offsetof(StateInfo, key) / sizeof(uint64_t) + 1;


/// The Position class stores the information regarding the board representation
/// like pieces, side to move, hash keys, castling info, etc. The most important
/// methods are do_move() and undo_move(), used by the search to update node info
/// when traversing the search tree.
/*
�ǖʂ�\������N���X
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
  fenStr�͋ǖʂ𕶎���ŕ\����������
	���၄
	"rnbqkbnr/pppppppp/8/8/8/8/PPPPPPPP/RNBQKBNR w KQkq - 0 1",
	R--rook,N--night,B--bishop,Q--queen,K--king,P--pawn,�啶����white�i���j
	r--rook,n--night,b--bishop,q--queen,k--king,p--pawn,��������black�i���j
	�����͋󔒂̐�,/�͍s�̏I��� �ǖʂ�\�����镶����̂��Ƌ󔒂����Ă��̋ǖʂŎ���
	�w���J���[��w/b�ŕ\���A���̎���KQkq�͕s�� -���s���@0 1���s��
	�i�ǋL�j
	KQkq�̓L���X�����O�Ɋ֌W����Ȃɂ�
  */
	/*
	set�֐���Position�N���X�̃R���X�g���N�^����Ă΂�
	fenStr����͂��ē����f�[�^���X�V���Ă���B
	�X�V���������f�[�^��
		board[s]
		byTypeBB[ALL_PIECES]
	  byColorBB[c]
		index[s]
		pieceList[c][pt][index[s]]
	���ɂ����낢��ǖʕێ��A�ǖʍX�V�ɕK�v�Ȃ��̂����������Ă��邪
	�ڍוs��
	*/
  void set(const std::string& fenStr, bool isChess960, Thread* th);
	/*
	���݂̋ǖʂ�fenStr������ɕϊ�����
	*/
  const std::string fen() const;
	/*
	�w������ƌ��݂̋ǖʂ𕶎���ɂ��ĕԂ�
	*/
  const std::string pretty(Move m = MOVE_NONE) const;

  // Position representation
	/*
	�p�r�s��
	�i�ǋL�j
	�S�Ă̋�̃r�b�g��������bitboard��Ԃ�
	�֐��̐錾�̂��Ƃɂ��Ă���const�͂���
	position�N���X�̃����o�[�ϐ��͕ύX�ł��Ȃ�
	���Ƃ������Ă���B�ibitboard��ǂ݂����ĕԂ�
	�֐��Ȃ̂ŕύX���邱�Ƃ��Ȃ��j
	*/
  Bitboard pieces() const;
	/*
	�w�肵����킪������bitboard��Ԃ�
	*/
  Bitboard pieces(PieceType pt) const;
	/*
	pt1,pt2�Ŏw�肵����킪������bitboard��Ԃ�
	*/
  Bitboard pieces(PieceType pt1, PieceType pt2) const;
	/*
	�w�肵���J���[�̋�����Ă���bitboard��Ԃ�
	*/
  Bitboard pieces(Color c) const;
	/*
	�w�肵�����A�w�肵���J���[�̋�����Ă���bitboard��Ԃ�
	*/
  Bitboard pieces(Color c, PieceType pt) const;
	/*
	�w�肵�����pt1,pt2�A�w�肵���J���[�̋�����Ă���bitboard��Ԃ�
	*/
  Bitboard pieces(Color c, PieceType pt1, PieceType pt2) const;
	/*
	�����o�[�ϐ�board[64]�̍��W���w�肵�ċ�R�[�h��Ԃ�
	*/
  Piece piece_on(Square s) const;
	/*
	�J���[���w�肵��pieceList[COLOR][PIECE_TYPE][16]�ϐ�����king�̍��W��Ԃ�
	*/
  Square king_square(Color c) const;
	/*
	�p�r�s��
	*/
  Square ep_square() const;
	/*
	�w�肵�����W�ɋ�Ȃ�������true��Ԃ�
	*/
  bool empty(Square s) const;
	/*
	�J���[�Ƌ����w�肷���pieceCount[COLOR][PIECE_TYPE]��
	�Ԃ��BpieceCount�z�񂪂Ȃ�Ȃ̂��������Ă��Ȃ���
	�������A�J���[���Ƃ̋���Ǝv����
	�Ƃ����template<PieceType Pt>�Ƃ͂Ȃ񂾂낤
	�����炭��킲�Ƃ�count�֐����e���v���[�g�����Ă���̂ł�
	���ꂾ�������ɓ��삷��H
	*/
  template<PieceType Pt> int count(Color c) const;
	/*
	�J���[�Ƌ����w�肷���pieceList[COLOR][PIECE_TYPE]��
	���W���X�g�i�z��j��Ԃ��Ă���
	���̊֐�����킲�ƂɃe���v���[�g������Ă���H
	*/
  template<PieceType Pt> const Square* list(Color c) const;

  // Castling
	/*
	�p�r�s��
	�����A�L���X�e�C���O�̂��Ƃ��Ǝv��
	*/
  int can_castle(Color c) const;
	/*
	�p�r�s��
	�����A�L���X�e�C���O�̂��Ƃ��Ǝv��
	*/
  int can_castle(CastlingRight cr) const;
	/*
	�p�r�s��
	�����A�L���X�e�C���O�̂��Ƃ��Ǝv��
	*/
  bool castling_impeded(CastlingRight cr) const;
	/*
	�p�r�s��
	�����A�L���X�e�C���O�̂��Ƃ��Ǝv��
	*/
  Square castling_rook_square(CastlingRight cr) const;

  // Checking
	/*
	�p�r�s��
	StateInfo�̃����o�[checkersBB��Ԃ��Ă���
	*/
  Bitboard checkers() const;
	/*
	�p�r�s��
	�����Aking�ւ̃`�G�b�N�Ɋւ���Ȃɂ�
	*/
  Bitboard discovered_check_candidates() const;
	/*
	�p�r�s��
	�����Aking�ւ̃`�G�b�N�Ɋւ���Ȃɂ�
	*/
  Bitboard pinned_pieces(Color c) const;

  // Attacks to/from a given square
	/*
	�p�r�s��
	�����v�Z���Ƃ�����
	*/
  Bitboard attackers_to(Square s) const;
  Bitboard attackers_to(Square s, Bitboard occ) const;
  Bitboard attacks_from(Piece pc, Square s) const;
  template<PieceType> Bitboard attacks_from(Square s) const;
  template<PieceType> Bitboard attacks_from(Square s, Color c) const;

  // Properties of moves
	/*
	�n���ꂽ�肪���@�肩�`�G�b�N����
	*/
  bool legal(Move m, Bitboard pinned) const;
	/*
	�p�r�s��
	*/
  bool pseudo_legal(const Move m) const;
	/*
	���肪�������Ȃ�true��Ԃ�
	*/
  bool capture(Move m) const;
	/*
	�p�r�s��
	����̎�ʂ��`�G�b�N���Ă��邪
	�֐������画�f����ɋ���Ƃ��
	�������͐����𔻒肵�Ă���悤����
	�H
	*/
  bool capture_or_promotion(Move m) const;
	/*
	�p�r�s��
	���@��̔���H
	*/
  bool gives_check(Move m, const CheckInfo& ci) const;
	/*
	�p�r�s��
	PAWN�Ɋւ��邱�Ƃ̂悤�ł��邪�s��
	*/
  bool advanced_pawn_push(Move m) const;
	/*
	����f�[�^�����R�[�h���擾����
	*/
  Piece moved_piece(Move m) const;
	/*
	�p�r�s��
	�����H
	*/
  PieceType captured_piece_type() const;

  // Piece specific
	/*
	�p�r�s��
	*/
  bool pawn_passed(Color c, Square s) const;
	/*
	PAWN������w�n�̍ŉ��i�ɂ��邩�`�G�b�N�A������
	true��Ԃ��APAWN�����邽�߂̃`�G�b�N�Ɏg�p���Ă���̂����H
	*/
  bool pawn_on_7th(Color c) const;
	/*
	�p�r�s��
	�����J���[��bishop���Q�ȏ゠���āi�܂�Ƃ��Ă��Ȃ��j
	board�̃J���[�i�s���͗l�̐F�j���قȂ��Ă����true
	�����������m��bishop�݂͌��ɈقȂ�board�J���[�ɔz�u�����̂�
	true���Ԃ��Ă���͓̂�����O�̂悤�ȋC�����邪�Achess960�ł�
	�Ⴄ�̂����H
	*/
  bool bishop_pair(Color c) const;
	/*
	WHITE�ABLACK����bishop���݂��Ⴂ��board�J���[�Ɉʒu�����
	true��Ԃ�
	�܂�Ⴄboard�J���[�̏ꍇ���݂��Ɏ�荇�����Ƃ͂Ȃ���
	���f�ł���
	*/
  bool opposite_bishops() const;

  // Doing and undoing moves
  /*
	�ǖʍX�V�A�ڍוs��
	�������炢�낢�������ǉ����ĉ���
	do_move���Ă�ł��郉�b�p�[
	*/
	void do_move(Move m, StateInfo& st);
	/*
	�ǖʍX�V�A�ڍוs��
	*/
  void do_move(Move m, StateInfo& st, const CheckInfo& ci, bool moveIsCheck);
	/*
	�ǖʕ����A�ڍוs��
	*/
  void undo_move(Move m);
	/*
	null move�𓮂���,�ڍוs��
	*/
  void do_null_move(StateInfo& st);
	/*
	do_null_move�̕���
	*/
  void undo_null_move();

  // Static exchange evaluation
	/*
	�p�r�s��
	SEE�̖��O���炵�ĐÎ~�T���H
	*/
  Value see(Move m) const;
	/*
	�p�r�s��
	*/
  Value see_sign(Move m) const;

  // Accessing hash keys
	/*
	�p�r�s��
	StateInfo�̃����o�[key��Ԃ�����
	*/
  Key key() const;
	/*
	�p�r�s��
	*/
  Key exclusion_key() const;
	/*
	�p�r�s��
	StateInfo�̃����o�[pawnKey��Ԃ�����
	*/
  Key pawn_key() const;
	/*
	�p�r�s��
	StateInfo�̃����o�[materialKey��Ԃ�����
	*/
  Key material_key() const;

  // Incremental piece-square evaluation
	/*
	�p�r�s��
	StateInfo�̃����o�[psq��Ԃ�
	*/
  Score psq_score() const;
	/*
	�p�r�s��
	StateInfo�̃����o�[npMaterial
	*/
  Value non_pawn_material(Color c) const;

  // Other properties of the position
	/*
	�����o�[�ϐ�sideToMove��Ԃ�
	*/
  Color side_to_move() const;
	/*
	�Q�[���̉���ڂ���Ԃ�
	*/
  int game_ply() const;
	/*
	�`�F�X960�iChess 960�j�́A�ϑ��`�F�X�̈��
	chess960���ǂ�����Ԃ�
	*/
  bool is_chess960() const;
	/*
	�p�r�s��
	�����o�[�ϐ�thisThread��Ԃ�
	*/
  Thread* this_thread() const;
	/*
	�p�r�s��
	�����o�[�ϐ�nodes��Ԃ�
	*/
  uint64_t nodes_searched() const;
	/*
	�p�r�s��
	*/
  void set_nodes_searched(uint64_t n);
	/*
	�p�r�s��
	���������̔�������Ă���H
	*/
  bool is_draw() const;

  // Position consistency check, for debugging
  bool pos_is_ok(int* step = NULL) const;
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
����f�[�^�����f�[�^���擾����
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
�Q�[���̉���ڂ���Ԃ�
*/
inline int Position::game_ply() const {
  return gamePly;
}

/*
WHITE�ABLACK����bishop���݂��Ⴂ��board�J���[�Ɉʒu�����
true��Ԃ�
�܂�Ⴄboard�J���[�̏ꍇ���݂��Ɏ�荇�����Ƃ͂Ȃ���
���f�ł���
*/
inline bool Position::opposite_bishops() const {

  return   pieceCount[WHITE][BISHOP] == 1
        && pieceCount[BLACK][BISHOP] == 1
        && opposite_colors(pieceList[WHITE][BISHOP][0], pieceList[BLACK][BISHOP][0]);
}

/*
�p�r�s��
�����J���[��bishop���Q�ȏ゠���āi�܂�Ƃ��Ă��Ȃ��j
board�̃J���[�i�s���͗l�̐F�j���قȂ��Ă����true
�����������m��bishop�݂͌��ɈقȂ�board�J���[�ɔz�u�����̂�
true���Ԃ��Ă���͓̂�����O�̂悤�ȋC�����邪�Achess960�ł�
�Ⴄ�̂����H
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
