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

#include <cassert>

#include "movepick.h"
#include "thread.h"

namespace {

  /*
  着手リストの生成場合分け
  */
  enum Stages {
    MAIN_SEARCH, CAPTURES_S1, KILLERS_S1, QUIETS_1_S1, QUIETS_2_S1, BAD_CAPTURES_S1,
    EVASION,     EVASIONS_S2,
    QSEARCH_0,   CAPTURES_S3, QUIET_CHECKS_S3,
    QSEARCH_1,   CAPTURES_S4,
    PROBCUT,     CAPTURES_S5,
    RECAPTURE,   CAPTURES_S6,
    STOP
  };

  // Our insertion sort, which is guaranteed (and also needed) to be stable
  /*
  挿入ソート
  */
  void insertion_sort(ExtMove* begin, ExtMove* end)
  {
    ExtMove tmp, *p, *q;

    for (p = begin + 1; p < end; ++p)
    {
        tmp = *p;
        for (q = p; q != begin && *(q-1) < tmp; --q)
            *q = *(q-1);
        *q = tmp;
    }
  }

  // Unary predicate used by std::partition to split positive values from remaining
  // ones so as to sort the two sets separately, with the second sort delayed.
  /*
  generate_next_stage関数のみから呼ばれる
  指し手のvalue（多分指し手の仮評価値）が0以上であればtrueを返す
  */
  inline bool has_positive_value(const ExtMove& ms) { return ms.value > 0; }

  // Picks the best move in the range (begin, end) and moves it to the front.
  // It's faster than sorting all the moves in advance when there are few
  // moves e.g. possible captures.
  /*
  std::swapはmax_elementが示す要素と*beginを交換する
  std::max_elementはbeginからendの間の要素で最大もののイーサレータを返す
  ここで不明なのがExtMoveはMove,valueの２つの要素を持っているが（構造体）
  どの要素で比較するのか不明,いろいろ試験してみたがvalueで最大を検出している
  */
  inline ExtMove* pick_best(ExtMove* begin, ExtMove* end)
  {
      std::swap(*begin, *std::max_element(begin, end));
      return begin;
  }
}


/// Constructors of the MovePicker class. As arguments we pass information
/// to help it to return the (presumably) good moves first, to decide which
/// moves to return (in the quiescence search, for instance, we only want to
/// search captures, promotions and some checks) and how important good move
/// ordering is at the current node.
/*
コストラクタ
seach.cppの呼び出し
MovePicker mp(pos, ttMove, depth, History, countermoves, followupmoves, ss);
のようにsearch関数から呼び出される

メインの探索ルーチンに使用する
王手がかかっているようであればstageをEVASIONにそうでなければMAIN_SEARCHに設定しておく
まずトランスポジションテーブルの手が有効であればそれを返す。そうでなければcurとendを同じにして
next_move関数が呼ばれた時に新しく着手リストを生成する
*/
MovePicker::MovePicker(const Position& p, Move ttm, Depth d, const HistoryStats& h,Move* cm, Move* fm, Search::Stack* s) 
					: pos(p), history(h), depth(d) {

  assert(d > DEPTH_ZERO);

  cur = end = moves;
  endBadCaptures = moves + MAX_MOVES - 1;
  countermoves = cm;
  followupmoves = fm;
  ss = s;

  if (pos.checkers())
      stage = EVASION;
  else
      stage = MAIN_SEARCH;

  ttMove = (ttm && pos.pseudo_legal(ttm) ? ttm : MOVE_NONE);
  end += (ttMove != MOVE_NONE);		//ttMoveがMOVE_NONEでないなら8プラスされる
}

/*
コンストラクタ
seach.cppの呼び出し
MovePicker mp(pos, ttMove, depth, History, to_sq((ss-1)->currentMove))
他のMovepickerがsearch関数からよばれているのに対してqsearch関数から呼ばれている


王手がかかっているようであればstageをEVASIONにそうでなければ探索深度に応じて
DEPTH_QS_NO_CHECKS(-2)とりdepthが大きければQSEARCH_0(8)に設定しておく
DEPTH_QS_RECAPTURES(-10)よりdepthが大きければQSEARCH_1(11)
(探索深度がマイナスというのはどういうことだろう）
いずれでもない場合はRECAPTURE(15)に設定、ttm（トランスポジションテーブル)からの手は無視して最初から生成させる
*/
MovePicker::MovePicker(const Position& p, Move ttm, Depth d, const HistoryStats& h,Square s) 
					: pos(p), history(h), cur(moves), end(moves) {

  assert(d <= DEPTH_ZERO);

  if (pos.checkers())
      stage = EVASION;

  else if (d > DEPTH_QS_NO_CHECKS)
      stage = QSEARCH_0;

  else if (d > DEPTH_QS_RECAPTURES)
  {
      stage = QSEARCH_1;

      // Skip TT move if is not a capture or a promotion. This avoids qsearch
      // tree explosion due to a possible perpetual check or similar rare cases
      // when TT table is full.
      if (ttm && !pos.capture_or_promotion(ttm))
          ttm = MOVE_NONE;
  }
  else
  {
      stage = RECAPTURE;
      recaptureSquare = s;
      ttm = MOVE_NONE;
  }

  ttMove = (ttm && pos.pseudo_legal(ttm) ? ttm : MOVE_NONE);
  end += (ttMove != MOVE_NONE);
}

/*
コンストラクタ
search関数のなかのstep9 ProbCutのところで使われた。
MovePicker mp(pos, ttMove, History, pos.captured_piece_type())
のように呼び出されている

stageはPROBCUTの１択
*/
MovePicker::MovePicker(const Position& p, Move ttm, const HistoryStats& h, PieceType pt)
                       : pos(p), history(h), cur(moves), end(moves) {

  assert(!pos.checkers());

  stage = PROBCUT;

  // In ProbCut we generate only captures that are better than the parent's
  // captured piece.
  /*
  取った駒の評価値をcaptureThresholdに入れておく
  トランスポジションテーブルからとってきた指し手が合法手であれば（pseudo_legal関数は合法手かどうかを判定する関数だと思う詳細不明）
  */
  captureThreshold = PieceValue[MG][pt];
  ttMove = (ttm && pos.pseudo_legal(ttm) ? ttm : MOVE_NONE);
  /*
  ttMoveが駒を取らない手 OR ttMoveを指した後の静止探索評価値がすでに取られている駒評価値より小さい場合
  ttMoveを無視する。概略するとttMoveが駒を取るような価値の高い手なら採用するがそうでなければ改めて
  着手リストを生成する
  */
  if (ttMove && (!pos.capture(ttMove) || pos.see(ttMove) <= captureThreshold))
      ttMove = MOVE_NONE;

  end += (ttMove != MOVE_NONE);
}


/// score() assign a numerical value to each move in a move list. The moves with
/// highest values will be picked first.
/*
取る手の着手リストの評価値を設定する

*/
template<>
void MovePicker::score<CAPTURES>() {
  // Winning and equal captures in the main search are ordered by MVV/LVA.
  // Suprisingly, this appears to perform slightly better than SEE based
  // move ordering. The reason is probably that in a position with a winning
  // capture, capturing a more valuable (but sufficiently defended) piece
  // first usually doesn't hurt. The opponent will have to recapture, and
  // the hanging piece will still be hanging (except in the unusual cases
  // where it is possible to recapture with the hanging piece). Exchanging
  // big pieces before capturing a hanging piece probably helps to reduce
  // the subtree size.
  // In main search we want to push captures with negative SEE values to the
  // badCaptures[] array, but instead of doing it now we delay until the move
  // has been picked up in pick_move_from_list(). This way we save some SEE
  // calls in case we get a cutoff.
  Move m;

  for (ExtMove* it = moves; it != end; ++it)
  {
      m = it->move;
	  /*
	  評価値は取る駒の駒評価値から駒種をValueにキャストした数を引いている
	  引いている意味は不明
	  */
      it->value =  PieceValue[MG][pos.piece_on(to_sq(m))]
                 - Value(type_of(pos.moved_piece(m)));
	  /*
	  もし指し手パターンがアンパッサンだったら
	  PAWNの評価値を追加する
	  （アンパッサンの時のPAWNの動きではto座標に敵の駒がいないので）
	  */
      if (type_of(m) == ENPASSANT)
          it->value += PieceValue[MG][PAWN];
      /*
	  指し手パターンが成りだったら、なった駒種の駒評価値からPAWNの駒評価値を引いたもの（昇格評価値）を追加している
	  */
      else if (type_of(m) == PROMOTION)
          it->value += PieceValue[MG][promotion_type(m)] - PieceValue[MG][PAWN];
  }
}

/*
穏やかな手の着手リストの評価値をhistory[][]変数で初期化している
history[][]は駒種別、升目別の駒訪問履歴のようなもので駒がたくさん
その升目に来ている升は得点が高い
*/
template<>
void MovePicker::score<QUIETS>() {

  Move m;

  for (ExtMove* it = moves; it != end; ++it)
  {
      m = it->move;
      it->value = history[pos.moved_piece(m)][to_sq(m)];
  }
}

/*
MVV-LVAはオーダリングの手法のようだ
https://chessprogramming.wikispaces.com/MVV-LVA
意味はわからない
この関数は王手を回避する手の評価を行っている
静止探索の結果がマイナスなら評価値を-2000をつける
そうではなく静止探索の結果がプラスなら最初の取った駒の駒評価値＋2000
駒を取らない手ならhistory[][]を評価値にする
*/
template<>
void MovePicker::score<EVASIONS>() {
  // Try good captures ordered by MVV/LVA, then non-captures if destination square
  // is not under attack, ordered by history value, then bad-captures and quiet
  // moves with a negative SEE. This last group is ordered by the SEE value.
  Move m;
  Value see;

  for (ExtMove* it = moves; it != end; ++it)
  {
      m = it->move;
      if ((see = pos.see_sign(m)) < VALUE_ZERO)
          it->value = see - HistoryStats::Max; // At the bottom

      else if (pos.capture(m))
          it->value =  PieceValue[MG][pos.piece_on(to_sq(m))]
                     - Value(type_of(pos.moved_piece(m))) + HistoryStats::Max;
      else
          it->value = history[pos.moved_piece(m)][to_sq(m)];
  }
}


/// generate_next_stage() generates, scores and sorts the next bunch of moves,
/// when there are no more moves to try for the current stage.
/*
next_move関数からよばれる各stageによって生成する指し手のパターンが異なる

*/
void MovePicker::generate_next_stage() {

  cur = moves;

  switch (++stage) {
  /*
  取る手をgenerate関数で生成して、score<CAPTURES>()関数で指し手の評価をしている
  */
  case CAPTURES_S1: case CAPTURES_S3: case CAPTURES_S4: case CAPTURES_S5: case CAPTURES_S6:
      end = generate<CAPTURES>(pos, moves);
      score<CAPTURES>();
      return;
  /*
  キラー手を生成する
  */
  case KILLERS_S1:
      cur = killers;
      end = cur + 2;
	  /*なぜkillers[0,3,4,5]にMOVE_NONEを入れる*/
      killers[0].move = ss->killers[0];
      killers[1].move = ss->killers[1];
      killers[2].move = killers[3].move = MOVE_NONE;
      killers[4].move = killers[5].move = MOVE_NONE;
	  /*
	  ここから先詳細不明
	  */
      // Please note that following code is racy and could yield to rare (less
      // than 1 out of a million) duplicated entries in SMP case. This is harmless.

      // Be sure countermoves are different from killers
      for (int i = 0; i < 2; ++i)
          if (   countermoves[i] != (cur+0)->move && countermoves[i] != (cur+1)->move)
              (end++)->move = countermoves[i];

      // Be sure followupmoves are different from killers and countermoves
      for (int i = 0; i < 2; ++i)
          if (   followupmoves[i] != (cur+0)->move
              && followupmoves[i] != (cur+1)->move
              && followupmoves[i] != (cur+2)->move
              && followupmoves[i] != (cur+3)->move)
              (end++)->move = followupmoves[i];
      return;
  /*

  */
  case QUIETS_1_S1:
	  /*
	  generate<QUIETS>(pos, moves)は駒を取ったりしない穏やかな指し手を生成する
	  */
      endQuiets = end = generate<QUIETS>(pos, moves);
	  /*
	  上のgenerate<QUIETS>で生成した着手リストはプライベート変数moves[]に入っている
	  その着手リストのvalueをhistory[][]変数で初期化しておく
	  history[][]変数は駒種ごと、升目ごとで駒が移動するたびに加算される（訪問履歴評価といったところかな）
	  */
      score<QUIETS>();
	  /*
	  std::partition関数はhas_positive_valueがtrueを返す要素とfalseを返す要素を分けて
	  trueとなる要素を最初に並べ、falseとなる要素をその後ろに並べる
	  falseとなる要素の先頭のイーサレータを返す
	  真偽を判定する関数has_positive_valueは指し手の仮評価値が0以上ならtrueを返す
	  insertion_sort関数はインサートソート
	  endは有効な着手数だけ後ろに伸びている
	  */
      end = std::partition(cur, end, has_positive_value);
      insertion_sort(cur, end);
      return;

  case QUIETS_2_S1:
	  /*
	  いままで試してきた手は全て無視してリセット？
	  */
      cur = end;
      end = endQuiets;
      if (depth >= 3 * ONE_PLY)
          insertion_sort(cur, end);
      return;

  case BAD_CAPTURES_S1:
      // Just pick them in reverse order to get MVV/LVA ordering
      cur = moves + MAX_MOVES - 1;
      end = endBadCaptures;
      return;
  /*
  王手がかかっている場合、生成される着手リスト
  */
  case EVASIONS_S2:
      end = generate<EVASIONS>(pos, moves);
      if (end > moves + 1)
          score<EVASIONS>();
      return;

  case QUIET_CHECKS_S3:
      end = generate<QUIET_CHECKS>(pos, moves);
      return;

  case EVASION: case QSEARCH_0: case QSEARCH_1: case PROBCUT: case RECAPTURE:
      stage = STOP;
      /* Fall through */

  case STOP:
      end = cur + 1; // Avoid another next_phase() call
      return;

  default:
      assert(false);
  }
}


/// next_move() is the most important method of the MovePicker class. It returns
/// a new pseudo legal move every time it is called, until there are no more moves
/// left. It picks the move with the biggest value from a list of generated moves
/// taking care not to return the ttMove if it has already been searched.
/*
search関数から最初に呼ばれる。
curがmoves配列の先頭を指しており、endがmoves配列の最後を指している
genmove.cppで手を生成するとendが生成したかずだけのびる
生成する手がないとcurとendは同じになり、その都度stageが上がっていき
生成するパターンが変化していく

生成する順番はMAIN_SEARCHからスタートしBAD_CAPTURES_S1->KILLERS_S1 ... BAD_CAPTURES_S1と生成するパターンを変えながら
手をかえしていく。next_move関数を呼び出しているsearch関数の終了条件はMOVE_NONEなのでnext_move関数が終了するのはnext_move関数の
case STOPにいくときになる。
BAD_CAPTURES_S1からの手が全て終了するとgenerate_next_stage関数でstageがあがりEVASIONになるとcase EVASIONに辿りつき
stageはSTOPになりそのままcase STOPに滑り落ちる(break文がないので）curとendをちかうようにして
next_move関数に返るcur == endが成立しないのでswitch文に移動しnext_move関数内のcase STOPに移動してMOVE_NONEを返し
search関数からの呼び出しが終了する

この１行づつが生成パターンの呼び出し順番になる
MAIN_SEARCHはラベルでCAPTURES_S1からBAD_CAPTURES_S1まで生成し終了する

MAIN_SEARCH, CAPTURES_S1, KILLERS_S1, QUIETS_1_S1, QUIETS_2_S1, BAD_CAPTURES_S1,

EVASIONの生成パターンで呼び出されたらEVASIONS_S2で生成する
EVASION,     EVASIONS_S2,
QSEARCH_0の生成パターンで呼び出されたらCAPTURES_S3, QUIET_CHECKS_S3で生成する
QSEARCH_0,   CAPTURES_S3, QUIET_CHECKS_S3,
QSEARCH_1の生成パターンで呼び出されたらCAPTURES_S4で生成する
QSEARCH_1,   CAPTURES_S4,
PROBCUTの生成パターンはProbCutのとき呼ばれる
PROBCUT,     CAPTURES_S5,
RECAPTUREの生成パターンで呼び出されたらCAPTURES_S6で生成する
RECAPTURE,   CAPTURES_S6,
MOVE_NONEを返す
STOP
};


*/
template<>
Move MovePicker::next_move<false>() {

  Move move;

  while (true)
  {
	  /*
	  cur == endが成立するのは指し手が0ということ
	  */
      while (cur == end)
          generate_next_stage();

      switch (stage) {
	  /*
	  ここにくるということは手を新たに生成していない、ttMoveで十分ということなのでttMoveをかえす
	  */
      case MAIN_SEARCH: case EVASION: case QSEARCH_0: case QSEARCH_1: case PROBCUT:
          ++cur;
          return ttMove;
	  /*
	  CAPTURES_S1で手が生成できたとき（他の指し手は生成していない）ここにくる
	  pick_best関数でもっとも点数のよかった手を返す、
	  静止探索して評価値が0を下回らない（駒の取り合いに勝つか引き分け）
	  静止探索してマイナスになるようなら取る手の着手リストはendBadCapturesの指す
	  位置に移動し、最初から着手リストがなくなるまでつづける（cur == end）が成立するまで
	  */
      case CAPTURES_S1:
          move = pick_best(cur++, end)->move;
          if (move != ttMove)
          {
              if (pos.see_sign(move) >= VALUE_ZERO)
                  return move;

              // Losing capture, move it to the tail of the array
              (endBadCaptures--)->move = move;
          }
          break;

      case KILLERS_S1:
          move = (cur++)->move;
          if (    move != MOVE_NONE
              &&  move != ttMove
              &&  pos.pseudo_legal(move)
              && !pos.capture(move))
              return move;
          break;

      case QUIETS_1_S1: case QUIETS_2_S1:
          move = (cur++)->move;
          if (   move != ttMove
              && move != killers[0].move
              && move != killers[1].move
              && move != killers[2].move
              && move != killers[3].move
              && move != killers[4].move
              && move != killers[5].move)
              return move;
          break;

      case BAD_CAPTURES_S1:
          return (cur--)->move;

      case EVASIONS_S2: case CAPTURES_S3: case CAPTURES_S4:
          move = pick_best(cur++, end)->move;
          if (move != ttMove)
              return move;
          break;

      case CAPTURES_S5:
           move = pick_best(cur++, end)->move;
           if (move != ttMove && pos.see(move) > captureThreshold)
               return move;
           break;

      case CAPTURES_S6:
          move = pick_best(cur++, end)->move;
          if (to_sq(move) == recaptureSquare)
              return move;
          break;

      case QUIET_CHECKS_S3:
          move = (cur++)->move;
          if (move != ttMove)
              return move;
          break;

      case STOP:
          return MOVE_NONE;

      default:
          assert(false);
      }
  }
}


/// Version of next_move() to use at split point nodes where the move is grabbed
/// from the split point's shared MovePicker object. This function is not thread
/// safe so must be lock protected by the caller.
template<>
Move MovePicker::next_move<true>() { return ss->splitPoint->movePicker->next_move<false>(); }
