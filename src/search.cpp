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

/*
�p��W
http://misakirara.s296.xrea.com/misaki/words.html
*/

#include <algorithm>
#include <cassert>
#include <cfloat>
#include <cmath>
#include <cstring>
#include <iostream>
#include <sstream>

#include "book.h"
#include "evaluate.h"
#include "movegen.h"
#include "movepick.h"
#include "notation.h"
#include "search.h"
#include "timeman.h"
#include "thread.h"
#include "tt.h"
#include "ucioption.h"

namespace Search {
  /*
  Signals.stop�͒T�����~�߂�t���O
  stopOnPonderhit�p�r�s��
  firstRootMove�p�r�s��
  failedLowAtRoot;
	stackfish�͔����[���{Window�{alpha-beta�T�������Ă���Ǝv����
	failedLowAtRoot��Winodw�T����Low���s�ɂȂ��true�ɂȂ�
	�����ݒ��start_thinking�֐�����false�ɐݒ�
  */
  volatile SignalsType Signals;
  LimitsType Limits;
  /*
  ���[�g
  */
  std::vector<RootMove> RootMoves;
  /*
  ���[�g�ǖ�
  */
  Position RootPos;
  /*
  ���[�g�ł̎��
  */
  Color RootColor;
  /*
  �T���ɗv�������Ԃ��~���Z�R���h�Ōv��
  start_thinking�֐����Ō��ݎ��Ԃ����Ă����A�T���I����̎��Ԃƍ����������邱�ƂŌo�ߎ��Ԃ��v������
  */
  Time::point SearchTime;
  /*
  �p�r�s��
  */
  StateStackPtr SetupStates;
}

using std::string;
using Eval::evaluate;
using namespace Search;

namespace {

  // Set to true to force running with one thread. Used for debugging
  const bool FakeSplit = false;

  // Different node types, used as template parameter
  /*
  search�֐��̃e���v���[�g����
  */
  enum NodeType { Root, PV, NonPV };

  // Dynamic razoring margin based on depth
  /*
  �p�r�s��
  */
  inline Value razor_margin(Depth d) { return Value(512 + 16 * d); }

  // Futility lookup tables (initialized at startup) and their access functions
  /*
  �p�r�s��
  */
  int FutilityMoveCounts[2][32]; // [improving][depth]
  /*
  �}����iFutility Pruning�j�̂Ƃ��̃}�[�W�������߂�
  */
  inline Value futility_margin(Depth d) {
    return Value(100 * d);
  }

  // Reduction lookup tables (initialized at startup) and their access function
  /*
  Reduction�i�k���H�j
  �p�r�s��
  */
  int8_t Reductions[2][2][64][64]; // [pv][improving][depth][moveNumber]
  
  /*
  �p�r�s��
  */
  template <bool PvNode> inline Depth reduction(bool i, Depth d, int mn) {

    return (Depth) Reductions[PvNode][i][std::min(int(d) / ONE_PLY, 63)][std::min(mn, 63)];
  }

  /*
  �p�r�s��
  */
  size_t MultiPV, PVIdx;
  /*
  ���Ԑ���H
  */
  TimeManager TimeMgr;
  /*
  �p�r�s��
  */
  double BestMoveChanges;
  /*
  �p�r�s��
  */
  Value DrawValue[COLOR_NB];
  /*
  History�̓N���X�Ńv���C�x�[�g�ϐ���Value table[pieceType][SQ]�������Ă���
  �ŏ���id_loop�֐�����0�N���A��update�֐��ōX�V����
  ��ړ�������̍��W�ɓ��_���^�������̈ʒu�]���ő����̋�ړ�����قǍ����_
  ��̈ړ������̂悤�Ȃ���
  */
  HistoryStats History;
  /*
  �p�r�s��
  */
  GainsStats Gains;
  /*
  �p�r�s��
  */
  MovesStats Countermoves, Followupmoves;
  
  /*
  ��ʒT���֐�
  */
  template <NodeType NT, bool SpNode>
  Value search(Position& pos, Stack* ss, Value alpha, Value beta, Depth depth, bool cutNode);

  /*
  �������[�p�T���֐�
  */
  template <NodeType NT, bool InCheck>
  Value qsearch(Position& pos, Stack* ss, Value alpha, Value beta, Depth depth);

  /*
  �������[�v�������Ă��肱��id_loop�֐�����search�֐����Ă�
  idle_loop�֐�->think�֐�->id_loop�֐�->search�֐��ƌĂ΂��悤�ɂȂ��Ă���
  main�֐���Threads.init()���Ă��new_thread�֐�->thread_create�֐�->start_routine�֐�->idle_loop�֐��ň�U
  sleep��ԂɑJ�ڂ���
  UCI����̃R�}���hgo�ɂ��start_thking�֐�����sleep��Ԃ�������idle_loop�֐�����T�����J�n�����
  */
  void id_loop(Position& pos);
  /*
  �p�r�s��
  */
  Value value_to_tt(Value v, int ply);
  /*
  �p�r�s��
  */
  Value value_from_tt(Value v, int ply);
  /*
  killer�Ƃ�History�z����X�V����
  �X�V�̃^�C�~���O�͂킩���Ă��Ȃ�
  */
  void update_stats(const Position& pos, Stack* ss, Move move, Depth depth, Move* quiets, int quietsCnt);
  /*
  �p�r�s��
  */
  string uci_pv(const Position& pos, int depth, Value alpha, Value beta);
  /*
  �p�r�s��
  */
  struct Skill {
    Skill(int l) : level(l), best(MOVE_NONE) {}
   ~Skill() {
      if (enabled()) // Swap best PV line with the sub-optimal one
          std::swap(RootMoves[0], *std::find(RootMoves.begin(),
                    RootMoves.end(), best ? best : pick_move()));
    }

    bool enabled() const { return level < 20; }
    bool time_to_pick(int depth) const { return depth == 1 + level; }
    Move pick_move();

    int level;
    Move best;
  };

} // namespace


/// Search::init() is called during startup to initialize various lookup tables
/*
main�֐�����Ă΂�Ă���
search�n�̏�����
*/
void Search::init() {

  int d;  // depth (ONE_PLY == 2)
  int hd; // half depth (ONE_PLY == 1)
  int mc; // moveCount

  // Init reductions array
  for (hd = 1; hd < 64; ++hd) for (mc = 1; mc < 64; ++mc)
  {
      double    pvRed = 0.00 + log(double(hd)) * log(double(mc)) / 3.00;
      double nonPVRed = 0.33 + log(double(hd)) * log(double(mc)) / 2.25;
      Reductions[1][1][hd][mc] = int8_t(   pvRed >= 1.0 ?    pvRed * int(ONE_PLY) : 0);
      Reductions[0][1][hd][mc] = int8_t(nonPVRed >= 1.0 ? nonPVRed * int(ONE_PLY) : 0);

      Reductions[1][0][hd][mc] = Reductions[1][1][hd][mc];
      Reductions[0][0][hd][mc] = Reductions[0][1][hd][mc];

      if (Reductions[0][0][hd][mc] > 2 * ONE_PLY)
          Reductions[0][0][hd][mc] += ONE_PLY;

      else if (Reductions[0][0][hd][mc] > 1 * ONE_PLY)
          Reductions[0][0][hd][mc] += ONE_PLY / 2;
  }

  // Init futility move count array
  for (d = 0; d < 32; ++d)
  {
      FutilityMoveCounts[0][d] = int(2.4 + 0.222 * pow(d + 0.00, 1.8));
      FutilityMoveCounts[1][d] = int(3.0 + 0.300 * pow(d + 0.98, 1.8));
  }
}


/// Search::perft() is our utility to verify move generation. All the leaf nodes
/// up to the given depth are generated and counted and the sum returned.
/*
�����̊֐�Search::perft����Ă΂��
�W�J�ł���m�[�h�̐���Ԃ��B
benchmark����g�p�����
*/
static uint64_t perft(Position& pos, Depth depth) {

  StateInfo st;
  uint64_t cnt = 0;
  CheckInfo ci(pos);
  const bool leaf = depth == 2 * ONE_PLY;

  for (MoveList<LEGAL> it(pos); *it; ++it)
  {
      pos.do_move(*it, st, ci, pos.gives_check(*it, ci));
      cnt += leaf ? MoveList<LEGAL>(pos).size() : ::perft(pos, depth - ONE_PLY);
      pos.undo_move(*it);
  }
  return cnt;
}

uint64_t Search::perft(Position& pos, Depth depth) {
  return depth > ONE_PLY ? ::perft(pos, depth) : MoveList<LEGAL>(pos).size();
}

/// Search::think() is the external interface to Stockfish's search, and is
/// called by the main thread when the program receives the UCI 'go' command. It
/// searches from RootPos and at the end prints the "bestmove" to output.

/*
idle_loop�֐�->think�֐�->id_loop�֐�->search�֐��ƌĂ΂��悤�ɂȂ��Ă���
�O���[�o���ϐ���RootPos�ϐ���wait_for_think_finished�֐��Ō��݂̋ǖʂ��R�s�[���Ă�����Ă���


*/
void Search::think() {

  static PolyglotBook book; // Defined static to initialize the PRNG only once

  RootColor = RootPos.side_to_move();
  /*
  ���Ԑ���H
  */
  TimeMgr.init(Limits, RootPos.game_ply(), RootColor);

  /*
  �����������̕]���l
  Options["Contempt Factor"]�̓f�t�H���g��0
  */
  int cf = Options["Contempt Factor"] * PawnValueEg / 100; // From centipawns
  DrawValue[ RootColor] = VALUE_DRAW - Value(cf);
  DrawValue[~RootColor] = VALUE_DRAW + Value(cf);
  /*
  ���[�g�ł̍��@��̎肪�Ȃ����
  UCI��info�R�}���h�ŒʒB����
  finalize:�ɔ��
  UCI�v���g�R���ɂ̓G���W�������畉����ʒm����R�}���h���Ȃ��悤�ł�
  http://www.geocities.jp/shogidokoro/usi.html
  */
  if (RootMoves.empty())
  {
      RootMoves.push_back(MOVE_NONE);
      sync_cout << "info depth 0 score "
                << score_to_uci(RootPos.checkers() ? -VALUE_MATE : VALUE_DRAW)
                << sync_endl;

      goto finalize;
  }

  /*
  Limits.infinite��go ponder�R�}���h�̂��Ƃ�infinite�I�v�V�����������
  stop���|������܂Ŗ������T���𑱂���B
  Limits.mate��go ponder�̃I�v�V�����A�����T�������� x move�Ŏ萔�𐧌��ɂ���
  Limits.mate�ɂ��̎萔�������Ă���

  ���Book���g�p����Ȃ�i�f�t�H���g��false�j���T���A
  ���̎肪RootMoves�z��ɂ���΂��̎��RootMoves�̐擪�ɍs����finalize���x����
  �ړ����邱�ƁB�܂�T��������Վ��D��̂���
  */
  if (Options["OwnBook"] && !Limits.infinite && !Limits.mate)
  {
      Move bookMove = book.probe(RootPos, Options["Book File"], Options["Best Book Move"]);

      if (bookMove && std::count(RootMoves.begin(), RootMoves.end(), bookMove))
      {
          std::swap(RootMoves[0], *std::find(RootMoves.begin(), RootMoves.end(), bookMove));
          goto finalize;
      }
  }

  /*
  Options["Write Search Log"]�̓f�t�H���g��false
  Limits.time White,Black���ꂼ��̎�������
  Limits.inc�@winc,binc(�P��m sec)�ڍוs��
  Limits.movestogo �T���ɐ�����݂�����̂̂悤�����ڍוs��
  */
  if (Options["Write Search Log"])
  {
      Log log(Options["Search Log Filename"]);
      log << "\nSearching: "  << RootPos.fen()
          << "\ninfinite: "   << Limits.infinite
          << " ponder: "      << Limits.ponder
          << " time: "        << Limits.time[RootColor]
          << " increment: "   << Limits.inc[RootColor]
          << " moves to go: " << Limits.movestogo
          << "\n" << std::endl;
  }

  // Reset the threads, still sleeping: will wake up at split time
  /*
  �p�r�s��
  thread�̓��C���X���b�h�ƒT���X���b�h(start_routine�X���b�h)��timer�X���b�h������͗l
  ���ƒT�����ɕ����̃X���b�h�ŒT���؂�T�������@����������Ă���Ǝv�����ڍוs��
  */
  for (size_t i = 0; i < Threads.size(); ++i)
      Threads[i]->maxPly = 0;

  Threads.timer->run = true;
  Threads.timer->notify_one(); // Wake up the recurring timer
  /*
  �T�����J�n
  */
  id_loop(RootPos); // Let's start searching !

  Threads.timer->run = false; // Stop the timer

  /*
  search�̃��O���L�^����I�v�V������true�ł���΃f�t�H���g�ł�false
  �t�@�C������SearchLog.txt�ɂȂ�
  */
  if (Options["Write Search Log"])
  {
      Time::point elapsed = Time::now() - SearchTime + 1;

      Log log(Options["Search Log Filename"]);
      log << "Nodes: "          << RootPos.nodes_searched()
          << "\nNodes/second: " << RootPos.nodes_searched() * 1000 / elapsed
          << "\nBest move: "    << move_to_san(RootPos, RootMoves[0].pv[0]);

      StateInfo st;
      RootPos.do_move(RootMoves[0].pv[0], st);
      log << "\nPonder move: " << move_to_san(RootPos, RootMoves[0].pv[1]) << std::endl;
      RootPos.undo_move(RootMoves[0].pv[0]);
  }

finalize:

  // When search is stopped this info is not printed
  /*
  UCI�v���g�R���Ńm�[�h���ƌo�ߎ��Ԃ�Ԃ��Ă���
  */
  sync_cout << "info nodes " << RootPos.nodes_searched()
            << " time " << Time::now() - SearchTime + 1 << sync_endl;

  // When we reach the maximum depth, we can arrive here without a raise of
  // Signals.stop. However, if we are pondering or in an infinite search,
  // the UCI protocol states that we shouldn't print the best move before the
  // GUI sends a "stop" or "ponderhit" command. We therefore simply wait here
  // until the GUI sends one of those commands (which also raises Signals.stop).
  /*
  �p�r�s��
  */
  if (!Signals.stop && (Limits.ponder || Limits.infinite))
  {
      Signals.stopOnPonderhit = true;
      RootPos.this_thread()->wait_for(Signals.stop);
  }

  // Best move could be MOVE_NONE when searching on a stalemate position
  /*
  UCI�ɒT�����ʂ�Ԃ��Ă���
  */
  sync_cout << "bestmove " << move_to_uci(RootMoves[0].pv[0], RootPos.is_chess960())
            << " ponder "  << move_to_uci(RootMoves[0].pv[1], RootPos.is_chess960())
            << sync_endl;
}


namespace {

  // id_loop() is the main iterative deepening loop. It calls search() repeatedly
  // with increasing depth until the allocated thinking time has been consumed,
  // user stops the search, or the maximum search depth is reached.
  /*
  think�֐�����Ăяo����Ă���
  ��������search�֐���NodeType(Root, PV, NonPV)��ݒ肵�ČĂяo��
  */
  void id_loop(Position& pos) {

    Stack stack[MAX_PLY_PLUS_6], *ss = stack+2; // To allow referencing (ss-2)
    int depth;
    Value bestValue, alpha, beta, delta;

    std::memset(ss-2, 0, 5 * sizeof(Stack));
    (ss-1)->currentMove = MOVE_NULL; // Hack to skip update gains

    depth = 0;
    BestMoveChanges = 0;
    bestValue = delta = alpha = -VALUE_INFINITE;
    beta = VALUE_INFINITE;

    TT.new_search();
    History.clear();
    Gains.clear();
    Countermoves.clear();
    Followupmoves.clear();

	/*
	�f�t�H���g�Ȃ�MultiPV�ɂP��Ԃ�
	�f�t�H���g�Ȃ�skill.level�ɂQ�O��Ԃ�
	*/
    MultiPV = Options["MultiPV"];
    Skill skill(Options["Skill Level"]);

    // Do we have to play with skill handicap? In this case enable MultiPV search
    // that we will use behind the scenes to retrieve a set of possible moves.
	/*
	�f�t�H���g�Ȃ�MultiPV��1�̂܂�
	*/
    if (skill.enabled() && MultiPV < 4)
        MultiPV = 4;

    MultiPV = std::min(MultiPV, RootMoves.size());

    // Iterative deepening loop until requested to stop or target depth reached
	/*
	�����[���@�iIterative deepening loop�j
	�w��̐[�x�iMAX_PLY�j�܂ŒB���邩�Astop���|��܂Ŕ����T�����s��
	Limits.depth��UCI�v���g�R������T���[�x���w�肵�Ă���΂�����ɏ]����������MAX_PLY���傫�Ȑ[�x��
	�Ӗ��Ȃ�
	depth=1����J�n�����MAX_PLY��120
	*/
    while (++depth <= MAX_PLY && !Signals.stop && (!Limits.depth || depth <= Limits.depth))
    {
        // Age out PV variability metric
		/*
		�ŏ���0.0�ɏ������i����id_loop�֐��̖`���Łj
		�p�r�s��
		*/
        BestMoveChanges *= 0.5;

        // Save the last iteration's scores before first PV line is searched and
        // all the move scores except the (new) PV are set to -VALUE_INFINITE.
		/*
		RootMoves�̓R���X�g���N�^�̂Ƃ�prevScore�ϐ�,score�ϐ��Ƃ�
		-VALUE_INFINITE(32001)�ɏ����ݒ肳��Ă���
		prevScore�͂��̂��ƂŕύX�����̂ł����ōď����������̂���
		*/
        for (size_t i = 0; i < RootMoves.size(); ++i)
            RootMoves[i].prevScore = RootMoves[i].score;

        // MultiPV loop. We perform a full root search for each PV line
		/*
		*/
        for (PVIdx = 0; PVIdx < MultiPV && !Signals.stop; ++PVIdx)
        {
            // Reset aspiration window starting size
			/*
			prevScore��-32001�Ȃ̂�
			alpha=max(-32001-16,-32001)= -32001
			beta=min(-32001+16,+32001)=  -31985
			window = 16
			*/
            if (depth >= 5)
            {
                delta = Value(16);
                alpha = std::max(RootMoves[PVIdx].prevScore - delta,-VALUE_INFINITE);
                beta  = std::min(RootMoves[PVIdx].prevScore + delta, VALUE_INFINITE);
            }

            // Start with a small aspiration window and, in the case of a fail
            // high/low, re-search with a bigger window until we're not failing
            // high/low anymore.
            while (true)
            {
                bestValue = search<Root, false>(pos, ss, alpha, beta, depth * ONE_PLY, false);

                // Bring the best move to the front. It is critical that sorting
                // is done with a stable algorithm because all the values but the
                // first and eventually the new best one are set to -VALUE_INFINITE
                // and we want to keep the same order for all the moves except the
                // new PV that goes to the front. Note that in case of MultiPV
                // search the already searched PV lines are preserved.
				/*
				RootMoves�z�������\�[�g���g���Ă���
				�\�[�g���ɂ��C���T�[�g�\�[�g�ƃ}�[�W�\�[�g���g�������Ă��邪
				��r�֐����w�肵�Ă��Ȃ��̂ŕW����less�֐��Ŕ�r���Ă��邪
				�W����less�֐�����RootMoves��bool operator<(const RootMove& m) const { return score > m.score; }
				���g����score���m���r���Ă���
				*/
                std::stable_sort(RootMoves.begin() + PVIdx, RootMoves.end());

                // Write PV back to transposition table in case the relevant
                // entries have been overwritten during the search.
				/*
				����ꂽ�����i�őP��菇�j��TT�ɓo�^���Ă���
				*/
                for (size_t i = 0; i <= PVIdx; ++i)
                    RootMoves[i].insert_pv_in_tt(pos);

                // If search has been stopped break immediately. Sorting and
                // writing PV back to TT is safe because RootMoves is still
                // valid, although it refers to previous iteration.
				/*
				stop��������΂��̉i�v���[�v����ł�
				*/
                if (Signals.stop)
                    break;

                // When failing high/low give some update (without cluttering
                // the UI) before a re-search.
				/*
				high/low���s�����Ƃ�
				uci_pv�֐��̓��e��W���o�͂ɏo��
				*/
                if (  (bestValue <= alpha || bestValue >= beta)
                    && Time::now() - SearchTime > 3000)
                    sync_cout << uci_pv(pos, depth, alpha, beta) << sync_endl;

                // In case of failing low/high increase aspiration window and
                // re-search, otherwise exit the loop.
				/*
				Low���s�����ꍇ�̍ĒT���̂��߂̕]���l�ݒ�
				*/
                if (bestValue <= alpha)
                {
					/*
					alpha�l��Ԃ��Ă����l���炳���delta(16)������A�܂�Window���L���čĒT������A������VALUE_INFINITE���͉����Ȃ�
					*/
                    alpha = std::max(bestValue - delta, -VALUE_INFINITE);
					/*
					failedLowAtRoot��Low���s�̃t���O
					check_time�֐����Ŏg�p����Ă���
					*/
                    Signals.failedLowAtRoot = true;
                    Signals.stopOnPonderhit = false;
                }
				/*
				High���s�����ꍇ�̍ĒT���̂��߂̕]���l�ݒ�
				alpha�̔��΂�beta�l��delta(16)�����グ��A�܂�Window���L���čĒT������A������VALUE_INFINITE���͏グ�Ȃ�
				*/
                else if (bestValue >= beta)
                    beta = std::min(bestValue + delta, VALUE_INFINITE);
				/*
				Low,High���s���Ȃ��̂Ő^�̕]���l��Window���ŕԂ��Ă����̂Ŏ��̔����[���Ɉڂ�
				*/
                else
                    break;
				/*
				Low,High���s�����ꍇ�͂P�U����
				16->24.0->36.0->54.0->81.0->121.5
				�Ə��X��delta���L���Ď��s���Ȃ��T�����s��
				*/
                delta += delta / 2;

                assert(alpha >= -VALUE_INFINITE && beta <= VALUE_INFINITE);
            }	//while(true)�I��

            // Sort the PV lines searched so far and update the GUI
			/*
			RootMove��RootMove�N���X��score�l�ň���\�[�g����
			*/
            std::stable_sort(RootMoves.begin(), RootMoves.begin() + PVIdx + 1);

            if (PVIdx + 1 == MultiPV || Time::now() - SearchTime > 3000)
                sync_cout << uci_pv(pos, depth, alpha, beta) << sync_endl;
        }	//MultiPV�I��

        // If skill levels are enabled and time is up, pick a sub-optimal best move
		/*
		���ƂŃR�����g������
		*/
        if (skill.enabled() && skill.time_to_pick(depth))
            skill.pick_move();

		/*
		search Log���ݒ肵�Ă����(�f�t�H���g�ł�false�jSearchLog.txt�Ƀ��O���c��
		*/
        if (Options["Write Search Log"])
        {
            RootMove& rm = RootMoves[0];
            if (skill.best != MOVE_NONE)
                rm = *std::find(RootMoves.begin(), RootMoves.end(), skill.best);

            Log log(Options["Search Log Filename"]);
            log << pretty_pv(pos, depth, rm.score, Time::now() - SearchTime, &rm.pv[0])
                << std::endl;
        }

        // Have we found a "mate in x"?
		/*
		Limits.mate�͉����T����𐧌�����I�v�V������
		�����ɂ����Ă������������������T�����~�ł��邪
		���̏����̈Ӗ����悭�킩���
		*/
        if (   Limits.mate
            && bestValue >= VALUE_MATE_IN_MAX_PLY
            && VALUE_MATE - bestValue <= 2 * Limits.mate)
            Signals.stop = true;

        // Do we have time for the next iteration? Can we stop searching now?
		/*
		�T����Limits�ɂ�鐧���A�T�����~����stop�t���O�Ȃǂ��|���Ă��Ȃ����
		*/
        if (Limits.use_time_management() && !Signals.stop && !Signals.stopOnPonderhit)
        {
            // Take some extra time if the best move has changed
			/*
			�p�r�s��
			*/
            if (depth > 4 && depth < 50 &&  MultiPV == 1)
                TimeMgr.pv_instability(BestMoveChanges);

            // Stop the search if only one legal move is available or all
            // of the available time has been used.
			/*
			���݂̏���Ԃ����Ԑ���̗L�����Ԃ��傫��������T�����I������iponder���|���Ă���ꍇ�͏����j
			*/
            if (   RootMoves.size() == 1
                || Time::now() - SearchTime > TimeMgr.available_time())
            {
                // If we are allowed to ponder do not stop the search now but
                // keep pondering until the GUI sends "ponderhit" or "stop".
                if (Limits.ponder)
                    Signals.stopOnPonderhit = true;
                else
                    Signals.stop = true;
            }
        }
    }	//�����[���I��
  }


  // search<>() is the main search function for both PV and non-PV nodes and for
  // normal and SplitPoint nodes. When called just after a split point the search
  // is simpler because we have already probed the hash table, done a null move
  // search, and searched the first move before splitting, so we don't have to
  // repeat all this work again. We also don't need to store anything to the hash
  // table here: This is taken care of after we return from the split point.
  /*
  NodeType�Ƃ́H�@
  SpNode�Ƃ́@split��Sp�����i�T������j
  id_loop����search�֐����ĂԎ���NodeType=Root,SpNode=false�ŌĂ΂��
  Depth depth�͔����[���̂��т�2->4->6�Ƒ����Ă���

  �T����search,qsearch�֐����ł�depth��ONE_PLY(=2)�Â����Ă���
  */
  template <NodeType NT, bool SpNode>
  Value search(Position& pos, Stack* ss, Value alpha, Value beta, Depth depth, bool cutNode) {

    const bool RootNode = NT == Root;
    const bool PvNode   = NT == PV || NT == Root;

    assert(-VALUE_INFINITE <= alpha && alpha < beta && beta <= VALUE_INFINITE);
    assert(PvNode || (alpha == beta - 1));
    assert(depth > DEPTH_ZERO);

    Move quietsSearched[64];
    StateInfo st;
    const TTEntry *tte;
    SplitPoint* splitPoint;
    Key posKey;
    Move ttMove, move, excludedMove, bestMove;
    Depth ext, newDepth, predictedDepth;
    Value bestValue, value, ttValue, eval, nullValue, futilityValue;
    bool inCheck, givesCheck, pvMove, singularExtensionNode, improving;
    bool captureOrPromotion, dangerous, doFullDepthSearch;
    int moveCount, quietCount;

    // Step 1. Initialize node
	/*
	�����̓m�[�h�̏������H
	*/
    Thread* thisThread = pos.this_thread();
    inCheck = pos.checkers();
    /*
	Root�̎��͂����͂Ƃ���Ȃ�
	�T�����򂷂�Ƃ��ɂƂ���Ǝv����
	*/
    if (SpNode)
    {
        splitPoint = ss->splitPoint;
        bestMove   = splitPoint->bestMove;
        bestValue  = splitPoint->bestValue;
        tte = NULL;
        ttMove = excludedMove = MOVE_NONE;
        ttValue = VALUE_NONE;

        assert(splitPoint->bestValue > -VALUE_INFINITE && splitPoint->moveCount > 0);

        goto moves_loop;
    }
	/*
	struct Stack {
		SplitPoint* splitPoint;
		int ply;
		Move currentMove;
		Move ttMove;
		Move excludedMove;
		Move killers[2];
		Depth reduction;
		Value staticEval;
		int skipNullMove;
	};
	ss�͍ŏ��̂Q�̓J�b�g����index�Q����search�֐��ɓn�����
	*/
    moveCount = quietCount = 0;
    bestValue = -VALUE_INFINITE;
    ss->currentMove = ss->ttMove = (ss+1)->excludedMove = bestMove = MOVE_NONE;
    ss->ply = (ss-1)->ply + 1;
    (ss+1)->skipNullMove = false; 
	(ss+1)->reduction = DEPTH_ZERO;
    (ss+2)->killers[0] = (ss+2)->killers[1] = MOVE_NONE;

    // Used to send selDepth info to GUI
	/*
	�p�r�s��
	*/
    if (PvNode && thisThread->maxPly < ss->ply)
        thisThread->maxPly = ss->ply;
	
	/*
	RootNode�ȊO�̂�
	���ƂŃR�����g����
	*/
    if (!RootNode)
    {
        // Step 2. Check for aborted search and immediate draw
        if (Signals.stop || pos.is_draw() || ss->ply > MAX_PLY)
            return ss->ply > MAX_PLY && !inCheck ? evaluate(pos) : DrawValue[pos.side_to_move()];

        // Step 3. Mate distance pruning. Even if we mate at the next move our score
        // would be at best mate_in(ss->ply+1), but if alpha is already bigger because
        // a shorter mate was found upward in the tree then there is no need to search
        // because we will never beat the current alpha. Same logic but with reversed
        // signs applies also in the opposite condition of being mated instead of giving
        // mate. In this case return a fail-high score.
        alpha = std::max(mated_in(ss->ply), alpha);
        beta = std::min(mate_in(ss->ply+1), beta);
        if (alpha >= beta)
            return alpha;
    }

    // Step 4. Transposition table lookup
    // We don't want the score of a partial search to overwrite a previous full search
    // TT value, so we use a different position key in case of an excluded move.
	/*
	ttMove�Ƃ̓g�����X�|�W�V�����e�[�u��������o�����w����
	�A�����[�g�m�[�h�ł�RootMoves[PVIdx].pv[0]������o�����w����
	RootMoves[PVIdx].pv[0]�̏ڍׂ͕s��
	*/
    excludedMove = ss->excludedMove;
    posKey = excludedMove ? pos.exclusion_key() : pos.key();
    tte = TT.probe(posKey);
    ss->ttMove = ttMove = RootNode ? RootMoves[PVIdx].pv[0] : tte ? tte->move() : MOVE_NONE;
    ttValue = tte ? value_from_tt(tte->value(), ss->ply) : VALUE_NONE;

    // At PV nodes we check for exact scores, whilst at non-PV nodes we check for
    // a fail high/low. The biggest advantage to probing at PV nodes is to have a
    // smooth experience in analysis mode. We don't probe at Root nodes otherwise
    // we should also update RootMoveList to avoid bogus output.
	/*
	�p�r�s��
	*/
    if (   !RootNode
        && tte
        && tte->depth() >= depth
        && ttValue != VALUE_NONE // Only in case of TT access race
        && (           PvNode ?  tte->bound() == BOUND_EXACT
            : ttValue >= beta ? (tte->bound() &  BOUND_LOWER)
                              : (tte->bound() &  BOUND_UPPER)))
    {
        ss->currentMove = ttMove; // Can be MOVE_NONE

        // If ttMove is quiet, update killers, history, counter move and followup move on TT hit
        if (ttValue >= beta && ttMove && !pos.capture_or_promotion(ttMove) && !inCheck)
            update_stats(pos, ss, ttMove, depth, NULL, 0);

        return ttValue;
    }

    // Step 5. Evaluate the position statically and update parent's gain statistics
	/*
	���肪�����Ă���Ȃ�moves_lopp���x���ɂƂ�ŒT�����n�߂�
	��������moves_lopp���x���܂ł͎}����̏����Ȃ̂ŉ��肪�������Ă���ꍇ�͖��Ӗ�
	*/
    if (inCheck)
    {
        ss->staticEval = eval = VALUE_NONE;
        goto moves_loop;
    }
	/*
	���肪�������Ă��Ȃ��Ē�Վ肪����Ȃ�
	*/
    else if (tte)
    {
        // Never assume anything on values stored in TT
        if ((ss->staticEval = eval = tte->eval_value()) == VALUE_NONE)
            eval = ss->staticEval = evaluate(pos);

        // Can ttValue be used as a better position evaluation?
        if (ttValue != VALUE_NONE)
            if (tte->bound() & (ttValue > eval ? BOUND_LOWER : BOUND_UPPER))
                eval = ttValue;
    }
    else
    {
        eval = ss->staticEval = evaluate(pos);
        TT.store(posKey, VALUE_NONE, BOUND_NONE, DEPTH_NONE, MOVE_NONE, ss->staticEval);
    }

    if (   !pos.captured_piece_type()
        &&  ss->staticEval != VALUE_NONE
        && (ss-1)->staticEval != VALUE_NONE
        && (move = (ss-1)->currentMove) != MOVE_NULL
        &&  type_of(move) == NORMAL)
    {
        Square to = to_sq(move);
        Gains.update(pos.piece_on(to), to, -(ss-1)->staticEval - ss->staticEval);
    }

    // Step 6. Razoring (skipped when in check)
	/*
	�Ȃɂ��̎}����H
	���������[�����ł�->qsearch�֐�����]���֐��ւ����̂ł�
	*/
    if (   !PvNode
        &&  depth < 4 * ONE_PLY
        &&  eval + razor_margin(depth) <= alpha
        &&  ttMove == MOVE_NONE
        &&  abs(beta) < VALUE_MATE_IN_MAX_PLY
        && !pos.pawn_on_7th(pos.side_to_move()))
    {
        if (   depth <= ONE_PLY
            && eval + razor_margin(3 * ONE_PLY) <= alpha)
            return qsearch<NonPV, false>(pos, ss, alpha, beta, DEPTH_ZERO);

        Value ralpha = alpha - razor_margin(depth);
        Value v = qsearch<NonPV, false>(pos, ss, ralpha, ralpha+1, DEPTH_ZERO);
        if (v <= ralpha)
            return v;
    }

    // Step 7. Futility pruning: child node (skipped when in check)
	/*
	Futility Pruning �́C�`�F�X�ōL���p�����Ă���
	�}�����@�ł���D�{�����[�ɂ����Ĕ��肳��郿��
	�@�̎}��������̔�������̐e�m�[�h (frontier node)
	�ŉ��̒l��p���čs�����Ƃɂ��C�s�v�ȐÎ~�T��
	�m�[�h�̓W�J�ƐÓI�]���֐��̌Ăяo�����팸����D

	���Ȃ킿�C�e�m�[�h P �ɂ�����]���l�ɁC�w���� m
	�ɑ΂��ĉςȃ}�[�W���l Vdiff(m) �������C�Ȃ���
	�̍ŏ��l�ɖ����Ȃ��ꍇ�� p �ȉ����}���肷�邱�Ƃ�
	�\�ƂȂ�D�œK�� Vdiff(m) �̒l�͕]���֐��ɂ��
	�ĈقȂ�C�������قǎ}���肪�L���ɓ����D
	http://www-als.ics.nitech.ac.jp/paper/H18-B/kanai.pdf

	eval�͌��ǖʂ̕]���l������futility_margin�֐����Ԃ��}�[�W�����������l
	���\�z�]���l�ł��̗\�z�]���l��beta�l�𒴂���̂ŃJ�b�g����

	���̎}����͌��[����2*7 = 14 ��菬�������Ɓi���[�ɋ߂����Ɓj�������̂ЂƂ�
	*/
    if (   !PvNode
        && !ss->skipNullMove
        &&  depth < 7 * ONE_PLY
        &&  eval - futility_margin(depth) >= beta
        &&  abs(beta) < VALUE_MATE_IN_MAX_PLY
        &&  abs(eval) < VALUE_KNOWN_WIN
        &&  pos.non_pawn_material(pos.side_to_move()))
        return eval - futility_margin(depth);

    // Step 8. Null move search with verification search (is omitted in PV nodes)
	/*
	�k�����[�u�i�}����j
	�k�����[�u�̏����A���[�x��2*ONE_PLY��葽�����ƁA�܂薖�[�ǖʈȊO�ł̓k�����[�uOK
	non_pawn_material��PAWN�ȊO�̋�]���l�̍��v��Ԃ�
	*/
    if (   !PvNode
        && !ss->skipNullMove
        &&  depth >= 2 * ONE_PLY
        &&  eval >= beta
        &&  abs(beta) < VALUE_MATE_IN_MAX_PLY
        &&  pos.non_pawn_material(pos.side_to_move()))
    {
        ss->currentMove = MOVE_NULL;

        assert(eval - beta >= 0);

        // Null move dynamic reduction based on depth and value
		/*
		Null�@Move�̒T���[�x�����߂�
		*/
        Depth R =  3 * ONE_PLY
                 + depth / 4
                 + int(eval - beta) / PawnValueMg * ONE_PLY;
		/*
		�p�X�̎�����s�i�ǖʂ͍X�V���Ȃ��j
		*/
        pos.do_null_move(st);
		/*
		�p�X�̎�̎��̓p�X���Ȃ�
		*/
        (ss+1)->skipNullMove = true;
		/*
		�k�����[�u�̐[�x�����݂̎c��̐[�x���[���Ȃ�qsearch�֐��i���[�T���j�����łȂ����search�֐��i��ʒT���j��
		�T������
		���ݎ��(null move)->������(�ʏ�move)->���ݎ��(null move)>>>
		*/
        nullValue = depth-R < ONE_PLY ? -qsearch<NonPV, false>(pos, ss+1, -beta, -beta+1, DEPTH_ZERO)
                                      : - search<NonPV, false>(pos, ss+1, -beta, -beta+1, depth-R, !cutNode);
		/*
		skipNullMove�����Ƃɖ߂�
		null move�����ɖ߂�
		*/
        (ss+1)->skipNullMove = false;
        pos.undo_null_move();
		
		/*
		����p�X���Ă����̕]���l��beta���傫���Ȃ�ʏ�Ɏ���w�����
		beta Cut���N�����Ɛ��������
		*/
        if (nullValue >= beta)
        {
            // Do not return unproven mate scores
            if (nullValue >= VALUE_MATE_IN_MAX_PLY)
                nullValue = beta;
			/*
			���[�x��12*ONE_PLAY��菬�����Ȃ牓���Ȃ��}�؂�i���[�ǖʂɋ߂��Ȃ�j
			*/
            if (depth < 12 * ONE_PLY)
                return nullValue;

            // Do verification search at high depths
			/*
			12*ONE_PLY>depth�iRoot�ǖʂɋ߂��ꍇ�j
			�ēx�T�����Ă��邪cutNode�̏������ǂ̂悤�ɒT���ɉe����^���Ă���̂��ڍוs��
			*/
            ss->skipNullMove = true;
            Value v = depth-R < ONE_PLY ? qsearch<NonPV, false>(pos, ss, beta-1, beta, DEPTH_ZERO)
                                        :  search<NonPV, false>(pos, ss, beta-1, beta, depth-R, false);
            ss->skipNullMove = false;
			/*
			������ς��ĒT�����Ă�beta�l�𒴂���悤�ł���Ή����Ȃ�Null Move Cut
			*/
            if (v >= beta)
                return nullValue;
        }
    }

    // Step 9. ProbCut (skipped when in check)
    // If we have a very good capture (i.e. SEE > seeValues[captured_piece_type])
    // and a reduced search returns a value much above beta, we can (almost) safely
    // prune the previous move.
	/*
	ProbCut �Ƃ̓I�Z���v���O����Logistello �̊J����M.Buro ���Ă̑O�����}����@�ŁA���W�ł�Multi-ProbCut �Ƃ�����@������B

	��{�A�C�f�A�́A�󂢒T���͐[���T���̋ߎ��ɂȂ�Ƃ������ƂŁA�[���T��������O�ɐ󂢒T�����s���A���ʂ��T��������
	�O�ꂽ��J�b�g���Ă��ǂ��񂶂�Ȃ��H�Ƃ�����@�B

	�c��[���������d �ɂȂ�����A�[��d' �̒T�����s��
	�]���l�Ƀ}�[�W��m �����A�T��������O��Ă��Ȃ����`�F�b�N
	�O�ꂽ=>�J�b�g
	�O��Ȃ�=>�[��d �̒T�����s��
	�Ȃ��A������O�ꂽ�O��Ȃ��̃`�F�b�N�́Anull-window search �����������B

	�ׂ������Ƃ��l�����ɁA���������ȃR�[�h������(�Ȃ��AD_probcut ��ProbCut 
	���s���Ƃ��炩���ߌ��߂Ă������[���ŁA��ł���d, D_shallow �͐󂢒T���̐[���ŁA��ł���d')�B
	http://d.hatena.ne.jp/tawake/20060710/1152520755 ���p�ӏ�
	�J�b�g�Ƃ����̂�
		if (value >= rbeta)
			return value;
	�̕������Ɣ��f�����
	*/
    if (   !PvNode
        &&  depth >= 5 * ONE_PLY
        && !ss->skipNullMove
        &&  abs(beta) < VALUE_MATE_IN_MAX_PLY)
    {
        Value rbeta = std::min(beta + 200, VALUE_INFINITE);
        Depth rdepth = depth - 4 * ONE_PLY;

        assert(rdepth >= ONE_PLY);
        assert((ss-1)->currentMove != MOVE_NONE);
        assert((ss-1)->currentMove != MOVE_NULL);
		/*
		pos.captured_piece_type()�͂Ƃ������ido_move�֐���st->capturedType�ɓo�^�����)
		*/
        MovePicker mp(pos, ttMove, History, pos.captured_piece_type());
        CheckInfo ci(pos);

        while ((move = mp.next_move<false>()) != MOVE_NONE)
            if (pos.legal(move, ci.pinned))
            {
                ss->currentMove = move;
                pos.do_move(move, st, ci, pos.gives_check(move, ci));
                value = -search<NonPV, false>(pos, ss+1, -rbeta, -rbeta+1, rdepth, !cutNode);
                pos.undo_move(move);
                if (value >= rbeta)
                    return value;
            }
    }

    // Step 10. Internal iterative deepening (skipped when in check)
	/*
	�p�r�s��
	*/
    if (    depth >= (PvNode ? 5 * ONE_PLY : 8 * ONE_PLY)
        && !ttMove
        && (PvNode || ss->staticEval + 256 >= beta))
    {
        Depth d = depth - 2 * ONE_PLY - (PvNode ? DEPTH_ZERO : depth / 4);

        ss->skipNullMove = true;
        search<PvNode ? PV : NonPV, false>(pos, ss, alpha, beta, d, true);
        ss->skipNullMove = false;

        tte = TT.probe(posKey);
        ttMove = tte ? tte->move() : MOVE_NONE;
    }

moves_loop: // When in check and at SpNode search starts from here
	/*
	countermoves�͂����ŏ���������Ă���

	*/
    Square prevMoveSq = to_sq((ss-1)->currentMove);
    Move countermoves[] = { Countermoves[pos.piece_on(prevMoveSq)][prevMoveSq].first,
                            Countermoves[pos.piece_on(prevMoveSq)][prevMoveSq].second };

    Square prevOwnMoveSq = to_sq((ss-2)->currentMove);
    Move followupmoves[] = { Followupmoves[pos.piece_on(prevOwnMoveSq)][prevOwnMoveSq].first,
                             Followupmoves[pos.piece_on(prevOwnMoveSq)][prevOwnMoveSq].second };
	/*
	���胊�X�g����
	*/
    MovePicker mp(pos, ttMove, depth, History, countermoves, followupmoves, ss);
    CheckInfo ci(pos);
    value = bestValue; // Workaround a bogus 'uninitialized' warning under gcc
    improving =   ss->staticEval >= (ss-2)->staticEval
               || ss->staticEval == VALUE_NONE
               ||(ss-2)->staticEval == VALUE_NONE;
	/*
	Singular Extension�Ƃ͂��̃m�[�h�̕]���l���Z��m�[�h�̕]���l���傫���ꍇ�n����
	���ʂȂǂ��^����̂ł��̃m�[�h�̒T�������������@
	�T���������L���ɂȂ������RootNode�ł͂Ȃ����Ɓ@�����@
	SpNode�ł͂Ȃ����ƁiSpNode�͒T������̂��ƁH�j�@����
	depth��8*ONE_PLY���傫�����Ɓi�܂�RootNode�ɋ߂��A���[�߂��ł͂Ȃ����Ɓ@����
	!excludedMove���ĂȂɁ@����
	�g�����X�|�W�V�����e�[�u���̕]���l�������l�@����
	�g�����X�|�W�V�����e�[�u���̎w����̐[�x�����݂̐[�x���3*ONE_PLAY���������̂�肨����������

	�Z��̃m�[�h�̕]���l�Ƃ��S�R�o�Ă��Ȃ��͉̂���
	*/
    singularExtensionNode =   !RootNode
                           && !SpNode
                           &&  depth >= 8 * ONE_PLY
                           &&  ttMove != MOVE_NONE
                           && !excludedMove // Recursive(�ċA�I) singular search is not allowed(�����ꂽ)->�ċA�I�ȃV���M���[�g���͋�����Ȃ��H
                           && (tte->bound() & BOUND_LOWER)
                           &&  tte->depth() >= depth - 3 * ONE_PLY;

    // Step 11. Loop through moves
    // Loop through all pseudo-legal moves until no moves remain or a beta cutoff occurs
	/*
	�������炪���C���̒T��
	*/
    while ((move = mp.next_move<SpNode>()) != MOVE_NONE)
    {
      assert(is_ok(move));

	  /*
	  �p�r�s��
	  */
      if (move == excludedMove)
          continue;

      // At root obey the "searchmoves" option and skip moves not listed in Root
      // Move List. As a consequence any illegal move is also skipped. In MultiPV
      // mode we also skip PV moves which have been already searched.
	  /*
	  RootNode���[�h��next_move�֐���RootMoves�ɂȂ���������Ă����炻��̓p�X����
	  RootMoves�ɂ���肵���ǂ܂Ȃ�
	  */
      if (RootNode && !std::count(RootMoves.begin() + PVIdx, RootMoves.end(), move))
          continue;

	  /*
	  RootNode�̂Ƃ���SpNode��false
	  SpNode�͒T������̂��ƁH
	  */
      if (SpNode)
      {
          // Shared counter cannot be decremented later if the move turns out to be illegal
          if (!pos.legal(move, ci.pinned))
              continue;

          moveCount = ++splitPoint->moveCount;
          splitPoint->mutex.unlock();
      }
      else
          ++moveCount;

	  /*
	  RootNode���[�h��p
	  RootNode�ő�1��ڂ̂Ƃ�Signals.firstRootMove��true�ɂ���
	  �p�r�s��
	  */
      if (RootNode)
      {
          Signals.firstRootMove = (moveCount == 1);

          if (thisThread == Threads.main() && Time::now() - SearchTime > 3000)
              sync_cout << "info depth " << depth / ONE_PLY
                        << " currmove " << move_to_uci(move, pos.is_chess960())
                        << " currmovenumber " << moveCount + PVIdx << sync_endl;
      }

      ext = DEPTH_ZERO;
      captureOrPromotion = pos.capture_or_promotion(move);

	  /*
	  �w����p�^�[�����m�[�}���œGKING�ւ̗������ז����Ă����Ȃ��ꍇ
		�ړ���ɓGKING�ɗ��������������bitboard��givesCheck�ɗ^����givesCheck��bool�^�Ȃ̂�bitboard�������true�ɂȂ�A�������Ȃ��Ƃ������Ƃ�����
	  �����łȂ���� gives_check�֐����Ăщ��肪�Ȃ������d�Ƀ`�G�b�N����
	  */
      givesCheck =  type_of(move) == NORMAL && !ci.dcCandidates
                  ? ci.checkSq[type_of(pos.piece_on(from_sq(move)))] & to_sq(move)
                  : pos.gives_check(move, ci);

	  /*
	  dangerous���̂́u�댯�ȁv�ƌ����Ӗ�
	  ���肪�\�Ȏ�@OR�@�w����p�^�[����NORMAL�ȊO(PROMOTION,ENPASSANT,CASTLING)
	  PAWN��RANK4�ȏ�̈ʒu�ɂ���.(RANK4��WHITE������݂��ʒu�ABLACK����݂��RANK5�ȉ��j
	  */
      dangerous =   givesCheck
                 || type_of(move) != NORMAL
                 || pos.advanced_pawn_push(move);

      // Step 12. Extend checks
	  /*
	  ���肪�\�ȏ�Ԃɂ���(givesCheck=true),�Î~�T�����Ă��]���l��0�ȏ�ł���
	  ���Ƃ�������ONE_PLY�����T����������
	  */
      if (givesCheck && pos.see_sign(move) >= VALUE_ZERO)
          ext = ONE_PLY;

      // Singular extension search. If all moves but one fail low on a search of
      // (alpha-s, beta-s), and just one fails high on (alpha, beta), then that move
      // is singular and should be extended. To verify this we do a reduced search
      // on all the other moves but the ttMove and if the result is lower than
      // ttValue minus a margin then we extend the ttMove.
	  /*
	  Singular extension search�Ƃ�
		���̃m�[�h�̕]���l���Z��m�[�h�̕]���l���傫���ꍇ�n����
		���ʂȂǂ��^����̂ł��̃m�[�h�̒T�������������@

		������if���ŕʒT�������[�x�̔����ōs���Ă���
		���̕����ł̕]���l�ɂ����ext�ϐ��i�[�x�����ϐ��j
		��ONE_PLY�����������Ă��镔���ł͂Ȃ���

		�ڍוs��
	  */
      if (    singularExtensionNode
          &&  move == ttMove
          && !ext
          &&  pos.legal(move, ci.pinned)
          &&  abs(ttValue) < VALUE_KNOWN_WIN)
      {
          assert(ttValue != VALUE_NONE);

          Value rBeta = ttValue - int(depth);
          ss->excludedMove = move;
          ss->skipNullMove = true;
          value = search<NonPV, false>(pos, ss, rBeta - 1, rBeta, depth / 2, cutNode);
          ss->skipNullMove = false;
          ss->excludedMove = MOVE_NONE;

          if (value < rBeta)
              ext = ONE_PLY;
      }

      // Update the current move (this must be done after singular extension search)
	  /*
	  �ŏI������Singular extension�̒T�������������Ă���Ƃ���
	  �����I��ONE_PLY�����Ă���͉̂���
	  */
      newDepth = depth - ONE_PLY + ext;

      // Step 13. Pruning at shallow depth (exclude PV nodes)
	  /*
	  �}����
	  P��Node�ł͂Ȃ��iPvNode���ĂȂɁj�@����
	  ���肪�������Ă��Ȃ��@�����i���肪�������Ă���悤�ȃm�[�h���}���肵�Ă͊댯�j
	  ���肪�\�Ȏ蓙�d�v�Ȏ�ł͂Ȃ��@����
	  bestValue > VALUE_MATED_IN_MAX_PLY(= �}�C�i�X31880)�ɒ[�ɕ]���l�������킯�ł͂Ȃ����悭���Ȃ��H
	  */
      if (   !PvNode
          && !captureOrPromotion
          && !inCheck
          && !dangerous
       /* &&  move != ttMove Already implicit in the next condition */
          &&  bestValue > VALUE_MATED_IN_MAX_PLY)
      {
          // Move count based pruning
		  /*
		  ���ݐ[�x�����[�ɋ߂��萔��FutilityMoveCounts[improving][depth]�i���ꂪ�����͕s���j
		  ��葽���@�܂茋�\�A�萔���[�x���[���ǂ񂾏o�̃p�X�Ƃ����}����H
		  */
          if (   depth < 16 * ONE_PLY
              && moveCount >= FutilityMoveCounts[improving][depth] )
          {
              if (SpNode)
                  splitPoint->mutex.lock();

              continue;
          }

		  /*
		  �p�r�s��
		  */
          predictedDepth = newDepth - reduction<PvNode>(improving, depth, moveCount);	//predictedDepth�@=�@�\�z�̐[�x�H

          // Futility pruning: parent node
		  /*
		  �p�r�s��
		  */
          if (predictedDepth < 7 * ONE_PLY)
          {
              futilityValue = ss->staticEval + futility_margin(predictedDepth)
                            + 128 + Gains[pos.moved_piece(move)][to_sq(move)];

              if (futilityValue <= alpha)
              {
                  bestValue = std::max(bestValue, futilityValue);

                  if (SpNode)
                  {
                      splitPoint->mutex.lock();
                      if (bestValue > splitPoint->bestValue)
                          splitPoint->bestValue = bestValue;
                  }
                  continue;
              }
          }

          // Prune moves with negative SEE at low depths
		  /*
		  �p�r�s��
		  */
          if (predictedDepth < 4 * ONE_PLY && pos.see_sign(move) < VALUE_ZERO)
          {
              if (SpNode)
                  splitPoint->mutex.lock();

              continue;
          }
      }

      // Check for legality just before making the move
	  /*
	  ���@��ł��邩�̃`�G�b�N�A���@��łȂ���΂��̃m�[�h�̓p�X
	  ���̂����Ń`�G�b�N�Ȃ̂������Ƒ����ł��Ȃ��̂���
	  */
      if (!RootNode && !SpNode && !pos.legal(move, ci.pinned))
      {
          moveCount--;
          continue;
      }

      pvMove = PvNode && moveCount == 1;
      ss->currentMove = move;
      if (!SpNode && !captureOrPromotion && quietCount < 64)
          quietsSearched[quietCount++] = move;

      // Step 14. Make the move
	  /*
	  �����ŋǖʍX�V
	  */
      pos.do_move(move, st, ci, givesCheck);

      // Step 15. Reduced depth search (LMR). If the move fails high it will be
      // re-searched at full depth.
	  /*
	  �}�؂�Ȃ̂��T�������Ȃ̂�����s���A�k���H
	  https://chessprogramming.wikispaces.com/Late+Move+Reductions
	  �ŏI�I�ɂ�doFullDepthSearch��ݒ肷��Ƃ��낪���̕����̌��_�H
	  */
      if (    depth >= 3 * ONE_PLY
          && !pvMove
          && !captureOrPromotion
          &&  move != ttMove
          &&  move != ss->killers[0]
          &&  move != ss->killers[1])
      {
          ss->reduction = reduction<PvNode>(improving, depth, moveCount);

          if (!PvNode && cutNode)
              ss->reduction += ONE_PLY;

          else if (History[pos.piece_on(to_sq(move))][to_sq(move)] < 0)
              ss->reduction += ONE_PLY / 2;

          if (move == countermoves[0] || move == countermoves[1])
              ss->reduction = std::max(DEPTH_ZERO, ss->reduction - ONE_PLY);

          Depth d = std::max(newDepth - ss->reduction, ONE_PLY);
          if (SpNode)
              alpha = splitPoint->alpha;

          value = -search<NonPV, false>(pos, ss+1, -(alpha+1), -alpha, d, true);

          // Re-search at intermediate depth if reduction is very high
		  /*
		  */
          if (value > alpha && ss->reduction >= 4 * ONE_PLY)
          {
              Depth d2 = std::max(newDepth - 2 * ONE_PLY, ONE_PLY);
              value = -search<NonPV, false>(pos, ss+1, -(alpha+1), -alpha, d2, true);
          }

          doFullDepthSearch = (value > alpha && ss->reduction != DEPTH_ZERO);
          ss->reduction = DEPTH_ZERO;
      }
      else
          doFullDepthSearch = !pvMove;

      // Step 16. Full depth search, when LMR is skipped or fails high
	  /*
	  ���������̊K�w�ɍ~��čs���Ƃ���
	  qsearch�֐���search�֐�����I�����Ă��邪�������G�ŏڍוs��
	  SpNode�i�T������j�AdoFullDepthSearch�͂Ȃ�
	  doFullDepthSearch��seach�֐��̎����ϐ���bool�^��step15�Őݒ肳��Ă���
	  �ڍוs��

	  */
      if (doFullDepthSearch)
      {
          if (SpNode)
              alpha = splitPoint->alpha;
		  /*
		  newDepth < ONE_PLY�����������
		  givesCheck ? -qsearch<NonPV,  true>(pos, ss+1, -(alpha+1), -alpha, DEPTH_ZERO) : -qsearch<NonPV, false>(pos, ss+1, -(alpha+1), -alpha, DEPTH_ZERO)
		  �����s����
		  newDepth < ONE_PLY���������Ȃ����
		  -search<NonPV, false>(pos, ss+1, -(alpha+1), -alpha, newDepth, !cutNode)
		  �����s����B
		  �܂�V�����ݒ肳�ꂽ�T���[����ONE_PLY(2)��菬�����ꍇ��qsearch�֐���
		  ONE_PLY���傫���ꍇ��search�֐����Ă�

		  */
          value = newDepth < ONE_PLY ?
                          givesCheck ? -qsearch<NonPV,  true>(pos, ss+1, -(alpha+1), -alpha, DEPTH_ZERO)
                                     : -qsearch<NonPV, false>(pos, ss+1, -(alpha+1), -alpha, DEPTH_ZERO)
                                     : - search<NonPV, false>(pos, ss+1, -(alpha+1), -alpha, newDepth, !cutNode);
      }

      // For PV nodes only, do a full PV search on the first move or after a fail
      // high (in the latter case search only if value < beta), otherwise let the
      // parent node fail low with value <= alpha and to try another move.
	  /*
	  newDepth < ONE_PLY����������΁i�܂肠�ƂP��Ŗ��[�ǖʂȂ�qsearch�֐��j
	  givesCheck ? -qsearch<PV,  true>(pos, ss+1, -beta, -alpha, DEPTH_ZERO) : -qsearch<PV, false>(pos, ss+1, -beta, -alpha, DEPTH_ZERO)
      �����s��
	  newDepth < ONE_PLY���������Ȃ���΁i�܂�܂����[�ǖʂł͂Ȃ��Ȃ�j
	  -search<PV, false>(pos, ss+1, -beta, -alpha, newDepth, false)
	  �����s����
	  �܂�V�����ݒ肳�ꂽ�T���[����ONE_PLY(2)��菬�����ꍇ��qsearch�֐���
	  ONE_PLY���傫���ꍇ��search�֐����Ă�
	  ��Ɠ����悤�ȍ\���ł��邪�A�Ⴄ�̂�NonPV��PV���̈Ⴂ���Ǝv��
	  */
      if (PvNode && (pvMove || (value > alpha && (RootNode || value < beta))))
          value = newDepth < ONE_PLY ?
                          givesCheck ? -qsearch<PV,  true>(pos, ss+1, -beta, -alpha, DEPTH_ZERO)
                                     : -qsearch<PV, false>(pos, ss+1, -beta, -alpha, DEPTH_ZERO)
                                     : - search<PV, false>(pos, ss+1, -beta, -alpha, newDepth, false);
      // Step 17. Undo move
	  /*
	  �����ŋǖʕ���
	  */
      pos.undo_move(move);

      assert(value > -VALUE_INFINITE && value < VALUE_INFINITE);

      // Step 18. Check for new best move
	  /*
	  �T�����򂵂Ă���Thread�Ȃ�
	  �܂��ڍוs��
	  */
      if (SpNode)
      {
          splitPoint->mutex.lock();
          bestValue = splitPoint->bestValue;
          alpha = splitPoint->alpha;
      }

      // Finished searching the move. If a stop or a cutoff occurred, the return
      // value of the search cannot be trusted, and we return immediately without
      // updating best move, PV and TT.
	  /*
	  �p�r�s��
	  */
      if (Signals.stop || thisThread->cutoff_occurred())
          return VALUE_ZERO;

      if (RootNode)
      {
          RootMove& rm = *std::find(RootMoves.begin(), RootMoves.end(), move);

          // PV move or new best move ?
          if (pvMove || value > alpha)
          {
              rm.score = value;
              rm.extract_pv_from_tt(pos);

              // We record how often the best move has been changed in each
              // iteration. This information is used for time management: When
              // the best move changes frequently, we allocate some more time.
              if (!pvMove)
                  ++BestMoveChanges;
          }
          else
              // All other moves but the PV are set to the lowest value: this is
              // not a problem when sorting because the sort is stable and the
              // move position in the list is preserved - just the PV is pushed up.
              rm.score = -VALUE_INFINITE;
      }

      if (value > bestValue)
      {
          bestValue = SpNode ? splitPoint->bestValue = value : value;

          if (value > alpha)
          {
              bestMove = SpNode ? splitPoint->bestMove = move : move;

              if (PvNode && value < beta) // Update alpha! Always alpha < beta
                  alpha = SpNode ? splitPoint->alpha = value : value;
              else
              {
                  assert(value >= beta); // Fail high

                  if (SpNode)
                      splitPoint->cutoff = true;

                  break;
              }
          }
      }

      // Step 19. Check for splitting the search
	  /*
	  SpNode�o�Ȃ����̕Ԃ���H
	  */
      if (   !SpNode
          &&  Threads.size() >= 2
          &&  depth >= Threads.minimumSplitDepth
          &&  (   !thisThread->activeSplitPoint
               || !thisThread->activeSplitPoint->allSlavesSearching)
          &&  thisThread->splitPointsSize < MAX_SPLITPOINTS_PER_THREAD)
      {
          assert(bestValue > -VALUE_INFINITE && bestValue < beta);

          thisThread->split<FakeSplit>(pos, ss, alpha, beta, &bestValue, &bestMove,
                                       depth, moveCount, &mp, NT, cutNode);

          if (Signals.stop || thisThread->cutoff_occurred())
              return VALUE_ZERO;

          if (bestValue >= beta)
              break;
      }
    }	//���������C���T���̏I��

	/*
	SpNode�̎��̕Ԃ���H
	*/
    if (SpNode)
        return bestValue;

    // Following condition would detect a stop or a cutoff set only after move
    // loop has been completed. But in this case bestValue is valid because we
    // have fully searched our subtree, and we can anyhow save the result in TT.
    /*
       if (Signals.stop || thisThread->cutoff_occurred())
        return VALUE_DRAW;
    */

    // Step 20. Check for mate and stalemate
    // All legal moves have been searched and if there are no legal moves, it
    // must be mate or stalemate. If we are in a singular extension search then
    // return a fail low score.
	/*
	�p�r�s��
	*/
    if (!moveCount)
        bestValue = excludedMove ? alpha
                   :     inCheck ? mated_in(ss->ply) : DrawValue[pos.side_to_move()];

    // Quiet best move: update killers, history, countermoves and followupmoves
	/*
	beta Cut�����������Ƃ���update_stats�֐�������
	*/
    else if (bestValue >= beta && !pos.capture_or_promotion(bestMove) && !inCheck)
        update_stats(pos, ss, bestMove, depth, quietsSearched, quietCount - 1);

	/*
	�g�����X�|�W�V�����e�[�u���ɓo�^
	*/
    TT.store(posKey, value_to_tt(bestValue, ss->ply),
             bestValue >= beta  ? BOUND_LOWER :
             PvNode && bestMove ? BOUND_EXACT : BOUND_UPPER,
             depth, bestMove, ss->staticEval);

    assert(bestValue > -VALUE_INFINITE && bestValue < VALUE_INFINITE);

	/*
	�Ԃ�
	*/
    return bestValue;
  }	//������search�֐��̏I��


  // qsearch() is the quiescence search function, which is called by the main
  // search function when the remaining depth is zero (or, to be more precise,
  // less than ONE_PLY).
  /*
  ���[�ǖʐ�p�T���֐��̂͂�
  search�֐��ɔ�ׂ�Ƃ����ԍs�������Ȃ�
  */
  template <NodeType NT, bool InCheck>
  Value qsearch(Position& pos, Stack* ss, Value alpha, Value beta, Depth depth) {

    const bool PvNode = NT == PV;

    assert(NT == PV || NT == NonPV);
    assert(InCheck == !!pos.checkers());
    assert(alpha >= -VALUE_INFINITE && alpha < beta && beta <= VALUE_INFINITE);
    assert(PvNode || (alpha == beta - 1));
    assert(depth <= DEPTH_ZERO);

    StateInfo st;
    const TTEntry* tte;
    Key posKey;
    Move ttMove, move, bestMove;
    Value bestValue, value, ttValue, futilityValue, futilityBase, oldAlpha;
    bool givesCheck, evasionPrunable;
    Depth ttDepth;

    // To flag BOUND_EXACT a node with eval above alpha and no available moves
	/*
	�T������Ȃ�alpha��oldalpha�ɑޔ��H
	*/
    if (PvNode)
        oldAlpha = alpha;

    ss->currentMove = bestMove = MOVE_NONE;
    ss->ply = (ss-1)->ply + 1;

    // Check for an instant draw or if the maximum ply has been reached
	/*
	�����������肵�Ĉ��������@OR �ő�[�x�܂ł����Ȃ�@�]���l��������ċA��
	*/
    if (pos.is_draw() || ss->ply > MAX_PLY)
        return ss->ply > MAX_PLY && !InCheck ? evaluate(pos) : DrawValue[pos.side_to_move()];

    // Decide whether or not to include checks: this fixes also the type of
    // TT entry depth that we are going to use. Note that in qsearch we use
    // only two types of depth in TT: DEPTH_QS_CHECKS or DEPTH_QS_NO_CHECKS.
	/*
	�p�r�s��
	*/
    ttDepth = InCheck || depth >= DEPTH_QS_CHECKS ? DEPTH_QS_CHECKS
                                                  : DEPTH_QS_NO_CHECKS;

    // Transposition table lookup
	/*
	�g�����X�|�W�V�����e�[�u���Ɏ肪���邩���ׂ�A�����ttMove�Ɏ�����Ă���
	*/
    posKey = pos.key();
    tte = TT.probe(posKey);
    ttMove = tte ? tte->move() : MOVE_NONE;
    ttValue = tte ? value_from_tt(tte->value(),ss->ply) : VALUE_NONE;

	/*
	�g�����X�|�W�V�����e�[�u���̎�̕]���l���^�l�Ȃ�g�����X�|�W�V�����e�[�u���̕]���l��Ԃ�
	*/
    if (   tte
        && tte->depth() >= ttDepth
        && ttValue != VALUE_NONE // Only in case of TT access race
        && (           PvNode ?  tte->bound() == BOUND_EXACT
            : ttValue >= beta ? (tte->bound() &  BOUND_LOWER)
                              : (tte->bound() &  BOUND_UPPER)))
    {
        ss->currentMove = ttMove; // Can be MOVE_NONE
        return ttValue;
    }

    // Evaluate the position statically
	/*
	���̏�Ԃŉ��肪�������Ă���ꍇ�̓}�C�i�X32001��Ԃ�
	*/
    if (InCheck)
    {
        ss->staticEval = VALUE_NONE;
        bestValue = futilityBase = -VALUE_INFINITE;
    }
    else
	/*
	���肪�������Ă��Ȃ��ꍇ
	*/
    {
        if (tte)
        {
            // Never assume anything on values stored in TT
			/*
			�p�r�s��
			*/
            if ((ss->staticEval = bestValue = tte->eval_value()) == VALUE_NONE)
                ss->staticEval = bestValue = evaluate(pos);

            // Can ttValue be used as a better position evaluation?
            if (ttValue != VALUE_NONE)
                if (tte->bound() & (ttValue > bestValue ? BOUND_LOWER : BOUND_UPPER))
                    bestValue = ttValue;
        }
		/*
		�p�r�s��
		*/
        else
            ss->staticEval = bestValue = evaluate(pos);

        // Stand pat. Return immediately if static value is at least beta
		/*
		�p�r�s��
		*/
        if (bestValue >= beta)
        {
            if (!tte)
                TT.store(pos.key(), value_to_tt(bestValue, ss->ply), BOUND_LOWER,
                         DEPTH_NONE, MOVE_NONE, ss->staticEval);

            return bestValue;
        }

        if (PvNode && bestValue > alpha)
            alpha = bestValue;

        futilityBase = bestValue + 128;
    }

    // Initialize a MovePicker object for the current position, and prepare
    // to search the moves. Because the depth is <= 0 here, only captures,
    // queen promotions and checks (only if depth >= DEPTH_QS_CHECKS) will
    // be generated.
	/*
	���̐�͂ȂɁH
	*/
    MovePicker mp(pos, ttMove, depth, History, to_sq((ss-1)->currentMove));
    CheckInfo ci(pos);

    // Loop through the moves until no moves remain or a beta cutoff occurs
    while ((move = mp.next_move<false>()) != MOVE_NONE)
    {
      assert(is_ok(move));

      givesCheck =  type_of(move) == NORMAL && !ci.dcCandidates
                  ? ci.checkSq[type_of(pos.piece_on(from_sq(move)))] & to_sq(move)
                  : pos.gives_check(move, ci);

      // Futility pruning
      if (   !PvNode
          && !InCheck
          && !givesCheck
          &&  move != ttMove
          &&  futilityBase > -VALUE_KNOWN_WIN
          && !pos.advanced_pawn_push(move))
      {
          assert(type_of(move) != ENPASSANT); // Due to !pos.advanced_pawn_push

          futilityValue = futilityBase + PieceValue[EG][pos.piece_on(to_sq(move))];

          if (futilityValue < beta)
          {
              bestValue = std::max(bestValue, futilityValue);
              continue;
          }

          if (futilityBase < beta && pos.see(move) <= VALUE_ZERO)
          {
              bestValue = std::max(bestValue, futilityBase);
              continue;
          }
      }

      // Detect non-capture evasions that are candidates to be pruned
      evasionPrunable =    InCheck
                       &&  bestValue > VALUE_MATED_IN_MAX_PLY
                       && !pos.capture(move)
                       && !pos.can_castle(pos.side_to_move());

      // Don't search moves with negative SEE values
      if (   !PvNode
          && (!InCheck || evasionPrunable)
          &&  move != ttMove
          &&  type_of(move) != PROMOTION
          &&  pos.see_sign(move) < VALUE_ZERO)
          continue;

      // Check for legality just before making the move
      if (!pos.legal(move, ci.pinned))
          continue;

      ss->currentMove = move;

      // Make and search the move
      pos.do_move(move, st, ci, givesCheck);
      value = givesCheck ? -qsearch<NT,  true>(pos, ss+1, -beta, -alpha, depth - ONE_PLY)
                         : -qsearch<NT, false>(pos, ss+1, -beta, -alpha, depth - ONE_PLY);
      pos.undo_move(move);

      assert(value > -VALUE_INFINITE && value < VALUE_INFINITE);

      // Check for new best move
      if (value > bestValue)
      {
          bestValue = value;

          if (value > alpha)
          {
              if (PvNode && value < beta) // Update alpha here! Always alpha < beta
              {
                  alpha = value;
                  bestMove = move;
              }
              else // Fail high
              {
                  TT.store(posKey, value_to_tt(value, ss->ply), BOUND_LOWER,
                           ttDepth, move, ss->staticEval);

                  return value;
              }
          }
       }
    }

    // All legal moves have been searched. A special case: If we're in check
    // and no legal moves were found, it is checkmate.
    if (InCheck && bestValue == -VALUE_INFINITE)
        return mated_in(ss->ply); // Plies to mate from the root

    TT.store(posKey, value_to_tt(bestValue, ss->ply),
             PvNode && bestValue > oldAlpha ? BOUND_EXACT : BOUND_UPPER,
             ttDepth, bestMove, ss->staticEval);

    assert(bestValue > -VALUE_INFINITE && bestValue < VALUE_INFINITE);

    return bestValue;
  }


  // value_to_tt() adjusts a mate score from "plies to mate from the root" to
  // "plies to mate from the current position". Non-mate scores are unchanged.
  // The function is called before storing a value in the transposition table.
  /*
  �p�r�s��
  */
  Value value_to_tt(Value v, int ply) {

    assert(v != VALUE_NONE);

    return  v >= VALUE_MATE_IN_MAX_PLY  ? v + ply
          : v <= VALUE_MATED_IN_MAX_PLY ? v - ply : v;
  }


  // value_from_tt() is the inverse of value_to_tt(): It adjusts a mate score
  // from the transposition table (which refers to the plies to mate/be mated
  // from current position) to "plies to mate/be mated from the root".
  /*
  �p�r�s��
  */
  Value value_from_tt(Value v, int ply) {

    return  v == VALUE_NONE             ? VALUE_NONE
          : v >= VALUE_MATE_IN_MAX_PLY  ? v - ply
          : v <= VALUE_MATED_IN_MAX_PLY ? v + ply : v;
  }


  // update_stats() updates killers, history, countermoves and followupmoves stats after a fail-high
  // of a quiet move.
  /*
  killers,history,countermoves,followupmoves���X�V����֐��̂悤�ł�
  beta Cut�����������Ƃ��Ă΂��
  */
  void update_stats(const Position& pos, Stack* ss, Move move, Depth depth, Move* quiets, int quietsCnt) {

    if (ss->killers[0] != move)
    {
        ss->killers[1] = ss->killers[0];
        ss->killers[0] = move;
    }

    // Increase history value of the cut-off move and decrease all the other
    // played quiet moves.
    Value bonus = Value(int(depth) * int(depth));
    History.update(pos.moved_piece(move), to_sq(move), bonus);
    for (int i = 0; i < quietsCnt; ++i)
    {
        Move m = quiets[i];
        History.update(pos.moved_piece(m), to_sq(m), -bonus);
    }

    if (is_ok((ss-1)->currentMove))
    {
        Square prevMoveSq = to_sq((ss-1)->currentMove);
        Countermoves.update(pos.piece_on(prevMoveSq), prevMoveSq, move);
    }

    if (is_ok((ss-2)->currentMove) && (ss-1)->currentMove == (ss-1)->ttMove)
    {
        Square prevOwnMoveSq = to_sq((ss-2)->currentMove);
        Followupmoves.update(pos.piece_on(prevOwnMoveSq), prevOwnMoveSq, move);
    }
  }


  // When playing with a strength handicap, choose best move among the MultiPV
  // set using a statistical rule dependent on 'level'. Idea by Heinz van Saanen.
  /*
  �p�r�s��	
  */
  Move Skill::pick_move() {

    static RKISS rk;

    // PRNG sequence should be not deterministic
    for (int i = Time::now() % 50; i > 0; --i)
        rk.rand<unsigned>();

    // RootMoves are already sorted by score in descending order
    int variance = std::min(RootMoves[0].score - RootMoves[MultiPV - 1].score, PawnValueMg);
    int weakness = 120 - 2 * level;
    int max_s = -VALUE_INFINITE;
    best = MOVE_NONE;

    // Choose best move. For each move score we add two terms both dependent on
    // weakness. One deterministic and bigger for weaker moves, and one random,
    // then we choose the move with the resulting highest score.
    for (size_t i = 0; i < MultiPV; ++i)
    {
        int s = RootMoves[i].score;

        // Don't allow crazy blunders even at very low skills
        if (i > 0 && RootMoves[i-1].score > s + 2 * PawnValueMg)
            break;

        // This is our magic formula
        s += (  weakness * int(RootMoves[0].score - s)
              + variance * (rk.rand<unsigned>() % weakness)) / 128;

        if (s > max_s)
        {
            max_s = s;
            best = RootMoves[i].pv[0];
        }
    }
    return best;
  }


  // uci_pv() formats PV information according to the UCI protocol. UCI
  // requires that all (if any) unsearched PV lines are sent using a previous
  // search score.
  /*
  �ڍוs���ł��邪,UCI�����Ɍ��i�K�̋ǖʏ��iscore,nodes,nps�Ȃǁj���o�͂���
  id_loop�֐�����̂݌Ăяo�����
  */
  string uci_pv(const Position& pos, int depth, Value alpha, Value beta) {

    std::stringstream ss;
    Time::point elapsed = Time::now() - SearchTime + 1;
    size_t uciPVSize = std::min((size_t)Options["MultiPV"], RootMoves.size());
    int selDepth = 0;

    for (size_t i = 0; i < Threads.size(); ++i)
        if (Threads[i]->maxPly > selDepth)
            selDepth = Threads[i]->maxPly;

    for (size_t i = 0; i < uciPVSize; ++i)
    {
        bool updated = (i <= PVIdx);

        if (depth == 1 && !updated)
            continue;

        int d   = updated ? depth : depth - 1;
        Value v = updated ? RootMoves[i].score : RootMoves[i].prevScore;

        if (ss.rdbuf()->in_avail()) // Not at first line
            ss << "\n";

        ss << "info depth " << d
           << " seldepth "  << selDepth
           << " score "     << (i == PVIdx ? score_to_uci(v, alpha, beta) : score_to_uci(v))
           << " nodes "     << pos.nodes_searched()
           << " nps "       << pos.nodes_searched() * 1000 / elapsed
           << " time "      << elapsed
           << " multipv "   << i + 1
           << " pv";

        for (size_t j = 0; RootMoves[i].pv[j] != MOVE_NONE; ++j)
            ss << " " << move_to_uci(RootMoves[i].pv[j], pos.is_chess960());
    }

    return ss.str();
  }

} // namespace


/// RootMove::extract_pv_from_tt() builds a PV by adding moves from the TT table.
/// We also consider both failing high nodes and BOUND_EXACT nodes here to
/// ensure that we have a ponder move even when we fail high at root. This
/// results in a long PV to print that is important for position analysis.
/*
�p�r�s��
*/
void RootMove::extract_pv_from_tt(Position& pos) {

  StateInfo state[MAX_PLY_PLUS_6], *st = state;
  const TTEntry* tte;
  int ply = 1;    // At root ply is 1...
  Move m = pv[0]; // ...instead pv[] array starts from 0
  Value expectedScore = score;

  pv.clear();

  do {
      pv.push_back(m);

      assert(MoveList<LEGAL>(pos).contains(pv[ply - 1]));

      pos.do_move(pv[ply++ - 1], *st++);
      tte = TT.probe(pos.key());
      expectedScore = -expectedScore;

  } while (   tte
           && expectedScore == value_from_tt(tte->value(), ply)
           && pos.pseudo_legal(m = tte->move()) // Local copy, TT could change
           && pos.legal(m, pos.pinned_pieces(pos.side_to_move()))
           && ply < MAX_PLY
           && (!pos.is_draw() || ply <= 2));

  pv.push_back(MOVE_NONE); // Must be zero-terminating

  while (--ply) pos.undo_move(pv[ply - 1]);
}


/// RootMove::insert_pv_in_tt() is called at the end of a search iteration, and
/// inserts the PV back into the TT. This makes sure the old PV moves are searched
/// first, even if the old TT entries have been overwritten.
/*
�����ɂ���^�C�~���O�͔����[����search�֐�����Ԃ��Ă����Ƃ���Ȃ̂�
pv�ɂ͂P�Â�RootMove�ɑ΂��čőP��菇���o�^����Ă���
�����ł͂����ɓ����Ă�����TT�i�g�����X�|�W�V�����e�[�u���j�ɓo�^���Ă���
���肾���ł͂Ȃ�PV�S��TT�ɓo�^���Ă���
*/
void RootMove::insert_pv_in_tt(Position& pos) {

  StateInfo state[MAX_PLY_PLUS_6], *st = state;
  const TTEntry* tte;
  int idx = 0; // Ply starts from 1, we need to start from 0
  /*
  TT�̓g�����X�|�W�V�����e�[�u����\���O���[�o���ϐ�
  pv:�i�őP��菇�j������ׂ̂P�����̉ϒ��z��łP��RootMove�N���X�Ɋi�[����Ă���
  */
  do {
      tte = TT.probe(pos.key());
	  /*
	  �n���ꂽ�ǖʂŁATT�ɓo�^���ꂽ�肪�Ȃ��܂��͂����Ă������̎�Ƃ͈قȂ�Ƃ�
	  �������TT�ɓo�^����
	  */
      if (!tte || tte->move() != pv[idx]) // Don't overwrite correct entries
          TT.store(pos.key(), VALUE_NONE, BOUND_NONE, DEPTH_NONE, pv[idx], VALUE_NONE);

      assert(MoveList<LEGAL>(pos).contains(pv[idx]));

      pos.do_move(pv[idx++], *st++);

  } while (pv[idx] != MOVE_NONE);

  while (idx) pos.undo_move(pv[--idx]);
}


/// Thread::idle_loop() is where the thread is parked when it has no work to do
/*
Thread::split�֐�����Ă΂��
option��o["Threads"]�𕡐��ɂ��Ȃ��ƌĂ΂�Ȃ��A�f�t�H���g��1

*/
void Thread::idle_loop() {

  // Pointer 'this_sp' is not null only if we are called from split(), and not
  // at the thread creation. This means we are the split point's master.
  SplitPoint* this_sp = splitPointsSize ? activeSplitPoint : NULL;

  assert(!this_sp || (this_sp->masterThread == this && searching));

  while (true)
  {
      // If we are not searching, wait for a condition to be signaled instead of
      // wasting CPU time polling for work.
      while (!searching || exit)
      {
          if (exit)
          {
              assert(!this_sp);
              return;
          }

          // Grab the lock to avoid races with Thread::notify_one()
          mutex.lock();

          // If we are master and all slaves have finished then exit idle_loop
          if (this_sp && this_sp->slavesMask.none())
          {
              mutex.unlock();
              break;
          }

          // Do sleep after retesting sleep conditions under lock protection. In
          // particular we need to avoid a deadlock in case a master thread has,
          // in the meanwhile, allocated us and sent the notify_one() call before
          // we had the chance to grab the lock.
          if (!searching && !exit)
              sleepCondition.wait(mutex);

          mutex.unlock();
      }

      // If this thread has been assigned work, launch a search
      if (searching)
      {
          assert(!exit);

          Threads.mutex.lock();

          assert(searching);
          assert(activeSplitPoint);
          SplitPoint* sp = activeSplitPoint;

          Threads.mutex.unlock();

          Stack stack[MAX_PLY_PLUS_6], *ss = stack+2; // To allow referencing (ss-2)
          Position pos(*sp->pos, this);

          std::memcpy(ss-2, sp->ss-2, 5 * sizeof(Stack));
          ss->splitPoint = sp;

          sp->mutex.lock();

          assert(activePosition == NULL);

          activePosition = &pos;

          if (sp->nodeType == NonPV)
              search<NonPV, true>(pos, ss, sp->alpha, sp->beta, sp->depth, sp->cutNode);

          else if (sp->nodeType == PV)
              search<PV, true>(pos, ss, sp->alpha, sp->beta, sp->depth, sp->cutNode);

          else if (sp->nodeType == Root)
              search<Root, true>(pos, ss, sp->alpha, sp->beta, sp->depth, sp->cutNode);

          else
              assert(false);

          assert(searching);

          searching = false;
          activePosition = NULL;
          sp->slavesMask.reset(idx);
          sp->allSlavesSearching = false;
          sp->nodes += pos.nodes_searched();

          // Wake up the master thread so to allow it to return from the idle
          // loop in case we are the last slave of the split point.
          if (    this != sp->masterThread
              &&  sp->slavesMask.none())
          {
              assert(!sp->masterThread->searching);
              sp->masterThread->notify_one();
          }

          // After releasing the lock we can't access any SplitPoint related data
          // in a safe way because it could have been released under our feet by
          // the sp master.
          sp->mutex.unlock();

          // Try to late join to another split point if none of its slaves has
          // already finished.
          if (Threads.size() > 2)
              for (size_t i = 0; i < Threads.size(); ++i)
              {
                  int size = Threads[i]->splitPointsSize; // Local copy
                  sp = size ? &Threads[i]->splitPoints[size - 1] : NULL;

                  if (   sp
                      && sp->allSlavesSearching
                      && available_to(Threads[i]))
                  {
                      // Recheck the conditions under lock protection
                      Threads.mutex.lock();
                      sp->mutex.lock();

                      if (   sp->allSlavesSearching
                          && available_to(Threads[i]))
                      {
                           sp->slavesMask.set(idx);
                           activeSplitPoint = sp;
                           searching = true;
                      }

                      sp->mutex.unlock();
                      Threads.mutex.unlock();

                      break; // Just a single attempt
                  }
              }
      }

      // If this thread is the master of a split point and all slaves have finished
      // their work at this split point, return from the idle loop.
      if (this_sp && this_sp->slavesMask.none())
      {
          this_sp->mutex.lock();
          bool finished = this_sp->slavesMask.none(); // Retest under lock protection
          this_sp->mutex.unlock();
          if (finished)
              return;
      }
  }
}


/// check_time() is called by the timer thread when the timer triggers. It is
/// used to print debug info and, more importantly, to detect when we are out of
/// available time and thus stop the search.

void check_time() {

  static Time::point lastInfoTime = Time::now();
  int64_t nodes = 0; // Workaround silly 'uninitialized' gcc warning

  if (Time::now() - lastInfoTime >= 1000)
  {
      lastInfoTime = Time::now();
      dbg_print();
  }

  if (Limits.ponder)
      return;

  if (Limits.nodes)
  {
      Threads.mutex.lock();

      nodes = RootPos.nodes_searched();

      // Loop across all split points and sum accumulated SplitPoint nodes plus
      // all the currently active positions nodes.
      for (size_t i = 0; i < Threads.size(); ++i)
          for (int j = 0; j < Threads[i]->splitPointsSize; ++j)
          {
              SplitPoint& sp = Threads[i]->splitPoints[j];

              sp.mutex.lock();

              nodes += sp.nodes;

              for (size_t idx = 0; idx < Threads.size(); ++idx)
                  if (sp.slavesMask.test(idx) && Threads[idx]->activePosition)
                      nodes += Threads[idx]->activePosition->nodes_searched();

              sp.mutex.unlock();
          }

      Threads.mutex.unlock();
  }

  Time::point elapsed = Time::now() - SearchTime;
  bool stillAtFirstMove =    Signals.firstRootMove
                         && !Signals.failedLowAtRoot
                         &&  elapsed > TimeMgr.available_time() * 75 / 100;

  bool noMoreTime =   elapsed > TimeMgr.maximum_time() - 2 * TimerThread::Resolution
                   || stillAtFirstMove;

  if (   (Limits.use_time_management() && noMoreTime)
      || (Limits.movetime && elapsed >= Limits.movetime)
      || (Limits.nodes && nodes >= Limits.nodes))
      Signals.stop = true;
}
