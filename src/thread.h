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
std::Thread�ŏ��������ƕ�����₷�����ȁiC++11�j
http://yohhoy.hatenablog.jp/entry/2013/12/15/204116
*/

#ifndef THREAD_H_INCLUDED
#define THREAD_H_INCLUDED

#include <bitset>
#include <vector>

#include "material.h"
#include "movepick.h"
#include "pawns.h"
#include "position.h"
#include "search.h"

const int MAX_THREADS = 128;
const int MAX_SPLITPOINTS_PER_THREAD = 8;

/*
Mutex��Windows�ł̓N���e�C�J���Z�N�V�����i�r������p�j
Windows�̃X���b�h�֌W���߂�Mutex�����邪stackfish�ł͎g�p���Ă��Ȃ�
Windows��Mutex�̓v���Z�X�Ԃŋ��L���铯���Ɏg�p����
�����ł�Mutex�̓N���e�C�J���Z�N�V������UNIX����Mutex�ƌĂ�ł��邾��
*/
struct Mutex {
  Mutex() { lock_init(l); }
 ~Mutex() { lock_destroy(l); }

  void lock() { lock_grab(l); }
  void unlock() { lock_release(l); }

private:
  friend struct ConditionVariable;

  Lock l;
};

/*
ConditionVariable�͏����ϐ����������}���`�X���b�h�̂P��@�̖��O
Windows�ł�Event�ɑ�������

UNIX�ɂ͏����t���ϐ��Ƃ��������Z�p�����邪Windows�ɂ͂Ȃ��悤��
C++11�Ɏ������ꂽstd::Thread�N���X�ɂ͂��邪�����ł͎g�p����Ă��Ȃ�

cond_init�֐���Event�ϐ��iWaitCondition c�ϐ��j����������
�����l��FALSE�Őݒ肳���B

wait�֐��͎��Ԗ�������Event�ϐ����V�O�i����Ԃɂ���
�܂肱��wait�֐����Ăяo�����X���b�h��Event�ϐ������L��
���̃X���b�h���T�X�y���h��Ԃɂ���
notify_one�֐����Ăяo����Event�ϐ��̓V�O�i����ԂƂȂ�
���̃X���b�h���V�O�i����Ԃɂ���

*/
struct ConditionVariable {
  ConditionVariable() { cond_init(c); }
 ~ConditionVariable() { cond_destroy(c); }

  void wait(Mutex& m) { cond_wait(c, m.l); }
  void wait_for(Mutex& m, int ms) { timed_wait(c, m.l, ms); }
  void notify_one() { cond_signal(c); }

private:
  WaitCondition c;
};

struct Thread;

struct SplitPoint {

  // Const data after split point has been setup
  const Position* pos;
  const Search::Stack* ss;
  Thread* masterThread;
  Depth depth;
  Value beta;
  int nodeType;
  bool cutNode;

  // Const pointers to shared data
  MovePicker* movePicker;
  SplitPoint* parentSplitPoint;

  // Shared data
  Mutex mutex;
  std::bitset<MAX_THREADS> slavesMask;
  volatile bool allSlavesSearching;
  volatile uint64_t nodes;
  volatile Value alpha;
  volatile Value bestValue;
  volatile Move bestMove;
  volatile int moveCount;
  volatile bool cutoff;
};


/// ThreadBase struct is the base of the hierarchy from where we derive all the
/// specialized thread classes.

struct ThreadBase {

  ThreadBase() : handle(NativeHandle()), exit(false) {}
  virtual ~ThreadBase() {}
  virtual void idle_loop() = 0;
  void notify_one();
  void wait_for(volatile const bool& b);

  Mutex mutex;
  ConditionVariable sleepCondition;
  NativeHandle handle;
  volatile bool exit;
};


/// Thread struct keeps together all the thread related stuff like locks, state
/// and especially split points. We also use per-thread pawn and material hash
/// tables so that once we get a pointer to an entry its life time is unlimited
/// and we don't have to care about someone changing the entry under our feet.

struct Thread : public ThreadBase {

  Thread();
  virtual void idle_loop();
  bool cutoff_occurred() const;
  bool available_to(const Thread* master) const;

  template <bool Fake>
  void split(Position& pos, const Search::Stack* ss, Value alpha, Value beta, Value* bestValue, Move* bestMove,
             Depth depth, int moveCount, MovePicker* movePicker, int nodeType, bool cutNode);

  SplitPoint splitPoints[MAX_SPLITPOINTS_PER_THREAD];
  Material::Table materialTable;
  Endgames endgames;
  Pawns::Table pawnsTable;
  Position* activePosition;
  size_t idx;
  int maxPly;
  SplitPoint* volatile activeSplitPoint;
  volatile int splitPointsSize;
  volatile bool searching;
};


/// MainThread and TimerThread are derived classes used to characterize the two
/// special threads: the main one and the recurring timer.

struct MainThread : public Thread {
  MainThread() : thinking(true) {} // Avoid a race with start_thinking()
  virtual void idle_loop();
  volatile bool thinking;
};

struct TimerThread : public ThreadBase {
  TimerThread() : run(false) {}
  virtual void idle_loop();
  bool run;
  static const int Resolution = 5; // msec between two check_time() calls
};


/// ThreadPool struct handles all the threads related stuff like init, starting,
/// parking and, most importantly, launching a slave thread at a split point.
/// All the access to shared thread data is done through this class.

struct ThreadPool : public std::vector<Thread*> {

  void init(); // No c'tor and d'tor, threads rely on globals that should
  void exit(); // be initialized and are valid during the whole thread lifetime.

  MainThread* main() { return static_cast<MainThread*>((*this)[0]); }
  void read_uci_options();
  Thread* available_slave(const Thread* master) const;
  void wait_for_think_finished();
  void start_thinking(const Position&, const Search::LimitsType&, Search::StateStackPtr&);

  Depth minimumSplitDepth;
  Mutex mutex;
  ConditionVariable sleepCondition;
  TimerThread* timer;
};

extern ThreadPool Threads;

#endif // #ifndef THREAD_H_INCLUDED
