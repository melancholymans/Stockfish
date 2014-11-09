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

#ifndef MISC_H_INCLUDED
#define MISC_H_INCLUDED

#include <fstream>
#include <string>
#include <vector>

#include "types.h"

/*
stackfish�̊ȒP�ȏЉ���\�z����
version�ԍ��A�@�g�pOS��bit���ABMI2���g���Ă��邩�ASIMD���߂̎g�p
�Ȃǂ�\������B
*/
extern const std::string engine_info(bool to_uci = false);
/*
�w�肳�ꂽ�~���b�����҂�
thread.h��wait_for�֐�����Ăяo�����
*/
extern void timed_wait(WaitCondition&, Lock&, int);
/*
�v���t�F�b�`�i��ǂ݁j
���O�Ƀv���Z�b�T�ɋ߂��L���b�V���K�w�Ƀf�[�^�����[�h���Ă��������ꍇ�Ɏg�p������@�ł��D
�Ƃ���ǂ���ɐ�ǂ݂����������������w�肵�Ă���悤�����A���Ȃ�ׂ������͂��K�v�ȋC������
*/
extern void prefetch(char* addr);
/*
�p�r�s��
*/
extern void start_logger(bool b);
/*
�p�r�s���A�g�p���Ă��鍭�ՂȂ�
*/
extern void dbg_hit_on(bool b);
extern void dbg_hit_on_c(bool c, bool b);
extern void dbg_mean_of(int v);
extern void dbg_print();


struct Log : public std::ofstream {
  Log(const std::string& f = "log.txt") : std::ofstream(f.c_str(), std::ios::out | std::ios::app) {}
 ~Log() { if (is_open()) close(); }
};


namespace Time {
  typedef int64_t point;
  inline point now() { return system_time_to_msec(); }
}


template<class Entry, int Size>
struct HashTable {
  HashTable() : table(Size, Entry()) {}
  Entry* operator[](Key k) { return &table[(uint32_t)k & (Size - 1)]; }

private:
  std::vector<Entry> table;
};


enum SyncCout { IO_LOCK, IO_UNLOCK };
std::ostream& operator<<(std::ostream&, SyncCout);

#define sync_cout std::cout << IO_LOCK
#define sync_endl std::endl << IO_UNLOCK

#endif // #ifndef MISC_H_INCLUDED
