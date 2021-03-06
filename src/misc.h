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
stackfishの簡単な紹介文を構築する
version番号、使用OSのbit数、BMI2を使っているか、SIMD命令の使用
などを表示する。
*/
extern const std::string engine_info(bool to_uci = false);
/*
指定されたミリ秒だけ待つ
thread.hのwait_for関数から呼び出される
*/
extern void timed_wait(WaitCondition&, Lock&, int);
/*
プリフェッチ（先読み）
事前にプロセッサに近いキャッシュ階層にデータをロードしておきたい場合に使用する方法です．
ところどころに先読みさせたいメモリを指定してあるようだが、かなり細かい分析が必要な気がする
*/
extern void prefetch(char* addr);

/*
用途不明
*/
extern void start_logger(bool b);

/*
用途不明、使用している痕跡なし
*/
extern void dbg_hit_on(bool b);
extern void dbg_hit_on_c(bool c, bool b);
extern void dbg_mean_of(int v);
extern void dbg_print();

/*
用途不明
*/

struct Log : public std::ofstream {
  Log(const std::string& f = "log.txt") : std::ofstream(f.c_str(), std::ios::out | std::ios::app) {}
 ~Log() { if (is_open()) close(); }
};

/*
用途不明
*/
namespace Time {
  typedef int64_t point;
  inline point now() { return system_time_to_msec(); }
}

/*
用途不明
ハッシュテーブルだと思うが詳細不明
*/
template<class Entry, int Size>
struct HashTable {
  HashTable() : table(Size, Entry()) {}
  Entry* operator[](Key k) { return &table[(uint32_t)k & (Size - 1)]; }

private:
  std::vector<Entry> table;
};

/*
用途不明
*/
enum SyncCout { IO_LOCK, IO_UNLOCK };
std::ostream& operator<<(std::ostream&, SyncCout);

#define sync_cout std::cout << IO_LOCK
#define sync_endl std::endl << IO_UNLOCK

#endif // #ifndef MISC_H_INCLUDED
