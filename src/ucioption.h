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

#ifndef UCIOPTION_H_INCLUDED
#define UCIOPTION_H_INCLUDED

#include <map>
#include <string>

namespace UCI {

class Option;

/*
用途不明、OptionMapに使われている
*/
/// Custom comparator because UCI options should be case insensitive
struct CaseInsensitiveLess {
  bool operator() (const std::string&, const std::string&) const;
};

/// Our options container is actually a std::map
/*
オプションを保持しているマップコンテナ
*/
typedef std::map<std::string, Option, CaseInsensitiveLess> OptionsMap;

/// Option class implements an option as defined by UCI protocol
class Option {
  /*
  この宣言の意味がよくわからん
  */
  typedef void (*OnChange)(const Option&);

public:
  /*


  オプションの取り方が４種類ある
  関数（省略可）buttn型　ー＞実装はないようである
  bool型＋関数（省略可）　check型
  char型＋関数（省略可）　string型
  int型,int型,int型,関数（省略可）　spin型

  button,check,string型,spin型はこのクラスのプライベート変数に
  文字列として保存されている
  */
  Option(OnChange = NULL);
  Option(bool v, OnChange = NULL);
  Option(const char* v, OnChange = NULL);
  Option(int v, int min, int max, OnChange = NULL);

  Option& operator=(const std::string& v);
  void operator<<(const Option& o);
  operator int() const;
  operator std::string() const;

private:
  friend std::ostream& operator<<(std::ostream&, const OptionsMap&);

  std::string defaultValue, currentValue, type;
  int min, max;
  size_t idx;
  OnChange on_change;
};
/*
この関数はOptionではなくUser interfaceでの関数
*/
void init(OptionsMap&);
void loop(int argc, char* argv[]);

} // namespace UCI

extern UCI::OptionsMap Options;

#endif // #ifndef UCIOPTION_H_INCLUDED
