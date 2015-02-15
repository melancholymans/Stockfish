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

#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "evaluate.h"
#include "notation.h"
#include "position.h"
#include "search.h"
#include "thread.h"
#include "tt.h"
#include "ucioption.h"

using namespace std;

extern void benchmark(const Position& pos, istream& is);

namespace {

  // FEN string of the initial position, normal chess
  const char* StartFEN = "rnbqkbnr/pppppppp/8/8/8/8/PPPPPPPP/RNBQKBNR w KQkq - 0 1";

  // Keep a track of the position keys along the setup moves (from the start position
  // to the position just before the search starts). This is needed by the repetition
  // draw detection code.
  Search::StateStackPtr SetupStates;


  // position() is called when engine receives the "position" UCI command.
  // The function sets up the position described in the given FEN string ("fen")
  // or the starting position ("startpos") and then makes the moves given in the
  // following move list ("moves").
  /*
  position startposまたは fen 局面を構成するfen文字列を入力することで
  局面を再設定できる
  */
  void position(Position& pos, istringstream& is) {

    Move m;
    string token, fen;

    is >> token;

    if (token == "startpos")
    {
        fen = StartFEN;
        is >> token; // Consume "moves" token if any
    }
    else if (token == "fen")
        while (is >> token && token != "moves")
            fen += token + " ";
    else
        return;

    pos.set(fen, Options["UCI_Chess960"], Threads.main());
    SetupStates = Search::StateStackPtr(new std::stack<StateInfo>());

    // Parse move list (if any)
    while (is >> token && (m = move_from_uci(pos, token)) != MOVE_NONE)
    {
        SetupStates->push(StateInfo());
        pos.do_move(m, SetupStates->top());
    }
  }


  // setoption() is called when engine receives the "setoption" UCI command. The
  // function updates the UCI option ("name") to the given value ("value").
  /*
  オプションを設定する
  */
  void setoption(istringstream& is) {

    string token, name, value;

    is >> token; // Consume "name" token

    // Read option name (can contain spaces)
	/*
	nameとvalueの間にある文字列をspaceを挟みながら連結してname変数に入れておく
	*/
    while (is >> token && token != "value")
        name += string(" ", !name.empty()) + token;

    // Read option value (can contain spaces)
    while (is >> token)
        value += string(" ", !value.empty()) + token;

    if (Options.count(name))
        Options[name] = value;
    else
        sync_cout << "No such option: " << name << sync_endl;
  }


  // go() is called when engine receives the "go" UCI command. The function sets
  // the thinking time and other parameters from the input string, and starts
  // the search.
  /*
  User　Interfaceからこのコマンドがきたらオプションを設定の上
  Threads.start_thinking関数を呼んで探索開始
  */
  void go(const Position& pos, istringstream& is) {

    Search::LimitsType limits;
    string token;

    while (is >> token)
    {
		/*
		go のあとにsearchmoves を続けて特定の指し手をa2a3,c2c4などと指定するとその手のみ探索する
		*/
        if (token == "searchmoves")
            while (is >> token)
                limits.searchmoves.push_back(move_from_uci(pos, token));

        else if (token == "wtime")     is >> limits.time[WHITE];
        else if (token == "btime")     is >> limits.time[BLACK];
        else if (token == "winc")      is >> limits.inc[WHITE];
        else if (token == "binc")      is >> limits.inc[BLACK];
        else if (token == "movestogo") is >> limits.movestogo;
        else if (token == "depth")     is >> limits.depth;
        else if (token == "nodes")     is >> limits.nodes;
        else if (token == "movetime")  is >> limits.movetime;
        else if (token == "mate")      is >> limits.mate;
        else if (token == "infinite")  limits.infinite = true;
        else if (token == "ponder")    limits.ponder = true;
    }

    Threads.start_thinking(pos, limits, SetupStates);
  }

} // namespace


/// Wait for a command from the user, parse this text string as an UCI command,
/// and call the appropriate functions. Also intercepts EOF from stdin to ensure
/// that we exit gracefully if the GUI dies unexpectedly. In addition to the UCI
/// commands, the function also supports a few debug commands.
/*
main関数から呼ばれる,"quit"と入力するまでここで永久ループ
*/
void UCI::loop(int argc, char* argv[]) {
  /*
  pos(StartFEN, false, Threads.main())
  局面をStartFEN文字で初期設定する
  falseはchess960ではなく正規のchessであることを指示している
  Threads.main()は詳細不明
  */
  Position pos(StartFEN, false, Threads.main()); // The root position
  string token, cmd;
  /*
  stackfish自身の引数をここで文字列にしている、セパレータはspace
  */
  for (int i = 1; i < argc; ++i)
      cmd += std::string(argv[i]) + " ";

  do {
	  /*
	  ここでユーザーの入力を待機している
	  */
      if (argc == 1 && !getline(cin, cmd)) // Block here waiting for input
          cmd = "quit";

      istringstream is(cmd);

      is >> skipws >> token;

      if (token == "quit" || token == "stop" || token == "ponderhit")
      {
          // The GUI sends 'ponderhit' to tell us to ponder on the same move the
          // opponent has played. In case Signals.stopOnPonderhit is set we are
          // waiting for 'ponderhit' to stop the search (for instance because we
          // already ran out of time), otherwise we should continue searching but
          // switch from pondering to normal search.
		  /*
		  tokenがquitかstopなら何らかのスレッド操作をしているようで
		  詳細不明
		  */
          if (token != "ponderhit" || Search::Signals.stopOnPonderhit)
          {
              Search::Signals.stop = true;
              Threads.main()->notify_one(); // Could be sleeping
          }
          else
              Search::Limits.ponder = false;
      }
	  /*
	  用途不明
	  */
      else if (token == "perft" || token == "divide")
      {
          int depth;
          stringstream ss;

          is >> depth;
          ss << Options["Hash"]    << " "
             << Options["Threads"] << " " << depth << " current " << token;

          benchmark(pos, ss);
      }
	  /*
	  position key,material key,pawn keyを表示するだけ
	  */
      else if (token == "key")
          sync_cout << hex << uppercase << setfill('0')
                    << "position key: "   << setw(16) << pos.key()
                    << "\nmaterial key: " << setw(16) << pos.material_key()
                    << "\npawn key:     " << setw(16) << pos.pawn_key()
                    << dec << nouppercase << setfill(' ') << sync_endl;
	  /*
	  chessエンジンの簡単な紹介(engine_info)を表示したあと
	  << Optionsは現在のオプションを表示する <<に対して演算子オーバーロードしてある
	  */
      else if (token == "uci")
          sync_cout << "id name " << engine_info(true)
                    << "\n"       << Options
                    << "\nuciok"  << sync_endl;
	  /*
	  Eval term |    White    |    Black    |    Total
	  |   MG    EG  |   MG    EG  |   MG    EG
	  ---------------------+-------------+-------------+-------------
	  Material, PST, Tempo |   ---   --- |   ---   --- |  0.00  0.00
	  Material imbalance |   ---   --- |   ---   --- |  0.00  0.00
	  Pawns |   ---   --- |   ---   --- |  0.00  0.00
	  Knights |  0.12  0.00 |  0.12  0.00 |  0.00  0.00
	  Bishops | -0.12 -0.37 | -0.12 -0.37 |  0.00  0.00
	  Rooks | -0.32  0.00 | -0.32  0.00 |  0.00  0.00
	  Queens |  0.00  0.00 |  0.00  0.00 |  0.00  0.00
	  Mobility | -0.93 -0.98 | -0.93 -0.98 |  0.00  0.00
	  King safety |  1.02 -0.06 |  1.02 -0.06 |  0.00  0.00
	  Threats |  0.00  0.00 |  0.00  0.00 |  0.00  0.00
	  Passed pawns |  0.00  0.00 |  0.00  0.00 |  0.00  0.00
	  Space |  0.36  0.00 |  0.36  0.00 |  0.00  0.00
	  ---------------------+-------------+-------------+-------------
	  Total |   ---   --- |   ---   --- |  0.09  0.04

	  Total Evaluation: 0.09 (white side)	
	  このような表示を返してくる
	  */
      else if (token == "eval")
      {
          Search::RootColor = pos.side_to_move(); // Ensure it is set
          sync_cout << Eval::trace(pos) << sync_endl;
      }
	  /*
	  このコマンドはなにもしない
	  */
      else if (token == "ucinewgame") { /* Avoid returning "Unknown command" */ }
	  /*
	  探索を開始する
	  */
      else if (token == "go")         go(pos, is);
	  /*
	  任意または標準局面にする
	  */
      else if (token == "position")   position(pos, is);
	  /*
	  setoption name <option name> value <x>のフォーマットで入力して
	  オプションを設定する
	  */
      else if (token == "setoption")  setoption(is);
	  /*
	  先手を後手に、後手を先手に変更
	  */
      else if (token == "flip")       pos.flip();
	  /*
	  ベンチマーク
	  */
      else if (token == "bench")      benchmark(pos, is);
	  /*
	  標準出力に盤を表示
	  */
      else if (token == "d")          sync_cout << pos.pretty() << sync_endl;
	  /*
	  User　Interfaceからisreadyとのコマンドがきたら無条件でreadyokと返す
	  */
      else if (token == "isready")    sync_cout << "readyok" << sync_endl;
      else
          sync_cout << "Unknown command: " << cmd << sync_endl;

  } while (token != "quit" && argc == 1); // Passed args have one-shot behaviour

  Threads.wait_for_think_finished(); // Cannot quit whilst the search is running
}

/*
http://www.glaurungchess.com/shogi/usi.html

USIのコマンド部分(User Interfaceからゲームエンジンへのコマンド）

These are all the command the engine gets from the interface.

usi
	Tell engine to use the USI (universal shogi interface). This will be sent once as a first command after program boot to tell the engine to switch to USI mode. After receiving the usi command the engine must identify itself with the id command and send the option commands to tell the GUI which engine settings the engine supports. After that, the engine should send usiok to acknowledge the USI mode. If no usiok is sent within a certain time period, the engine task will be killed by the GUI.

debug [ on | off ]
	Switch the debug mode of the engine on and off. In debug mode the engine should send additional infos to the GUI, e.g. with the info string command, to help debugging, e.g. the commands that the engine has received etc. This mode should be switched off by default and this command can be sent any time, also when the engine is thinking.

isready
	This is used to synchronize the engine with the GUI. When the GUI has sent a command or multiple commands that can take some time to complete, this command can be used to wait for the engine to be ready again or to ping the engine to find out if it is still alive. This command is also required once before the engine is asked to do any search to wait for the engine to finish initializing. This command must always be answered with readyok and can be sent also when the engine is calculating in which case the engine should also immediately answer with readyok without stopping the search.

setoption name <id> [value <x>]
	This is sent to the engine when the user wants to change the internal parameters of the engine. For the button type no value is needed. One string will be sent for each parameter and this will only be sent when the engine is waiting. The name and value of the option in <id> should not be case sensitive and can not include spaces.

register
	This is the command to try to register an engine or to tell the engine that registration will be done later. This command should always be sent if the engine has sent registration error at program startup.

	The following tokens are allowed:

	later
	The user doesn't want to register the engine now.
	name <x>
	The engine should be registered with the name <x>
	code <y>
	The engine should be registered with the code <y>
	Example:


	"register later"
	"register name Stefan MK code 4359874324"

usinewgame
	This is sent to the engine when the next search (started with position and go) will be from a different game. This can be a new game the engine should play or a new game it should analyse but also the next position from a testsuite with positions only. As the engine's reaction to usinewgame can take some time the GUI should always send isready after usinewgame to wait for the engine to finish its operation.

position [sfen <sfenstring> | startpos ] moves <move1> ... <movei>
	Set up the position described in sfenstring on the internal board and play the moves on the internal board. If the game was played from the start position, the string startpos will be sent.

	Note: If this position is from a different game than the last position sent to the engine, the GUI should have sent a usinewgame inbetween.

go
	Start calculating on the current position set up with the position command. There are a number of commands that can follow this command, all will be sent in the same string.

	If one command is not sent its value should be interpreted as if it would not influence the search.

	searchmoves <move1> ... <movei>
	Restrict search to this moves only

	Example: After position startpos and go infinite searchmoves 7g7f 2g2f, the engine should only search the two moves P-7f and P-2f in the initial position.

ponder
	Start searching in pondering mode. Do not exit the search in ponder mode, even if it's mate! This means that the last move sent in in the position string is the ponder move. The engine can do what it wants to do, but after a ponderhit command it should execute the suggested move to ponder on. This means that the ponder move sent by the GUI can be interpreted as a recommendation about which move to ponder. However, if the engine decides to ponder on a different move, it should not display any mainlines as they are likely to be misinterpreted by the GUI because the GUI expects the engine to ponder on the suggested move.

	btime <x>
		Black has x milliseconds left on the clock
	wtime <x>
		White has x milliseconds left on the clock
	binc <x>
		Black increment per move i milliseconds if x > 0.
	winc <x>
		White increment per move i milliseconds if x > 0.
	movestogo <x>
		There are x moves to the next time control. This will only be sent if x > 0. If you don't get this and get the wtime and btime, it's sudden death.
	depth <x>
		Search x plies only.
	nodes <x>
		Search x nodes only.
	mate <x>
		Search for a mate in x moves. Programmers coming from computer chess should note that we follow the shogi convention for counting moves, i.e. we count what chess players usually describe as "plies" or "half moves". For instance, mate 5 means mate in 5 plies.
	movetime <x>
		Search exactly x milliseconds.
	infinite
		Search until the stop command is received. Do not exit the search without being told so in this mode!

stop
	Stop calculating as soon as possible. Don't forget the bestmove and possibly the ponder token when finishing the search.

ponderhit
	The user has played the expected move. This will be sent if the engine was told to ponder on the same move the user has played. The engine should continue searching but switch from pondering to normal search.

quit
	Quit the program as soon as possible.

5.3. Engine to GUI（ゲームエンジンからUser Interfaceへのコマンド)

id
	name <x>
	This must be sent after receiving the usi command to identify the engine, e.g. id name Shredder X.Y\n
	author <x>
	This must be sent after receiving the usi command to identify the engine, e.g. id author Stefan MK\n

usiok
	Must be sent after the id and optional options to tell the GUI that the engine has sent all infos and is ready in usi mode.

readyok
	This must be sent when the engine has received an isready command and has processed all input and is ready to accept new commands now. It is usually sent after a command that can take some time to be able to wait for the engine, but it can be used anytime, even when the engine is searching, and must always be answered with readyok.

bestmove <move1> [ponder <move2>]
	The engine has stopped searching and found the move <move> best in this position. The engine can send the move it likes to ponder on. The engine must not start pondering automatically. this command must always be sent if the engine stops searching, also in pondering mode if there is a stop command, so for every go command a bestmove command is needed! Directly before that the engine should send a final info command with the final search information, the the GUI has the complete statistics about the last search.

copyprotection
	This is needed for copyprotected engines. After the usiok command the engine can tell the GUI, that it will check the copy protection now. This is done by copyprotection checking. If the check is ok the engine should send copyprotection ok, otherwise copyprotection error.

	If there is an error the engine should not function properly but should not quit alone. If the engine reports copyprotection error the GUI should not use this engine and display an error message instead!

	The code in the engine can look like this:


	TellGUI("copyprotection checking\n");
	// ... check the copy protection here ...
	if(ok)
	TellGUI("copyprotection ok\n");
	else
	TellGUI("copyprotection error\n");

registration
	This is needed for engines that need a username and/or a code to function with all features. Analogously to the copyprotection command the engine can send registration checking after the usiok command followed by either registration ok or registration error. Also after every attempt to register the engine it should answer with registration checking and then either registration ok or registration error.

	In contrast to the copyprotection command, the GUI can use the engine after the engine has reported an error, but should inform the user that the engine is not properly registered and might not use all its features.

	In addition the GUI should offer to open a dialog to enable registration of the engine. To try to register an engine the GUI can send the register command. The GUI has to always answer with the register command if the engine sends registration error at engine startup (this can also be done with register later) and tell the user somehow that the engine is not registered. This way the engine knows that the GUI can deal with the registration procedure and the user will be informed that the engine is not properly registered.

info
	The engine wants to send information to the GUI. This should be done whenever one of the info has changed.

	The engine can send only selected infos or multiple infos with one info command, e.g. info currmove 2g2f currmovenumber 1 or info depth 12 nodes 123456 nps 100000.

	Also all infos belonging to the pv should be sent together, e.g.


	info depth 2 score cp 214 time 1242 nodes 2124 nps 34928 pv 2g2f 8c8d 2f2e

	I suggest to start sending currmove, currmovenumber, currline and refutation only after one second in order to avoid too much traffic.

	Additional info:

	depth <x>
		Search depth in plies.
	seldepth <x>
		Selective search depth in plies. If the engine sends seldepth there must also be a depth present in the same string.
	time <x>
		The time searched in ms. This should be sent together with the pv.
	nodes <x>
		x nodes searched. The engine should send this info regularly.
	pv <move1> ... <movei>
		The best line found.
	multipv <num>
		This for the multi pv mode. For the best move/pv add multipv 1 in the string when you send the pv. In k-best mode, always send all k variants in k strings together.
		score
	cp <x>
		The score from the engine's point of view, in centipawns.
	mate <y>
		Mate in y plies. If the engine is getting mated, use negative values for y.
	lowerbound
		The score is just a lower bound.
	upperbound
		The score is just an upper bound.
	currmove <move>
		Currently searching this move.
	currmovenumber <x>
		Currently searching move number x, for the first move x should be 1, not 0.
	hashfull <x>
		The hash is x permill full. The engine should send this info regularly.
	nps <x>
		x nodes per second searched. the engine should send this info regularly.
	cpuload <x>
		The cpu usage of the engine is x permill.
	string <str>
		Any string str which will be displayed be the engine. if there is a string command the rest of the line will be interpreted as <str>.
	refutation <move1> <move2> ... <movei>
		Move <move1> is refuted by the line <move2> ... <movei>, where i can be any number >= 1. Example: after move 8h2b+ is searched, the engine can send info refutation 8h2b+ 1c2b if 1c2b is the best answer after 8h2b+ or if 1c2b refutes the move 8h2b+. If there is no refutation for 8h2b+ found, the engine should just send info refutation 8h2b+. The engine should only send this if the option USI_ShowRefutations is set to true.
	currline <cpunr> <move1> ... <movei>
		This is the current line the engine is calculating. <cpunr> is the number of the cpu if the engine is running on more than one cpu. <cpunr> = 1,2,3.... If the engine is just using one cpu, <cpunr> can be omitted. If <cpunr> is greater than 1, always send all k lines in k strings together. The engine should only send this if the option USI_ShowCurrLine is set to true.

option
	This command tells the GUI which parameters can be changed in the engine. This should be sent once at engine startup after the usi and the id commands if any parameter can be changed in the engine. The GUI should parse this and build a dialog for the user to change the settings. Note that not every option should appear in this dialog, as some options like USI_Ponder, USI_AnalyseMode, etc. are better handled elsewhere or are set automatically.

	If the user wants to change some settings, the GUI will send a setoption command to the engine.

	Note that the GUI need not send the setoption command when starting the engine for every option if it doesn't want to change the default value. For all allowed combinations see the examples below, as some combinations of this tokens don't make sense.

	One string will be sent for each parameter.

	name <id>
	The option has the name <id>. Whitespace is not allowed in an option name. Note that the name should normally not be displayed directly in the GUI: The GUI should look up the option name in the translation file, and present the translation into the users preferred language in the engine's option dialog.

	Certain options have a fixed value for <id>, which means that the semantics of this option is fixed. Usually those options should not be displayed in the normal engine options window of the GUI but get a special treatment. USI_Pondering for example should be set automatically when pondering is enabled or disabled in the GUI options. The same for USI_AnalyseMode which should also be set automatically by the GUI. All those certain options have the prefix USI_. If the GUI gets an unknown option with the prefix USI_, it should just ignore it and not display it in the engine's options dialog.

	The options with fixed semantics are:

	<id> = USI_Hash, type spin
	The value in MB for memory for hash tables can be changed, this should be answered with the first setoptions command at program boot if the engine has sent the appropriate option name Hash command, which should be supported by all engines! So the engine should use a very small hash first as default.
	<id> = USI_Ponder, type check
	This means that the engine is able to ponder (i.e. think during the opponent's time). The GUI will send this whenever pondering is possible or not. Note: The engine should not start pondering on its own if this is enabled, this option is only needed because the engine might change its time management algorithm when pondering is allowed.
	<id> = USI_OwnBook, type check
	This means that the engine has its own opening book which is accessed by the engine itself. If this is set, the engine takes care of the opening book and the GUI will never execute a move out of its book for the engine. If this is set to false by the GUI, the engine should not access its own book.
	<id> = USI_MultiPV, type spin
	The engine supports multi best line or k-best mode. The default value is 1.
	<id> = USI_ShowCurrLine, type check
	The engine can show the current line it is calculating. See info currline above. This option should be false by default.
	<id> = USI_ShowRefutations, type check
	The engine can show a move and its refutation in a line. See info refutations above. This option should be false by default.
	<id> = USI_LimitStrength, type check
	The engine is able to limit its strength to a specific dan/kyu number. This should always be implemented together with USI_Strength. This option should be false by default.
	<id> = USI_Strength, type spin
	The engine can limit its strength within the given interval. Negative numbers are kyu levels, while positive numbers are amateur dan levels. If USI_LimitStrength is set to false, this value should be ignored. If USI_LimitStrength is set to true, the engine should play with this specific strength. This option should always be implemented together with USI_LimitStrength.
	<id> = USI_AnalyseMode, type check
	The engine wants to behave differently when analysing or playing a game. For example when playing it can use some kind of learning, or an asymetric evaluation function. The GUI should set this option to false if the engine is playing a game, and to true if the engine is analysing.
	type <t>
	The option has type t. There are 5 different types of options the engine can send:

	check
	A checkbox that can either be true or false.
	spin
	A spin wheel or slider that can be an integer in a certain range.
	combo
	A combo box that can have different predefined strings as a value.
	button
	A button that can be pressed to send a command to the engine
	string
	A text field that has a string as a value, an empty string has the value <empty>.
	filename
	Similar to string, but is presented as a file browser instead of a text field in the GUI.
	default <x>
	The default value of this parameter is x.

	min <x>
	The minimum value of this parameter is x.

	max <x>
	The maximum value of this parameter is x.

	Here are some examples illustrating the different types of options:


	"option name Nullmove type check default true\n"
	"option name Selectivity type spin default 2 min 0 max 4\n"
	"option name Style type combo default Normal var Solid var Normal var Risky\n"
	"option name LearningFile type filename default /shogi/my-shogi-engine/learn.bin"
	"option name ResetLearning type button\n"

	
	*/