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
#include <iostream>

#include "bitboard.h"
#include "evaluate.h"
#include "position.h"
#include "search.h"
#include "thread.h"
#include "tt.h"
#include "ucioption.h"

/*test start*/
using namespace Bitboards;
enum LLType{
	python,ruby,perl
};
void test(void);
template<LLType LT> void print(void);
/*test end*/

int main(int argc, char* argv[]) {

  std::cout << engine_info() << std::endl;

  UCI::init(Options);
  Bitboards::init();
  Position::init();
  Bitbases::init_kpk();
  Search::init();
  Pawns::init();
  Eval::init();
  Threads.init();
  TT.resize(Options["Hash"]);
  test();
  UCI::loop(argc, argv);

  Threads.exit();
}

//codeをみただけではわからない
//いろいろ試して理解を促進する
void test(void)
{
	/*
    Score temp = make_score(24,11);
    printf("A1 %d\n",(~SQ_A1));
    printf("A2 %d\n",(~SQ_A2));
    printf("A3 %d\n",(~SQ_A3));
    printf("A4 %d\n",(~SQ_A4));
    printf("A5 %d\n",(~SQ_A5));
    printf("A6 %d\n",(~SQ_A6));
    printf("A7 %d\n",(~SQ_A7));
    printf("A8 %d\n",(~SQ_A8));
    printf("B1 %d\n",(~SQ_B1));
	*/
	/*
    for(File f = FILE_A;f < FILE_NB;++f){
        for(Rank r =  RANK_1;r < RANK_NB;++r){
           printf("square = %d\n",make_square(f,r));
        }
    }
    printf("%d\n",SQ_A1 & 7);
		*/
    //relative_square
		/*
    Square sq;
    for(sq = SQ_A1;sq < SQ_NONE;++sq){
        printf(" %2d",relative_square(WHITE,sq));
        if((sq % 8)==7){
            printf("\n");
        }
    }
    printf("\n");
    for(sq = SQ_A1;sq < SQ_NONE;++sq){
        printf(" %2d",relative_square(BLACK,sq));
        if((sq % 8)==7){
            printf("\n");
        }
    }
	*/
    //relative_rank
	/*
    printf("relative_rank\n");
    for(sq = SQ_A1;sq < SQ_NONE;++sq){
        printf(" %2d",relative_rank(WHITE,sq));
        if((sq % 8)==7){
            printf("\n");
        }
    }
    printf("\n");
    for(sq = SQ_A1;sq < SQ_NONE;++sq){
        printf(" %2d",relative_rank(BLACK,sq));
        if((sq % 8)==7){
            printf("\n");
        }
    }
	*/
    //opposite_colors
	/*
    printf("\n");
    Square s2 = SQ_A1;
    for(sq = SQ_A1;sq < SQ_NONE;++sq){
        printf(" %2d",opposite_colors(sq,SQ_A2));
        if((sq % 8)==7){
            printf("\n");
        }
    }
	
    File f = FILE_B;
    printf("%c\n",to_char(f,true));
    printf("%c\n",to_char(f,false));

    Rank r = RANK_1;
    printf("%c\n",to_char(r));
	*/
	//templateの実験
	/*
	print<python>();
	*/
	//DistanceRingsBB変数の初期値確認
	/*
	for (Square s1 = SQ_A1; s1 <= SQ_H8; ++s1){
			for (Square s2 = SQ_A1; s2 <= SQ_H8; ++s2){
					if (s1 != s2)
					{
							printf("s1=%d,s2=%d\n",s1,s2);
							printf("%s",Bitboards::pretty(DistanceRingsBB[s1][s2]).c_str());
					}
			}
	}
	*/
	//StepAttacksBBのなかみ確認
	/*
	Piece piece_code = W_KING;
		
	sq = SQ_A1;
	printf("piece=%d,sq=%d\n",piece_code,sq);
	printf("%s",Bitboards::pretty(StepAttacksBB[piece_code][sq]).c_str());

	sq = SQ_B1;
	printf("piece=%d,sq=%d\n",piece_code,sq);
	printf("%s",Bitboards::pretty(StepAttacksBB[piece_code][sq]).c_str());

	sq = SQ_C1;
	printf("piece=%d,sq=%d\n",piece_code,sq);
	printf("%s",Bitboards::pretty(StepAttacksBB[piece_code][sq]).c_str());
		
	sq = SQ_A2;
	printf("piece=%d,sq=%d\n",piece_code,sq);
	printf("%s",Bitboards::pretty(StepAttacksBB[piece_code][sq]).c_str());
		
	sq = SQ_B2;
	printf("piece=%d,sq=%d\n",piece_code,sq);
	printf("%s",Bitboards::pretty(StepAttacksBB[piece_code][sq]).c_str());
	*/
	//RAttacksのなかみ確認
	/*
	sq = SQ_A1;
	printf("piece=%d,sq=%d\n",ROOK,sq);
	printf("%s",Bitboards::pretty(*RAttacks[sq]).c_str());
	sq = SQ_B2;
	printf("piece=%d,sq=%d\n",ROOK,sq);
	printf("%s",Bitboards::pretty(*RAttacks[sq]).c_str());
	sq = SQ_C3;
	printf("piece=%d,sq=%d\n",ROOK,sq);
	printf("%s",Bitboards::pretty(*RAttacks[sq]).c_str());
	sq = SQ_D4;
	printf("piece=%d,sq=%d\n",ROOK,sq);
	printf("%s",Bitboards::pretty(*RAttacks[sq]).c_str());
	sq = SQ_E5;
	printf("piece=%d,sq=%d\n",ROOK,sq);
	printf("%s",Bitboards::pretty(*RAttacks[sq]).c_str());
	sq = SQ_F6;
	printf("piece=%d,sq=%d\n",ROOK,sq);
	printf("%s",Bitboards::pretty(*RAttacks[sq]).c_str());
	sq = SQ_G7;
	printf("piece=%d,sq=%d\n",ROOK,sq);
	printf("%s",Bitboards::pretty(*RAttacks[sq]).c_str());
	sq = SQ_H8;
	printf("piece=%d,sq=%d\n",ROOK,sq);
	printf("%s",Bitboards::pretty(*RAttacks[sq]).c_str());

	sq = SQ_A1;
	printf("piece=%d,sq=%d\n",BISHOP,sq);
	printf("%s",Bitboards::pretty(*BAttacks[sq]).c_str());
	sq = SQ_B2;
	printf("piece=%d,sq=%d\n",BISHOP,sq);
	printf("%s",Bitboards::pretty(*BAttacks[sq]).c_str());
	sq = SQ_C3;
	printf("piece=%d,sq=%d\n",BISHOP,sq);
	printf("%s",Bitboards::pretty(*BAttacks[sq]).c_str());
	sq = SQ_D4;
	printf("piece=%d,sq=%d\n",BISHOP,sq);
	printf("%s",Bitboards::pretty(*BAttacks[sq]).c_str());
	sq = SQ_E5;
	printf("piece=%d,sq=%d\n",BISHOP,sq);
	printf("%s",Bitboards::pretty(*BAttacks[sq]).c_str());
	sq = SQ_F6;
	printf("piece=%d,sq=%d\n",BISHOP,sq);
	printf("%s",Bitboards::pretty(*BAttacks[sq]).c_str());
	sq = SQ_G7;
	printf("piece=%d,sq=%d\n",BISHOP,sq);
	printf("%s",Bitboards::pretty(*BAttacks[sq]).c_str());
	sq = SQ_H8;
	printf("piece=%d,sq=%d\n",BISHOP,sq);
	printf("%s",Bitboards::pretty(*BAttacks[sq]).c_str());
	*/
	//flipの確認
	/*
	Position pos("rnbqkbnr/pppppppp/8/8/8/8/PPPPPPPP/RNBQKBNR w KQkq - 0 1", false, Threads.main()); 
	pos.flip();
	*/
	//SquareBBの確認
	/*
	printf("SquareBB\n");
	for (Square s = SQ_A1; s <= SQ_H8; ++s){
			printf("%s",Bitboards::pretty(SquareBB[s]).c_str());
	}
	*/
	//PseudoAttacksの確認
	/*
	printf("PseudoAttacks Rook\n");
	for (Square s = SQ_A1; s <= SQ_H8; ++s){
			printf("%s",Bitboards::pretty(PseudoAttacks[ROOK][s]).c_str());
	}
		
	printf("PseudoAttacks Bishop\n");
	for (Square s = SQ_A1; s <= SQ_H8; ++s){
			printf("%s",Bitboards::pretty(PseudoAttacks[BISHOP][s]).c_str());
	}
	printf("PseudoAttacks Queen\n");
	for (Square s = SQ_A1; s <= SQ_H8; ++s){
			printf("%s",Bitboards::pretty(PseudoAttacks[QUEEN][s]).c_str());
	}
	*/
	//LineBB[sq1][sq2]の確認
	/*
	printf("LineBB[sq1][sq2]\n");
	for (Square s1 = SQ_A1; s1 <= SQ_H8; ++s1){
		for(Square s2 = SQ_A1; s2 <= SQ_H8; ++s2){
			printf("%d %d \n %s",s1,s2,Bitboards::pretty(LineBB[s1][s2]).c_str());
		}
	}
	*/
	//BetweenBB[s1][s2]の確認
	/*
	printf("BetweenBB[sq1][sq2]\n");
	for (Square s1 = SQ_A1; s1 <= SQ_H8; ++s1){
		for(Square s2 = SQ_A1; s2 <= SQ_H8; ++s2){
			printf("%d %d \n %s",s1,s2,Bitboards::pretty(BetweenBB[s1][s2]).c_str());
		}
	}
	*/
	//pop_lsbの確認、コードが理解できないので、挙動を把握
	//pop_lsbは直接呼び出せないのでlsb関数で呼ぶ
	/*
	Bitboard bb;
	bb = 0xDD9;	//b1101 1101 1001
	printf("sq=%d\n",lsb(bb));
	printf("sq=%d\n",lsb(bb-1));
	*/
	//ForwardBBの確認
	/*
	printf("ForwardBB[color][sq]\n");
	for (Color c = WHITE; c <= BLACK; ++c){
		for(Square s = SQ_A1; s <= SQ_H8; ++s){
			printf("color=%d sq=%d \n %s",c,s,Bitboards::pretty(ForwardBB[c][s]).c_str());
		}
	}
	*/
	//PawnAttackSpanの確認
	/*
	printf("PawnAttackSpan[color][sq]\n");
	for (Color c = WHITE; c <= BLACK; ++c){
		for(Square s = SQ_A1; s <= SQ_H8; ++s){
			printf("color=%d sq=%d \n %s",c,s,Bitboards::pretty(PawnAttackSpan[c][s]).c_str());
		}
	}
	*/
	/*
	printf("PassedPawnMask[color][sq]\n");
	for (Color c = WHITE; c <= BLACK; ++c){
		for(Square s = SQ_A1; s <= SQ_H8; ++s){
			printf("color=%d sq=%d \n %s",c,s,Bitboards::pretty(PassedPawnMask[c][s]).c_str());
		}
	}
	*/
	/*
	printf("FileABB\n");
	printf(" %s", Bitboards::pretty(FileABB).c_str());
	*/
	/*
	printf("SquareBB\n");
	printf(" %s", Bitboards::pretty(SquareBB[SQ_A1]).c_str());
	printf(" %s", Bitboards::pretty(SquareBB[SQ_A2]).c_str());
	printf(" %s", Bitboards::pretty(SquareBB[SQ_A3]).c_str());
	printf(" %s", Bitboards::pretty(SquareBB[SQ_A4]).c_str());
	printf(" %s", Bitboards::pretty(SquareBB[SQ_A5]).c_str());
	printf(" %s", Bitboards::pretty(SquareBB[SQ_A6]).c_str());
	printf(" %s", Bitboards::pretty(SquareBB[SQ_A7]).c_str());
	printf(" %s", Bitboards::pretty(SquareBB[SQ_A8]).c_str());

	printf(" %s", Bitboards::pretty(SquareBB[SQ_B1]).c_str());
	printf(" %s", Bitboards::pretty(SquareBB[SQ_B2]).c_str());
	printf(" %s", Bitboards::pretty(SquareBB[SQ_B3]).c_str());
	printf(" %s", Bitboards::pretty(SquareBB[SQ_B4]).c_str());
	printf(" %s", Bitboards::pretty(SquareBB[SQ_B5]).c_str());
	printf(" %s", Bitboards::pretty(SquareBB[SQ_B6]).c_str());
	printf(" %s", Bitboards::pretty(SquareBB[SQ_B7]).c_str());
	printf(" %s", Bitboards::pretty(SquareBB[SQ_B8]).c_str());
	*/
	printf("Rank1BB\n");
	printf(" %s", Bitboards::pretty(Rank1BB).c_str());
return;
}


template<LLType LT> void print(void)
{
	switch(LT){
	case python:
		printf("python\n");
		break;
	case ruby:
		printf("ruby\n");
		break;
	case perl:
		printf("perl\n");
		break;
	}
}