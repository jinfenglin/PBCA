#pragma once
#include <list>
#include <iostream>
#include <fstream>
#include <ctime>
#include <cmath>
#include<cstdlib>
#include<algorithm>
#include <malloc.h>
#include <stdlib.h>
#include <string>
#include <deque>
#include <io.h>

using namespace std;
#define MaxSize 1000 
#define MaxLength 300
#define PATH "../input file.txt"
#define LOG "../log.txt"
enum raw_map_type {EUC2D,MATRIX,UNKNOWN};
struct TSP_input
{
	int L;            //L is the length of sequence,the number of cities
	int start_point;  //where the merchant start
	int round;         //rounds need to run	
	int times;        //How many times this parameter will be used

	int pz;           //pz is the size of safe zone
	int pmax;         //pmax is the max population
};
struct sequence      //struct of the phage/bacteria
{
	int seq[MaxSize]; 
	int length;
	int fitness;
};
struct map
{
	int matrix[MaxLength][MaxLength];
	int length;
	int number;
};
struct EUC_2D_Node
{
	int node;
	int x,y;
};
struct HGT_Segment
{
	int start_point;
	int segments[MaxLength];
	int L;
	int value;
};
class PBCA
{
	friend class HGT_Class;
public:
	PBCA(void);
	~PBCA(void);
	int valid_time;//allow others to use valid_time from tsp_in
	bool input(ifstream &ifs );//Get parameter from input to fullfill tsp_in, have 2 ways to read 1.consoler 2.text file read. This is the text version
	bool input();         //consoler version reading
	bool random_produce(int amount,list<sequence> &list); //randomly produce genes
	bool random_produce();           //if amount is miss will produce the same number of pz
	void Run();

	void Mutate_Selection();        //let the phage evolve
	int fitness(sequence&);          //calculate the fitness
	void Make_Map_EUD();             //change the content of map.txt according to raw_map.txt
	void Read_Map();                //only 1 map pertime
	void Select_Map();              //select map for the parameter
	void mutate(sequence);          //mutate on a single sequence
	void mutate_standard(sequence); //mutate use the same operator with GA example
	void mutate_multi(sequence);  //reverse a segement
	void All_fitness(list<sequence>&);
	void replace(sequence b,sequence p);
	void replace_EPR(sequence b,sequence p);
	void select(int number);       //amount want to delete
	void fight();                  //competition between phage and bacteria
	void Test_Run();
	void population_regulation();
	/*Test tools*/
	void Test_showparameters();//output parameter
	void Test_showSequenceArray();//output bacteria and phage array
	void Test_showMapVector();//output map vecto
	int Output();
	bool Is_SequenceInList(sequence,list<sequence>);
	friend bool operator ==(sequence s1,sequence s2);
	/*HGT*/
	void HGT();
	void HGT_Crossover();
	void Sample();
	void Distribution();
	void HGT_Replace(deque<HGT_Segment>::iterator it,list<sequence>::iterator it_b);

private:
	int replace_count,mutate_count;
	TSP_input tsp_in;
	int PBCA_Run_count;
	map current_map;
	list<map> map_vector;
	list<sequence> bacteria,phage;
	int* auto_select();
	void make_e2_map(list<EUC_2D_Node> elist,int);
	int largest_amount();
	float Survive_Rate(sequence,sequence);
	void Load_To_Log(ofstream &ofs);
	bool is_InTheGene(int numb,sequence S,int length);//Check out whether number int is in the sequence 
	raw_map_type map_type;
	int HGT_Count;
	deque<HGT_Segment> Hlist;
};                           
                            
