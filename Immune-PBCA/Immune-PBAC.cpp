// Immune-PBAC.cpp : 定义控制台应用程序的入口点。
//


#include "PBCA.h"
using namespace std;
 
int main(void)
{
	/***********/

  PBCA *a=new PBCA();
  char num=0;
  //a->Run();
  //a->Test_Run();
  cout<<"1.Read raw map 2.Read map directly"<<endl;
  num=getchar();
  while(num!='1'&&num!='2')
  {
	  cout<<"fail to chose,chose again."<<endl;
	  num=getchar();
  }
  if(num='1')
	  a->Make_Map_EUD();

	a->Run();
	getchar();

	return 0;
}

