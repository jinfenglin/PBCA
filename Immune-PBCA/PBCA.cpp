//#include "PBCA.h"
#include <unistd.h>
#include <stdio.h>
#include "HGT.h"
#include <limits.h>
#include "ga/ga.h"


#define intact -1 
bool CompareRule_max(sequence s1,sequence s2)
{
	if(s1.fitness>s2.fitness)
		return true;
	else
		return false;
}
bool CompareRule_min(sequence s1,sequence s2) //s1 is better than s2 then return true
{
	if(s1.fitness<s2.fitness)
		return true;
	else
		return false;
}
PBCA::PBCA(void)
{
	//initialize the parameter
	tsp_in.L=intact;
	tsp_in.pmax=intact;
	tsp_in.pz=intact;
	tsp_in.round=intact;
	tsp_in.start_point=intact;
	tsp_in.times=intact;
	current_map.length=0;
	replace_count=0;
	mutate_count=0;
	srand (time(NULL)); // This will ensure a really randomized number by help of valid_time.
}
PBCA::~PBCA(void)
{
}
bool PBCA::input() //each round get 5 parameter always true
{
	cout<<"Please input parameters in order:"<<endl;
	cout<<"1.Number of cities 2.Start point 3.Rounds to run"<<endl;
	cout<<"4.Max population 5.Size of the safe zone"<<endl;
	while(1)
	{
		cin>>tsp_in.L>>tsp_in.start_point>>tsp_in.round>>tsp_in.pmax>>tsp_in.pz;
		if(tsp_in.L<=intact||tsp_in.start_point<=intact||tsp_in.round<=intact||tsp_in.pmax<=intact||tsp_in.pz<=intact)//0 are not included here...but it should be...
		cout<<"Error detected,please input again"<<endl;
		else
			{
				tsp_in.times=1;
				cin.ignore(INT_MAX,'\n');                      //everyround have 5 paremeter if not clean may cause bug
				cout<<"Input completed"<<endl;
				return true;
		}
	}
}
bool PBCA::input(ifstream &ifs) //each valid_time get 1 line and kept it for tsp.times rounds when file is over return false
{
	cout<<"Begin to read file"<<endl;
	if(!ifs.eof())
	{
		ifs>>tsp_in.L>>tsp_in.start_point>>tsp_in.round>>tsp_in.pmax>>tsp_in.pz>>tsp_in.times;
		cout<<"reading finished"<<endl;
		valid_time=tsp_in.times;
		return true;
	}
	else
		return false;
}
void PBCA::Test_showparameters()
{
	cout<<"tsp_in.L="<<tsp_in.L<<" "<<"tsp_in.start_point="<<tsp_in.start_point<<" "<<"tsp_in.round="<<tsp_in.round<<" ";
	cout<<"tsp_in.pmax="<<tsp_in.pmax<<" "<<"tsp_in.pz="<<tsp_in.pz<<" "<<"tsp_in.times="<<tsp_in.times<<" "<<endl;
}
int PBCA::largest_amount()
{
	int amount=1;
	for(int i=1;i<=tsp_in.L;i++)
			amount=amount*i;
	return amount;
}
int PBCA::fitness(sequence &S)
{ 
	int fitness=0;
	int now=0,next=0;
	int now_location,next_location;
	for(int i=0;i<S.length;i++)
	{

		now_location=i;
		next_location=(now_location+1)%S.length;			
		now=S.seq[now_location]-1;
		next=S.seq[next_location]-1;                 //sequence start form 1 matrix start from 0 thus -1
		fitness+=current_map.matrix[now][next];
	}
	if(fitness<0)
		cout<<"error in fitness"<<endl;
	S.fitness=fitness;
	return fitness;
}
bool PBCA::random_produce(int amount,list<sequence> &list)
{
	int dice_roll;
	/*unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	default_random_engine generator(seed);
	uniform_int_distribution<int> distribution(1,tsp_in.L); //produce random number in the interval [1,tsp_in.L]*/	
	if(amount>largest_amount()&&tsp_in.L<=15)//prevent the sequence is too short thus make the overlaping of exited individuals
	{
		amount=largest_amount();
		cout<<"amount overflow,after adjusting produce "<<amount<<"individuals"<<endl;
	}

	for(int i=0;i<amount;i++)//amount stand for the quantity of sequence
	{
		sequence S;          /*initialize the sequence*/
		S.fitness=intact;
		S.length=tsp_in.L;
		for(int n=0;n<tsp_in.L;n++)
		{
			while(1)// To generate leagal genom to put into sequence
			{
				dice_roll=rand()%tsp_in.L+1; // Randomizing the number between 1-tsp_in.
				if(!is_InTheGene(dice_roll,S,n))
					break;
			}
			S.seq[n]=dice_roll;
		}
		if(find(list.begin(),list.end(),S)==list.end())
			list.push_front(S);
		else
			i--;

	}
	return true; //tp
}
bool PBCA::is_InTheGene(int numb,sequence S,int length) //numb in the S return true 
{
	for(int i=0;i<length;i++)
	{
		if(numb==S.seq[i])
			return true;
	}
	return false;
}
void PBCA::Test_showSequenceArray()
{
	list<sequence>::iterator it;
	if(!bacteria.empty())
	{
		cout<<"bacteria sequence:"<<endl;
		for(it=bacteria.begin();it!=bacteria.end();it++)
			{
				for(int i=0;i<tsp_in.L;i++)
					cout<<it->seq[i];
				cout<<"   fitness:"<<it->fitness<<endl;
		}
	}
	else
		cout<<"Bacteria list is empty!"<<endl;

	if(!phage.empty())
	{
		cout<<"phage sequence:"<<endl;
		for(it=phage.begin();it!=phage.end();it++)
			{
				for(int i=0;i<tsp_in.L;i++)
					cout<<it->seq[i];
				cout<<"   fitness:"<<it->fitness<<endl;
		}
	}
	else
		cout<<"phage list is empty!"<<endl;
}
void PBCA::mutate(sequence S)  //add a new array into the phage array
{
	sequence copy_s=S;
	int r1=rand()%S.length;
	int r2=rand()%S.length;
	int temp=copy_s.seq[r2];
	copy_s.seq[r2]=copy_s.seq[r1];
	copy_s.seq[r1]=temp;//exchange
	/*if(find(phage.begin(),phage.end(),copy_s)==phage.end())*/
	phage.push_back(copy_s);
}
void PBCA::mutate_standard(sequence S)
{
	sequence copy_s=S;
	if(GARandomFloat()<0.5)
	{
		int r1=rand()%S.length;
		int r2=rand()%S.length;
		int temp=copy_s.seq[r2];
		copy_s.seq[r2]=copy_s.seq[r1];
		copy_s.seq[r1]=temp;
	}
	else
	{
		list<int> copy_s_list;
		for(int i=0;i<S.length;i++)//transfer into list
		{
			copy_s_list.push_back(copy_s.seq[i]);
		}
		int nNodes=GARandomInt(1,((int)(S.length/2-1)));
		int start_place=GARandomInt(0,S.length-nNodes-1);
		list<int> temp;
		list<int>::iterator it=copy_s_list.begin();
		for(int i=0;i<start_place;i++)
			it++;
		for(int i=0;i<nNodes;i++)//store the array want to displace
		{
			if(*it<0)
				cout<<"error"<<endl;
			temp.push_back(*it);
			it=copy_s_list.erase(it);
			
		}
		int insert_place=GARandomInt(0,S.length-nNodes-1);
		it=copy_s_list.begin();
		for(int i=0;i<insert_place;i++)
			it++;
		if(GARandomFloat()<0.5) //reverse nNodes
		{
			temp.reverse();	
		}
		copy_s_list.insert(it,temp.begin(),temp.end());	
		for(int i=0;i<S.length;i++)//transfer list to array
		{
			if(copy_s_list.front()<0)
				cout<<"error"<<endl;
			copy_s.seq[i]=copy_s_list.front();
			copy_s_list.pop_front();
		}
	}
	phage.push_back(copy_s);
}
void PBCA::mutate_multi(sequence S)
{
	sequence copy_s=S;
	int round=GARandomInt(1,S.length/2);
	for(int i=0;i<round;i++)
	{
		int r1=rand()%S.length;
		int r2=rand()%S.length;
		
		int temp=copy_s.seq[r2];
		copy_s.seq[r2]=copy_s.seq[r1];
		copy_s.seq[r1]=temp;//exchange
	}		
	/*if(find(phage.begin(),phage.end(),copy_s)==phage.end())*/
		phage.push_back(copy_s);

}
bool operator ==(sequence s1,sequence s2)
{
	if(s1.length!=s2.length)
		return false;
	else
	{
		for(int i=0;i<s1.length;i++)
		{
			if(s1.seq[i]!=s2.seq[i])
				return false;
		}
		return true;
	}
}
void PBCA::All_fitness(list<sequence> &lis)
{
	list<sequence>::iterator it;
	for( it=lis.begin();it!=lis.end();it++)
	{
		fitness(*it);
	}
}
void PBCA::select(int number) //chop the last numbers items form phage array
{
	for(int i=0;i<number;i++)
	{
		if(phage.size()<tsp_in.pz) //chop the worst part above safe zone
			break;
		else
		{
			phage.pop_back();
		}
	}	
}
void PBCA::Mutate_Selection()//mutate operateor works here
{
	float mp;
	int population=phage.size();
	float dicer; 
	list<sequence>::iterator it;

	mp=1.0-(float)phage.size()/(phage.size()+bacteria.size());
	for(it=phage.begin();it!=phage.end();it++)
	{
		dicer=(float)rand()/(float)RAND_MAX;//product random number between 0-1
		if(dicer<mp)                        //another way to mutate
			mutate_standard(*it);
			//mutate_multi(*it);
			//mutate(*it);
	}
	All_fitness(phage);
	phage.sort(CompareRule_min);
	//select(phage.size()-population);  //mutation will not expand the population
	//Test_showSequenceArray();

	
}
int* PBCA::auto_select()
{
	int ptr[MaxLength];
	int count=0;
	//ptr=(int *)malloc(sizeof(int));
	//Test_showMapVector();
	for(list<map>::iterator it=map_vector.begin();it!=map_vector.end();it++)
	{
		if(it->length==tsp_in.L)
		{
			ptr[count]=it->number;
			count++;
		}
	}	
	ptr[count]=-1;
	return ptr;
}
void PBCA::Select_Map()
{
	int record;
	if(current_map.length==0)
		cout<<"Have no map selected!"<<endl;
	else
		{
			cout<<"Current map is:No."<<current_map.number<<endl;
			/*cout<<"Show current map?Y/N"<<endl;
			
			char temp=' ';
			while(temp!='y'&&temp!='Y'&&temp!='n'&&temp!='N')
				temp=getchar();
			if(temp=='y'||temp=='Y')
			{
				for(int i=0;i<current_map.length;i++)
				{
					for(int j=0;j<current_map.length;j++)
					{
						cout<<current_map.matrix[i][j]<<" ";
					}
					cout<<endl;
				}
			}*/
			cout<<endl;
	}//show ccurrent map-end of if(have a current map)

			cout<<"Available maps:"<<endl; 
			int *pointer=auto_select();  //pick out the map whose length match with parameter
			int temp_array[MaxLength];
			for(int i=0;pointer[i]>=0;i++)  //copy picked map(their number) to the temp_array
			{
				temp_array[i]=pointer[i];
			}
			for(int i=0;temp_array[i]>=0;i++)   //out put the choice we have
			{                
				cout<<"No."<<temp_array[i]<<" ";//Abnormally change->pointer maybe caused by the attrubute of pointer.So copy the number to array,should use cpy
				record=i;                      //record the available map array length        
			}
			cout<<endl;
			/*end of building temp_array*/

		    cout<<"Auto select?Y/N?"<<endl; 
			char temp_char=' ';
			while(temp_char!='y'&&temp_char!='Y'&&temp_char!='n'&&temp_char!='N')
				temp_char=getchar();
			int temp_location;
			if(temp_char=='y'||temp_char=='Y')  //randomly pick one
			{
				                       
				cout<<"record="<<record<<endl;
				temp_location=rand()%(record+1);
				/*copy pick matrix to the current one*/
				current_map.number=temp_array[temp_location];
				for(list<map>::iterator it=map_vector.begin();it!=map_vector.end();it++) //Find No.temp MAP in the array and copy it to current map
				{
					if(it->number==temp_array[temp_location])
					{
						current_map.length=it->length;
						//current_map.matrix=it->matrix;  Be careful here
						for(int i=0;i<current_map.length;i++)
							for(int n=0;n<current_map.length;n++)
							{
								current_map.matrix[i][n]=it->matrix[i][n];
							}
					}
				}//copy-over
			}
			else //pick one by hand
			{
				while(1) //check out the input is legle or not
				{
					cout<<"Input the number of map:"<<endl;
					int temp_num=0;
					cin>>temp_num;
					for(list<map>::iterator it=map_vector.begin();it!=map_vector.end();it++)
					{
						if(it->number==temp_num&&it->length==tsp_in.L)
						{
							current_map.length=it->length;
							for(int i=0;i<current_map.length;i++)
								for(int n=0;n<current_map.length;n++)
								{
									current_map.matrix[i][n]=it->matrix[i][n];
								}
								current_map.number=it->number;
								goto out;
						}
						else
							cout<<"Invalid map:1.Dismatch in length/ 2.Number is invalid."<<endl;
					}
				}
			}
			out:cout<<"current map is No."<<current_map.number<<endl;

}
void PBCA::Read_Map()
{
	int number=0;
	ifstream ifs;
	ifs.open("../map.txt");
	
	if(!ifs)
		cout<<"Fail to open map file"<<endl;
	else
		cout<<"Success to open map file"<<endl;
	while(!ifs.eof())                   //each time read a block of date
	{  
		map M;
		M.number=number;
		number++;
		ifs>>M.length;
		for(int i=0;i<M.length;i++)
			for(int j=0;j<M.length;j++)
			{
				ifs>>M.matrix[i][j];
			}
			map_vector.push_back(M);
	}
	map_vector.pop_back();//Because read block, the last block will be counted for twice 
}
void PBCA::Test_showMapVector()
{
	list<map>::iterator it;
	for(it=map_vector.begin();it!=map_vector.end();it++)
		{
			for(int i=0;i<it->length;i++)
			{
				for(int j=0;j<it->length;j++)
				{
					cout<<it->matrix[i][j]<<" ";
				}
				cout<<endl;
			}
			cout<<endl;
	}
}
void PBCA::HGT()//function are stored in HGT.h
{
	Sample();
	Distribution();

}

float PBCA::Survive_Rate(sequence p,sequence b)
{
	float S,prate;
	int maxf,de;

	if(p.fitness>b.fitness)
	{
		maxf=p.fitness;
	}
	else
	{
		maxf=b.fitness;
		
	}
	de=b.fitness-p.fitness;
	prate=phage.size()/(phage.size()+bacteria.size());
	S=0.5+de/maxf+(1.0-prate)/4;
	return S;

}
void PBCA::fight()
{
	All_fitness(bacteria);
	//bacteria.sort(CompareRule);
	list<sequence>::iterator it_b,it_p;
	int p_length,b_length;
	p_length=phage.size()*0.5;
	b_length=bacteria.size()*0.3;
	if(phage.size()>tsp_in.pz&&p_length<tsp_in.pz)
		p_length=tsp_in.pz;
	if(bacteria.size()>tsp_in.pz&&b_length<tsp_in.pz)
		b_length=tsp_in.pz;
	int i,j;
	for(i=0,it_p=phage.begin();i<p_length;it_p++,i++)//this is the phages which take part in the competition
	{
		/*if(i%10==0)
			cout<<i<<endl;*/
		for(j=0,it_b=bacteria.begin();j<b_length;it_b++,j++) //this is the bacteria which are challenged by phages
		{
			if(CompareRule_min(*it_p,*it_b))
			{
				float dicer=(float)rand()/(float)RAND_MAX;
				float S=Survive_Rate(*it_p,*it_b);
				if(dicer<S*/*0.2*/0.3)//bacteria win
				{
					//replace(*it_b,*it_p);
					replace_EPR(*it_b,*it_p);
					replace_count++;
				
				}
				else//phage win
				{
					//mutate_multi(*it_p);//multi swap
					//mutate(*it_p);
					mutate_standard(*it_p);
					mutate_count++;
				}
			}
		}

	}
	All_fitness(bacteria);
	All_fitness(phage);


}
void PBCA::replace(sequence b,sequence p)  //more item could be added here
{
	int replace_length=tsp_in.L/15+1;
	list<int> temp_list_b,temp_list_p;
	for(int i=0;i<replace_length;i++)
	{
		temp_list_b.push_front(b.seq[i]);//preserve the numbers will be insert back to the bacteria sequence
		temp_list_p.push_front(p.seq[i]);//preserve the numbers transmitted into bacteria
	}
	for(int i=0;i<replace_length;i++)
	{
		b.seq[i]=p.seq[i];//copy the segment
		list<int>::iterator it=find(temp_list_b.begin(),temp_list_b.end(),p.seq[i]);
		if(it!=temp_list_b.end())
			temp_list_b.erase(it);//delete the overlap ones
	}
	/*reorgnize the sequence*/
	for(int i=replace_length;i<tsp_in.L;i++)
	{
		/*if(b.seq[i]in the temp_list_p )
			replace this and pop list_b*/
		if(find(temp_list_p.begin(),temp_list_p.end(),b.seq[i])!=temp_list_p.end())
		{
			b.seq[i]=temp_list_b.front();
			temp_list_b.pop_front();
		}
	}
	//if(!Is_SequenceInList(b,bacteria))
		bacteria.push_back(b);

}
void PBCA::replace_EPR(sequence b,sequence p)
{
  sequence sis;
  int sis_iter=0;
  sis.length=b.length;
  sis.fitness=0;
  int i,j,k,t1,t2,town;
  int ntowns=b.length;
  static char CM[MaxLength][MaxLength],visit[MaxLength];
  memset(CM, 0, MaxLength*MaxLength*sizeof(char));
  memset(visit, 0, MaxLength*sizeof(char));
  // create connection matrix

  for(j=0; j<ntowns; j++) {
	  t1 = b.seq[j]; t2 = b.seq[(j+1)%ntowns];
    CM[t1][t2]=1; CM[t2][t1]=1;
  }
 // mate2.head();
  for(j=0; j<ntowns; j++) {
	  t1 = p.seq[j]; t2 = p.seq[(j+1)%ntowns];
    CM[t1][t2]=1; CM[t2][t1]=1;
  }
  // select 1st town randomly
  town=GARandomInt(0,ntowns-1);
  visit[town]=1; memset(CM[town], 0, MaxLength*sizeof(char));
  //sis.insert(town);
  sis.seq[sis_iter]=town+1; // the head node 
  sis_iter++;
  
  GAList<int> PossFollowList;
  GAList<int> FollowersList[5];
  while (PossFollowList.head()) PossFollowList.destroy();
  for(k=0; k<5; k++) {
    while (FollowersList[k].head()) FollowersList[k].destroy(); 
  }
  
  // select the following town with the minimal no of next folling towns
  int nPoss,nFollow;
  for(i=1; i<ntowns; i++) {           
    nPoss = 0;
    for(j=0; j<ntowns; j++) {          // no of poss. following towns
      if (CM[j][town]) {
	nPoss += 1;
	PossFollowList.insert(j);}
    }
    // nPoss = 0;
    if (nPoss == 0) {
      do {town=GARandomInt(0,ntowns-1);} while (visit[town]); // no follower
      visit[town]=1; memset(CM[town], 0, MaxLength*sizeof(char));
	  //sis insert
	  sis.seq[sis_iter]=town+1; 
      sis_iter++;
    }
    else {
      PossFollowList.head();
      for(j=0; j<nPoss; j++) {
	nFollow = 0; 
	town = (*PossFollowList.current());
	for(k=0; k<ntowns; k++) {
	  if (CM[k][town]) nFollow++; 
	}
	FollowersList[nFollow].insert(town);
	PossFollowList.next();
      }
      k=0;
      while (FollowersList[k].size() == 0) k++;
      FollowersList[k].warp(GARandomInt(0,FollowersList[k].size()));
      town = (*FollowersList[k].current());
      visit[town]=1; memset(CM[town], 0, MaxLength*sizeof(char));
      //sis.insert(town);
	  sis.seq[sis_iter]=town+1; 
      sis_iter++;
    }
    while (PossFollowList.head()) PossFollowList.destroy();
    for(k=0; k<5; k++) {
      while (FollowersList[k].head()) FollowersList[k].destroy(); 
    }
  }
  bacteria.push_back(sis);
  return ;
  //sis.head(); 
}
bool PBCA::Is_SequenceInList(sequence S,list<sequence> L)
{
	if(find(L.begin(),L.end(),S)!=L.end())
		return true;
	else
		return false;

}
void PBCA::population_regulation()
{
	bacteria.sort(CompareRule_min);
	phage.sort(CompareRule_min);
	int population_now=bacteria.size()+phage.size();
	while(population_now>tsp_in.pmax)
	{
		
		if(bacteria.size()>tsp_in.pz&&phage.size()>tsp_in.pz)
		{
			if(CompareRule_min(bacteria.back(),phage.back()))
				phage.pop_back();
			else
 				bacteria.pop_back();
		}
		else if(bacteria.size()>tsp_in.pz&&phage.size()<=tsp_in.pz)
			bacteria.pop_back();
		else if(bacteria.size()<=tsp_in.pz&&phage.size()>tsp_in.pz)
			phage.pop_back();
		population_now=bacteria.size()+phage.size();
	}
	
}
void PBCA::Test_Run()
{	tsp_in.L=10;
	tsp_in.pmax=100;
	tsp_in.pz=10;
	tsp_in.round=10;
	tsp_in.start_point=1;
	tsp_in.times=1;

	sequence s1;
	s1.length=10;
	s1.seq[0]=1;
	s1.seq[1]=2;
	s1.seq[2]=3;
	s1.seq[3]=4;
	s1.seq[4]=5;
	s1.seq[5]=6;
	s1.seq[6]=7;
	s1.seq[7]=8;
	s1.seq[8]=9;
	s1.seq[9]=10;

	sequence s2;
	s2.length=10;	
	s2.seq[0]=1;
	s2.seq[1]=5;
	s2.seq[2]=3;
	s2.seq[3]=4;
	s2.seq[4]=2;
	s2.seq[5]=6;
	s2.seq[6]=7;
	s2.seq[7]=8;
	s2.seq[8]=10;
	s2.seq[9]=9;

	sequence s3;
	s3.length=10;	
	s3.seq[0]=10;
	s3.seq[1]=9;
	s3.seq[2]=8;
	s3.seq[3]=7;
	s3.seq[4]=6;
	s3.seq[5]=5;
	s3.seq[6]=4;
	s3.seq[7]=3;
	s3.seq[8]=2;
	s3.seq[9]=1;
	bacteria.push_back(s1);
	bacteria.push_back(s2);
	bacteria.push_back(s3);
	Sample();
	Hlist.pop_front();
	Hlist.front().value=1;
	Distribution();
	
}
int PBCA::Output()
{
	//All_fitness(bacteria);
	int record=0;
	bacteria.sort(CompareRule_min);
	phage.sort(CompareRule_min);
	list<sequence>::iterator it;
	cout<<"The best solution is:"<<endl;
	if(CompareRule_min(bacteria.front(),phage.front()))
	{
		it=bacteria.begin();
		record=it->fitness;
		for(int i=0;i<5;i++,it++)
		{
			cout<<"Fitness="<<it->fitness<<endl;
		}
	}
	else
	{
		it=phage.begin();
		record=it->fitness;
		for(int i=0;i<5;i++,it++)
		{
			cout<<"Fitness="<<it->fitness<<endl;
		}
	}
	return record;
}
void PBCA::Make_Map_EUD()
{
	ifstream ifs;
	int length;
	int flag=0;//sign the reading process has gone to date section
	list<EUC_2D_Node> elist;
	//int count=0;
	ifs.open("../raw_map.txt");
	if(!ifs)
		cout<<"fail to open raw_map.txt"<<endl;
	while(!ifs.eof())
	{
		string temp;
		int i=0;
		if(flag==0)
			ifs>>temp; //stop reading here in the date section
		if(temp=="DIMENSION")
		{
				cout<<"In the line of DIMENSION"<<endl;
				ifs>>temp;
				if(temp!=":")
					cout<<"dismatch with ':'"<<endl;
				ifs>>length;
				temp.clear();
		}
		else if(temp=="EDGE_WEIGHT_TYPE")//expand map_type here
		{
			ifs>>temp;
			if(temp!=":")
					cout<<"dismatch with ':'"<<endl;
			ifs>>temp;
			if(temp=="EUC_2D") 
				map_type=EUC2D;
			else 
				map_type=UNKNOWN;
			temp.clear();
		}
		else if(temp=="NODE_COORD_SECTION")
			flag=1;

		if(flag==1)
		{
			if(map_type==EUC2D)
			{	
				//euc2d raw map make
				EUC_2D_Node temp_node;
				ifs>>temp_node.node>>temp_node.x>>temp_node.y;
				elist.push_back(temp_node);
				/*cout<<node<<" "<<x<<" "<<y<<"   count="<<count<<endl;
				count++;*/
			
			}
			else if(map_type=UNKNOWN)
			{
				cout<<"Unknown raw map type"<<endl;
				getchar();
				exit(0);
			}
		}//if-flag end
	}//read end-while
	ifs.close();
	//add making process here
	if(map_type==EUC2D)
		make_e2_map(elist,length);
	else if(map_type=UNKNOWN)
		cout<<"error"<<endl;//in fact it would not occur.Use it as a landmark for coding
}
void PBCA::make_e2_map(list<EUC_2D_Node> elist,int length)//make map for e2 map
{
	ofstream ofs;
	ofs.open("../map.txt");
	ofs<<length<<endl;
	list<EUC_2D_Node>::iterator self,other;
	for(self=elist.begin();self!=elist.end();self++)
		{
			for(other=elist.begin();other!=elist.end();other++)
			{
				int xd=self->x-other->x;
				int yd=self->y-other->y;
				int distance=sqrt(float(xd*xd+yd*yd));
				ofs<<distance<<" ";
			}
			ofs<<endl;
	}
	ofs.close();

}
void PBCA::Load_To_Log(ofstream &ofs)
{
	ofs<<"Dimension:"<<tsp_in.L<<" Total round:"<<tsp_in.round<<" population:"<<tsp_in.pmax<<" safe zone:"<<tsp_in.pz<<endl;
}

void strreverse(char* begin, char* end) {
	
	char aux;
	
	while(end>begin)
	
		aux=*end, *end--=*begin, *begin++=aux;
	
}	
void itoa(int value, char* str, int base) {
	
	static char num[] = "0123456789abcdefghijklmnopqrstuvwxyz";
	
	char* wstr=str;
	
	int sign;
	

	
	// Validate base
	
	if (base<2 || base>35){ *wstr='\0'; return; }
	

	
	// Take care of sign
	
	if ((sign=value) < 0) value = -value;
	

	
	// Conversion. Number is reversed.
	
	do *wstr++ = num[value%base]; while(value/=base);
	
	if(sign<0) *wstr++='-';
	
	*wstr='\0';
	

	
	// Reverse string
	
	strreverse(str,wstr-1);
	
}
void Log_Move()
{
	char log_path[]="../testing log/version2.3a_log1.txt";//cut from 16
	int i=1;
	while(access(log_path,F_OK)==0)//have confilct will never stop
	{
		i++;
		char num_array[4];
		itoa(i,num_array,10);
		if(i<10)
		{
			log_path[30]=num_array[0];
		}
		else if(i<100)
		{
			
			log_path[29]=num_array[0];
			log_path[30]=num_array[1];
			log_path[31]='.';
			log_path[32]='t';
			log_path[33]='x';
			log_path[34]='t';
			log_path[35]='\0';
		}
		else
		cout<<"Too many file."<<endl;\
	}

	rename(LOG,log_path);
}
void PBCA::Run()
{
	/*reading*/
	ifstream ifs;
	ofstream ofs;
	ifs.open(PATH);
	ofs.open(LOG);
	Read_Map();
	HGT_Count=0;
	if(!ifs)
		cout<<"Fail to open input file"<<endl;
	else 
		cout<<"sucess to open input file"<<endl;
	input(ifs);	
	Select_Map();
	for(int i=0;;i++)
	{
		PBCA_Run_count=1;
		if(valid_time<1)                            //if the last input parameter is becoming invalid
			if(!input(ifs)&&valid_time==0)          //read in this step.only when ifs reach the end of file and the last set of parameter become invalid then jump out
					break;
	/*MainBody*/
		//Select_Map();	 each round could change the map move it out of loop,1map for all parameter
		if(!ofs)
			cout<<"fail to creat log.txt"<<endl;
		else
			Load_To_Log(ofs);
		random_produce(tsp_in.pz/2,bacteria);
		random_produce(tsp_in.pz/2,phage);
		time_t time_begin,time_end;
		time_begin=time(NULL);
		replace_count=0;
		mutate_count=0;
		for(int j=0;j<tsp_in.round;j++)
		{
			/*Test_showparameters();*/
			Mutate_Selection();
			fight();
			population_regulation();
			PBCA_Run_count++;
			//cout<<"round="<<count<<" "<<"fitness="<<Output()<<endl;
			if(PBCA_Run_count%50==0)
				{
					cout<<"round:"<<PBCA_Run_count<<endl;
					
					ofs<<"------"<<PBCA_Run_count<<"round-----";
					ofs<<"Fitness="<<Output()<<endl;
					time_end=time(NULL);
					cout<<"Time="<<time_end-time_begin<<endl;
					ofs<<"Time="<<time_end-time_begin<<endl;
			}
			
			if(PBCA_Run_count%5==0)
			{
				//HGT();
				HGT_Crossover();
				cout<<PBCA_Run_count<<endl;
			}
			//Test_showSequenceArray();
		}
		time_end=time(NULL);
		ofs<<"Time="<<time_end-time_begin<<endl;
		ofs<<"Replace_count:"<<replace_count<<"Mutate_count:"<<mutate_count<<"HGT_Count"<<HGT_Count<<endl;
		Output();
	/*Rest area*/

		valid_time--;
		bacteria.clear();
		phage.clear();
		
		
		cout<<i<<endl;
	} 

	ifs.close();
    ofs.close();
	Log_Move();
	getchar();

}
