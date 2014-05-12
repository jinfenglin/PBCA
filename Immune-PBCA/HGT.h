#include "PBCA.h"
int find_in_array(int *a,int number,int start,int length);
bool sortrule(sequence s1,sequence s2)
{
	
	if(s1.fitness>s2.fitness)
		return true;
	else
		return false;
}
void PBCA::Sample()
{
	int projection[MaxLength];
	int sample_rate=3;
	All_fitness(bacteria);
	bacteria.sort(sortrule);
	list<sequence>::iterator it=bacteria.begin();
	for(int i=0;i<tsp_in.L;i++)
	{
		projection[i]=it->seq[i];
	}
	
	for(int i=0;i<sample_rate-1;i++)
	{
		it++;
		for(int j=0;j<tsp_in.L;j++)
		{
			if(projection[j]!=it->seq[j])
				projection[j]=-1;
		}	
	}
	int count=0;
	for(int i=0;i<tsp_in.L;i++)
	{
		if(projection[i]!=-1)
		{
			count++;
			if(projection[i+1]<0||i+1==tsp_in.L) //it's the end of a consequence array
			{
				float progress=PBCA_Run_count/tsp_in.round;
				int block_length=progress*tsp_in.L/2+2;
				block_length=1;
				if(count<block_length)
				{
					count=0;
					continue;
				}
				HGT_Segment hcell;
				int i_back=i;
				hcell.L=count;
				hcell.value=0;
				hcell.start_point=i-count+1;
				for(int z=count-1;z>=0;z--,i_back--)
				{
					hcell.segments[z]=projection[i_back];
				}
				Hlist.push_back(hcell);
				count=0;
			}
		}
	}
}
struct Score_Compare
{
  bool operator ()(const HGT_Segment &a, const HGT_Segment &b) {
	  return a.value>b.value;
  }
}SC;
void PBCA::Distribution()//Spread the segments to population.
{
	sort(Hlist.begin(),Hlist.end(),SC);
	int sum_score=0;
	deque<HGT_Segment>::iterator it=Hlist.begin();
	list<sequence>::iterator it_b=bacteria.begin();
	int length=Hlist.size();
	while(it!=Hlist.end())
	{
		sum_score+=it->value;
		it++;
	}
	it=Hlist.begin();
	for(int i=0;i<length;i++)
	{
		float percent;
		int valid_round;
		if(sum_score==0)
		{
			percent=0;
			valid_round=1;
		}
		else
		{
			percent=it->value/sum_score;
			valid_round=bacteria.size()*percent;
		}
		for(int j=0;j<valid_round;j++)
		{
			HGT_Replace(it,it_b);
			it_b++;
			if(it_b==bacteria.end())
				return;
		}
		it++;
	}
}
void PBCA::HGT_Replace(deque<HGT_Segment>::iterator it,list<sequence>::iterator it_b)//replace the segment in it_b with segments from it                                                                                 //and find out if it improve the performance
{
	sequence copy_b;
	for(int i=0;i<it_b->length;i++)
	{
		copy_b.seq[i]=it_b->seq[i];
	}
	copy_b.length=it_b->length;
	int temp;
	int j=0;
	for(int i=it->start_point;j<it->L;i++,j++)//replcaed one segment
	{
		//copy the replaced segment to temp 
		//replace the copy_b
		//place the supplented ones
		if(copy_b.seq[i]!=it->segments[j])
		{
			temp=copy_b.seq[i];
			copy_b.seq[i]=it->segments[j];
			int location=find_in_array(copy_b.seq,it->segments[j],it->start_point,it->L);
			if(location<0)
				cout<<"error"<<endl;
			else
			{
				copy_b.seq[location]=temp;
			}
		}	
	}
	fitness(copy_b);
	if(copy_b.fitness<it_b->fitness)
	{
		//update both deque
		it->value++;
		bacteria.push_back(copy_b);
	}
	else if(copy_b.fitness>it_b->fitness)
	{
		it->value--;
	}

}
void PBCA::HGT_Crossover()
{
	list<sequence>::iterator it=bacteria.begin();
	for(int i=0;i<bacteria.size()*0.3;i+=2)
	{
		replace_EPR(*it,*(++it));
		HGT_Count++;
		it++;
	}
	
}
int find_in_array(int *a,int number,int start,int length)//1.array 2.the number want to find 3.the start point of zone 4.length of the zone 
	                                                     //because we have double number in a single array,the one in zone is invalid
{
	
	int front_line=start,back_line=start+length-1;
	for(int i=0;a[i]>0;i++)
	{
		if(a[i]==number)
		{
			if(i>=front_line&&i<=back_line)
			{
				continue;//invalid one found
			}
			else
			{
				return i;
			}
		}
	}
	return -1;
}
