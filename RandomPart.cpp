#include <bitset>
#include <map>
#include <random>

using namespace std;

#include "data.h"

bool check_partition(map<unsigned int, uint64_t> Partition)
{
  map<unsigned int, uint64_t>::iterator Part;
  uint64_t sum = 0;
  unsigned int rank = 0; 

  for (Part = Partition.begin(); Part != Partition.end(); Part++)
  {
    sum |= (*Part).second;
    rank += bitset<n>((*Part).second).count();
    //cout << bitset<n>( (*Part).second ) << " \t";
  }
  //cout << bitset<n>(sum) << endl;

  return (bitset<n>(sum).count() == rank);
}


/*map<unsigned int, uint64_t> Random_Partition(unsigned int A, mt19937 gen) // A = number of parts
{
	uniform_int_distribution<int> uni(0.,A-1);
	map<unsigned int, uint64_t> Partition_rand;
	uint64_t si = 0;
	unsigned int part = 0;

	for(unsigned int a=0; a<A; a++)
	{
		Partition_rand[a]=si;
	}

	si = 1;

	for(int i=0; i<n; i++)
	{
		part = uni(gen);
		Partition_rand[part] += si;
		si = (si << 1);
	}

	for(map<unsigned int, uint64_t>::iterator it = Partition_rand.begin(); it != Partition_rand.end(); it++)
	{
		cout << (*it).second << "\t " << bitset<n>((*it).second) << endl;
	}

	cout << "Is this partition correct?  " << ((check_partition(Partition_rand))?"Yes":"No") << endl;

	
	return Partition_rand;
}*/
