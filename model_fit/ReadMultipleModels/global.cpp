#include "global.h"
#include "memorypool.h"




/*
boost::array<double,nob> generate_dirichlet_distributed_value(boost::array<double,nob> belief_vector)
{
	boost::array<double,nob> dirichlet_values;
	boost::array<double,nob> gamma_samples;

	double sum = 0.0;
	for(unsigned int belief_values =0; belief_values < nob; ++belief_values)
	{
		boost::random::gamma_distribution<double> gamrand(belief_vector[belief_values], 1.0);
		gamma_samples[belief_values] = gamrand(RandomSingleton::Instance());
		sum += gamma_samples[belief_values];
	}

	for(unsigned int belief_values =0; belief_values < nob; ++belief_values)
	{
		dirichlet_values[belief_values] = gamma_samples[belief_values]/sum;
	}

	return dirichlet_values;
}*/


int pmax(int value)
{
	int maxp = (value < 1 ? 0: value);

	return maxp;
}

void print_matrix(Matrix const & m)
{
	for(unsigned int i = 0; i < m.size1(); ++i)
	{
		cout << i << ": ";
		for(unsigned int j = 0; j < m.size2(); ++j)
		{
			cout << m(i,j) << "\t";
		}
		cout << "\n";
	}
	cout << endl;
}



std::vector<double> greedy_probabilities(unsigned int preferred_action, unsigned int actions)
{
	if(actions > 1)
	{
		std::vector<double> greedy_probabilities(actions, eps/(actions-1));

		greedy_probabilities[preferred_action] = 1-eps;

		return greedy_probabilities;
	}
	else
	{
		return std::vector<double>(1, 1.0);
	}
}
