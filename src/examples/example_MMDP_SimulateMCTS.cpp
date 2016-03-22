/* This file is part of the Multiagent Decision Process (MADP) Toolbox v0.3.
*
* The majority of MADP is free software released under GNUP GPL v.3. However,
* some of the included libraries are released under a different license. For
* more information, see the included COPYING file. For other information,
* please refer to the included README file.
*
* This file has been written and/or modified by the following people:
*
* Frans Oliehoek
* Matthijs Spaan
*
* For contact information please see the included AUTHORS file.
*/

#include <iostream>

#include "DecPOMDPDiscrete.h"
#include "TreeNode.h"

#include "argumentHandlers.h"
#include "argumentUtils.h"
#include <list>
#include <algorithm>
using namespace std;
using namespace ArgumentUtils;

const char *argp_program_version = "example_MMDP_SolveAndSimulate";

// Program documentation
static char doc[] =
"example_MMDP_SolveAndSimulate - loads an (MMDP) problem and simulates Monte Carlo Tree Search (MCTS). \
\vFor more information please consult the MADP documentation.";

//NOTE: make sure that the below value (nrChildParsers) is correct!
const int nrChildParsers = 5;
const struct argp_child childVector[] = {
    ArgumentHandlers::problemFile_child,
    ArgumentHandlers::globalOptions_child,
    ArgumentHandlers::modelOptions_child,
    ArgumentHandlers::solutionMethodOptions_child,
    ArgumentHandlers::simulation_child,
    { 0 }
};
#include "argumentHandlersPostChild.h"

struct tree_node
{
	int state;
	int iterations;
	double average;
	map<int, list<tree_node> > children; //Map of action to list of children
	//bool operator()(const tree_node & node) {
	//	return node.state == state && node.iterations == iterations && node.average == average;
	//}
	string toString() {
		ostringstream os;
		os << "Node is in state " << state << " and has average " << average << " from " << iterations << " iterations.";
		return os.str();
	}
};

tree_node newNode(int state)
{
	tree_node node;
	node.state = state;
	node.iterations = 0;
	node.average = 0.0;
	return node;
}

struct find_by_state {
	find_by_state(const int state) : state(state) {}
	bool operator()(const tree_node & node) {
		return node.state == state;
	}
private:
	int state;
};

double simulation(DecPOMDPDiscreteInterface* decpomdp, tree_node* startNode, int horizon) {
	int currentState = startNode->state;
	int nrActions = decpomdp->GetNrJointActions();
	double discount = decpomdp->GetDiscount();
	double reward = 0;

	for (int i = 0; i < horizon; i++) {
		int action = rand() % nrActions;
		currentState = decpomdp->SampleSuccessorState(currentState, action);
		reward += pow(discount, i) * decpomdp->GetReward(currentState, action);
	}

	return reward;
}

double MCTS(DecPOMDPDiscreteInterface* decpomdp, tree_node* currentNode, int horizon)
{
	if (horizon <= 0) {
		return 0;
	}
	else {
		double discount = decpomdp->GetDiscount();
		int nrActions = decpomdp->GetNrJointActions();
		int action = rand() % nrActions;
		int nextState = decpomdp->SampleSuccessorState(currentNode->state, action);
		double actionReward = decpomdp->GetReward(currentNode->state, action);
		double reward;
		list<tree_node>* stateList = &(currentNode->children[action]); // Will instantiate new object if not done this action before
		list<tree_node>::iterator result = find_if(stateList->begin(), stateList->end(), find_by_state(nextState));
		if (result != stateList->end()) { // Gotten here before, go deeper
			tree_node* state = &(*result); // Magic to get rid of the iterator
			double deepReward = MCTS(decpomdp, state, horizon - 1) * discount;
			reward = actionReward + deepReward;
		}
		else { // Never gotten here before
			tree_node state = newNode(nextState);
			double simulationReward = simulation(decpomdp, &state, horizon - 1) * discount;
			reward = actionReward + simulationReward;
			stateList->push_back(state);
		}

		currentNode->average = (currentNode->iterations * currentNode->average + reward) / (currentNode->iterations + 1);
		currentNode->iterations++;
		return currentNode->average;
	}
}

//BFS currently does not take into account the stochastic nature of the environment.
double BFS(DecPOMDPDiscreteInterface* decpomdp, long steps, int state) {
    if (steps <= 0)
    return 0;

    int nrActions = decpomdp->GetNrJointActions();
    long stepsToDivide = floor(double(steps) / nrActions);
    double discount = decpomdp->GetDiscount();
    double maxReward = -std::numeric_limits<double>::max();

    for (int i = 0; i < nrActions; i++) {
        int action = i;
        int newState = decpomdp->SampleSuccessorState(state, action);
        double actionReward = decpomdp->GetReward(state, action);
        double bfsReward = BFS(decpomdp, stepsToDivide, newState);
        double reward = actionReward + bfsReward*discount;
        if (steps == 10000) {
            cout << "Choosing action " << action << " in state " << state << " gives " << actionReward << " immediate pay-off. Combined with the rest, the pay-off becomes " << reward
            << endl;
        }
        if (reward > maxReward) maxReward = reward;
    }

    return maxReward;
}

int main(int argc, char **argv)
{
    ArgumentHandlers::Arguments args;
    argp_parse (&ArgumentHandlers::theArgpStruc, argc, argv, 0, 0, &args);

    try
    {
		cout << "Instantiating the problem..."<<endl;
		DecPOMDPDiscreteInterface* decpomdp = GetDecPOMDPDiscreteInterfaceFromArgs(args);
		cout << "...done."<<endl;

		//srand(time(NULL));
		srand(42);

		//Build transition tree (I don't think the domains are large enough that this won't work)
		int nrStates = decpomdp->GetNrStates();
		int nrActions = decpomdp->GetNrJointActions();

		size_t initialState = decpomdp->SampleInitialState();

		// Pure BFS. Takes a long time, answer is only marginally better than pure random.
		long steps = pow(nrActions, args.horizon);
		cout << "Maximum BFS reward found (" << steps << " steps) is " << BFS(decpomdp, steps, initialState) << endl;

		// Pure random
		double discount = decpomdp->GetDiscount();
		double maxReward = -std::numeric_limits<double>::max();

		for (int i = 0; i < 10000; i++) {
			size_t currentState = initialState;
			double sumReward = 0;
			for (int j = 0; j < args.horizon; j++) {
				int action = rand() % nrActions;
				currentState = decpomdp->SampleSuccessorState(currentState, action);
				sumReward += decpomdp->GetReward(currentState, action)*pow(discount, j);
			}
			if (sumReward > maxReward) maxReward = sumReward;
			//cout << "Run " << i << " with reward " << sumReward << endl;
		}
		cout << "Maximum random search reward found is " << maxReward << endl;
		
		maxReward = -std::numeric_limits<double>::max();
		tree_node r = newNode(initialState);
		tree_node* root = &r;
		// MCTS
		for (int i = 0; i < 100; i++) {
			double reward = MCTS(decpomdp, root, args.horizon);
			if (reward > maxReward) maxReward = reward;
		}
		cout << "Root: " << root->toString() << endl;
		for (map<int, list<tree_node> >::iterator ii = root->children.begin(); ii != root->children.begin(); ++ii)
		{
			cout << (*ii).first << ": " << endl;

			for (list<tree_node>::iterator jj = (*ii).second.begin(); jj != (*ii).second.end(); ++jj)
			{
				cout << "\t" << root->toString() << endl;				
			}
		}

		cout << "Maximum MCTS search reward found is " << maxReward << endl;
	}
	catch(E& e){ e.Print(); }

	return(0);
}