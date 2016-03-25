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
	double min;
	double max;
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

double simulation(DecPOMDPDiscreteInterface* decpomdp, int winningState, int losingState, tree_node* startNode, int horizon) {
	int currentState = startNode->state;
	int nrActions = decpomdp->GetNrJointActions();

	for (int i = 0; i < horizon; i++) {
		int action = rand() % nrActions;
		currentState = decpomdp->SampleSuccessorState(currentState, action);
		if (currentState == winningState) return 1.0;
		if (currentState == losingState) return 0.0;
	}

	return 0.0;
}

// Select an action using the UCT bestchild algorithm
int select_action(tree_node* currentNode, DecPOMDPDiscreteInterface* decpomdp)
{
	double exploration = sqrt(2);
	int nrActions = decpomdp->GetNrJointActions();
	double maxUct = 0;
	int selectedAction = 0;

	for(int action = 0; action<nrActions; action++)
	{
		int iterations = 0;
		double sumValue = 0.0;
		for (list<tree_node>::iterator child=currentNode->children[action].begin(); child != currentNode->children[action].end(); ++child)
		{
			iterations += child->iterations;
			sumValue += child->iterations * child->average;
		}
		double averageValue = sumValue / iterations;
		double uctValue = averageValue + 2*exploration*sqrt(2*log(currentNode->iterations)/iterations);
		if(iterations==0)
		{
			uctValue = 10000+rand()%1000; // explore new actions first.
		}
		if(uctValue > maxUct)
		{
			selectedAction = action;
			maxUct = uctValue;
		}
	}
	return selectedAction;
}

double MCTS(DecPOMDPDiscreteInterface* decpomdp, int winningState, int losingState, tree_node* currentNode, int horizon)
{
	if (horizon <= 0) {
		return 0;
	}
	else {
		int nrActions = decpomdp->GetNrJointActions();
		int action = select_action(currentNode, decpomdp);
		int nextState = decpomdp->SampleSuccessorState(currentNode->state, action);
		double reward;
		if (nextState == winningState) {
			reward = 1.0;
		}
		else if (nextState == losingState) {
			reward = 0.0;
		}
		else {
			list<tree_node>* stateList = &(currentNode->children[action]); // Will instantiate new object if not done this action before
			list<tree_node>::iterator result = find_if(stateList->begin(), stateList->end(), find_by_state(nextState));
			if (result != stateList->end()) { // Gotten here before, go deeper
				tree_node* state = &(*result); // Magic to get rid of the iterator
				reward = MCTS(decpomdp, winningState, losingState, state, horizon - 1);
			}
			else { // Never gotten here before
				tree_node state = newNode(nextState);
				reward = simulation(decpomdp, winningState, losingState, &state, horizon - 1);
				stateList->push_back(state);
			}
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

void printTree(tree_node* node, int depth, int maxDepth) {
	cout << string(depth, '\t') << node->toString() << endl;
	if (maxDepth == 0 || depth < maxDepth) {
		for (map<int, list<tree_node> >::iterator ii = node->children.begin(); ii != node->children.end(); ++ii)
		{
			cout << string(depth, '\t') << ii->first << ": " << endl;

			for (list<tree_node>::iterator jj = ii->second.begin(); jj != ii->second.end(); ++jj)
			{
				printTree(&*jj, depth + 1, maxDepth);
			}
		}
	}
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

		int nrStates = decpomdp->GetNrStates();
		int nrActions = decpomdp->GetNrJointActions();

		double maxStateReward = -std::numeric_limits<double>::max();
		double minStateReward = std::numeric_limits<double>::max();
		int winningState = 0;
		int losingState = 0;
		for (int i = 0; i < nrStates; i++) {
			double reward = 0;
			for (int j = 0; j < nrActions; j++) {
				reward += decpomdp->GetReward(i, j);
			}
			reward = reward / nrStates;
			if (reward > maxStateReward) {
				maxStateReward = reward;
				winningState = i;
			}
			if (reward < minStateReward) {
				minStateReward = reward;
				losingState = i;
			}
		}
		cout << "State " << winningState << " is the winning state and " << losingState << " is the losing state" << endl;

		size_t initialState = decpomdp->SampleInitialState();

		// Pure BFS. Takes a long time, answer is only marginally better than pure random.
		//long steps = pow(nrActions, args.horizon);
		//cout << "Maximum BFS reward found (" << steps << " steps) is " << BFS(decpomdp, steps, initialState) << endl;

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
		for (int i = 0; i < 10000; i++) {
			double reward = MCTS(decpomdp, winningState, losingState, root, args.horizon);
			if (reward > maxReward) maxReward = reward;
		}
		printTree(root, 0, 2);

		cout << "Maximum MCTS search reward found is " << maxReward << endl;

		int maxIterations = 0;
		int maxAction = 0;
		for(int action = 0; action<nrActions; action++)
		{
			int iterations = 0;
			for (list<tree_node>::iterator child=root->children[action].begin(); child != root->children[action].end(); ++child)
			{
				iterations += child->iterations;
			}

			if (iterations > maxIterations) {
				maxIterations = iterations;
				maxAction = action;
			}
			cout << action << ": " << iterations << endl;
		}
		cout << "Best action is " << maxAction << endl;
	}
	catch(E& e){ e.Print(); }

	return(0);
}
