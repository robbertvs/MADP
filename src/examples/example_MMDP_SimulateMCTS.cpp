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
struct action_node;
struct state_node
{
	int state;
	int iterations;
	double average;
	bool isWinning;
	bool isLosing;
	set<action_node*> children;

	bool operator==(const state_node &other) const {
		return state == other.state;
	}

	string toString() {
		ostringstream os;
		os << "Node is in state " << state << " and has average " << average << " from " << iterations << " iterations.";
		return os.str();
	}

	state_node(int s) {
		state = s;
		iterations = 0;
		average = 0.0;
		isWinning = false;
		isLosing = false;
	}
};

struct action_node {
	state_node* state;
	int action;
	int iterations;
	double average;
	set<state_node*> children;

	bool operator==(const action_node &other) const {
		return state == other.state && action == other.action;
	}
	action_node() {}
	action_node(state_node* s, int a) {
		state = s;
		action = a;
		iterations = 0;
		average = 0.0;
	}
};

state_node* getNode(map<int, state_node*> *states, int state) {
	map<int, state_node*>::iterator tryState = states->find(state);
	if (tryState == states->end()) {
		//cout << "Making new node with state " << state << endl;
		state_node* node = new state_node(state);
		(*states)[state] = node;
		return node;
	}
	else {
		return tryState->second;
	}
}

void addWinningAndLosingStates(DecPOMDPDiscreteInterface* decpomdp, map<int, state_node*> *states) {
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
	state_node* winningNode = getNode(states, winningState);
	winningNode->isWinning = true;

	state_node* losingNode = getNode(states, losingState);
	losingNode->isLosing = true;
}

double simulation(DecPOMDPDiscreteInterface* decpomdp, map<int, state_node*> *states, state_node* startNode, int horizon) {
	int currentState = startNode->state;
	int nrActions = decpomdp->GetNrJointActions();

	for (int i = 0; i < horizon; i++) {
		int action = rand() % nrActions;
		currentState = decpomdp->SampleSuccessorState(currentState, action);
		state_node* stateNode = getNode(states, currentState);
		if (stateNode->isWinning) return 1.0;
		if (stateNode->isLosing) return 0.0;
	}

	return 0.0; // This also is a loss
}

// Select an action using the UCT bestchild algorithm
action_node* select_action(state_node* currentNode, DecPOMDPDiscreteInterface* decpomdp)
{
	double exploration = sqrt(2);
	int nrActions = decpomdp->GetNrJointActions();
	double maxUct = 0;
	action_node* selectedAction;
	set<action_node*> *actions = &(currentNode->children);

	for (int action = 0; action < nrActions; action++) {
		double uctValue;
		action_node* dummy = new action_node(currentNode, action);
		action_node* actionNode;
		set<action_node*>::iterator result = actions->find(dummy);
		if (result != actions->end() && (*result)->iterations > 0) {
			actionNode = (*result);
			uctValue = actionNode->average + 2 * exploration * sqrt(2 * log(currentNode->iterations) / actionNode->iterations);
		}
		else {
			uctValue = 10000 + rand() % 1000;
			actionNode = dummy;
			actions->insert(actionNode);
		}
		
		if (uctValue > maxUct)
		{
			selectedAction = actionNode;
			maxUct = uctValue;
		}
	}

	cout << "Choosing action " << selectedAction->action << " from state " << currentNode->state << endl;
	return selectedAction;
}

double MCTS(DecPOMDPDiscreteInterface* decpomdp, map<int, state_node*> *states, state_node* currentNode, int horizon)
{
	if (horizon <= 0) {
		return 0.0;
	}
	else {
		action_node* action = select_action(currentNode, decpomdp);
		int nextState = decpomdp->SampleSuccessorState(currentNode->state, action->action);
		cout << "Going from state " << currentNode->state << " to state " << nextState << " with action " << action->action << endl;
		double reward;
		
		state_node* stateNode = getNode(states, nextState);
		if (stateNode->isWinning) {
			reward = 1.0;
		}
		else if (stateNode->isLosing) {
			reward = 0.0;
		}
		else {
			set<state_node*>::iterator result = action->children.find(stateNode);
			if (result != action->children.end()) { // Gotten here before, go deeper
				reward = MCTS(decpomdp, states, stateNode, horizon - 1);
			}
			else { // Never gotten here before
				reward = simulation(decpomdp, states, stateNode, horizon - 1);
				cout << "Simulation for node " << nextState << " gives reward " << reward << endl;;
				//(&(action->children))->insert(stateNode);
			}
		}

		currentNode->average = (currentNode->iterations * currentNode->average + reward) / (currentNode->iterations + 1);
		currentNode->iterations++;
		action->average = (action->iterations * action->average + reward) / (action->iterations + 1);
		action->iterations++;
		
		return reward;
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

//void printTree(state_node* node, int depth, int maxDepth) {
//	cout << string(depth, '\t') << node->toString() << endl;
//	if (maxDepth == 0 || depth < maxDepth) {
//		for (map<int, list<state_node*> >::iterator ii = node->children.begin(); ii != node->children.end(); ++ii)
//		{
//			cout << string(depth, '\t') << ii->first << ": " << endl;
//
//			for (list<state_node*>::iterator jj = ii->second.begin(); jj != ii->second.end(); ++jj)
//			{
//				printTree(*jj, depth + 1, maxDepth);
//			}
//		}
//	}
//}

int main(int argc, char **argv)
{
    ArgumentHandlers::Arguments args;
    argp_parse (&ArgumentHandlers::theArgpStruc, argc, argv, 0, 0, &args);

    try
    {
		cout << "Instantiating the problem..."<<endl;
		DecPOMDPDiscreteInterface* decpomdp = GetDecPOMDPDiscreteInterfaceFromArgs(args);
		cout << "...done."<<endl;

		srand(42);

		int nrStates = decpomdp->GetNrStates();
		int nrActions = decpomdp->GetNrJointActions();

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
		map<int, state_node*> states;
		state_node* root = getNode(&states, initialState);
		addWinningAndLosingStates(decpomdp, &states);

		// MCTS
		for (int i = 0; i < 10; i++) {
			cout << "-----------------------" << endl;
			cout << "Simulation " << i << endl;
			double reward = MCTS(decpomdp, &states, root, args.horizon);
			if (reward > maxReward) maxReward = reward;
		}
		//printTree(root, 0, 2);
		for (int i = 0; i < nrStates; i++) {
			cout << states[i]->toString() << endl;
		}
		
		cout << "Maximum MCTS search reward found is " << maxReward << endl;

		//int maxIterations = 0;
		//int maxAction = 0;
		//for (map<int, list<state_node*> >::iterator ii = root->children.begin(); ii != root->children.end(); ++ii)
		//{
		//	int iterations = 0;
		//	for (list<state_node*>::iterator jj = ii->second.begin(); jj != ii->second.end(); ++jj)
		//	{
		//		iterations += (*jj)->iterations;
		//	}

		//	if (iterations > maxIterations) {
		//		maxIterations = iterations;
		//		maxAction = ii->first;
		//	}
		//	cout << ii->first << ": " << iterations << endl;
		//}
		//cout << "Best action is " << maxAction << endl;
	}
	catch(E& e){ e.Print(); }

	return(0);
}
