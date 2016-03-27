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
#include "SimulationDecPOMDPDiscrete.h"
#include "NullPlanner.h"
#include "directories.h"

#include "MDPValueIteration.h"
#include "QTable.h"
#include "AgentMDP.h"

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
struct state_node;
struct action_node;
template<class T> struct comparator
{
	bool operator()(const T* lhs, const T* rhs) const { return *lhs < *rhs; }
};

struct state_node
{
	int state;
	int iterations;
	double sum;
	bool isWinning;
	bool isLosing;
	set<action_node*, comparator<action_node> > children;

	bool operator==(const state_node &other) const {
		return state == other.state;
	}

	bool operator<(const state_node &other) const {
		return (state < other.state);
	}

	double average() {
		return (iterations == 0) ? 0 : sum / iterations;
	}

	string toString() {
		ostringstream os;
		os << "Node is in state " << state << " and has average " << average() << " from " << iterations << " iterations.";
		return os.str();
	}

	state_node(int s) {
		state = s;
		iterations = 0;
		sum = 0.0;
		isWinning = false;
		isLosing = false;
	}
};

struct action_node {
	state_node* state;
	int action;
	int iterations;
	double sum;
	set<state_node*, comparator<state_node> > children;

	bool operator==(const action_node &other) const {
		return state == other.state && action == other.action;
	}

	bool operator<(const action_node &other) const {
		if (*state == *other.state)
			return action < other.action;
		return *state < *other.state;
	}

	double average() {
		return (iterations == 0) ? 0 : sum / iterations;
	}

	action_node() {}
	action_node(state_node* s, int a) {
		state = s;
		action = a;
		iterations = 0;
		sum = 0.0;
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
	double maxStateReward = -DBL_MAX;
	double minStateReward = DBL_MAX;
	int winningState = 0;
	int losingState = 0;
	bool losingDuplicate = false;
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
		if (reward == minStateReward) {
			losingDuplicate = true;
		}
		if (reward < minStateReward) {
			minStateReward = reward;
			losingState = i;
			losingDuplicate = false;
		}
	}
	cout << "State " << winningState << " is the winning state and " << losingState << " is the losing state" << endl;
	state_node* winningNode = getNode(states, winningState);
	winningNode->isWinning = true;

	if (!losingDuplicate) {
		state_node* losingNode = getNode(states, losingState);
		losingNode->isLosing = true;
	}
}

double simulation(DecPOMDPDiscreteInterface* decpomdp, map<int, state_node*> *states, state_node* startNode, int horizon) {
	int currentState = startNode->state;
	int nrActions = decpomdp->GetNrJointActions();

	for (int i = 0; i < horizon; i++) {
		int action = rand() % nrActions;
		currentState = decpomdp->SampleSuccessorState(currentState, action);
		state_node* stateNode = getNode(states, currentState);
		if (stateNode->isWinning) return 1.0 * pow(decpomdp->GetDiscount(), i);
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
	set<action_node*, comparator<action_node> > *actions = &(currentNode->children);

	for (int action = 0; action < nrActions; action++) {
		double uctValue;
		action_node* dummy = new action_node(currentNode, action);
		action_node* actionNode;
		set<action_node*>::iterator result = actions->find(dummy);
		if (result != actions->end()) {
			actionNode = (*result);
		}
		else {
			actionNode = dummy;
			actions->insert(actionNode);
		}

		if (actionNode->iterations > 0) {
			uctValue = actionNode->average() + 2 * exploration * sqrt(2 * log(currentNode->iterations) / actionNode->iterations);
		}
		else {
			uctValue = 10000 + rand() % 1000;
		}

		if (uctValue > maxUct)
		{
			selectedAction = actionNode;
			maxUct = uctValue;
		}
	}

	return selectedAction;
}

double MCTS(DecPOMDPDiscreteInterface* decpomdp, map<int, state_node*> *states, state_node* currentNode, int horizon)
{
	if (horizon <= 0) {
		return 0.0;
	}
	else {
		double discount = decpomdp->GetDiscount();
		action_node* action = select_action(currentNode, decpomdp);
		int nextState = decpomdp->SampleSuccessorState(currentNode->state, action->action);
		//cout << "Going from state " << currentNode->state << " to state " << nextState << " with action " << action->action << endl;
		double reward;

		currentNode->iterations++;
		action->iterations++;

		state_node* stateNode = getNode(states, nextState);
		if (stateNode->isWinning) {
			reward = 1.0;
		}
		else if (stateNode->isLosing) {
			reward = 0.0;
		}
		else {
			set<state_node*>::iterator result = action->children.find(stateNode);
			if (result == action->children.end()) {
				action->children.insert(stateNode);
			}

			if (stateNode->iterations > 0) { // Gotten here before, go deeper
				reward = MCTS(decpomdp, states, stateNode, horizon - 1) * discount;
			}
			else { // Never gotten here before
				reward = simulation(decpomdp, states, stateNode, horizon - 1) * discount;
				//cout << "Simulation for node " << nextState << " gives reward " << reward << endl;
				stateNode->iterations++;
				stateNode->sum += reward;
			}
		}
		action->sum += reward;
		currentNode->sum += reward;

		return reward;
	}
}

action_node* select_final_action(state_node* currentNode, DecPOMDPDiscreteInterface* decpomdp)
{
	int nrActions = decpomdp->GetNrJointActions();
	int maxIterations = 0;
	action_node* selectedAction;
	set<action_node*, comparator<action_node> > *actions = &(currentNode->children);
	for (set<action_node*, comparator<action_node> >::iterator ii = actions->begin(); ii != actions->end(); ++ii)
	{
		if ((*ii)->iterations > maxIterations) {
			maxIterations = (*ii)->iterations;
			selectedAction = (*ii);
		}
	}

	if (maxIterations == 0) { // No idea what to do, pick randomly
		selectedAction = new action_node(currentNode, rand() % nrActions);
	}

	return selectedAction;
}

double final_simulation(DecPOMDPDiscreteInterface* decpomdp, map<int, state_node*> *states, state_node* currentNode, int* horizon)
{
	if (*horizon <= 0) {
		return 0.0;
	}
	
	(*horizon)--;
	action_node* action = select_final_action(currentNode, decpomdp);
	int nextState = decpomdp->SampleSuccessorState(currentNode->state, action->action);
	//cout << "Choosing action " << action->action << " in state " << currentNode->state << " takes us into state " << nextState << endl;

	state_node* stateNode = getNode(states, nextState);
	if (stateNode->isWinning) {
		return 1.0;
	} 
	if (stateNode->isLosing) {
		return 0.0;
	}				

	return final_simulation(decpomdp, states, stateNode, horizon);
	
}

double vi_simulation(DecPOMDPDiscreteInterface* decpomdp, QTable q, map<int, state_node*> *states, state_node* currentNode, int* horizon)
{
	if (*horizon <= 0) {
		return 0.0;
	}

	(*horizon)--;
	int nrActions = decpomdp->GetNrJointActions();
	int actionMax = 0;
	double qMax = -DBL_MAX;
	for (int i = 0; i < nrActions; ++i)
	{
		double qVal = q(currentNode->state, i);
		if (qVal > qMax)
		{
			qMax = qVal;
			actionMax = i;
		}
	}

	int nextState = decpomdp->SampleSuccessorState(currentNode->state, actionMax);
	state_node* stateNode = getNode(states, nextState);
	if (stateNode->isWinning) {
		return 1.0;
	}
	if (stateNode->isLosing) {
		return 0.0;
	}

	return vi_simulation(decpomdp, q, states, stateNode, horizon);
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

		srand(10); //Seed to start in state 7 for 4x3

		int nrStates = decpomdp->GetNrStates();
		int nrActions = decpomdp->GetNrJointActions();

		int initialState = decpomdp->SampleInitialState();
		cout << "Starting in state " << initialState << endl;

		// Pure BFS. Takes a long time, answer is only marginally better than pure random.
		//long steps = pow(nrActions, args.horizon);
		//cout << "Maximum BFS reward found (" << steps << " steps) is " << BFS(decpomdp, steps, initialState) << endl;

		// Pure random
		//double discount = decpomdp->GetDiscount();
		//double maxReward = -std::numeric_limits<double>::max();

		//for (int i = 0; i < 10000; i++) {
		//	size_t currentState = initialState;
		//	double sumReward = 0;
		//	for (int j = 0; j < args.horizon; j++) {
		//		int action = rand() % nrActions;
		//		sumReward += decpomdp->GetReward(currentState, action)*pow(discount, j);
		//		currentState = decpomdp->SampleSuccessorState(currentState, action);
		//	}
		//	if (sumReward > maxReward) maxReward = sumReward;
		//}
		//cout << "Maximum random search reward found is " << maxReward << endl;

		map<int, state_node*> states;
		state_node* root = getNode(&states, initialState);
		addWinningAndLosingStates(decpomdp, &states);


		clock_t startTime = clock();
		// MCTS
		for (int i = 0; i < 10000; i++) {
			MCTS(decpomdp, &states, root, args.horizon);
		}
		cout << "MCTS took " << double(clock() - startTime) / (double)CLOCKS_PER_SEC << " seconds." << endl;

		//// Print Nodes
		//for (int i = 0; i < nrStates; i++) {
		//	cout << getNode(&states, i)->toString() << endl;
		//}

		//// Print best actions
		//for (map<int, state_node*>::iterator ii = states.begin(); ii != states.end(); ++ii)
		//{
		//	int maxIterations = 0;
		//	int maxAction = -1;
		//	state_node* state = ii->second;
		//	for (set<action_node*, comparator<action_node> >::iterator jj = state->children.begin(); jj != state->children.end(); ++jj)
		//	{
		//		if ((*jj)->iterations > maxIterations) {
		//			maxIterations = (*jj)->iterations;
		//			maxAction = (*jj)->action;
		//		}
		//	}
		//	cout << "Best action for state " << ii->first << " is " << maxAction << endl;
		//}

		// Simulate results
		double finalSum = 0.0;
		int finalStepsSum = 0;
		int sims = 1000;
		//cout << "strtoi(unlist(strsplit(\"";
		for(int i = 0; i<sims; i++)
		{
			int stepsLeft = args.horizon;
			double reward = final_simulation(decpomdp, &states, root, &stepsLeft);
			finalSum += reward;
			finalStepsSum += args.horizon - stepsLeft;
			if (reward > 0.0)
				cout << (args.horizon - stepsLeft) << ",";
		}
		//cout << "\", split=\",\")))" << endl;
		cout << endl;
		cout << "MCTS has a win chance of " << setprecision(3) << (finalSum / sims) << " in an average of " << ((double)finalStepsSum / sims) << " steps" << endl;

		srand(10); // Reset seed so VI is not affected by MCTS

		startTime = clock();
		PlanningUnitDecPOMDPDiscrete *np = new NullPlanner(args.horizon, decpomdp);
		MDPValueIteration vi(*np);
		vi.Plan();
		QTable q = vi.GetQTable(0);
		cout << "VI took " << double(clock() - startTime) / (double)CLOCKS_PER_SEC << " seconds." << endl;

		finalSum = 0.0;
		finalStepsSum = 0;
		//cout << "strtoi(unlist(strsplit(\"";
		for (int i = 0; i<sims; i++)
		{
			int stepsLeft = args.horizon;
			double reward =	vi_simulation(decpomdp, q, &states, root, &stepsLeft);
			finalSum += reward;
			finalStepsSum += args.horizon - stepsLeft;
			if(reward > 0.0)
				cout << (args.horizon - stepsLeft) << ",";
		}
		//cout << "\", split=\",\")))" << endl;
		cout << endl;
		cout << "Value iteration has a win chance of " << setprecision(3) << (finalSum / sims) << " in an average of " << ((double)finalStepsSum / sims) << " steps" << endl;

	}
	catch(E& e){ e.Print(); }

	return(0);
}
