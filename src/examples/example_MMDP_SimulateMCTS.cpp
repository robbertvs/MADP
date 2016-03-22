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

double MCTS(DecPOMDPDiscreteInterface* decpomdp, long steps, int state) {
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

    map<int, map<int, map<int, double> > > tree;
    for (int i = 0; i < nrStates; i++) {
      map<int, map<int, double> > stateMap;
      for (int j = 0; j < nrActions; j++) {
        map<int, double> stateActionMap;
        for (int k = 0; k < nrStates; k++) {
          double p = decpomdp->GetTransitionProbability(i, j, k);
          if (p > 0) stateActionMap[k] = p;
        }
        stateMap[j] = stateActionMap;
      }
      tree[i] = stateMap;
    }

    for (map<int, map<int, map<int, double> > >::iterator ii = tree.begin(); ii != tree.end(); ++ii)
    {
<<<<<<< Updated upstream
        cout << "Instantiating the problem..."<<endl;
        DecPOMDPDiscreteInterface* decpomdp = GetDecPOMDPDiscreteInterfaceFromArgs(args);
        cout << "...done."<<endl;

        //srand(time(NULL));
		srand(42);

		//Build transition tree (I don't think the domains are large enough that this won't work)
		int nrStates = decpomdp->GetNrStates();
		int nrActions = decpomdp->GetNrJointActions();

		size_t initialState = decpomdp->SampleInitialState();
=======
      cout << (*ii).first << ": " << endl;

      for (map<int, map<int, double> >::iterator jj = (*ii).second.begin(); jj != (*ii).second.end(); ++jj)
      {
        cout << "\t" << (*jj).first << ": " << endl;
        for (map<int, double>::iterator kk = (*jj).second.begin(); kk != (*jj).second.end(); ++kk)
        {
          cout << "\t\t" << (*kk).first << ": " << (*kk).second << endl;
        }
      }
    }

    size_t initialState = decpomdp->SampleInitialState();
>>>>>>> Stashed changes

    // Pure BFS. Takes a long time, answer is only marginally better than pure random.
    long steps = pow(nrActions, args.horizon);
    cout << "Maximum BFS reward found (" << steps << " steps) is " << BFS(decpomdp, steps, initialState) << endl;

    // Pure random
    double discount = decpomdp->GetDiscount();
    double maxReward = -std::numeric_limits<double>::max();
    tree_node root = newNode(0.0);

    for(int i = 0; i < 10000; i++) {
      size_t currentState = initialState;
      tree_node currentNode = root;
      double sumReward = 0;
<<<<<<< Updated upstream
      double sumReward = doActions(decpomdp, currentState, currentNode, horizon);
=======
      double sumReward = doActions(currentState, currentNode, horizon);
>>>>>>> Stashed changes
      if (sumReward > maxReward) maxReward = sumReward;
      //cout << "Run " << i << " with reward " << sumReward << endl;
    }
    cout << "Maximum random search reward found is " << maxReward << endl;
  }
  catch(E& e){ e.Print(); }
<<<<<<< Updated upstream

  return(0);
}

double doActions(tree_node currentNode, int horizon)
{
  // should do random move
  if(horizon == 0)
    return 0;

  currentNode.iterations++;

  int action = rand() % nrActions;
  if(!currentNode.children.contains(action))
  {
    currentNode.children.insert(action, newNode(0));
  }
  currentState = decpomdp->SampleSuccessorState(currentState, action);
  tree_node nextNode = currentNode.children.get(action);
  sumReward = decpomdp->GetReward(currentState, action)*pow(discount, j);

}

double simulation(DecPOMDPDiscreteInterface* decpomdp, tree_node startNode, int horizon) {
	int currentState = startNode.state;
	int nrActions = decpomdp->GetNrJointActions();
	double reward = 0;

	for (int i = 0; i < horizon; i++) {
		int action = rand() % nrActions;
		currentState = decpomdp->SampleSuccessorState(currentState, action);
		reward += pow(discount, i) * decpomdp->GetReward(currentState, action);
	}

	return reward;
}

=======

  return(0);
}

double doActions(tree_node currentNode, int horizon)
{
  // should do random move
  if(horizon == 0)
    return 0;

  currentNode.iterations++;

  int action = rand() % nrActions;
  if(!currentNode.children.contains(action))
  {
    currentNode.children.insert(action, newNode(0));
  }
  currentState = decpomdp->SampleSuccessorState(currentState, action);
  sumReward += decpomdp->GetReward(currentState, action)*pow(discount, j);
}

>>>>>>> Stashed changes
tree_node newNode(double value)
{
  tree_node node->average = value;
  return node;
}

struct tree_node
{
<<<<<<< Updated upstream
  int state;
=======
>>>>>>> Stashed changes
  int iterations;
  double average;
  std::map<int, tree_node> children;
};
