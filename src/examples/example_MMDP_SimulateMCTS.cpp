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

int main(int argc, char **argv)
{
    ArgumentHandlers::Arguments args;
    argp_parse (&ArgumentHandlers::theArgpStruc, argc, argv, 0, 0, &args);

    try
    {
        cout << "Instantiating the problem..."<<endl;
        DecPOMDPDiscreteInterface* decpomdp = GetDecPOMDPDiscreteInterfaceFromArgs(args);
        cout << "...done."<<endl;

        srand(time(NULL));

        size_t nrStates = decpomdp->GetNrStates();
        size_t nrActions = decpomdp->GetNrJointActions();
        size_t initialState = decpomdp->SampleInitialState();

		double discount = decpomdp->GetDiscount();
		double maxReward = -1000; // Needs to be Double min value
		
		cout << "States " << nrStates <<
			", Actions " << nrActions <<
			", Initial State " << initialState << endl;
		
        for(int i = 0; i < 500; i++) {
            size_t currentState = initialState;
			double sumReward = 0;
            for(int j = 0; j < 100; j++) {
				int action = rand() % nrActions;
                currentState = decpomdp->SampleSuccessorState(currentState, action);
				sumReward += decpomdp->GetReward(currentState, action)*pow(discount, j);
            }
			if (sumReward > maxReward) maxReward = sumReward;
			cout << "Run " << i << " with reward " << sumReward << endl;
        }
		cout << "Maximum reward found is " << maxReward << endl;
    }
    catch(E& e){ e.Print(); }

    return(0);
}
