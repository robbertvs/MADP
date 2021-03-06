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

/* Only include this header file once. */
#ifndef _BELIEFSET_H_
#define _BELIEFSET_H_ 1

#include <vector>
#include "JointBeliefInterface.h"

/// Represents a belief set.
typedef std::vector<JointBeliefInterface*> BeliefSet;

#endif /* !_BELIEFSET_H_ */

// Local Variables: ***
// mode:c++ ***
// End: ***
