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
#include "Timing.h"
#include "directories.h"
#include "DecPOMDPDiscrete.h"
#include "BGBackupType.h"

using namespace std;

int main(int argc, char **argv)
{
    if(argc!=2)
    {
        cout << "Use as follows: analyzeTimings "
             << "<timingsFile>" << endl;
        return(1);
    }

    string timingsFile=string(argv[1]);

    try {
    Timing time;
    time.Load(timingsFile);
    time.PrintSummary();
    }
    catch(E& e){
        e.Print();
        return(1);
    }

    return(0);
}
