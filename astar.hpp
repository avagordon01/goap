#pragma once

#include <cstdio>

#include "goap.hpp"

//!< A node in our network of world states.
struct astarnode {
    worldstate_t ws;        //!< The state of the world at this node.
    int g;                  //!< The cost so far.
    int h;                  //!< The heuristic for remaining cost (don't overestimate!)
    int f;                  //!< g+h combined.
    const char *actionname; //!< How did we get to this node?
    worldstate_t parentws;  //!< Where did we come from?
};

//! Make a plan of actions that will reach desired world state. Returns total cost of the plan.
int astar_plan(
    actionplanner_t *ap,       //!< the goap action planner that holds atoms and action repertoire
    worldstate_t start,        //!< the current world state
    worldstate_t goal,         //!< the desired world state
    const char **plan,         //!< for returning all actions that make up plan
    worldstate_t *worldstates, //!< for returning intermediate world states
    int *plansize              //!< in: size of plan buffer, out: size of plan (in nr of steps)
);
