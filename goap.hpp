#pragma once

#include <cstdbool>
#include <cstdint>
#include <bitset>

#define MAXATOMS 256
#define MAXACTIONS 64

typedef std::bitset<MAXATOMS> bfield;

//!< Describes the world state by listing values (t/f) for all known atoms.
typedef struct {
    bfield values;   //!< Values for atoms.
    bfield dontcare; //!< Mask for atoms that do not matter.
} worldstate;

//!< Action planner that keeps track of world state atoms and its action repertoire.
typedef struct {
    const char *atm_names[MAXATOMS]; //!< Names associated with all world state atoms.
    int numatoms;                    //!< Number of world state atoms.

    const char *act_names[MAXACTIONS]; //!< Names of all actions in repertoire.
    worldstate act_pre[MAXACTIONS];  //!< Preconditions for all actions.
    worldstate act_pst[MAXACTIONS];  //!< Postconditions for all actions (action effects).
    int act_costs[MAXACTIONS];         //!< Cost for all actions.
    int numactions;                    //!< The number of actions in out repertoire.

} actionplanner;

//!< Initialize an action planner. It will clear all information on actions and state.
extern void goap_actionplanner_clear(actionplanner *ap);

//!< Initialize a worldstate to 'dontcare' for all state atoms.
extern void goap_worldstate_clear(worldstate *ws);

//!< Set an atom of worldstate to specified value.
extern bool goap_worldstate_set(actionplanner *ap, worldstate *ws, const char *atomname, bool value);

//!< Add a precondition for named action.
extern bool goap_set_pre(actionplanner *ap, const char *actionname, const char *atomname, bool value);

//!< Add a postcondition for named action.
extern bool goap_set_pst(actionplanner *ap, const char *actionname, const char *atomname, bool value);

//!< Set the cost for named action.
extern bool goap_set_cost(actionplanner *ap, const char *actionname, int cost);

//!< Describe the action planner by listing all actions with pre and post conditions. For debugging purpose.
extern void goap_description(actionplanner *ap, char *buf, int sz);

//!< Describe the worldstate by listing atoms that matter, in lowercase for false-valued, and uppercase for true-valued atoms.
extern void goap_worldstate_description(const actionplanner *ap, const worldstate *ws, char *buf, int sz);

//!< Given the specified 'from' state, list all possible 'to' states along with the action required, and the action cost. For internal use.
extern int goap_get_possible_state_transitions(
    actionplanner *ap, worldstate fr, worldstate *to, const char **actionnames, int *actioncosts, int cnt);
