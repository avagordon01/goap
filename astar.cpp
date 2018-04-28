#include "astar.hpp"
#include "goap.hpp"

#include <limits.h>

#define MAXOPEN 1024 //!< The maximum number of nodes we can store in the opened set.
#define MAXCLOS 1024 //!< The maximum number of nodes we can store in the closed set.

static astarnode_t opened[MAXOPEN]; //!< The set of nodes we should consider.
static astarnode_t closed[MAXCLOS]; //!< The set of nodes we already visited.

static int numOpened = 0; //!< The nr of nodes in our opened set.
static int numClosed = 0; //!< The nr of nodes in our closed set.

//!< This is our heuristic: estimate for remaining distance is the nr of mismatched atoms that matter.
static int calc_h(worldstate_t fr, worldstate_t to) {
    const bfield_t care = (to.dontcare ^ -1LL);
    const bfield_t diff = ((fr.values & care) ^ (to.values & care));
    int dist = 0;
    for (int i = 0; i < MAXATOMS; ++i)
        if ((diff & (1LL << i)) != 0)
            dist++;
    return dist;
}

//!< Internal function to look up a world state in our opened set.
static int idx_in(worldstate_t ws, astarnode_t* set, size_t size) {
    for (int i = 0; i < size; ++i)
        if (set[i].ws.values == ws.values)
            return i;
    return -1;
}

//!< Internal function to reconstruct the plan by tracing from last node to initial node.
static void reconstruct_plan(
    actionplanner_t *ap, astarnode_t *goalnode, const char **plan, worldstate_t *worldstates, int *plansize) {
    astarnode_t *curnode = goalnode;
    int idx = *plansize - 1;
    int numsteps = 0;
    while (curnode && curnode->actionname) {
        if (idx >= 0) {
            plan[idx] = curnode->actionname;
            worldstates[idx] = curnode->ws;
            const int i = idx_in(curnode->parentws, closed, numClosed);
            curnode = (i == -1) ? 0 : closed + i;
        }
        --idx;
        numsteps++;
    }
    idx++; // point to last filled

    if (idx > 0)
        for (int i = 0; i < numsteps; ++i) {
            plan[i] = plan[i + idx];
            worldstates[i] = worldstates[i + idx];
        }
    if (idx < 0)
        fprintf(stderr, "error Plan of size %d cannot be returned in buffer of size %d\n", numsteps, *plansize);

    *plansize = numsteps;
}

/* from: http://theory.stanford.edu/~amitp/GameProgramming/ImplementationNotes.html
OPEN = priority queue containing START
CLOSED = empty set
while lowest rank in OPEN is not the GOAL:
  current = remove lowest rank item from OPEN
  add current to CLOSED
  for neighbors of current:
    cost = g(current) + movementcost(current, neighbor)
    if neighbor in OPEN and cost less than g(neighbor):
      remove neighbor from OPEN, because new path is better
    if neighbor in CLOSED and cost less than g(neighbor): **
      remove neighbor from CLOSED
    if neighbor not in OPEN and neighbor not in CLOSED:
      set g(neighbor) to cost
      add neighbor to OPEN
      set priority queue rank to g(neighbor) + h(neighbor)
      set neighbor's parent to current
 */

int astar_plan(
    actionplanner_t *ap,
    worldstate_t start,
    worldstate_t goal,
    const char **plan,
    worldstate_t *worldstates,
    int *plansize) {
    // put start in opened list
    numOpened = 0;
    astarnode_t n0;
    n0.ws = start;
    n0.parentws = start;
    n0.g = 0;
    n0.h = calc_h(start, goal);
    n0.f = n0.g + n0.h;
    n0.actionname = 0;
    opened[numOpened++] = n0;
    // empty closed list
    numClosed = 0;

    do {
        if (numOpened == 0) {
            printf("Did not find a path.\n");
            return -1;
        }
        // find the node with lowest rank
        int lowestIdx = -1;
        int lowestVal = INT_MAX;
        for (int i = 0; i < numOpened; ++i) {
            if (opened[i].f < lowestVal) {
                lowestVal = opened[i].f;
                lowestIdx = i;
            }
        }
        // remove the node with the lowest rank
        astarnode_t cur = opened[lowestIdx];
        if (numOpened)
            opened[lowestIdx] = opened[numOpened - 1];
        numOpened--;
        //static char dsc[2048];
        //goap_worldstate_description( ap, &cur.ws, dsc, sizeof(dsc) );
        //printf( dsc );
        // if it matches the goal, we are done!
        const bfield_t care = (goal.dontcare ^ -1LL);
        const bool match = ((cur.ws.values & care) == (goal.values & care));
        if (match) {
            reconstruct_plan(ap, &cur, plan, worldstates, plansize);
            return cur.f;
        }
        // add it to closed
        closed[numClosed++] = cur;
        if (numClosed == MAXCLOS) {
            printf("Closed set overflow\n");
            return -1;
        } // ran out of storage for closed set
        // iterate over neighbours
        const char *actionnames[MAXACTIONS];
        int actioncosts[MAXACTIONS];
        worldstate_t to[MAXACTIONS];
        const int numtransitions =
            goap_get_possible_state_transitions(ap, cur.ws, to, actionnames, actioncosts, MAXACTIONS);
        //printf( "%d neighbours\n", numtransitions );
        for (int i = 0; i < numtransitions; ++i) {
            astarnode_t nb;
            const int cost = cur.g + actioncosts[i];
            int idx_o = idx_in(to[i], opened, numOpened);
            int idx_c = idx_in(to[i], closed, numClosed);
            // if neighbor in OPEN and cost less than g(neighbor):
            if (idx_o >= 0 && cost < opened[idx_o].g) {
                // remove neighbor from OPEN, because new path is better
                if (numOpened)
                    opened[idx_o] = opened[numOpened - 1];
                numOpened--;
                idx_o = -1; // BUGFIX: neighbor is no longer in OPEN, signal this so that we can re-add it.
            }
            // if neighbor in CLOSED and cost less than g(neighbor):
            if (idx_c >= 0 && cost < closed[idx_c].g) {
                // remove neighbor from CLOSED
                if (numClosed)
                    closed[idx_c] = closed[numClosed - 1];
                numClosed--;
                idx_c = -1; // BUGFIX: neighbour is no longer in CLOSED< signal this so that we can re-add it.
            }
            // if neighbor not in OPEN and neighbor not in CLOSED:
            if (idx_c == -1 && idx_o == -1) {
                nb.ws = to[i];
                nb.g = cost;
                nb.h = calc_h(nb.ws, goal);
                nb.f = nb.g + nb.h;
                nb.actionname = actionnames[i];
                nb.parentws = cur.ws;
                opened[numOpened++] = nb;
            }
            if (numOpened == MAXOPEN) {
                printf("Opened set overflow\n");
                return -1;
            } // ran out of storage for opened set
        }
    } while (true);

    return -1;
}
