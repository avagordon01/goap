#include "astar.hpp"
#include "goap.hpp"

#include <limits.h>
#include <vector>

//!< This is our heuristic: estimate for remaining distance is the nr of mismatched atoms that matter.
static int calc_h(worldstate fr, worldstate to) {
    const bfield care = (to.dontcare ^ -1LL);
    const bfield diff = ((fr.values & care) ^ (to.values & care));
    int dist = 0;
    for (int i = 0; i < MAXATOMS; ++i)
        if ((diff & (1LL << i)) != 0)
            dist++;
    return dist;
}

//!< Internal function to look up a world state in our opened set.
static int idx_in(worldstate ws, std::vector<astarnode> set) {
    for (int i = 0; i < set.size(); ++i)
        if (set[i].ws.values == ws.values)
            return i;
    return -1;
}

//!< Internal function to reconstruct the plan by tracing from last node to initial node.
static void reconstruct_plan(
    actionplanner *ap, astarnode *goalnode, const char **plan, worldstate *worldstates, int *plansize, std::vector<astarnode> closed) {
    astarnode *curnode = goalnode;
    int idx = *plansize - 1;
    int numsteps = 0;
    while (curnode && curnode->actionname) {
        if (idx >= 0) {
            plan[idx] = curnode->actionname;
            worldstates[idx] = curnode->ws;
            const int i = idx_in(curnode->parentws, closed);
            curnode = (i == -1) ? 0 : &closed[i];
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
    actionplanner *ap,
    worldstate start,
    worldstate goal,
    const char **plan,
    worldstate *worldstates,
    int *plansize) {
    std::vector<astarnode> opened;
    std::vector<astarnode> closed;
    // put start in opened list
    astarnode n0;
    n0.ws = start;
    n0.parentws = start;
    n0.g = 0;
    n0.h = calc_h(start, goal);
    n0.f = n0.g + n0.h;
    n0.actionname = 0;
    opened.push_back(n0);

    while (!opened.empty()) {
        // find the node with lowest rank
        int lowestIdx = -1;
        int lowestVal = INT_MAX;
        for (int i = 0; i < opened.size(); ++i) {
            if (opened[i].f < lowestVal) {
                lowestVal = opened[i].f;
                lowestIdx = i;
            }
        }
        // remove the node with the lowest rank
        astarnode cur = opened[lowestIdx];
        opened[lowestIdx] = opened.back();
        opened.pop_back();
        // if it matches the goal, we are done!
        const bfield care = (goal.dontcare ^ -1LL);
        const bool match = ((cur.ws.values & care) == (goal.values & care));
        if (match) {
            reconstruct_plan(ap, &cur, plan, worldstates, plansize, closed);
            return cur.f;
        }
        // add it to closed
        closed.push_back(cur);
        // iterate over neighbours
        const char *actionnames[MAXACTIONS];
        int actioncosts[MAXACTIONS];
        worldstate to[MAXACTIONS];
        const int numtransitions =
            goap_get_possible_state_transitions(ap, cur.ws, to, actionnames, actioncosts, MAXACTIONS);
        for (int i = 0; i < numtransitions; ++i) {
            const int cost = cur.g + actioncosts[i];
            int idx_o = idx_in(to[i], opened);
            int idx_c = idx_in(to[i], closed);
            // if neighbor in OPEN and cost less than g(neighbor):
            if (idx_o != -1 && cost < opened[idx_o].g) {
                // remove neighbor from OPEN, because new path is better
                opened[idx_o] = opened.back();
                opened.pop_back();
                idx_o = -1; // BUGFIX: neighbor is no longer in OPEN, signal this so that we can re-add it.
            }
            // if neighbor not in OPEN and neighbor not in CLOSED:
            if (idx_c == -1 && idx_o == -1) {
                astarnode nb;
                nb.ws = to[i];
                nb.g = cost;
                nb.h = calc_h(nb.ws, goal);
                nb.f = nb.g + nb.h;
                nb.actionname = actionnames[i];
                nb.parentws = cur.ws;
                opened.push_back(nb);
            }
        }
    }

    printf("Did not find a path.\n");
    return -1;
}
