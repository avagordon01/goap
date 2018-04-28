#include "goap.hpp"
#include "astar.hpp"

#include <string.h>
#include <stdio.h>

static int idx_for_atomname(actionplanner *ap, const char *atomname) {
    int idx;
    for (idx = 0; idx < ap->numatoms; ++idx)
        if (!strcmp(ap->atm_names[idx], atomname))
            return idx;

    if (idx < MAXATOMS) {
        ap->atm_names[idx] = atomname;
        ap->numatoms++;
        return idx;
    }

    return -1;
}

static int idx_for_actionname(actionplanner *ap, const char *actionname) {
    int idx;
    for (idx = 0; idx < ap->numactions; ++idx)
        if (!strcmp(ap->act_names[idx], actionname))
            return idx;

    if (idx < MAXACTIONS) {
        ap->act_names[idx] = actionname;
        ap->act_costs[idx] = 1; // default cost is 1
        ap->numactions++;
        return idx;
    }

    return -1;
}

void goap_actionplanner_clear(actionplanner *ap) {
    ap->numatoms = 0;
    ap->numactions = 0;
    for (int i = 0; i < MAXATOMS; ++i) {
        ap->atm_names[i] = 0;
    }
    for (int i = 0; i < MAXACTIONS; ++i) {
        ap->act_names[i] = 0;
        ap->act_costs[i] = 0;
        goap_worldstate_clear(ap->act_pre + i);
        goap_worldstate_clear(ap->act_pst + i);
    }
}

void goap_worldstate_clear(worldstate *ws) {
    ws->values = bfield();
    ws->dontcare = bfield().set();
}

bool goap_worldstate_set(actionplanner *ap, worldstate *ws, const char *atomname, bool value) {
    const int idx = idx_for_atomname(ap, atomname);
    if (idx == -1)
        return false;
    ws->values.set(idx, value);
    ws->dontcare.set(idx, false);
    return true;
}

extern bool goap_set_pre(actionplanner *ap, const char *actionname, const char *atomname, bool value) {
    const int actidx = idx_for_actionname(ap, actionname);
    const int atmidx = idx_for_atomname(ap, atomname);
    if (actidx == -1 || atmidx == -1)
        return false;
    goap_worldstate_set(ap, ap->act_pre + actidx, atomname, value);
    return true;
}

bool goap_set_pst(actionplanner *ap, const char *actionname, const char *atomname, bool value) {
    const int actidx = idx_for_actionname(ap, actionname);
    const int atmidx = idx_for_atomname(ap, atomname);
    if (actidx == -1 || atmidx == -1)
        return false;
    goap_worldstate_set(ap, ap->act_pst + actidx, atomname, value);
    return true;
}

bool goap_set_cost(actionplanner *ap, const char *actionname, int cost) {
    const int actidx = idx_for_actionname(ap, actionname);
    if (actidx == -1)
        return false;
    ap->act_costs[actidx] = cost;
    return true;
}

void goap_worldstate_description(const actionplanner *ap, const worldstate *ws, char *buf, int sz) {
    int added = 0;
    for (int i = 0; i < MAXATOMS; ++i) {
        if (!ws->dontcare.test(i)) {
            const char *val = ap->atm_names[i];
            char upval[128];
            size_t j;
            for (j = 0; j < strlen(val); ++j)
                upval[j] = (val[j] - 32);
            upval[j++] = 0;
            const bool set = ws->values.test(i);
            added = snprintf(buf, sz, "%s,", set ? upval : val);
            buf += added;
            sz -= added;
        }
    }
}

void goap_description(actionplanner *ap, char *buf, int sz) {
    int added = 0;
    for (int a = 0; a < ap->numactions; ++a) {
        added = snprintf(buf, sz, "%s:\n", ap->act_names[a]);
        sz -= added;
        buf += added;

        worldstate pre = ap->act_pre[a];
        worldstate pst = ap->act_pst[a];
        for (int i = 0; i < MAXATOMS; ++i)
            if (!pre.dontcare.test(i)) {
                bool v = pre.values.test(i);
                added = snprintf(buf, sz, "  %s==%d\n", ap->atm_names[i], v);
                sz -= added;
                buf += added;
            }
        for (int i = 0; i < MAXATOMS; ++i)
            if (!pst.dontcare.test(i)) {
                bool v = pst.values.test(i);
                added = snprintf(buf, sz, "  %s:=%d\n", ap->atm_names[i], v);
                sz -= added;
                buf += added;
            }
    }
}

static worldstate goap_do_action(actionplanner *ap, int actionnr, worldstate fr) {
    const worldstate pst = ap->act_pst[actionnr];
    const bfield unaffected = pst.dontcare;
    const bfield affected = ~unaffected;

    fr.values = (fr.values & unaffected) | (pst.values & affected);
    fr.dontcare &= pst.dontcare;
    return fr;
}

int goap_get_possible_state_transitions(
    actionplanner *ap, worldstate fr, worldstate *to, const char **actionnames, int *actioncosts, int cnt) {
    int writer = 0;
    for (int i = 0; i < ap->numactions && writer < cnt; ++i) {
        // see if precondition is met
        const worldstate pre = ap->act_pre[i];
        const bfield care = ~pre.dontcare;
        const bool met = ((pre.values & care) == (fr.values & care));
        if (met) {
            actionnames[writer] = ap->act_names[i];
            actioncosts[writer] = ap->act_costs[i];
            to[writer] = goap_do_action(ap, i, fr);
            ++writer;
        }
    }
    return writer;
}
