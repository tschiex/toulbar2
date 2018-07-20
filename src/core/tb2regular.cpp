/*----------------------------------------------------------------------
 *
 * Incremental regular constraint
 *
 */

#include "tb2regular.hpp"
#include "search/tb2clusters.hpp"

//--------------------------------------------------
// output operator for vector, pair
namespace std {

template <typename T>
ostream& operator<<(ostream& os, vector<T> const& v)
{
    os << "v(sz=" << v.size() << ")[";
    bool first = true;
    for (auto&& t : v) {
        if (first)
            first = false;
        else
            os << ",";
        os << t;
    }
    os << "]";
    return os;
}

template <typename U, typename T>
ostream& operator<<(ostream& os, pair<U, T> const& p)
{
    return os << "p{" << p.first << "," << p.second << "}";
}
} // namespace std

bool DACCompare(const EnumeratedVariable* v1, const EnumeratedVariable* v2)
{
    return (v1->getDACOrder() < v2->getDACOrder());
}

// ---------------------------------------------------
// ------------ The WRegular class -------------------
// ---------------------------------------------------

// A counting one, reads a vector of values and a min/max distance
WRegular::WRegular(WCSP* wcsp, EnumeratedVariable** scope_in, int arity_in, istream& file)
    : AbstractNaryConstraint(wcsp, scope_in, arity_in)
    , intDLinkStore(arity_in * 10) // TODO something less naive would be good
    , lb(MIN_COST)
{
    static const bool debug{ true };

    Cost weight;
    file >> weight; // violation cost
    bool boundByAbove;
    file >> boundByAbove; // max dist or min dist
    int distBound;
    file >> distBound; // the actual bound

    // Copy the scope and DAC order it
    Scope.assign(scope_in, scope_in + arity_in);
    DACScope.assign(scope_in, scope_in + arity_in);
    sort(DACScope.begin(), DACScope.end(), DACCompare);

    allArcs.resize(get_layer_num());
    arcsAtLayerValue.resize(get_layer_num());
    delta.resize(get_layer_num());
    //    alpha.resize(get_layer_num()+1);
    //    beta.resize(get_layer_num()+1);
    alphap.resize(get_layer_num() + 1);
    //    betap.resize(get_layer_num()+1);
    layerWidth.push_back(1); // one node to start

    // Create counting arcs
    // Would need to be minimized on the final layers
    for (int layer = 0; layer < get_layer_num(); layer++) {
        conflictWeights.push_back(0);
        delta[layer].resize(DACScope[layer]->getDomainInitSize(), MIN_COST);
        //    alpha[layer].resize(layerWidth[layer],MAX_COST);
        //    beta[layer].resize(layerWidth[layer],MAX_COST);
        alphap[layer].resize(layerWidth[layer], (layer ? MAX_COST : MIN_COST));
        //    betap[layer].resize(layerWidth[layer],MAX_COSTq);

        int corrNext = max(0, distBound - get_layer_num() + layer + 1);

        if (layer < get_layer_num() - 1) {
            // the count in next layer cannot exceed the layer number
            // all counts above the bound are equivalent (monotonic)
            int maxw = min(layer + 2, distBound + 2);
            corrNext = min(corrNext, maxw - 1); // one path at least
            layerWidth.push_back(maxw - corrNext);
        } else
            layerWidth.push_back(1); // last layer has one exit node
        int corr = max(0, corrNext - 1);

        int arcRef = 0;
        arcsAtLayerValue[layer].resize(DACScope[layer]->getDomainInitSize(), &intDLinkStore);
        unsigned forbidValue;
        file >> forbidValue;
        for (int node = 0; node < layerWidth[layer]; node++) {
            for (unsigned val = 0; val < DACScope[layer]->getDomainInitSize(); val++) {
                int nextCount = node + corr + (val != forbidValue);
                int nextNode = min(layerWidth[layer + 1] - 1, max(0, nextCount - corrNext));
                bool needToPay = (layer == get_layer_num() - 1) && ((nextCount > distBound) ^ boundByAbove);
                allArcs[layer].push_back(Arc(node, val, needToPay ? weight : MIN_COST, nextNode));
                arcsAtLayerValue[layer][val].push_back(arcRef);
                arcRef++;
            }
        }
    }

    if (debug) {
        ofstream os("wregular.dot");
        printLayers(os);
        os.close();
        String sol = L"0000";
        cout << eval(sol) << endl;
    }
}

// This one needs to be finished. may be a cleanDanglingNodes method outside of it would be better.
WRegular::WRegular(WCSP* wcsp, EnumeratedVariable** scope_in, int arity_in, WFA& automaton)
    : AbstractNaryConstraint(wcsp, scope_in, arity_in)
    , intDLinkStore(arity_in * 10) // TODO something less naive would be good
    , lb(MIN_COST)
{
    // Useless, should be in super class ?
    Scope.assign(scope_in, scope_in + arity_in);
    // Copy the scope and DAC order it
    DACScope.assign(scope_in, scope_in + arity_in);
    sort(DACScope.begin(), DACScope.end(), DACCompare);

    // Unroll the automata
    // Only DFA for now: one initial state, one transition per value
    if (automaton.getAcceptingStates().size() != 1) {
        cerr << "Only single accepting states automata can be used with regular.\n";
        exit(EXIT_FAILURE);
    }
    map<unsigned int, unsigned int> stateTranslate;
    stateTranslate[automaton.getAcceptingStates()[0].first] = 0;
    vector<unsigned int> stateRevTranslate;
    stateRevTranslate.push_back(automaton.getAcceptingStates()[0].first);

    Cost initialWeight = automaton.getAcceptingStates()[0].second;
    vector<map<int, pair<int, Cost> > > transitions; //per state and value: target and weight
    transitions.resize(automaton.getNbStates());
    for (const auto& transition : automaton.getTransitions()) {
        if (transitions[transition->start].find(transition->symbol) == transitions[transition->start].end())
            transitions[transition->start][transition->symbol] = make_pair(transition->end, transition->weight);
        else {
            cerr << "Only deterministic automata can be used with regular.";
            exit(EXIT_FAILURE);
        }
    }

    // map of final states to prune / weight last layer
    map<int, Cost> finalStates;
    for (auto accepting : automaton.getAcceptingStates())
        finalStates[accepting.first] = accepting.second;

    //Start layering
    map<unsigned int, unsigned int> nextTranslate; // maps automaton states to nodeIdx in layer
    vector<unsigned int> nextRevTranslate; // the other way around
    // Refs to swap at each layer
    map<unsigned int, unsigned int>& currentTranslateRef = stateTranslate;
    map<unsigned int, unsigned int>& nextTranslateRef = nextTranslate;
    vector<unsigned int>& currentRevTranslateRef = stateRevTranslate;
    vector<unsigned int>& nextRevTranslateRef = nextRevTranslate;
    // variables starting with "a/u" are resp. automaton/unrolled automaton states.
    // second letter may be c or n for current/next layer
    // Translate got from a to u, RevTranslate from u to a

    // allocate for all layers
    allArcs.resize(DACScope.size());
    degrees.resize(DACScope.size() + 1);

    for (unsigned int layer = 0; layer < DACScope.size(); layer++) {
        unsigned int unStateId = 0;
        ArcRef arcIdx = 0;
        conflictWeights.push_back(0);
        // allocate for all layer states
        degrees[layer].resize(stateRevTranslate.size());
        for (unsigned int ucState = 0; ucState < stateRevTranslate.size(); ucState++) {
            unsigned int acState = currentRevTranslateRef[ucState];

            for (const auto& transition : transitions[acState]) {
                unsigned int anTarget = transition.second.first;
                if (layer < DACScope.size() - 1 || finalStates.count(anTarget)) { // last layer: prune with accepting states
                    if (!nextTranslate.count(anTarget)) // the target is unknown and must be mapped
                        nextTranslateRef[anTarget] = unStateId++;
                    // let's create arcs
                    unsigned int unTarget = nextTranslateRef[anTarget]; // get the target's id in the new layer
                    Cost arcWeight = transition.second.second;
                    if (layer == 0)
                        arcWeight += initialWeight;
                    if (layer == DACScope.size() - 1 && finalStates.count(anTarget))
                        arcWeight += finalStates[anTarget];
                    Arc myArc(ucState, transition.first, transition.second.second, unTarget);
                    allArcs[layer].push_back(myArc);
                    // TODO push in the active arcs at layer/value
                    arcsAtLayerValue[layer][transition.first].push_back(arcIdx);
                    arcIdx++;
                }
            }
        }
        // all states done, lets update translators: clear current and swap
        currentTranslateRef.clear();
        currentRevTranslateRef.clear();
        swap(currentTranslateRef, nextTranslateRef);
        swap(currentRevTranslateRef, nextRevTranslateRef);
    }

    // Now we should clean states with no successors
    // to be put in a separate method later (useful for various creators)
    for (unsigned int layer = DACScope.size() - 1; layer >= 0; layer--) {
    }
    // Init conflict weights
    for (int i = 0; i != arity_; ++i) {
        conflictWeights.push_back(0);
    }
}

WRegular::~WRegular()
{
}

// Export the unrolled automata for debugging purposes //
// Could be nice to see All arcs and separate deleted from non deleted
// but the deletion status is in the DLink, not in the arc.
std::ostream& WRegular::printLayers(std::ostream& os)
{

    os << "digraph \"wregular\" {" << endl;
    os << "\tgraph [hierarchic=1];" << endl;
    // Draw vertices
    int nodeShift = 0;
    for (int layer = 0; layer <= get_layer_num(); layer++) {
        for (int node = 0; node < layerWidth[layer]; node++) {
            os << "\t" << nodeShift + node << " [name=\"" << layer << "," << node << "\"];" << endl;
        }
        nodeShift += layerWidth[layer];
    }
    // and Arcs
    nodeShift = 0;
    for (int layer = 0; layer < get_layer_num(); layer++) {
        for (unsigned val = 0; val < DACScope[layer]->getDomainInitSize(); val++) {
            for (ArcRef arc : arcsAtLayerValue[layer][val]) {
                os << "\t" << nodeShift + allArcs[layer][arc].source << " -> " << nodeShift + layerWidth[layer] + allArcs[layer][arc].target << " [label=\"";
                os << allArcs[layer][arc].value << "," << allArcs[layer][arc].weight << "\"];" << endl;
            }
        }
        nodeShift += layerWidth[layer];
    }
    os << "}";
    return os;
}

std::ostream& WRegular::printstate(std::ostream& os)
{
    os << "lb = " << lb << " connected = " << connected()
       << " depth = " << Store::getDepth() << "\n";
    for (int i = 0; i != arity_; ++i) {
        auto* x = scope[i];
        if (connected(i))
            os << " * ";
        else
            os << "   ";
        os << "var " << i << " ";
        x->print(os);
        if (x->assigned() && connected(i))
            os << "*****";
        os << "\n";
    }
    return os;
}

void WRegular::projectLB(Cost c)
{
    lb += c;
    assert(lb <= c);
    Constraint::projectLB(c);
    for (int layer = 0; layer < get_layer_num() + 1; layer++) {
        for (auto& t : alphap[layer]) {
            t -= c;
        }
    }
}

void WRegular::extend(int idx, unsigned val, Cost c)
{
    EnumeratedVariable* x = DACScope[idx];
    if (x->unassigned()) {
        delta[idx][val] += c;
        x->extend(x->toValue(val), c);
    }
}

void WRegular::forwardoic()
{
    static const bool debug{ true };

    vector<vector<Cost> > unaryCostExtension(get_layer_num());
    for (int layer = 0; layer < get_layer_num(); layer++) {
        alphap[layer + 1].assign(layerWidth[layer + 1], wcsp->getUb());
        EnumeratedVariable* x = DACScope[layer];
        for (unsigned val = 0; val < DACScope[layer]->getDomainInitSize(); val++) {
            for (auto arc : arcsAtLayerValue[layer][val]) {
                if (debug) {
                    cout << "alphap[" << layer + 1 << "][" << allArcs[layer][arc].get_target() << "] is " << alphap[layer + 1][allArcs[layer][arc].get_target()];
                    cout << " compared to " << alphap[layer][allArcs[layer][arc].get_source()] << "+" << allArcs[layer][arc].get_weight() << "+" << x->getCost(x->toValue(val)) << endl;
                }
                if (alphap[layer + 1][allArcs[layer][arc].get_target()] > alphap[layer][allArcs[layer][arc].get_source()] + allArcs[layer][arc].get_weight() + x->getCost(x->toValue(val))) {
                    alphap[layer + 1][allArcs[layer][arc].get_target()] = alphap[layer][allArcs[layer][arc].get_source()] + allArcs[layer][arc].get_weight() + x->getCost(x->toValue(val));

                    if (debug) {
                        cout << "\t=> set to " << alphap[layer + 1][allArcs[layer][arc].get_target()] << endl;
                    }
                }
            }
        }

        unaryCostExtension[layer].resize(DACScope[layer]->getDomainInitSize(), MIN_COST);
        for (unsigned val = 0; val < DACScope[layer]->getDomainInitSize(); val++) {
            for (auto it = arcsAtLayerValue[layer][val].begin(); it != arcsAtLayerValue[layer][val].end(); ++it) {
                if (alphap[layer + 1][allArcs[layer][*it].get_target()] < wcsp->getUb()) {
                    unaryCostExtension[layer][val] = max(
                        unaryCostExtension[layer][val],
                        alphap[layer + 1][allArcs[layer][*it].get_target()] - alphap[layer][allArcs[layer][*it].get_source()] - allArcs[layer][*it].get_weight());
                } else {
                    arcsAtLayerValue[layer][val].erase(it);
                    if (debug) {
                        cout << "Erasing arc " << *it << " at layer " << layer << endl;
                    }
                }
            }
            extend(layer, val, unaryCostExtension[layer][val]);
        }
    }
    Cost cmin = wcsp->getUb();
    for (const auto& t : alphap[get_layer_num()]) {
        if (cmin > t) {
            cmin = t;
        }
    }
    if (debug)
        cout << "regular S0IC initial bound: " << cmin << endl;

    projectLB(cmin);
}

void WRegular::propagate()
{
    forwardoic();
}

void WRegular::updateap(int layer, vector<int> states)
{
    // Updates alphap[idx][state] for state in states
    // always called with layer >= 1
    assert(layer > 0);
    vector<Cost> alphap_old(states.size());
    EnumeratedVariable* x = DACScope[layer - 1];
    for (const auto& state : states) {
        alphap_old[state] = alphap[layer][state];
        alphap[layer][state] = MAX_COST;
    }
    vector<bool> toExtend(DACScope[layer - 1]->getDomainInitSize(), false);
    for (unsigned val = 0; val < DACScope[layer - 1]->getDomainInitSize(); val++) {
        for (auto it = arcsAtLayerValue[layer - 1][val].begin(); it != arcsAtLayerValue[layer - 1][val].end(); ++it) {
            if (count(states.begin(), states.end(), allArcs[layer - 1][*it].get_target())) {
                toExtend[val] = true;
                if (alphap[layer][allArcs[layer - 1][*it].get_target()] > alphap[layer - 1][allArcs[layer - 1][*it].get_source()]
                        + delta[layer - 1][allArcs[layer - 1][*it].get_source()] + allArcs[layer - 1][*it].get_weight()
                        + x->getCost(x->toValue(val))) {
                    alphap[layer][allArcs[layer - 1][*it].get_target()] = alphap[layer - 1][allArcs[layer - 1][*it].get_source()]
                        + delta[layer - 1][allArcs[layer - 1][*it].get_source()] + allArcs[layer - 1][*it].get_weight()
                        + x->getCost(x->toValue(val));
                }
            }
        }
    }
    Cost ext;
    for (unsigned val = 0; val < DACScope[layer - 1]->getDomainInitSize(); val++) {
        if (toExtend[val]) {
            ext = MIN_COST;
            for (auto it = arcsAtLayerValue[layer - 1][val].begin(); it != arcsAtLayerValue[layer - 1][val].end(); ++it) {
                if (alphap[layer][allArcs[layer][*it].get_target()] < wcsp->getUb()) {
                    ext = max(ext,
                        alphap[layer][allArcs[layer][*it].get_target()] - alphap[layer - 1][allArcs[layer][*it].get_source()] - allArcs[layer - 1][*it].get_weight() - delta[layer - 1][val]);
                } else {
                    arcsAtLayerValue[layer - 1][val].erase(it);
                }
            }
            extend(layer - 1, val, ext);
        }
    }
    if (layer < get_layer_num()) {
        vector<int> toUpdate;
        for (unsigned val = 0; val < DACScope[layer]->getDomainInitSize(); val++) {
            for (auto it = arcsAtLayerValue[layer][val].begin(); it != arcsAtLayerValue[layer][val].end(); ++it) {
                if (count(states.begin(), states.end(), allArcs[layer][*it].get_source()) && alphap_old[allArcs[layer][*it].get_source()] != alphap[layer][allArcs[layer][*it].get_source()]) {
                    toUpdate.push_back(allArcs[layer][*it].get_target());
                }
            }
        }
        updateap(layer + 1, toUpdate);
    }
}

Cost WRegular::eval(const String& s)
{
    static const bool debug{ false };

    Cost res = -lb;
    int q = 0;
    for (int i = 0; i < arity_; i++) {
        BTListWrapper<ArcRef>::iterator it;

        for (it = arcsAtLayerValue[i][s[i] - CHAR_FIRST].begin(); it != arcsAtLayerValue[i][s[i] - CHAR_FIRST].end(); ++it) {
            if (debug)
                cout << "Arc " << allArcs[i][*it].source << "-" << allArcs[i][*it].value << "-" << allArcs[i][*it].target << " " << allArcs[i][*it].weight << endl;
            if (allArcs[i][*it].get_source() == q)
                break;
        }

        if (it == arcsAtLayerValue[i][s[i] - CHAR_FIRST].end()) {
            return wcsp->getUb();
        } else {
            res += allArcs[i][*it].get_weight() + delta[i][s[i] - CHAR_FIRST];
        }
    }
    assert(res >= MIN_COST);
    return res;
}

void WRegular::assign(int idx)
{
    static const bool debug{ false };

    auto* x = scope[idx];

    if (debug)
        cout << "In assign " << idx << "=" << x->getValue() << "\n";

    if (!connected(idx))
        return;
    deconnect(idx);

    vector<int> toUpdate(get_layer_width(idx + 1));
    iota(begin(toUpdate), end(toUpdate), 0);
    updateap(idx + 1, toUpdate);
}

void WRegular::remove(int idx)
{
    static const bool debug{ false };

    if (debug)
        cout << "In remove " << idx << endl;

    vector<int> toUpdate(get_layer_width(idx + 1));
    iota(begin(toUpdate), end(toUpdate), 0);
    updateap(idx + 1, toUpdate);
}

void WRegular::increase(int idx)
{
    remove(idx);
}

void WRegular::decrease(int idx)
{
    remove(idx);
}

void WRegular::projectFromZero(int idx)
{
    static const bool debug{ false };

    if (!connected(idx))
        return;

    if (debug)
        cout << "In projectFromZero "
             << "\n";

    vector<int> toUpdate(get_layer_width(idx + 1));
    iota(begin(toUpdate), end(toUpdate), 0);
    updateap(idx + 1, toUpdate);
}

double WRegular::computeTightness()
{
    return 1.0;
}

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */
