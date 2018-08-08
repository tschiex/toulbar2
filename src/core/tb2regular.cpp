// TODO check lb is used when "true" costs are computed
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
    DACScope.assign(scope_in, scope_in + arity_in);
    //   sort(DACScope.begin(), DACScope.end(), DACCompare);

    allArcs.resize(get_layer_num());
    arcsAtLayerValue.resize(get_layer_num());
    Supp.resize(get_layer_num());
    delta.resize(get_layer_num());
    alpha.resize(get_layer_num() + 1);
    beta.resize(get_layer_num() + 1);
    alphap.resize(get_layer_num() + 1);
    //    betap.resize(get_layer_num()+1);
    layerWidth.push_back(1); // one node to start

    int nbSols;
    file >> nbSols;

    if (nbSols != 1) {
        cerr << "Error: regular can ony handle one solution." << endl;
        exit(EXIT_FAILURE);
    }
    // Create counting arcs
    // Would need to be minimized on the final layers
    for (int layer = 0; layer < get_layer_num(); layer++) {
        conflictWeights.push_back(0);
        delta[layer].resize(DACScope[layer]->getDomainInitSize(), MIN_COST);
        alpha[layer].resize(layerWidth[layer], (layer ? MAX_COST : MIN_COST));
        beta[layer].resize(layerWidth[layer], MAX_COST);
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
        Supp[layer].resize(DACScope[layer]->getDomainInitSize(), -1);
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
    alpha[get_layer_num()].resize(layerWidth[get_layer_num()], MAX_COST);
    beta[get_layer_num()].resize(layerWidth[get_layer_num()], MIN_COST);
    alphap[get_layer_num()].resize(layerWidth[get_layer_num()], MAX_COST);

    if (debug) {
        ofstream os(to_string(this) + "wregular.dot");
        printLayers(os);
        os.close();
        /*String sol = L"1302";
        cout << "Cost(1302) = " << eval(sol) << endl;
        sol = L"2031";
        cout << "Cost(2031) = " << eval(sol) << endl;*/
    }
}

// This one needs to be finished. may be a cleanDanglingNodes method outside of it would be better.
WRegular::WRegular(WCSP* wcsp, EnumeratedVariable** scope_in, int arity_in, WFA& automaton)
    : AbstractNaryConstraint(wcsp, scope_in, arity_in)
    , intDLinkStore(arity_in * 10) // TODO something less naive would be good
    , lb(MIN_COST)
{
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
            os
                << "*****";
        os << "\n";
    }
    return os;
}

void WRegular::projectLB(Cost c)
{
    static const bool debug{ true };

    if (debug)
        cout << "Project " << c << " on lb" << endl;

    lb += c;
    Constraint::projectLB(c);
    for (int layer = 0; layer < get_layer_num() + 1; layer++) {
        for (auto& t : alphap[layer]) {
            t -= c;
        }
        for (auto& t : alpha[layer]) {
            t -= c;
        }
    }
    if (debug) {
        cout << this->getName();
        /*String sol = L"1302";
        cout << "ProjectLB - Cost(1302) = " << eval(sol) << endl;
        sol = L"2031";
        cout << "ProjectLB - Cost(2031) = " << eval(sol) << endl;*/
    }
}

void WRegular::extend(int idx, unsigned val, Cost c)
{
    static const bool debug{ true };

    EnumeratedVariable* x = DACScope[idx];

    if (debug)
        cout << "Extend cost " << c << " from var " << idx << " value " << val << " (unary cost is " << x->getCost(x->toValue(val)) << ")" << endl;

    if (x->unassigned()) {
        delta[idx][val] += c;
        x->extend(x->toValue(val), c);
    }

    if (debug) {
        /*String sol = L"1302";
        cout << "Extend - Cost(1302) = " << eval(sol) << endl;

        sol = L"2031";
        cout << "Extend - Cost(2031) = " << eval(sol) << endl;*/
    }
}

void WRegular::project1(int idx, unsigned val, Cost c)
{
    static const bool debug{ true };

    EnumeratedVariable* x = DACScope[idx];

    if (debug)
        cout << "Project cost " << c << " to var " << idx << " value " << val << endl;
    cout << *x << endl;

    if (x->unassigned()) {
        delta[idx][val] -= c;
        x->project(x->toValue(val), c);
    }
    cout << *x << endl;
}

void WRegular::forwardoic() // TODO rendre idempotent // initialisation des alpha/beta à Top dans ce cas: ici ou créateur ?
{
    static const bool debug{ true };

    if (debug) {
        cout << "Entering forwardoic" << endl;
    }

    for (int layer = 0; layer < get_layer_num(); layer++) {

        alphap[layer + 1].assign(layerWidth[layer + 1], wcsp->getUb());
        EnumeratedVariable* x = DACScope[layer];
        for (auto itval = x->begin(); itval != x->end(); ++itval) {
            for (auto arc : arcsAtLayerValue[layer][*itval]) {
                /*if (debug) {
                    cout << "alphap[" << layer + 1 << "][" << allArcs[layer][arc].get_target() << "] is " << alphap[layer + 1][allArcs[layer][arc].get_target()];
                    cout << " compared to " << alphap[layer][allArcs[layer][arc].get_source()] << "+" << allArcs[layer][arc].get_weight() << "+" << delta[layer][*itval] << "+" << x->getCost(x->toValue(*itval)) << endl;
                }*/
                if (alphap[layer + 1][allArcs[layer][arc].get_target()] > alphap[layer][allArcs[layer][arc].get_source()] + allArcs[layer][arc].get_weight() + delta[layer][*itval] + x->getCost(x->toValue(*itval))) {
                    alphap[layer + 1][allArcs[layer][arc].get_target()] = alphap[layer][allArcs[layer][arc].get_source()] + allArcs[layer][arc].get_weight() + delta[layer][*itval] + x->getCost(x->toValue(*itval));

                    if (debug) {
                        //cout << "\t=> set to " << alphap[layer + 1][allArcs[layer][arc].get_target()] << endl;
                    }
                }
            }
        }
        Cost unaryCostExtension;
        for (auto itval = x->begin(); itval != x->end(); ++itval) {
            unaryCostExtension = MIN_COST;
            for (auto it = arcsAtLayerValue[layer][*itval].begin(); it != arcsAtLayerValue[layer][*itval].end(); ++it) {
                if (alphap[layer + 1][allArcs[layer][*it].get_target()] < wcsp->getUb() && alphap[layer][allArcs[layer][*it].get_source()] < wcsp->getUb()) { //
                    unaryCostExtension = max(
                        unaryCostExtension,
                        alphap[layer + 1][allArcs[layer][*it].get_target()] - alphap[layer][allArcs[layer][*it].get_source()] - allArcs[layer][*it].get_weight() - delta[layer][*itval]);
                } else {
                    arcsAtLayerValue[layer][*itval].erase(it);
                    if (debug) {
                        cout << "Erasing arc " << *it << " at layer " << layer << endl;
                    }
                }
            }
            if (unaryCostExtension > MIN_COST)
                extend(layer, *itval, unaryCostExtension);
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

// TODO A t'on vraiment besoin de 3 passes pour AC ? (beta/alpha/beta)
void WRegular::backwardac()
{
    static const bool debug{ true };

    if (debug)
        cout << "Backward" << endl;

    for (int layer = get_layer_num() - 1; layer >= 0; --layer) {
        beta[layer].assign(layerWidth[layer], wcsp->getUb());
        EnumeratedVariable* x = DACScope[layer];
        for (auto itval = x->begin(); itval != x->end(); ++itval) {
            for (auto arc : arcsAtLayerValue[layer][*itval]) {
                if (debug) {
                    cout << "beta[" << layer << "][" << allArcs[layer][arc].get_source() << "] is " << beta[layer][allArcs[layer][arc].get_source()];
                    cout << " compared to " << beta[layer + 1][allArcs[layer][arc].get_target()] << "+" << allArcs[layer][arc].get_weight() << "+" << delta[layer][*itval] << endl;
                }
                if (beta[layer][allArcs[layer][arc].get_source()] > beta[layer + 1][allArcs[layer][arc].get_target()] + allArcs[layer][arc].get_weight() + delta[layer][*itval]) {
                    beta[layer][allArcs[layer][arc].get_source()] = beta[layer + 1][allArcs[layer][arc].get_target()] + allArcs[layer][arc].get_weight() + delta[layer][*itval];

                    if (debug) {
                        cout << "\t=> set to " << beta[layer][allArcs[layer][arc].get_source()] << endl;
                    }
                }
            }
        }
        for (auto itval = x->begin(); itval != x->end(); ++itval) {
            for (auto it = arcsAtLayerValue[layer][*itval].begin(); it != arcsAtLayerValue[layer][*itval].end(); ++it) {
                if (beta[layer][allArcs[layer][*it].get_source()] > wcsp->getUb() || beta[layer + 1][allArcs[layer][*it].get_target()] > wcsp->getUb()) { //
                    arcsAtLayerValue[layer][*itval].erase(it);
                    if (debug) {
                        cout << "Erasing arc " << *it << " at layer " << layer << endl;
                    }
                }
            }
        }
    }
}

void WRegular::forwardac()
{
    static const bool debug{ true };

    if (debug)
        cout << "ForwardAC" << endl;

    for (int layer = 0; layer < get_layer_num(); layer++) {
        EnumeratedVariable* x = DACScope[layer];
        for (auto itval = x->begin(); itval != x->end(); ++itval) {
            Cost cmin = wcsp->getUb();
            for (auto arc : arcsAtLayerValue[layer][*itval]) {
                if (cmin > alpha[layer][allArcs[layer][arc].get_source()] + allArcs[layer][arc].get_weight() + delta[layer][*itval] + beta[layer + 1][allArcs[layer][arc].get_target()]) {
                    cmin = alpha[layer][allArcs[layer][arc].get_source()] + allArcs[layer][arc].get_weight() + delta[layer][*itval] + beta[layer + 1][allArcs[layer][arc].get_target()];
                    Supp[layer][*itval] = arc;
                }
            }
            if (cmin > MIN_COST)
                project1(layer, *itval, cmin);
        }
        alpha[layer + 1].assign(layerWidth[layer + 1], wcsp->getUb());
        for (auto itval = x->begin(); itval != x->end(); ++itval) {
            for (auto arc : arcsAtLayerValue[layer][*itval]) {
                if (debug) {
                    cout << "alpha[" << layer + 1 << "][" << allArcs[layer][arc].get_target() << "] is " << alpha[layer + 1][allArcs[layer][arc].get_target()];
                    cout << " compared to " << alpha[layer][allArcs[layer][arc].get_source()] << "+" << allArcs[layer][arc].get_weight() << "+" << delta[layer][*itval] << endl;
                }
                if (alpha[layer + 1][allArcs[layer][arc].get_target()] > alpha[layer][allArcs[layer][arc].get_source()] + allArcs[layer][arc].get_weight() + delta[layer][*itval]) {
                    alpha[layer + 1][allArcs[layer][arc].get_target()] = alpha[layer][allArcs[layer][arc].get_source()] + allArcs[layer][arc].get_weight() + delta[layer][*itval];

                    if (debug) {
                        cout << "\t=> set to " << alpha[layer + 1][allArcs[layer][arc].get_target()] << endl;
                    }
                }
            }
        }
    }
}

void WRegular::propagate()
{
    forwardoic();
    backwardac();
    forwardac();
    backwardac();
    /*String sol = L"1302";
    cout << "Propagate - Cost(1302) = " << eval(sol) << endl;
    sol = L"2031";
    cout << "Propagate - Cost(2031) = " << eval(sol) << endl;*/
}

// Called when a change occurs at layer on states
void WRegular::updatea(int layer, vector<int> states)
{
    static const bool debug{ true };

    if (debug)
        cout << "updatea on layer " << layer << endl;

    EnumeratedVariable* x = DACScope[layer - 1];

    assert(layer > 0);
    vector<Cost> alpha_old(get_layer_width(layer));
    for (const auto& state : states) {
        alpha_old[state] = alpha[layer][state];
        alpha[layer][state] = wcsp->getUb();
    }
    for (auto itval = x->begin(); itval != x->end(); ++itval) {
        for (auto arc : arcsAtLayerValue[layer - 1][*itval]) {
            if (count(states.begin(), states.end(), allArcs[layer - 1][arc].get_target())) {
                /*if (debug) {
                    cout << "alpha[" << layer << "][" << allArcs[layer - 1][arc].get_target() << "] is " << alpha[layer][allArcs[layer - 1][arc].get_target()];
                    cout << " compared to " << alpha[layer - 1][allArcs[layer - 1][arc].get_source()] << "+" << allArcs[layer - 1][arc].get_weight() << "+" << delta[layer - 1][*itval] << endl;
                }*/
                if (alpha[layer][allArcs[layer - 1][arc].get_target()] > alpha[layer - 1][allArcs[layer - 1][arc].get_source()]
                        + delta[layer - 1][*itval] + allArcs[layer - 1][arc].get_weight()) {
                    alpha[layer][allArcs[layer - 1][arc].get_target()] = alpha[layer - 1][allArcs[layer - 1][arc].get_source()]
                        + delta[layer - 1][*itval] + allArcs[layer - 1][arc].get_weight();
                    if (debug) {
                        //cout << "\t=> set to " << alpha[layer][allArcs[layer - 1][arc].get_target()] << endl;
                    }
                }
            }
        }
    }
    // TODO do we really need to remeber the old alpha?
    if (layer < get_layer_num()) {
        vector<int> toUpdate;
        for (auto itval = DACScope[layer]->begin(); itval != DACScope[layer]->end(); ++itval) {
            for (auto arc : arcsAtLayerValue[layer][*itval]) {
                if (count(states.begin(), states.end(), allArcs[layer][arc].get_source()) && alpha_old[allArcs[layer][arc].get_source()] != alpha[layer][allArcs[layer][arc].get_source()] && !count(toUpdate.begin(), toUpdate.end(), allArcs[layer][arc].get_target())) {
                    toUpdate.push_back(allArcs[layer][arc].get_target());
                }
            }
        }
        if (!toUpdate.empty())
            updatea(layer + 1, toUpdate);
    }
}

void WRegular::updateb(int layer, vector<int> states)
{
    static const bool debug{ true };

    if (debug)
        cout << "updateb on layer " << layer << endl;

    EnumeratedVariable* x = DACScope[layer];

    if (debug)
        cout << *x << endl;

    assert(layer < get_layer_num());
    vector<Cost> beta_old(get_layer_width(layer));
    for (const auto& state : states) {
        beta_old[state] = beta[layer][state];
        beta[layer][state] = wcsp->getUb();
    }
    for (auto itval = x->begin(); itval != x->end(); ++itval) {
        for (auto arc : arcsAtLayerValue[layer][*itval]) {
            if (count(states.begin(), states.end(), allArcs[layer][arc].get_source())) {
                /*if (debug) {
                    cout << "beta[" << layer << "][" << allArcs[layer][arc].get_source() << "] is " << beta[layer][allArcs[layer][arc].get_source()];
                    cout << " compared to " << beta[layer + 1][allArcs[layer][arc].get_target()] << "+" << allArcs[layer][arc].get_weight() << "+" << delta[layer][*itval] << endl;
                }*/
                if (beta[layer][allArcs[layer][arc].get_source()] > beta[layer + 1][allArcs[layer][arc].get_target()]
                        + delta[layer][*itval] + allArcs[layer][arc].get_weight()) {
                    beta[layer][allArcs[layer][arc].get_source()] = beta[layer + 1][allArcs[layer][arc].get_target()]
                        + delta[layer][*itval] + allArcs[layer][arc].get_weight();
                    if (debug) {
                        //cout << "\t=> set to " << beta[layer][allArcs[layer][arc].get_source()] << endl;
                    }
                }
            }
        }
    }
    for (auto itval = x->begin(); itval != x->end(); ++itval) {
        for (auto it = arcsAtLayerValue[layer][*itval].begin(); it != arcsAtLayerValue[layer][*itval].end(); ++it) {
            if (beta[layer][allArcs[layer][*it].get_source()] > wcsp->getUb() || beta[layer + 1][allArcs[layer][*it].get_target()] > wcsp->getUb()) { //
                arcsAtLayerValue[layer][*itval].erase(it);
                if (debug) {
                    cout << "Erasing arc " << *it << " at layer " << layer << endl;
                }
            }
        }
    }
    if (layer > 0) {
        vector<int> toUpdate;
        for (auto itval = DACScope[layer - 1]->begin(); itval != DACScope[layer - 1]->end(); ++itval) {
            for (auto arc : arcsAtLayerValue[layer - 1][*itval]) {
                if (count(states.begin(), states.end(), allArcs[layer - 1][arc].get_target()) && beta_old[allArcs[layer - 1][arc].get_target()] != beta[layer][allArcs[layer - 1][arc].get_target()] && !count(toUpdate.begin(), toUpdate.end(), allArcs[layer - 1][arc].get_source())) {
                    toUpdate.push_back(allArcs[layer - 1][arc].get_source());
                }
            }
        }
        if (!toUpdate.empty())
            updateb(layer - 1, toUpdate);
    }
}

bool WRegular::updateap(int layer, vector<int> states)
{
    // Updates alphap[idx][state] for state in states
    // always called with layer >= 1

    static const bool debug{ true };

    if (debug)
        cout << "updateap on layer " << layer << endl;

    assert(layer > 0);
    vector<Cost> alphap_old(get_layer_width(layer));
    EnumeratedVariable* x = DACScope[layer - 1];
    for (const auto& state : states) {
        alphap_old[state] = alphap[layer][state];
        alphap[layer][state] = wcsp->getUb();
    }
    for (auto itval = x->begin(); itval != x->end(); ++itval) {
        for (auto arc : arcsAtLayerValue[layer - 1][*itval]) {
            if (count(states.begin(), states.end(), allArcs[layer - 1][arc].get_target())) {
                /*                if (debug) {
                    cout << "alphap[" << layer << "][" << allArcs[layer - 1][arc].get_target() << "] is " << alphap[layer][allArcs[layer - 1][arc].get_target()];
                    cout << " compared to " << alphap[layer - 1][allArcs[layer - 1][arc].get_source()] << "+" << allArcs[layer - 1][arc].get_weight() << "+" << delta[layer - 1][val] << "+" << x->getCost(x->toValue(val)) << endl;
                }*/
                if (alphap[layer][allArcs[layer - 1][arc].get_target()] > alphap[layer - 1][allArcs[layer - 1][arc].get_source()]
                        + delta[layer - 1][*itval] + allArcs[layer - 1][arc].get_weight()
                        + x->getCost(x->toValue(*itval))) {
                    alphap[layer][allArcs[layer - 1][arc].get_target()] = alphap[layer - 1][allArcs[layer - 1][arc].get_source()]
                        + delta[layer - 1][*itval] + allArcs[layer - 1][arc].get_weight()
                        + x->getCost(x->toValue(*itval));
                    if (debug) {
                        //cout << "\t=> set to " << alphap[layer][allArcs[layer - 1][arc].get_target()] << endl;
                    }
                }
            }
        }
    }
    for (auto itval = x->begin(); itval != x->end(); ++itval) {
        for (auto it = arcsAtLayerValue[layer - 1][*itval].begin(); it != arcsAtLayerValue[layer - 1][*itval].end(); ++it) {
            if (alphap[layer][allArcs[layer - 1][*it].get_target()] > wcsp->getUb() || alphap[layer - 1][allArcs[layer - 1][*it].get_source()] > wcsp->getUb()) { //
                arcsAtLayerValue[layer - 1][*itval].erase(it);
            }
        }
    }
    if (layer < get_layer_num()) {
        vector<int> toUpdate;
        for (auto itval = DACScope[layer]->begin(); itval != DACScope[layer]->end(); ++itval) {
            for (auto it = arcsAtLayerValue[layer][*itval].begin(); it != arcsAtLayerValue[layer][*itval].end(); ++it) {
                if (count(states.begin(), states.end(), allArcs[layer][*it].get_source()) && alphap_old[allArcs[layer][*it].get_source()] != alphap[layer][allArcs[layer][*it].get_source()] && !count(toUpdate.begin(), toUpdate.end(), allArcs[layer][*it].get_target())) {
                    toUpdate.push_back(allArcs[layer][*it].get_target());
                }
            }
        }
        if (toUpdate.empty()) {
            return true;
        } else {
            return updateap(layer + 1, toUpdate);
        }
    } else {
        bool is_oic = false;

        for (const auto& t : alphap[get_layer_num()]) {
            if (t == MIN_COST) {
                is_oic = true;
            }
        }

        return is_oic;
    }
}

// TODO states could be a set ?
// TODO Too much job ?
void WRegular::updateoic(int layer)
{
    // Updates alphap[idx][state] for state in states
    // always called with layer >= 1

    static const bool debug{ true };

    if (debug)
        cout << "updateoic on layer " << layer << endl;

    EnumeratedVariable* x = DACScope[layer - 1];

    assert(layer > 0);

    vector<int> toUpdatea; // TODO set ?
    vector<int> toUpdateb; // TODO set ?
    Cost ext;
    for (auto itval = x->begin(); itval != x->end(); ++itval) {
        //if (toExtend[*itval]) {
        ext = MIN_COST;
        for (auto it = arcsAtLayerValue[layer - 1][*itval].begin(); it != arcsAtLayerValue[layer - 1][*itval].end(); ++it) {
            if (alphap[layer][allArcs[layer - 1][*it].get_target()] < wcsp->getUb() && alphap[layer - 1][allArcs[layer - 1][*it].get_source()] < wcsp->getUb()) { //
                // TODO the condition in the if was already tested in updateap, and infeasible arcs have alerady been removed
                ext = max(ext, alphap[layer][allArcs[layer - 1][*it].get_target()] - alphap[layer - 1][allArcs[layer - 1][*it].get_source()] - allArcs[layer - 1][*it].get_weight() - delta[layer - 1][*itval]);
            } else {
                arcsAtLayerValue[layer - 1][*itval].erase(it);
            }
            if (!count(toUpdatea.begin(), toUpdatea.end(), allArcs[layer - 1][*it].get_target()))
                toUpdatea.push_back(allArcs[layer - 1][*it].get_target());
            if (!count(toUpdateb.begin(), toUpdateb.end(), allArcs[layer - 1][*it].get_source()))
                toUpdateb.push_back(allArcs[layer - 1][*it].get_source());
        }
        if (ext > MIN_COST)
            extend(layer - 1, *itval, ext);
        //}
    }
    if (layer < get_layer_num() && !toUpdatea.empty())
        updatea(layer, toUpdatea);

    if (!toUpdateb.empty())
        updateb(layer - 1, toUpdateb);

    if (layer < get_layer_num()) {
        updateoic(layer + 1);
    } else {
        Cost cmin = wcsp->getUb();
        for (const auto& t : alphap[get_layer_num()]) {
            if (cmin > t) {
                cmin = t;
            }
        }
        if (debug)
            cout << "regular S0IC updated bound: " << cmin << endl;

        projectLB(cmin);
    }
}

void WRegular::checkSupport(int layer)
{
    static const bool debug{ true };

    EnumeratedVariable* x = DACScope[layer];

    if (debug)
        cout << "Checking support of variable " << *x << endl;

    assert(layer < get_layer_num());

    vector<int> toUpdatea;
    vector<int> toUpdateb;
    for (auto itval = x->begin(); itval != x->end(); ++itval) {
        if (alpha[layer][allArcs[layer][Supp[layer][*itval]].get_source()] + allArcs[layer][Supp[layer][*itval]].get_weight() + delta[layer][*itval] + beta[layer + 1][allArcs[layer][Supp[layer][*itval]].get_target()] != MIN_COST) {
            //look for or create new support
            Cost cmin = wcsp->getUb();
            for (auto arc : arcsAtLayerValue[layer][*itval]) {
                if (cmin > alpha[layer][allArcs[layer][arc].get_source()] + allArcs[layer][arc].get_weight() + delta[layer][*itval] + beta[layer + 1][allArcs[layer][arc].get_target()]) {
                    cmin = alpha[layer][allArcs[layer][arc].get_source()] + allArcs[layer][arc].get_weight() + delta[layer][*itval] + beta[layer + 1][allArcs[layer][arc].get_target()];
                    Supp[layer][*itval] = arc;
                }
                if (!count(toUpdatea.begin(), toUpdatea.end(), allArcs[layer][arc].get_target()))
                    toUpdatea.push_back(allArcs[layer][arc].get_target());

                if (!count(toUpdateb.begin(), toUpdateb.end(), allArcs[layer][arc].get_source()))
                    toUpdateb.push_back(allArcs[layer][arc].get_source());
            }
            if (cmin > MIN_COST)
                project1(layer, *itval, cmin);
        }
    }
    if (toUpdatea.size())
        updatea(layer + 1, toUpdatea);
    if (toUpdateb.size())
        updateb(layer, toUpdateb);
}

Cost WRegular::eval(const String& s)
{
    static const bool debug{ false };

    Cost res = -lb;
    int q = 0;
    for (int i = 0; i < arity_; i++) {
        BTListWrapper<ArcRef>::iterator it;

        for (it = arcsAtLayerValue[i][s[i] - CHAR_FIRST].begin(); it != arcsAtLayerValue[i][s[i] - CHAR_FIRST].end(); ++it) {
            if (allArcs[i][*it].get_source() == q)
                break;
        }

        if (it == arcsAtLayerValue[i][s[i] - CHAR_FIRST].end()) {
            return wcsp->getUb();
        } else {
            if (debug)
                cout << "Arc " << allArcs[i][*it].source << "-" << allArcs[i][*it].value << "-" << allArcs[i][*it].target << " " << allArcs[i][*it].weight << endl;
            res += allArcs[i][*it].get_weight() + delta[i][s[i] - CHAR_FIRST];
            q = allArcs[i][*it].target;
        }
    }
    assert(res >= MIN_COST);
    return res;
}

void WRegular::assign(int idx)
{
    static const bool debug{ true };

    auto* x = scope[idx];

    if (debug)
        cout << "In assign " << scope[idx]->getName() << "=" << x->getValue() << "\n";

    if (!connected(idx))
        return;
    deconnect(idx);

    vector<int> toUpdatea(get_layer_width(idx + 1));
    iota(begin(toUpdatea), end(toUpdatea), 0);
    updatea(idx + 1, toUpdatea);

    vector<int> toUpdateb(get_layer_width(idx));
    iota(begin(toUpdateb), end(toUpdateb), 0);
    updateb(idx, toUpdateb);

    if (!updateap(idx + 1, toUpdatea))
        updateoic(idx + 1);

    for (int layer = 0; layer < get_layer_num(); layer++) {
        checkSupport(layer);
    }
    if (debug) {
        /*String sol = L"1302";
        cout << "Assign - Cost(1302) = " << eval(sol) << endl;
        sol = L"2031";
        cout << "Assign - Cost(2031) = " << eval(sol) << endl;*/
        //for (const auto var : DACScope) {
        //cout << *var << endl;
        //}
    }
}

void WRegular::remove(int idx)
{
    static const bool debug{ true };

    if (debug) {
        cout << "In remove " << idx << endl;
    }

    vector<int> toUpdatea(get_layer_width(idx + 1));
    iota(begin(toUpdatea), end(toUpdatea), 0);
    updatea(idx + 1, toUpdatea);

    vector<int> toUpdateb(get_layer_width(idx));
    iota(begin(toUpdateb), end(toUpdateb), 0);
    updateb(idx, toUpdateb);

    if (!updateap(idx + 1, toUpdatea))
        updateoic(idx + 1);

    for (int layer = 0; layer < get_layer_num(); layer++) { // TODO All needed ?
        checkSupport(layer);
    }
    // TODO Should we delete arcs in Qia ?

    if (debug) {
        /*String sol = L"1302";
        cout << "Remove - Cost(1302) = " << eval(sol) << endl;
        sol = L"2031";
        cout << "Remove - Cost(2031) = " << eval(sol) << endl;*/
        //for (const auto var : DACScope) {
        //cout << *var << endl;
        //}
    }
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
    static const bool debug{ true };

    if (!connected(idx))
        return;

    if (debug)
        cout << "In projectFromZero "
             << "\n";

    vector<int> toUpdatea(get_layer_width(idx + 1));
    iota(begin(toUpdatea), end(toUpdatea), 0);

    if (!updateap(idx + 1, toUpdatea))
        updateoic(idx + 1);

    for (int layer = 0; layer < get_layer_num(); layer++) {
        checkSupport(layer);
    }
    if (debug) {
        /*String sol = L"1302";
        cout << "ProjectF0 - Cost(1302) = " << eval(sol) << endl;
        sol = L"2031";
        cout << "ProjectF0 - Cost(2031) = " << eval(sol) << endl;*/
        //for (const auto var : DACScope) {
        //cout << *var << endl;
        //}
    }
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
