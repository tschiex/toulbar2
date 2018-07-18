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

template <typename T> ostream& operator<<(ostream& os, vector<T> const& v)
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

}

bool DACCompare(const EnumeratedVariable* v1, const EnumeratedVariable* v2)
{
    return (v1->getDACOrder() < v2->getDACOrder());
}


// ---------------------------------------------------
// ------------ The WRegular class -------------------
// ---------------------------------------------------

// A first naive creator that allows everything with constant transition cost 
WRegular::WRegular(WCSP* wcsp, EnumeratedVariable** scope_in, int arity_in, Cost weight)
    : AbstractNaryConstraint(wcsp, scope_in, arity_in)
    , lb(MIN_COST)
{
   // Copy the scope and DAC order it
    Scope.assign(scope_in,scope_in+arity_in);
    DACScope.assign(scope_in,scope_in+arity_in);
    sort(DACScope.begin(),DACScope.end(),DACCompare);

    // Create a simple set of arcs that allows all values at each layer
    for (int layer = 0; layer < get_layer_num(); layer++){
        layerWidth.push_back(1); // one node
        int arcRef = 0;
        for (unsigned val = 0; val < DACScope[layer]->getDomainInitSize(); val++) {
            Arc newArc(0,val,weight,0);
            allArcs[layer].push_back(newArc);
            arcsAtLayerValue[layer][val].push_back(arcRef);
            arcRef++;
        }
    }
    layerWidth.push_back(1); // one node on the last layer too
}

WRegular::WRegular(WCSP* wcsp, EnumeratedVariable** scope_in, int arity_in, istream& file)
    : AbstractNaryConstraint(wcsp, scope_in, arity_in)
    , lb(MIN_COST)
{
    Cost weight;
    file >> weight;

   // Copy the scope and DAC order it
    Scope.assign(scope_in,scope_in+arity_in);
    DACScope.assign(scope_in,scope_in+arity_in);
    sort(DACScope.begin(),DACScope.end(),DACCompare);

    // Create a simple set of arcs that allows all values at each layer
    for (int layer = 0; layer < get_layer_num(); layer++){
        layerWidth.push_back(1); // one node
        int arcRef = 0;
        for (unsigned val = 0; val < DACScope[layer]->getDomainInitSize(); val++) {
            Arc newArc(0,val,weight,0);
            allArcs[layer].push_back(newArc);
            arcsAtLayerValue[layer][val].push_back(arcRef);
            arcRef++;
        }
    }
    layerWidth.push_back(1); // one node on the last layer too
}

    WRegular::WRegular(WCSP* wcsp, EnumeratedVariable** scope_in, int arity_in, WFA& automaton)
    : AbstractNaryConstraint(wcsp, scope_in, arity_in)
    , lb(MIN_COST)
{
    // Useless, should be in super class ?
    Scope.assign(scope_in,scope_in+arity_in);
    // Copy the scope and DAC order it
    DACScope.assign(scope_in,scope_in+arity_in);
    sort(DACScope.begin(),DACScope.end(),DACCompare);

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
    vector<map<int,pair<int,Cost>>> transitions; //per state and value: target and weight
    transitions.resize(automaton.getNbStates());
    for (const auto &transition : automaton.getTransitions()) {
        if (transitions[transition->start].find(transition->symbol) ==  transitions[transition->start].end())
            transitions[transition->start][transition->symbol]= make_pair(transition->end,transition->weight);
        else {
            cerr << "Only deterministic automata can be used with regular.";
            exit(EXIT_FAILURE);
        }
    }

    // map of final states to prune / weight last layer 
    map <int,Cost> finalStates;
    for (auto accepting : automaton.getAcceptingStates())
        finalStates[accepting.first] = accepting.second;
    
    //Start layering
    map<unsigned int, unsigned int> nextTranslate; // maps automaton states to nodeIdx in layer
    vector<unsigned int> nextRevTranslate; // the other way around
    // Refs to swap at each layer
    map<unsigned int, unsigned int> &currentTranslateRef = stateTranslate;
    map<unsigned int, unsigned int> &nextTranslateRef = nextTranslate;
    vector<unsigned int> &currentRevTranslateRef = stateRevTranslate;
    vector<unsigned int> &nextRevTranslateRef = nextRevTranslate;
    // variables starting with "a/u" are resp. automaton/unrolled automaton states.
    // second letter may be c or n for current/next layer
    // Translate got from a to u, RevTranslate from u to a

    // allocate for all layers
    allArcs.resize(DACScope.size());
    degrees.resize(DACScope.size()+1);

    for (unsigned int layer = 0; layer < DACScope.size(); layer ++){
        unsigned int unStateId = 0;
        ArcRef arcIdx = 0;
        // allocate for all layer states
        degrees[layer].resize(stateRevTranslate.size());
        for (unsigned int ucState = 0; ucState < stateRevTranslate.size();  ucState++) {
            unsigned int acState = currentRevTranslateRef[ucState];

            for (const auto& transition : transitions[acState]) {
                unsigned int anTarget = transition.second.first;
                if (layer < DACScope.size()-1 || finalStates.count(anTarget)) { // last layer: prune with accepting states
                    if (!nextTranslate.count(anTarget)) // the target is unknown and must be mapped
                        nextTranslateRef[anTarget] = unStateId++; 
                    // let's create arcs
                    unsigned int unTarget =  nextTranslateRef[anTarget]; // get the target's id in the new layer
                    Cost arcWeight = transition.second.second;
                    if (layer == 0) arcWeight += initialWeight;
                    if (layer == DACScope.size()-1 && finalStates.count(anTarget))
                        arcWeight += finalStates[anTarget];
                    Arc myArc(ucState,transition.first,transition.second.second,unTarget);
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
    for (unsigned int layer = DACScope.size()-1; layer >= 0; layer--) {

    }
    // Init conflict weights
    for (int i = 0; i != arity_; ++i) {
        conflictWeights.push_back(0);
    }
}

WRegular::~WRegular()
{
}

std::ostream& WRegular::printstate(std::ostream& os)
{
    os << "lb = " << lb << " connected = " << connected()
       << " depth = " << Store::getDepth() << "\n";
    for (int i = 0; i != arity_; ++i) {
        auto *x = scope[i];
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
    return  os;
}

void WRegular::propagate()
{
 
}

void WRegular::assign(int idx)
{
    static const bool debug{false};

    if (debug)
        cout << "After assign of " << idx << "\n";
}


void WRegular::remove(int idx)
{
    static const bool debug{false};

    if (debug)
        cout << "In remove " << idx << endl;

}

void WRegular::increase(int idx)
{
}

void WRegular::decrease(int idx)
{
}

void WRegular::projectFromZero(int idx)
{
    static const bool debug{false};

    if (!connected(idx))
        return;

    if (debug)
        cout << "In projectFromZero " << "\n";

}


double WRegular::computeTightness()
{
    return 1.0;
}

void WRegular::read(istream& is)
{

}

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

