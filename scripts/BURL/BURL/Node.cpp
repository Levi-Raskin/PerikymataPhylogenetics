#include "Node.hpp"
#include "RandomVariable.hpp"


Node::Node(void) : ancestor(nullptr), index(0), offset(0), isTip(false), name(""), flag(false), scratchInt(0) {
    
}


void Node::addNeighbor(Node* p){
    neighbors.insert(p);
    
    bool found = false;
    for(Node* n : p->getNeighbors())
        if(n == this)
            found = true;
    if(found == false)
        p->addNeighbor(this);
}
 
Node* Node::chooseNeighborAtRandom(RandomVariable* rng, Node* excludingNode) {

    std::vector<Node*> neighborVector;
    for (Node* p : neighbors)
        {
        if (p != excludingNode)
            neighborVector.push_back(p);
        }
    if (neighborVector.size() == 0)
        return nullptr;
    Node* p = neighborVector[(int)(rng->uniformRv()*neighborVector.size())];
    return p;
}

Node* Node::chooseNeighborAtRandom(RandomVariable* rng, std::set<Node*> excludedNodes) {

    std::vector<Node*> neighborVector;
    for (Node* p : neighbors)
        {
        std::set<Node*>::iterator it = excludedNodes.find(p);
        if (it == excludedNodes.end())
            neighborVector.push_back(p);
        }
    if (neighborVector.size() == 0)
        return nullptr;
    Node* p = neighborVector[(int)(rng->uniformRv()*neighborVector.size())];
    return p;

}

std::vector<Node*>& Node::getDescendants(void){
    descendantsVector.clear();
    for (Node* p : neighbors)
        {
            if (p != ancestor)
                descendantsVector.push_back(p);
        }
    return descendantsVector;
}

void Node::removeNeighbor(Node* p){
    neighbors.erase(p);
    bool found = false;
    for(Node* n : p->getNeighbors())
        if(n == this)
            found = true;
    if(found == true)
        p->removeNeighbor(this);
}
