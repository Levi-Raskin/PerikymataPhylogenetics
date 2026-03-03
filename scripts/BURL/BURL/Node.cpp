#include "Node.hpp"
#include "RandomVariable.hpp"


Node::Node(void) : ancestor(nullptr), index(0), offset(0), isTip(false), name(""), flag(false){
    
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
