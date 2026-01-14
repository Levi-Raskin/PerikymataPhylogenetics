#include "Msg.hpp"
#include "Node.hpp"

Node::Node(void) : ancestor(nullptr), name(""), branchLength(-1.0), index(-1), isTip(false){

}

void Node::addNeighbor(Node* p){
    auto it = std::find(neighbors.begin(), neighbors.end(), p);
    if (it == neighbors.end()) {
        neighbors.push_back(p);
    }else{
        Msg::error("tried to add neigboring node that is already a neighbor");
    }
}

std::vector<Node*> Node::getDescendants(void){
    std::vector<Node*> tmp;
    for(Node* n : neighbors)
        if(n != ancestor)
            tmp.push_back(n);
    return tmp;
}

void Node::removeNeighbor(Node* p){
    auto it = std::find(neighbors.begin(), neighbors.end(), p);
    if (it != neighbors.end()) {
        neighbors.erase(it);
    }
}
