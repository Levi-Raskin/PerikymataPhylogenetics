#include "Msg.hpp"
#include "Node.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "Tree.hpp"

#include <iomanip>
#include <iostream>

# pragma mark Constructors

Tree::Tree(int nt, double lambda) : numTaxa(nt){
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    
    std::vector<std::string> taxonNames;
    for(int i = 0; i < numTaxa; i++)
        taxonNames.push_back("t"+ std::to_string(i));
    root = addNode();
    
    for (int i=0; i<2; i++){
        Node* p = addNode();
        p->setIsTip(true);
        p->setIndex(i);
        p->setName(taxonNames[i]);
        setNeighbors(p, root);
        p->setAncestor(root);
    }
        
    // randomly add the remaining taxa to the branches of the tree
    for (int i=2; i<numTaxa; i++){
        Node* p = nullptr;
        do{
            double u = rng.uniformRv();
            int whichNode = (int)(u*nodes.size());
            p = nodes[whichNode];
        }while (p == root);
            
        Node* pAnc = p->getAncestor();
        removeAsNeighbors(p, pAnc);

        Node* newTip = addNode();
        newTip->setIndex(i);
        newTip->setIsTip(true);
        newTip->setName(taxonNames[i]);
        Node* newInt = addNode();

        setNeighbors(p, newInt);
        p->setAncestor(newInt);
        
        setNeighbors(pAnc, newInt);
        newInt->setAncestor(pAnc);
        
        setNeighbors(newInt, newTip);
        newTip->setAncestor(newInt);
    }
        
    // initialize the post-order traversal sequence
    initializeDownPassSequence();
    
    // index interior nodes
    int idx = numTaxa;
    for (Node* p : downPassSequence)
        {
        if (p->getIsTip() == false)
            p->setIndex(idx++);
        }
    
    // add the branch lengths to the tree
    for (Node* p : downPassSequence){
        if (p != root)
            p->setBranchLength(Probability::Exponential::rv(&rng, lambda));
    }
}

Tree::~Tree(void){
    for (int i=0; i<nodes.size(); i++)
        delete nodes[i];
    nodes.clear();
}

# pragma mark Methods


Node* Tree::addNode(void){
    Node* newNode = new Node;
    nodes.push_back(newNode);
    return newNode;
}

void Tree::initializeDownPassSequence(void){
    if(root == nullptr)
        Msg::error("root is nullptr");
    downPassSequence.clear();
    passDown(root, root);
}

void Tree::passDown(Node* p, Node* from) {

    if (p != nullptr){
        std::vector<Node*> pNeighbors = p->getNeighbors();
        for (Node* d : pNeighbors)
            {
            if (d != from)
                passDown(d, p);
            }
        p->setAncestor(from);
        downPassSequence.push_back(p);
    }
}

void Tree::print(void){
    showNode(root, 0);
    std::cout << "Postorder sequence: ";
    for (int i=0; i<downPassSequence.size(); i++)
        std::cout << downPassSequence[i]->getIndex() << " ";
    std::cout << std::endl;
}

void Tree::removeAsNeighbors(Node* p, Node* q){
    q->removeNeighbor(p);
    p->removeNeighbor(q);
}

void Tree::setNeighbors(Node* p, Node* q){
    p->addNeighbor(q);
    q->addNeighbor(p);
}

void Tree::showNode(Node* p, int indent) {

    if (p != nullptr)
        {
        std::vector<Node*> pNeighbors = p->getNeighbors();
        for (int i=0; i<indent; i++)
            std::cout << " ";
        std::cout << p->getIndex() << " [" << p << "] ( ";
        for (Node* n : pNeighbors)
            {
            if (n == p->getAncestor())
                std::cout << "a_";
            std::cout << n->getIndex() << " ";
            }
        std::cout << ") ";
        std::cout << p->getName() << " " << p->getIsTip() << " ";
        std::cout << std::fixed << std::setprecision(6);
        if (p != root)
            std::cout << p->getBranchLength() << " ";
        else
            std::cout << "--- ";
        if (p == root)
            std::cout << "<- Root ";
        std::cout << std::endl;

        for (Node* n : pNeighbors)
            {
            if (n != p->getAncestor())
                showNode(n, indent + 3);
            }
            
        }
}
