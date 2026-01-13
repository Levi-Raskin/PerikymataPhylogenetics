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
    setOffsets();
}

Tree::Tree(std::vector<std::string> taxonNames, double lambda) : Tree(taxonNames.size(), lambda) {
    int i = 0;
    for (Node* p : downPassSequence) {
        if (p->getIsTip() && i < taxonNames.size()) {
            p->setName(taxonNames[i++]);
        }
    }
}

Tree::Tree(std::string newick){
    numTaxa = 0;
    std::vector<std::string> newickTokens = parseNewickString(newick);
    Node* p = nullptr;
    bool readingBl = false;
    for(int i = 0; i < newickTokens.size(); i++){
        std::string token = newickTokens[i];
        if(token == "("){
            Node* newNode = addNode();
            if(p == nullptr){
                root = newNode;
            }else{
                setNeighbors(p, newNode);
                newNode->setAncestor(p);
            }
            p = newNode;
        }else if (token == ")" || token == ","){
            if(p->getAncestor() == nullptr)
                Msg::error("no anc found for p");
            p = p->getAncestor();
        }else if (token == ";"){
            if(p != root)
                Msg::error("expecting to be at root");
        }else if (token == ":"){
            readingBl = true;
        }else{
            if(readingBl == false){
                Node* newNode = addNode();
                
                setNeighbors(newNode, p);
                newNode->setAncestor(p);
                newNode->setName(token);
                newNode->setIsTip(true);
                numTaxa++;
                p = newNode;
            }else{
                double x = stod(token);
                if(p->getAncestor() != nullptr)
                    p->setBranchLength(x);
                readingBl = false;
            }
        }
    }
    initializeDownPassSequence();
    int idx = numTaxa;
    int tIdx = 0;
    for (Node* p : downPassSequence)
        {
            if(p->getIsTip() == true)
                p->setIndex(tIdx++);
            if (p->getIsTip() == false)
                p->setIndex(idx++);
        }
    setOffsets();
}

Tree::Tree(const Tree& t){
    clone(t);
}

Tree::~Tree(void){
    deleteNodes();
}

# pragma mark Methods
Tree& Tree::operator=(const Tree& t){
    if (this != &t)
        clone(t);
    return *this;
}

Node* Tree::addNode(void){
    Node* newNode = new Node;
    nodes.push_back(newNode);
    return newNode;
}

void Tree::clone(const Tree& t) {
    
    if (this->nodes.size() != t.nodes.size())
        {
        deleteNodes();
        for (int i=0; i<t.nodes.size(); i++)
            addNode();
        }
        
    this->numTaxa = t.numTaxa;
    this->root = this->nodes[t.root->getOffset()];
    
    for (int i=0; i<t.nodes.size(); i++)
        {
        Node* q = t.nodes[i];
        Node* p = this->nodes[i];
        p->setIndex(q->getIndex());
        p->setIsTip(q->getIsTip());
        p->setName(q->getName());
        p->setBranchLength(q->getBranchLength());
        if (q->getAncestor() != nullptr)
            p->setAncestor( this->nodes[q->getAncestor()->getOffset()] );
        else
            p->setAncestor(nullptr);
        std::vector<Node*> qNeighbors = q->getNeighbors();
        p->removeAllNeighbors();
        for (Node* qn : qNeighbors)
            p->addNeighbor( this->nodes[qn->getOffset()] );
        }
        
    initializeDownPassSequence();
}

void Tree::deleteNodes(void){
    for (int i=0; i<nodes.size(); i++)
        delete nodes[i];
    nodes.clear();
}

void Tree::initializeDownPassSequence(void){
    if(root == nullptr)
        Msg::error("root is nullptr");
    downPassSequence.clear();
    passDown(root, root);
}

std::vector<std::string>  Tree::parseNewickString(std::string newickStr){
    
    std::vector<std::string> tokens;
    
    std::string token = "";
    for(int i = 0; i < newickStr.length(); i++){
        char c = newickStr[i];
        if(c == '(' || c == ')' || c==',' || c==':' ||c == ';'){
            if(token != ""){
                tokens.push_back(token);
                token = "";
            }
            tokens.push_back(std::string(1,c));
        }else{
            token += std::string(1,c);
        }
    }
    
    if(token != "")
        tokens.push_back(token);
    
//    for(int i = 0; i < tokens.size(); i++){
//        std::cout << i << " " << tokens[i] << std::endl;
//    }
//        
    
    return tokens;
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

void Tree::setOffsets(void){
    for(int i = 0; i < nodes.size(); i++)
        nodes[i]->setOffset(i);
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
