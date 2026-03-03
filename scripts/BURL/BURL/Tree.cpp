#include "Eigen/Dense"
#include "Msg.hpp"
#include "Node.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "Tree.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

Tree::Tree(const Tree& t) {

    clone(t);
}

Tree::Tree(std::vector<std::string> taxonNames, double lambda) {

    numTaxa = (int)taxonNames.size();

    RandomVariable& rng = RandomVariable::randomVariableInstance();
    
    // start with three-species tree
    root = addNode();
    for (int i=0; i<3; i++)
        {
        Node* p = addNode();
        p->setIsTip(true);
        p->setIndex(i);
        p->setName(taxonNames[i]);
        p->addNeighbor(root);
        p->setAncestor(root);
        root->addNeighbor(p);
        }
        
    // randomly add the remaining taxa to the branches of the tree
    for (int i=3; i<numTaxa; i++)
        {
        Node* p = nullptr;
        do {
            double u = rng.uniformRv();
            int whichNode = (int)(u*nodes.size());
            p = nodes[whichNode];
            } while (p == root);
        Node* pAnc = p->getAncestor();
        p->removeNeighbor(pAnc);
        pAnc->removeNeighbor(p);

        Node* newTip = addNode();
        newTip->setIndex(i);
        newTip->setIsTip(true);
        newTip->setName(taxonNames[i]);
        Node* newInt = addNode();

        p->addNeighbor(newInt);
        p->setAncestor(newInt);
        pAnc->addNeighbor(newInt);
        newInt->addNeighbor(pAnc);
        newInt->addNeighbor(p);
        newInt->addNeighbor(newTip);
        newInt->setAncestor(pAnc);
        newTip->addNeighbor(newInt);
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
    for (Node* p : downPassSequence)
        {
        if (p != root)
            setBranch(p, p->getAncestor(), Probability::Exponential::rv(&rng, lambda));
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
                p->addNeighbor(newNode);
                newNode->addNeighbor(p);
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
                newNode->addNeighbor(p);
                newNode->setAncestor(p);
                p->addNeighbor(newNode);
                newNode->setName(token);
                newNode->setIsTip(true);
                tips.push_back(newNode);
                numTaxa++;
                p = newNode;
            }else{
                double x = stod(token);
                if(p->getAncestor() != nullptr)
                    setBranch(p, p->getAncestor(), x);
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
}

Tree::~Tree(void) {

    deleteNodes();
}

Tree& Tree::operator=(const Tree& t) {

    if (this != &t)
        clone(t);
    return *this;
}

Node* Tree::addNode(void) {

    Node* newNode = new Node;
    newNode->setOffset((int)nodes.size());
    nodes.push_back(newNode);
    return newNode;
}
 
void Tree::checkBranchLengthsNeg(void){
    for(auto const& x : branchLengths)
        if(x.second < 0)
            Msg::error("negative branch lengths");
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
        p->setFlag(q->getFlag());
        if (q->getAncestor() != nullptr)
            p->setAncestor( this->nodes[q->getAncestor()->getOffset()] );
        else
            p->setAncestor(nullptr);
        std::set<Node*>& qNeighbors = q->getNeighbors();
        p->removeAllNeighbors();
        for (Node* qn : qNeighbors)
            p->addNeighbor( this->nodes[qn->getOffset()] );
        }
        
    initializeDownPassSequence();
    this->branchLengths.clear();
    for (BranchLengths::const_iterator it = t.branchLengths.begin(); it != t.branchLengths.end(); it++)
        {
        if(it->first.first == nullptr)
            Msg::error("stop");
        Node* e1 = this->nodes[it->first.first->getOffset()];
        Node* e2 = this->nodes[it->first.second->getOffset()];
        setBranch(e1, e2, it->second);
        }
}

void Tree::collapseNode(Node* n){
    if(n->getNeighbors().size() != 2)
        Msg::error("cannot collapse this node because it does not have two neighbors");
        
//    std::cout << "deleting " << n->getIndex() << std::endl;
    
    std::vector<Node*> neighbors(n->getNeighbors().begin(), n->getNeighbors().end());
    double bl = getBranchLength(neighbors[0], n) + getBranchLength(neighbors[1], n);
    if(n == root){
        //need to select new root
        Node* foundNewRoot = nullptr;
        for(int i = 0; i < neighbors.size(); i++){
            if(neighbors[i]->getNeighbors().size() == 3){
                foundNewRoot = neighbors[i];
                neighbors.erase(neighbors.begin() + i);
            } //beacuse of how we set up the rot class, there should only be one node wiht this condition
        }
        
        if(neighbors.size() != 1)
            Msg::error("wonky neighbor sizing");
        if(foundNewRoot == nullptr)
            Msg::error("could not find suitable new root");
        
//        std::cout << "New root = " << foundNewRoot->getIndex() << std::endl;

        root = foundNewRoot;
        root->setAncestor(nullptr);
        Node* b = neighbors[0];
        b->addNeighbor(root);
        setBranch(root, b, bl);
        
        b->removeNeighbor(n);
        root->removeNeighbor(n);
        removeBranch(n, b);
        removeBranch(n, root);
        
    }else{
        neighbors[0]->addNeighbor(neighbors[1]);
        setBranch(neighbors[0], neighbors[1], bl);
        
        neighbors[0]->removeNeighbor(n);
        neighbors[1]->removeNeighbor(n);
        removeBranch(neighbors[0], n);
        removeBranch(neighbors[1], n);
    }
    
    n->removeAllNeighbors();
    n->setAncestor(nullptr);
    freeNode = n;
    freeNode->setIndex(-1);
//    for(int i = 0; i < nodes.size(); i++){
//        if(nodes[i] == n){
//            nodes.erase(nodes.begin() + i);
//            delete n;
//        }
//    }
    
    initializeDownPassSequence(); // resolves issues with DownPassSeq
}

void Tree::dropTip(std::string tip){
    Node* p = nullptr;
    for(Node* n : downPassSequence)
        if(n->getName() == tip)
            if(n->getIsTip() == true)
                p = n;
    if(p == nullptr)
        Msg::error("could not find tip in tree");
        
    // starting: (((p, q), a),b) p is tip to drop, a is int node to delete
    // ending:   ((q, some other tip that is already b descendent), b)
    Node* a = p->getAncestor();
    Node* b = a->getAncestor();
    Node* q = nullptr;
    std::vector<Node*> aDesc = a->getDescendants();
    if(aDesc[0] == p){
        q = aDesc[1];
    }else if (aDesc[1] == p){
        q = aDesc[0];
    }else{
        Msg::error("could not find tip in descendants");
    }
       
    if(a == root){
        //std::cout << "pongo pyg case" << std::endl;
        //q becomes new root
        removeBranch(q, a);
        q->removeNeighbor(a);
        root = q;
        initializeDownPassSequence();
    }else{
        double newBl = getBranchLength(q, a) + getBranchLength(a, b);
        //remove branches from branch map
        removeBranch(p, a);
        removeBranch(a, b);
        p->removeAllNeighbors();
        a->removeAllNeighbors();
        q->removeNeighbor(a);
        b->removeNeighbor(a);
        q->addNeighbor(b);
        q->setAncestor(b);
        b->addNeighbor(q);
        setBranch(q, b, newBl);
        initializeDownPassSequence();
    }
    
    
}

void Tree::dropTip(int idx){
    Node* p = nullptr;
    for(Node* n : downPassSequence)
        if(n->getIndex() == idx)
            if(n->getIsTip() == true)
                p = n;
    if(p == nullptr)
        Msg::error("could not find tip in tree");
        
    // starting: (((p, q), a),b) p is tip to drop, a is int node to delete
    // ending:   ((q, some other tip that is already b descendent), b)
    Node* a = p->getAncestor();
    Node* b = a->getAncestor();
    Node* q = nullptr;
    std::vector<Node*> aDesc = a->getDescendants();
    if(aDesc[0] == p){
        q = aDesc[1];
    }else if (aDesc[1] == p){
        q = aDesc[0];
    }else{
        Msg::error("could not find tip in descendants");
    }
       
    if(a == root){
        //std::cout << "pongo pyg case" << std::endl;
        //q becomes new root
        removeBranch(q, a);
        q->removeNeighbor(a);
        root = q;
        initializeDownPassSequence();
    }else{
        double newBl = getBranchLength(q, a) + getBranchLength(a, b);
        //remove branches from branch map
        removeBranch(p, a);
        removeBranch(a, b);
        p->removeAllNeighbors();
        a->removeAllNeighbors();
        q->removeNeighbor(a);
        b->removeNeighbor(a);
        q->addNeighbor(b);
        q->setAncestor(b);
        b->addNeighbor(q);
        setBranch(q, b, newBl);
        initializeDownPassSequence();
    }
}

void Tree::deleteNodes(void) {

    for (int i=0; i<nodes.size(); i++)
        delete nodes[i];
    nodes.clear();
}
    
std::vector<Node*> Tree::getAllDescendants(Node* n){
    std::vector<Node*> descendants;
    for (Node* child : n->getDescendants()) {
        descendants.push_back(child);
        if (!child->getIsTip()) {
            std::vector<Node*> childDescendants = getAllDescendants(child);
            descendants.insert(descendants.end(), childDescendants.begin(), childDescendants.end());
        }
    }
    return descendants;
}

void Tree::forceBinary(void){
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    std::vector<Node*> rootDesc = root->getDescendants();
    if(rootDesc.size() != 2){
        Node* newOutgroup = rootDesc[(int)(rng.uniformRv() * rootDesc.size())];
        Node* newRoot = addNode();
        newRoot->setIndex((int)nodes.size() + 1);
        double blHalf = getBranchLength(newOutgroup, root) / 2;
        
        removeBranch(newOutgroup, root);
        newOutgroup->removeNeighbor(root);
        
        newOutgroup->addNeighbor(newRoot);
        root->addNeighbor(newRoot);
        setBranch(newRoot, root, blHalf);
        setBranch(newOutgroup, newRoot, blHalf);
        reroot(newRoot);
    }
    
    for(Node* n : downPassSequence)
        if(n->getDescendants().size() != 2 && n->getIsTip() == false){
            print();
            Msg::error("Tree is not binary");
        }
        
}

double Tree::getBranchLength(Node* e1, Node* e2) {

    std::pair<Node*,Node*> key;
    initializeBranchLengthKey(key, e1, e2);
    BranchLengths::iterator it = branchLengths.find(key);
    if (it == branchLengths.end()){
        print("Tree with issue");
        std::cout << "for nodes: " << e1->getIndex() << " " << e2->getIndex() << std::endl;
        std::cout << "for nodes: " << e1 << " " << e2 << std::endl;
        std::cout << "e1 anc: " << e1->getAncestor()->getIndex() << std::endl;
        Msg::error("Could not find branch to return");
    }
    return it->second;
}

std::string Tree::getNewickString(void) {

    std::stringstream strm;
    writeTree(root, strm);
    strm << ";";
    return strm.str();
}

int Tree::getNumTaxa(void){
    initializeDownPassSequence();
    int idx = 0;
    for(Node* n : downPassSequence)
        if(n->getIsTip() == true)
            idx++;
    return idx;
}

Node* Tree::getTaxonNode(std::string name) {

    for (Node* p : downPassSequence)
        {
        if (p->getName() == name)
            return p;
        }
    return nullptr;
}

void Tree::initializeBranchLengthKey(std::pair<Node*,Node*>& key, Node* e1, Node* e2) {

    if (e1 < e2)
        key = std::make_pair(e1,e2);
    else
        key = std::make_pair(e2,e1);
}

void Tree::initializeDownPassSequence(void) {
    if(root == nullptr)
        Msg::error("root is nullptr");
    downPassSequence.clear();
    passDown(root, root);
    checkBranchLengthsNeg();
}

void Tree::keepTips(std::vector<std::string> t){
    //function drops all tips but those in the vector tips
    initializeDownPassSequence();
    
    for(Node* n : tips){
        //std::cout << n->getName() << std::endl;
        
        int cnt = (int)count(t.begin(), t.end(), n->getName());

        // Check if the target value was found
        if (cnt == 0)
            dropTip(n->getName());
        
    }
    
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

    if (p != nullptr)
        {
        std::set<Node*>& pNeighbors = p->getNeighbors();
        for (Node* d : pNeighbors)
            {
            if (d != from)
                passDown(d, p);
            }
        p->setAncestor(from);
        downPassSequence.push_back(p);
        }
}

void Tree::print(void) {
    
    showNode(root, 0);
    std::cout << "Postorder sequence: ";
    for (int i=0; i<downPassSequence.size(); i++)
        std::cout << downPassSequence[i]->getIndex() << " ";
    std::cout << std::endl;
    //printVCV();
}

void Tree::print(std::string header) {

    std::cout << header << std::endl;
    print();
}

void Tree::reindexNodes(void){
    int idx = 0;
    for (Node* p : downPassSequence)
        {
            if(p->getIsTip() == true)
                p->setIndex(idx++);
        }
    for (Node* p : downPassSequence)
        {
            if(p->getIsTip() == false)
                p->setIndex(idx++);
        }
}

void Tree::removeAllBranches(void) {

    branchLengths.clear();
}

void Tree::removeBranch(Node* e1, Node* e2) {

    std::pair<Node*,Node*> key;
    initializeBranchLengthKey(key, e1, e2);
    BranchLengths::iterator it = branchLengths.find(key);
    if (it != branchLengths.end())
        branchLengths.erase(it);
    else
        Msg::error("Could not find branch to remove");
}

void Tree::reroot(Node* r) {
    if(r->getIsTip() == true)
        Msg::error("why are you rerooting on a tip?");
    root = r;
    downPassSequence.clear();
    passDown(root, root);
}

void Tree::setBranch(Node* e1, Node* e2, double x) {
    std::pair<Node*,Node*> key;
    initializeBranchLengthKey(key, e1, e2);
    BranchLengths::iterator it = branchLengths.find(key);
    if (it != branchLengths.end())
        {
        it->second = x;
        }
    else
        {
        branchLengths.insert( std::make_pair(key, x) );
        }
}

void Tree::showNode(Node* p, int indent) {

    if (p != nullptr)
        {
        std::set<Node*>& pNeighbors = p->getNeighbors();
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
            std::cout << getBranchLength(p, p->getAncestor()) << " ";
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

void Tree::writeTree(Node* p, std::stringstream& strm) {

    if (p != nullptr)
        {
        std::set<Node*>& pNeighbors = p->getNeighbors();
        if (p->getIsTip() == false)
            strm << "(";
        else
            {
            std::string tName = p->getName();
            std::replace( tName.begin(), tName.end(), ' ', '_');
            strm << tName;
            }
        bool firstNode = true;
        for (Node* n : pNeighbors)
            {
            if (n != p->getAncestor())
                {
                if (firstNode == false)
                    strm << ",";
                writeTree(n, strm);
                strm << ":" << getBranchLength(p, n);
                firstNode = false;
                }
            }
        if (p->getIsTip() == false)
            strm << ")";

        }
}
