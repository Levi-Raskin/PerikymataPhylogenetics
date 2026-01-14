#include "MultivariateBrownianMotion.hpp"
#include "ReadCSV.hpp"
#include "TicToc.hpp"
#include "Tree.hpp"
#include <iostream>

int main(int argc, const char * argv[]) {
    ReadCSV r = ReadCSV("/Users/levir/Documents/GitHub/Raskin_et_al_perikymata_hsp/LC.csv", true, true);
    
    std::vector<std::string> rawReadDatNames = r.getRownames();
    for(int i = 0 ; i < rawReadDatNames.size(); i++){
        std::string trimmed_name = rawReadDatNames[i];
        trimmed_name.erase(0, trimmed_name.find_first_not_of('"'));
        trimmed_name.erase(trimmed_name.find_last_not_of('"') + 1);
        rawReadDatNames[i] =trimmed_name;
    }
    Eigen::MatrixXd readDat = r.getEigenMat();
        
    Tree gatree = Tree("(((Gorilla_beringei:2.558516,Gorilla_gorilla:2.558516):6.093717,((Homo_sapiens:0.568721,Neanderthal:0.538721):5.607159,(Pan_paniscus:2.333553,Pan_troglodytes:2.333553):3.842326):2.476353):6.480222,(Pongo_abelii:3.825854,Pongo_pygmaeus:3.825854):11.306601);");
    
    gatree.print();
    
    MultivariateBrownianMotion ctmc = MultivariateBrownianMotion(rawReadDatNames, &readDat);
    ctmc.setTree(&gatree);

    TicToc t("likelihood");
    double lnl = ctmc.lnLikelihood();
    t.end();
    return 0;
}
