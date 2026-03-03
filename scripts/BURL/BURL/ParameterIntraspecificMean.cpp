#include "ParameterIntraspecificMean.hpp"
#include "PhylogeneticModel.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include <iostream>
#include <set>

ParameterIntraspecificMean::ParameterIntraspecificMean(double prob, std::string n, Eigen::MatrixXd* data, PhylogeneticModel* p) : Parameter(prob, n), nTraits(data->cols()), nObs(data->rows()), tipData(data), windowSize(0.5), stepSize(1e-3), model(p), targetAcceptanceRate(0.43), upperAcceptanceRate(targetAcceptanceRate + 0.01),lowerAcceptanceRate(targetAcceptanceRate-0.01), ngAdaptive(50000), covarianceUpdateFreq(100), useEmpiricalCovariance(false), usePosteriorDraw(false), numAcceptances(0), numRejections(0){

    RandomVariable& rng = RandomVariable::randomVariableInstance();

    //mean is unconstrained pos/neg
    //prior
    mu0 = Eigen::VectorXd::Ones(nTraits);
    for(int i = 0; i < nTraits; i++)
        mu0(i) = Probability::Normal::rv(&rng, 1.0, 1.0);

    mean.push_back(tipData->colwise().mean());
    mean.push_back(tipData->colwise().mean());
    adaptiveProposalActive = true;
    
    dof = nTraits + 2; //such that the mean is the scale matrix
    lambda = dof - nTraits - 1;

    mu0 = Eigen::VectorXd::Zero(nTraits);
    
    psi = Eigen::MatrixXd::Constant(nTraits, nTraits, 1e-6);
    psi.diagonal().array() = 1.0;
    psi = Eigen::VectorXd::Ones(nTraits).asDiagonal() * psi * Eigen::VectorXd::Ones(nTraits).asDiagonal();
    
    
    empiricalMean = Eigen::VectorXd::Zero(nTraits);
    empiricalCovariance = Eigen::MatrixXd::Identity(nTraits, nTraits);
    
    proposalCov.resize(nTraits, nTraits);
    
    // Preconditioned RW initialization
    cholSigma = Eigen::MatrixXd::Identity(nTraits, nTraits);
    zDraw     = Eigen::VectorXd::Zero(nTraits);
    // Roberts-Gelman-Gilks optimal scaling: 2.38/sqrt(q), further scaled by 1/sqrt(nObs)
    rwScaleFactor = 2.38 / std::sqrt((double)nTraits * (double)nObs);
    
    // Data-likelihood MALA initialization
    dataMean   = tipData->colwise().mean();  // x̄ᵢ
    driftOld   = Eigen::VectorXd::Zero(nTraits);
    driftNew   = Eigen::VectorXd::Zero(nTraits);
    residFwd   = Eigen::VectorXd::Zero(nTraits);
    residRev   = Eigen::VectorXd::Zero(nTraits);
    // Optimal MALA scaling: ε ~ q^{-1/3} (Roberts & Rosenthal 1998)
    // Further scaled by 1/sqrt(nObs) to match posterior width
    malaEpsilon    = 1.0 / (std::pow((double)nTraits, 1.0/3.0) * std::sqrt((double)nObs));
    malaDriftCoeff = 0.5 * malaEpsilon * malaEpsilon * (double)nObs;
}

double ParameterIntraspecificMean::lnProbability(void){
    //not needed; accounted for by multivariateBrownain Motion
    return 0.0;
}

void ParameterIntraspecificMean::print(void){
    std::cout << " " << std::endl;
}

double ParameterIntraspecificMean::update(void){
    return updateMHmvNDraw();
}

double ParameterIntraspecificMean::updateConditionalDraw(void){
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    const size_t numAccepted = std::count(recentAcceptRej.begin(),
                                      recentAcceptRej.end(), true);
    const double acceptanceRate = static_cast<double>(numAccepted) / recentAcceptRej.size();

    // Adaptive step size adjustment
    const size_t totalSamples = numRejections + numAcceptances;
    if (totalSamples % 50 == 0 && totalSamples < ngAdaptive) {
        if (acceptanceRate < lowerAcceptanceRate)
            windowSize /= 1.05;
        else if (acceptanceRate > upperAcceptanceRate)
            windowSize *= 1.05;
            
//        std::cout << windowSize << "\n";
    }

    Eigen::VectorXd tipMean = tipData->colwise().mean();
    
    static double lambdaN = lambda + nObs;
    Eigen::VectorXd uN = (mu0 * lambda + nObs * (tipMean) )/ ( lambdaN) ;
    
    Eigen::MatrixXd vcv = varianceCovariance->getValue();
    Eigen::MatrixXd vcvMod = (1/lambdaN)* vcv;
    
    static int elementsUpdated = 1;
    static int elementsNotUpdated = nTraits - elementsUpdated;
    
    // Resize matrices only when necessary (reuse memory)
    Eigen::MatrixXd sigma11;
    Eigen::MatrixXd sigma12;
    Eigen::MatrixXd sigma21;
    Eigen::MatrixXd sigma22;
    sigma11.resize(elementsUpdated, elementsUpdated);
    sigma12.resize(elementsUpdated, elementsNotUpdated);
    sigma21.resize(elementsNotUpdated, elementsUpdated);
    sigma22.resize(elementsNotUpdated, elementsNotUpdated);
    
    Eigen::VectorXd u1;
    Eigen::VectorXd u2;
    Eigen::VectorXd x2;
    u1.resize(elementsUpdated);
    u2.resize(elementsNotUpdated);
    x2.resize(elementsNotUpdated);
    
    // Randomly select which elements to update
    std::set<int> updatingElements;
    do{
        updatingElements.insert((int)(rng.uniformRv() * nTraits));
    } while(updatingElements.size() < elementsUpdated);
    
    // Build vector of non-updating elements
    std::vector<int> updatingVec(updatingElements.begin(), updatingElements.end());
    std::vector<int> nonUpdatingElements;
    for(int i = 0; i < nTraits; i++) {
        if(updatingElements.find(i) == updatingElements.end()) {
            nonUpdatingElements.push_back(i);
        }
    }
    
    // Build conditional covariance matrices using block operations
    for(int i = 0; i < elementsUpdated; i++) {
        int updateIdx = updatingVec[i];
        u1(i) = uN(updateIdx);
        x2(i) = mean[0](updateIdx);
        
        for(int j = 0; j < elementsUpdated; j++) {
            sigma11(i, j) = vcvMod(updateIdx, updatingVec[j]);
        }
        for(int j = 0; j < elementsNotUpdated; j++) {
            sigma12(i, j) = vcvMod(updateIdx, nonUpdatingElements[j]);
        }
    }
    
    for(int i = 0; i < elementsNotUpdated; i++) {
        int nonUpdateIdx = nonUpdatingElements[i];
        u2(i) = uN(nonUpdateIdx);
        
        for(int j = 0; j < elementsUpdated; j++) {
            sigma21(i, j) = vcvMod(nonUpdateIdx, updatingVec[j]);
        }
        for(int j = 0; j < elementsNotUpdated; j++) {
            sigma22(i, j) = vcvMod(nonUpdateIdx, nonUpdatingElements[j]);
        }
    }
    
    // Use LLT decomposition instead of inverse (faster and more stable)
    Eigen::LLT<Eigen::MatrixXd> llt(sigma22);
    
    // Compute conditional distribution parameters
    Eigen::MatrixXd sigma22Inv_sigma21 = llt.solve(sigma21);
    Eigen::MatrixXd sigmaCond = sigma11 - sigma12 * sigma22Inv_sigma21;
    
    // Conditional mean: u1 + sigma12 * sigma22^(-1) * (x2_cond - u2)
    Eigen::VectorXd x2_cond;
    x2_cond.resize(elementsNotUpdated);
    for(int i = 0; i < elementsNotUpdated; i++) {
        x2_cond(i) = mean[0](nonUpdatingElements[i]);
    }
    
    Eigen::VectorXd x2_minus_u2 = x2_cond - u2;
    Eigen::VectorXd sigma22Inv_diff = llt.solve(x2_minus_u2);
    Eigen::VectorXd uCond = u1 + sigma12 * sigma22Inv_diff;
    
    sigmaCond *= windowSize;
    
    // Sample new values from conditional distribution
    Eigen::VectorXd newVals = Probability::MultivariateNormal::rv(&rng, uCond, &sigmaCond);
    
    // Update the parameter values in mean[0]
    for(int idx = 0; idx < elementsUpdated; idx++) {
        mean[0](updatingVec[idx]) = newVals(idx);
    }
    
    return 0.0;
}

double ParameterIntraspecificMean::updatePosteriorDraw(void){
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    Eigen::VectorXd tipMean = tipData->colwise().mean();
    
    static double lambdaN = lambda + nObs;
    Eigen::VectorXd uN = (mu0 * lambda + nObs * (tipMean) )/ ( lambdaN) ;
    
    Eigen::MatrixXd vcv = varianceCovariance->getValue();
    
    //draw directly from the joint posterior; not a true gibbs sampler but a high-acceptance probability move
    Eigen::MatrixXd vcvScaled = vcv / lambdaN;
    mean[0] = Probability::MultivariateNormal::rv(&rng, uN, &vcvScaled);
    return 0.0;
}

double ParameterIntraspecificMean::updateHMC(void){
    RandomVariable& rng = RandomVariable::randomVariableInstance();
//    t.end();
    double a = 0.0;
    for(bool b : recentAcceptRej)
        if(b == true)
            a++;
    a /= recentAcceptRej.size();
    
    if((numRejections + numAcceptances) % 50 ==0 && ((numRejections + numAcceptances) <  ngAdaptive)){
            if(a < 0.79)
                stepSize /= 1.05;
            else if (a > 0.81)
                stepSize *= 1.05;
    }else if ((numRejections + numAcceptances) == ngAdaptive){
        std::cout << parmName << " done adaptive sampling | final acceptRej: " << a << std::endl;
        adaptiveProposalActive = false;
    }
    // ---- 1. Flatten parameters: rates + Cholesky off-diagonals ----

    // ---- 2. Define log probability ----
    auto logProb = [&]() -> double {
        double lp = model->lnPriorProbability();
        lp += model->lnLikelihood();
        return lp;
    };
    
    auto numericalGradient = [&]() {
        constexpr double jiggle = 1e-8;
        constexpr double inv_eps = 1.0 / jiggle;
        double f0 = logProb();
        Eigen::VectorXd grad(nTraits);
        Eigen::VectorXd orig = mean[0];
        for (int i = 0; i < nTraits; ++i) {
            mean[0](i) += jiggle;
            double f1 = logProb();
            grad(i) = (f1 - f0) * inv_eps;
        }
        mean[0] = orig;
        return grad;
    };

    // ---- 4. Draw momenta ----
    Eigen::VectorXd momentum = Eigen::VectorXd::NullaryExpr(nTraits, [&]() {
        return Probability::Normal::rv(&rng);
    });
    
    

    // ---- 5. Compute initial energy ----
    Eigen::VectorXd gradU0 = -numericalGradient();

    // ---- 6. Leapfrog integration ----
    double eps = stepSize;
        
    constexpr int L = 10;
    Eigen::VectorXd p = momentum;

    p -= 0.5 * eps * gradU0;
    for (int i = 0; i < L; ++i) {
        mean[0]  += eps * p;
        Eigen::VectorXd gradU = -numericalGradient();
        if (i != L - 1)
            p -= eps * gradU;
        else
            p -= 0.5 * eps * gradU;
    }
    return 0.0;
}

double ParameterIntraspecificMean::updateHMCSingleElement(void){
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    double a = 0.0;
    for(bool b : recentAcceptRej)
        if(b == true)
            a++;
    a /= recentAcceptRej.size();
    
    if((numRejections + numAcceptances) % 50 ==0 && ((numRejections + numAcceptances) <  ngAdaptive)){
            if(a < 0.75)
                stepSize /= 1.05;
            else if (a > 0.85)
                stepSize *= 1.05;
    }else if ((numRejections + numAcceptances) == ngAdaptive){
        std::cout << parmName << " done adaptive sampling | final acceptRej: " << a << std::endl;
        adaptiveProposalActive = false;
    }


    // choose idx to update
    
    int elementToUpdate = (int)(rng.uniformRv() * nTraits);
    
    double momentum = Probability::Normal::rv(&rng);

    //initial energy
    double grad = -numericalGradientSingleElement(elementToUpdate);

    double eps = stepSize;
        
    constexpr int L = 20;
    double p = momentum;

    p -= 0.5 * eps * grad;
    for (int i = 0; i < L; ++i) {
        mean[0](elementToUpdate)  += eps * p;
        double gradU = -numericalGradientSingleElement(elementToUpdate);
        if (i != L - 1)
            p -= eps * gradU;
        else
            p -= 0.5 * eps * gradU;
    }
    return 0.0;
}

double ParameterIntraspecificMean::calculatePosterior(void){
    return model->lnPriorProbability() + model->lnLikelihood();
}

double ParameterIntraspecificMean::numericalGradientSingleElement(int idx) {
    constexpr double jiggle = 1e-3;
    constexpr double inv_eps = 1.0 / jiggle;
    
    double f0 = calculatePosterior();
    mean[0](idx) += jiggle;
    double f1 = calculatePosterior();
    mean[0](idx) -= jiggle;
        
    return (f1 - f0) * inv_eps;
};

double ParameterIntraspecificMean::updateMHmvNDraw(void){
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    const size_t totalSamples = numRejections + numAcceptances;
    
    if (totalSamples % 50 == 0 && totalSamples < ngAdaptive) {
    const size_t numAccepted = std::count(recentAcceptRej.begin(), recentAcceptRej.end(), true);
        const double acceptanceRate = static_cast<double>(numAccepted) / recentAcceptRej.size();
        if (acceptanceRate < lowerAcceptanceRate)
            windowSize /= 1.05;
        else if (acceptanceRate > upperAcceptanceRate)
            windowSize *= 1.05;
//        std::cout << windowSize << "\n";
    }else if ((numRejections + numAcceptances) == ngAdaptive){
        adaptiveProposalActive = false;
    }
    
     if (totalSamples % covarianceUpdateFreq == 0 && totalSamples < ngAdaptive) {
        updateEmpiricalCovariance();
    }
    
    if (totalSamples == 5000 && !useEmpiricalCovariance) {
        useEmpiricalCovariance = true;
    }
    
    // Choose proposal covariance
    if(totalSamples < ngAdaptive){
        if (useEmpiricalCovariance && recentSamples.size() > nTraits * 2) {
            proposalCov = windowSize * empiricalCovariance;
            proposalCholLower = proposalCov.llt().matrixL().toDenseMatrix();
        } else {
            proposalCov = windowSize * Eigen::MatrixXd::Identity(nTraits, nTraits);
            proposalCholLower = proposalCov.llt().matrixL().toDenseMatrix();
        }
    }
    
    // Propose
    Probability::MultivariateNormal::rv(&rng, mean[0], mean[1], proposalCholLower);
    return 0.0;

}

double ParameterIntraspecificMean::updateMHSingleElement(void){
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    const size_t totalSamples = numRejections + numAcceptances;
    
    if (totalSamples % 50 == 0 && totalSamples < ngAdaptive) {
    const size_t numAccepted = std::count(recentAcceptRej.begin(), recentAcceptRej.end(), true);
        const double acceptanceRate = static_cast<double>(numAccepted) / recentAcceptRej.size();
        if (acceptanceRate < lowerAcceptanceRate)
            windowSize /= 1.05;
        else if (acceptanceRate > upperAcceptanceRate)
            windowSize *= 1.05;
//        std::cout << windowSize << "\n";
    }else if ((numRejections + numAcceptances) == ngAdaptive){
        adaptiveProposalActive = false;
    }
    
    int idx = (int)(rng.uniformRv() * nTraits);
    mean[0](idx) = Probability::Uniform::rv(&rng, mean[0](idx) - windowSize, mean[0](idx) + windowSize);
    return 0.0;
}


double ParameterIntraspecificMean::updatePreconditionedRW(void){
    RandomVariable& rng = RandomVariable::randomVariableInstance();

    // ---- Adaptive tuning (same pattern as updateMHSingleElement) ----
    const size_t totalSamples = numRejections + numAcceptances;

    if (totalSamples % 50 == 0 && totalSamples < ngAdaptive) {
        const size_t numAccepted = std::count(recentAcceptRej.begin(),
                                              recentAcceptRej.end(), true);
        const double acceptanceRate =
            recentAcceptRej.empty()
                ? targetAcceptanceRate
                : static_cast<double>(numAccepted) / recentAcceptRej.size();

        if (acceptanceRate < lowerAcceptanceRate)
            rwScaleFactor /= 1.05;
        else if (acceptanceRate > upperAcceptanceRate)
            rwScaleFactor *= 1.05;
    } else if (totalSamples == ngAdaptive) {
        adaptiveProposalActive = false;
    }

    // ---- Refresh Cholesky of Σᵢ ----
    // Σᵢ is updated by ParameterIntraspecificVariance, so we recompute each
    // time. Cost is O(q³) = O(1000) for q=10, negligible vs. pruning.
    {
        const Eigen::MatrixXd& sigma = varianceCovariance->getValue();
        Eigen::LLT<Eigen::MatrixXd> llt(sigma);
        if (llt.info() == Eigen::Success) {
            cholSigma = llt.matrixL();
        }
        // On decomposition failure (shouldn't happen), keep previous cholSigma
    }

    // ---- Draw z ~ N(0, I_q) ----
    for (int k = 0; k < nTraits; ++k)
        zDraw(k) = Probability::Normal::rv(&rng);

    // ---- Propose: µ* = µ + λ · L · z  where L Lᵀ = Σᵢ ----
    // Equivalent to proposing from N(µ, λ² Σᵢ): symmetric, so HR = 0
    mean[0] += rwScaleFactor * (cholSigma * zDraw);

    return 0.0;  // symmetric proposal
}

double ParameterIntraspecificMean::updateDataLikelihoodMALA(void){
    /*
     * Data-Likelihood-Gradient MALA (Metropolis-Adjusted Langevin Algorithm)
     *
     * Uses the analytically available gradient of the data likelihood:
     *   ∇_µᵢ log Pr(Xᵢ | µᵢ, Σᵢ) = nᵢ · Σᵢ⁻¹ · (x̄ᵢ - µᵢ)
     *
     * Preconditioned by Σᵢ (simplified manifold MALA), the drift simplifies:
     *   h(µ) = (ε²/2) · Σᵢ · [nᵢ Σᵢ⁻¹ (x̄ᵢ - µ)] = (ε²/2) · nᵢ · (x̄ᵢ - µ)
     *
     * Proposal:  µ* = µ + h(µ) + ε · Lᵢ · z,   z ~ N(0,I)
     * Density:   q(µ* | µ) = N(µ*;  µ + h(µ),  ε² Σᵢ)
     *
     * Asymmetric → must return log q(µ_old | µ_new) - log q(µ_new | µ_old).
     */
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    
    // ---- Adaptive tuning ----
    const size_t totalSamples = numRejections + numAcceptances;
    
    if (totalSamples % 50 == 0 && totalSamples < ngAdaptive) {
        const size_t numAccepted = std::count(recentAcceptRej.begin(),
                                              recentAcceptRej.end(), true);
        const double acceptanceRate =
            recentAcceptRej.empty()
                ? targetAcceptanceRate
                : static_cast<double>(numAccepted) / recentAcceptRej.size();
        
        // MALA optimal acceptance ~57.4% (Roberts & Rosenthal 1998)
        // We use the class target rate (0.43) which is a reasonable compromise
        // since this is mixed with the preconditioned RW
        if (acceptanceRate < lowerAcceptanceRate)
            malaEpsilon /= 1.05;
        else if (acceptanceRate > upperAcceptanceRate)
            malaEpsilon *= 1.05;
        
        malaDriftCoeff = 0.5 * malaEpsilon * malaEpsilon * (double)nObs;
    } else if (totalSamples == ngAdaptive) {
        adaptiveProposalActive = false;
    }
    
    // ---- Refresh Cholesky of Σᵢ and data mean ----
    {
        const Eigen::MatrixXd& sigma = varianceCovariance->getValue();
        Eigen::LLT<Eigen::MatrixXd> llt(sigma);
        if (llt.info() == Eigen::Success) {
            cholSigma = llt.matrixL();
        }
    }
    dataMean = tipData->colwise().mean();  // x̄ᵢ (recompute in case of imputation)
    
    // ---- Save current state ----
    const Eigen::VectorXd muOld = mean[0];  // copy for Hastings ratio
    
    // ---- Compute drift at current state: h(µ_old) = malaDriftCoeff · (x̄ᵢ - µ_old) ----
    driftOld = malaDriftCoeff * (dataMean - muOld);
    
    // ---- Draw z ~ N(0, I_q) ----
    for (int k = 0; k < nTraits; ++k)
        zDraw(k) = Probability::Normal::rv(&rng);
    
    // ---- Propose: µ* = µ_old + h(µ_old) + ε · L · z ----
    mean[0] = muOld + driftOld + malaEpsilon * (cholSigma * zDraw);
    
    const Eigen::VectorXd& muNew = mean[0];
    
    // ---- Compute drift at proposed state: h(µ_new) = malaDriftCoeff · (x̄ᵢ - µ_new) ----
    driftNew = malaDriftCoeff * (dataMean - muNew);
    
    // ---- Log Hastings ratio ----
    //
    // q(µ* | µ) = N(µ*; µ + h(µ), ε²Σᵢ)
    //
    // log q(a | b) = -1/(2ε²) · (a - b - h(b))ᵀ Σᵢ⁻¹ (a - b - h(b))  + const
    //
    // log HR = log q(µ_old | µ_new) - log q(µ_new | µ_old)
    //        = -1/(2ε²) · [ |µ_old - µ_new - h(µ_new)|²_{Σ⁻¹}
    //                       -|µ_new - µ_old - h(µ_old)|²_{Σ⁻¹} ]
    //
    // The forward residual (µ_new - µ_old - h(µ_old)) is just ε·L·z, but
    // we compute both explicitly for clarity and numerical safety.
    
    double invEps2x2 = 1.0 / (2.0 * malaEpsilon * malaEpsilon);
    
    // Forward residual: µ_new - µ_old - h(µ_old)
    residFwd = muNew - muOld - driftOld;
    
    // Reverse residual: µ_old - µ_new - h(µ_new)
    residRev = muOld - muNew - driftNew;
    
    // Mahalanobis distances using the Cholesky solve: xᵀ Σ⁻¹ x = (L⁻¹x)ᵀ(L⁻¹x)
    Eigen::VectorXd solvedFwd = cholSigma.triangularView<Eigen::Lower>().solve(residFwd);
    Eigen::VectorXd solvedRev = cholSigma.triangularView<Eigen::Lower>().solve(residRev);
    
    double mahalFwd = solvedFwd.squaredNorm();
    double mahalRev = solvedRev.squaredNorm();
    
    double logHR = invEps2x2 * (mahalFwd - mahalRev);
    
    return logHR;
}

void ParameterIntraspecificMean::updateEmpiricalCovariance(void){
    if (recentSamples.size() < nTraits + 1) return;
    
    // Compute empirical mean
    empiricalMean = Eigen::VectorXd::Zero(nTraits);
    for (const auto& sample : recentSamples) {
        empiricalMean += sample;
    }
    empiricalMean /= recentSamples.size();
    
    // Compute empirical covariance
    empiricalCovariance = Eigen::MatrixXd::Zero(nTraits, nTraits);
    for (const auto& sample : recentSamples) {
        Eigen::VectorXd centered = sample - empiricalMean;
        empiricalCovariance += centered * centered.transpose();
    }
    empiricalCovariance /= (recentSamples.size() - 1);
    
    // Add small diagonal for numerical stability
    empiricalCovariance += 1e-6 * Eigen::MatrixXd::Identity(nTraits, nTraits);
}

void ParameterIntraspecificMean::updateForAcceptance(void){
    mean[1] = mean[0];
    numAcceptances++;
    recentAcceptRej.push_back(true);
    if(recentAcceptRej.size() > 1000)
        recentAcceptRej.pop_front();
    recentSamples.push_back(mean[0]);
    if (recentSamples.size() > 1000)  // Keep last 500 samples
        recentSamples.pop_front();
}

void ParameterIntraspecificMean::updateForRejection(void){
    mean[0] = mean[1];
    numRejections++;
    recentAcceptRej.push_back(false);
    if(recentAcceptRej.size() > 1000)
        recentAcceptRej.pop_front();
}
