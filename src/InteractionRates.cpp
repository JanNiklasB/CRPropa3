#include "crpropa/InteractionRates.h"
#include "crpropa/Referenced.h"
#include "crpropa/Vector3.h"
#include "crpropa/Units.h"
#include "crpropa/Random.h"

#include <nanoflann.hpp>

#include <vector>
#include <string>
#include <unordered_map>

namespace crpropa {

InteractionRatesHomogeneous::InteractionRatesHomogeneous(std::string ratesName, bool isPositionDependent) : InteractionRates() {
  this->ratesName = ratesName;
  this->isPositionDependent = isPositionDependent;
    
    // getTables
    // check for table consistency?
    
}

std::vector<double> InteractionRatesHomogeneous::getTabulatedEnergy() const {
    return tabEnergy;
}

std::vector<double> InteractionRatesHomogeneous::getTabulatedRate() const {
    return tabRate;
}

std::vector<double> InteractionRatesHomogeneous::getTabulatedE() const {
    return tabE;
}

std::vector<double> InteractionRatesHomogeneous::getTabulateds() const {
    return tabs;
}

std::vector<std::vector<double>> InteractionRatesHomogeneous::getTabulatedCDF() const {
    return tabCDF;
}

double InteractionRatesHomogeneous::getProcessRate(const double E, const Vector3d &position) const {
    if (!this->isPositionDependent) {
        
        // check if in tabulated energy range
        if ((E < this->tabEnergy.front()) or (E > this->tabEnergy.back())) {
            // throw std::runtime_error("Candidate energy out of tables!");
            
            return -1;
        }
        
        // compute the interaction rate for the given candidate energy, E
        double rate = interpolate(E, this->tabEnergy, this->tabRate);
        
        return rate;
        
    } else {
        throw std::runtime_error("Error in boolean isPositionDependent!");
    }
}

void InteractionRatesHomogeneous::loadPerformInteractionTabs(const Vector3d &position, std::vector<double> &tabE, std::vector<double> &tabs, std::vector<std::vector<double>> &tabCDF) const {
    tabE = this->getTabulatedE();
    tabs = this->getTabulateds();
    tabCDF = this->getTabulatedCDF();
}

void InteractionRatesHomogeneous::setTabulatedEnergy(std::vector<double>& tabEnergy) {
    this->tabEnergy = tabEnergy;
}

void InteractionRatesHomogeneous::setTabulatedRate(std::vector<double>& tabRate) {
    this->tabRate = tabRate;
}

void InteractionRatesHomogeneous::setTabulatedE(std::vector<double>& tabE) {
    this->tabE = tabE;
}

void InteractionRatesHomogeneous::setTabulateds(std::vector<double>& tabs) {
    this->tabs = tabs;
}

void InteractionRatesHomogeneous::setTabulatedCDF(std::vector<std::vector<double>>& tabCDF) {
    this->tabCDF = tabCDF;
}

InteractionRatesPositionDependent::InteractionRatesPositionDependent(std::string ratesName, bool isPositionDependent) : InteractionRates() {
    this->ratesName = ratesName;
    this->isPositionDependent = isPositionDependent;

}

int InteractionRatesPositionDependent::findClosestGridPoint(const Vector3d &position) const {
    
    if (!tree) {
        throw std::runtime_error("KD-Tree not initialized!");
    }
    
    unsigned int closestIndex;
    double closestDistSquared;
    double queryPoint[3] = { position.x, position.y, position.z };
    
    this->tree->knnSearch(queryPoint, 1, &closestIndex, &closestDistSquared);
    
    /**
    std::cout << "iMin: " << cloud.ids[closestIndex] << std::endl;
    
    std::cout << "position: " << position.x / kpc << " " << position.y / kpc << " " << position.z / kpc << std::endl;
    
    std::cout << "node: " << cloud.kdtree_get_pt(cloud.ids[closestIndex], 0) / kpc << " " << cloud.kdtree_get_pt(cloud.ids[closestIndex], 1) / kpc << " " << cloud.kdtree_get_pt(cloud.ids[closestIndex], 2) / kpc << std::endl;
     */
    
    return this->cloud.ids[closestIndex];
}

std::vector<double> InteractionRatesPositionDependent::getTabulatedEnergy() const {
    return tabEnergy;
}

std::vector<std::vector<double>> InteractionRatesPositionDependent::getTabulatedRate() const {
    return tabRate;
}

std::vector<double> InteractionRatesPositionDependent::getTabulatedE() const {
    return tabE;
}

std::vector<std::vector<double>> InteractionRatesPositionDependent::getTabulateds() const {
    return tabs;
}

std::vector<std::vector<std::vector<double>>> InteractionRatesPositionDependent::getTabulatedCDF() const {
    return tabCDF;
}

std::unordered_map<int, Vector3d> InteractionRatesPositionDependent::getPhotonDict() const {
    return photonDict;
}

std::vector<double> InteractionRatesPositionDependent::getClosestRate(const Vector3d &position) const {
    int iMin = findClosestGridPoint(position);
    return tabRate[iMin];
}

std::vector<double> InteractionRatesPositionDependent::getClosests(const Vector3d &position) const {
    int iMin = findClosestGridPoint(position);
    return tabs[iMin];
}

std::vector<std::vector<double>> InteractionRatesPositionDependent::getClosestCDF(const Vector3d &position) const {
    int iMin = findClosestGridPoint(position);
    return tabCDF[iMin];
}

double InteractionRatesPositionDependent::getProcessRate(const double E, const Vector3d &position) const {
    if (!this->isPositionDependent) {
    
        throw std::runtime_error("Error in boolean isPositionDependent!");
        
    } else {
        
        std::vector<double> tabRate = this->getClosestRate(position);
        
        if ((E < this->tabEnergy.front()) or (E > this->tabEnergy.back())) {
        //    throw std::runtime_error("Candidate energy out of tables!");
            
            return -1;
        }
        
        // compute the interaction rate for the given candidate energy, E
        double rate = interpolate(E, this->tabEnergy, tabRate);
        
        return rate;
    }
}

void InteractionRatesPositionDependent::loadPerformInteractionTabs(const Vector3d &position, std::vector<double> &tabE, std::vector<double> &tabs, std::vector<std::vector<double>> &tabCDF) const {
    
    std::vector<double> E = this->getTabulatedE();
    std::vector<double> s = this->getClosests(position);
    std::vector<std::vector<double>> CDF = this->getClosestCDF(position);
    
    tabE = E;
    tabs = s;
    tabCDF = CDF;
}

void InteractionRatesPositionDependent::setTabulatedEnergy(std::vector<double>& tabEnergy) {
    this->tabEnergy = tabEnergy;
}

void InteractionRatesPositionDependent::setTabulatedRate(std::vector<std::vector<double>>& tabRate) {
    this->tabRate = tabRate;
}

void InteractionRatesPositionDependent::setTabulatedE(std::vector<double>& tabE) {
    this->tabE = tabE;
}

void InteractionRatesPositionDependent::setTabulateds(std::vector<std::vector<double>>& tabs) {
    this->tabs = tabs;
}

void InteractionRatesPositionDependent::setTabulatedCDF(std::vector<std::vector<std::vector<double>>>& tabCDF) {
    this->tabCDF = tabCDF;
}

void InteractionRatesPositionDependent::setPhotonDict(std::unordered_map<int, Vector3d>& photonDict) {
    this->photonDict = photonDict;
    
    // delete old clouds
    this->cloud.points.clear();
    this->cloud.ids.clear();
    
    for (const auto& el : this->photonDict) {
        this->cloud.ids.push_back(el.first);
        this->cloud.points.push_back(el.second);
    }
    
    // delete old tree
    if (this->tree) {
        delete this->tree;
    }
    
    this->tree = new KDTree(3, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10));
    this->tree->buildIndex();
    
}

} //namespace crpropa
