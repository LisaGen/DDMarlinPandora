/**
 *  @file   DDMarlinPandora/src/DDPfoCreator.cc
 * 
 *  @brief  Implementation of the pfo creator class.
 * 
 *  $Log: $
 */

#include "DDPfoCreator.h"

#include "marlin/Global.h"
#include "marlin/Processor.h"

#include "EVENT/LCCollection.h"

#include "IMPL/ClusterImpl.h"
#include "IMPL/LCCollectionVec.h"
#include "IMPL/LCFlagImpl.h"
#include "IMPL/LCGenericObjectImpl.h"
#include "IMPL/LCRelationImpl.h"
#include "IMPL/ReconstructedParticleImpl.h"
#include "IMPL/VertexImpl.h"

#include "CalorimeterHitType.h"

#include "Api/PandoraApi.h"

#include "Objects/Cluster.h"
#include "Objects/ParticleFlowObject.h"
#include "Objects/Track.h"

#include "Pandora/PdgTable.h"

#include <algorithm>
#include <cmath>

DDPfoCreator::DDPfoCreator(const Settings &settings, const pandora::Pandora *const pPandora) :
    m_settings(settings),
    m_pandora(*pPandora)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

DDPfoCreator::~DDPfoCreator()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode DDPfoCreator::CreateParticleFlowObjects(EVENT::LCEvent *pLCEvent)
{
  std::cout<<"************************PFO CREATOR**********************"<<std::endl;
    const pandora::PfoList *pPandoraPfoList = NULL;
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::GetCurrentPfoList(m_pandora, pPandoraPfoList));

    IMPL::LCCollectionVec *pClusterCollection = new IMPL::LCCollectionVec(LCIO::CLUSTER);
    IMPL::LCCollectionVec *pReconstructedParticleCollection = new IMPL::LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
    IMPL::LCCollectionVec *pStartVertexCollection = new IMPL::LCCollectionVec(LCIO::VERTEX);

    IMPL::LCFlagImpl lcFlagImpl(pClusterCollection->getFlag());
    lcFlagImpl.setBit(LCIO::CLBIT_HITS);
    pClusterCollection->setFlag(lcFlagImpl.getFlag());

    pandora::StringVector subDetectorNames;
    this->InitialiseSubDetectorNames(subDetectorNames);
    pClusterCollection->parameters().setValues("ClusterSubdetectorNames", subDetectorNames);

    std::cout<<" Particle creator for N PFO= "<<pPandoraPfoList->size()<<std::endl;
    // Create lcio "reconstructed particles" from the pandora "particle flow objects"
   
    for (pandora::PfoList::const_iterator pIter = pPandoraPfoList->begin(), pIterEnd = pPandoraPfoList->end(); pIter != pIterEnd; ++pIter)
    {

        const pandora::ParticleFlowObject *const pPandoraPfo(*pIter);
        IMPL::ReconstructedParticleImpl *const pReconstructedParticle(new ReconstructedParticleImpl());

	std::cout<<"pfo energy="<<pPandoraPfo->GetEnergy()<<std::endl;
        const bool hasTrack(!pPandoraPfo->GetTrackList().empty());
        const pandora::ClusterList &clusterList(pPandoraPfo->GetClusterList());
	std::cout<<" hasTrack: "<< hasTrack<<std::endl;

        float clustersTotalEnergy(0.f);
        pandora::CartesianVector referencePoint(0.f, 0.f, 0.f), clustersWeightedPosition(0.f, 0.f, 0.f);
	std::cout<<"--------Loop on cluster list size="<<clusterList.size()<<std::endl;

        for (pandora::ClusterList::const_iterator cIter = clusterList.begin(), cIterEnd = clusterList.end(); cIter != cIterEnd; ++cIter)
        {
            const pandora::Cluster *const pPandoraCluster(*cIter);
            pandora::CaloHitList pandoraCaloHitList;
            pPandoraCluster->GetOrderedCaloHitList().FillCaloHitList(pandoraCaloHitList);
            
            pandoraCaloHitList.insert(pandoraCaloHitList.end(), pPandoraCluster->GetIsolatedCaloHitList().begin(), pPandoraCluster->GetIsolatedCaloHitList().end());
	    std::cout << "PandoraCAloHitList sixe " << pandoraCaloHitList.size() << std::endl;

	    std::cout<<" call SetClusterSubDetectorEnergies"<<std::endl;
            pandora::FloatVector hitE, hitX, hitY, hitZ;
            IMPL::ClusterImpl *const pLcioCluster(new ClusterImpl());
            this->SetClusterSubDetectorEnergies(subDetectorNames, pLcioCluster, pandoraCaloHitList, hitE, hitX, hitY, hitZ);
            std::cout<<" cluster energy before= "<<pLcioCluster->getEnergy()<<std::endl;
	    std::cout<<" call SetClusterEnergyAndError"<<std::endl;
            float clusterCorrectEnergy(0.f);
	    if(m_settings._digitalCalo==1)
	      {
		std::cout<<" is digital calo"<<std::endl;
		this->SetDigiCalClusterEnergyAndError(pPandoraPfo, pPandoraCluster, pLcioCluster, clusterCorrectEnergy, pandoraCaloHitList);
	      }
	     if(m_settings._digitalCalo==0)
              {
		std::cout<<" is analogic calo"<<std::endl;
		this->SetClusterEnergyAndError(pPandoraPfo, pPandoraCluster, pLcioCluster, clusterCorrectEnergy);
              }
            //this->SetClusterEnergyAndError(pPandoraPfo, pPandoraCluster, pLcioCluster, clusterCorrectEnergy, pandoraCaloHitList);
	     std::cout<<" cluster energy after="<<clusterCorrectEnergy<<" pLcioCluster ene="<<pLcioCluster->getEnergy()<<std::endl;

            pandora::CartesianVector clusterPosition(0.f, 0.f, 0.f);
            const unsigned int nHitsInCluster(pandoraCaloHitList.size());
	    //std::cout<<"before: position x="<<clusterPosition.GetX()<<" y="<<clusterPosition.GetY()<<" z="<<clusterPosition.GetZ()<<std::endl;
            this->SetClusterPositionAndError(nHitsInCluster, hitE, hitX, hitY, hitZ, pLcioCluster, clusterPosition);
	    std::cout<<"after: position x="<<clusterPosition.GetX()<<" y="<<clusterPosition.GetY()<<" z="<<clusterPosition.GetZ()<<std::endl;
            if (!hasTrack)
            {
                clustersWeightedPosition += clusterPosition * clusterCorrectEnergy;
                clustersTotalEnergy += clusterCorrectEnergy;
            }

            pClusterCollection->addElement(pLcioCluster);
            pReconstructedParticle->addCluster(pLcioCluster);
	    std::cout<<" cluster energy ="<<clusterCorrectEnergy<<" pLcioCluster ene="<<pLcioCluster->getEnergy()<<std::endl;
        }
	std::cout<<" ---endl loop on clusters"<<std::endl;

	    if (!hasTrack)
        {
            if (clustersTotalEnergy < std::numeric_limits<float>::epsilon())
            {
                streamlog_out(WARNING) << "DDPfoCreator::CreateParticleFlowObjects: invalid cluster energy " << clustersTotalEnergy << std::endl;
                throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);
            }
            else
            {
                referencePoint = clustersWeightedPosition * (1.f / clustersTotalEnergy);
            }
        }
        else
        {
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CalculateTrackBasedReferencePoint(pPandoraPfo, referencePoint));
        }

        this->SetRecoParticleReferencePoint(referencePoint, pReconstructedParticle);
        this->AddTracksToRecoParticle(pPandoraPfo, pReconstructedParticle);
        this->SetRecoParticlePropertiesFromPFO(pPandoraPfo, pReconstructedParticle);
        pReconstructedParticleCollection->addElement(pReconstructedParticle);
	std::cout<<" particle energy="<<pReconstructedParticle->getEnergy()<<std::endl;

        IMPL::VertexImpl *const pStartVertex(new VertexImpl());
        pStartVertex->setAlgorithmType(m_settings.m_startVertexAlgName.c_str());
        pStartVertex->setPosition(referencePoint.GetX(),referencePoint.GetY(),referencePoint.GetZ());
        pStartVertex->setAssociatedParticle(pReconstructedParticle);
        pStartVertexCollection->addElement(pStartVertex);
    }

    pLCEvent->addCollection(pClusterCollection, m_settings.m_clusterCollectionName.c_str());
    pLCEvent->addCollection(pReconstructedParticleCollection, m_settings.m_pfoCollectionName.c_str());
    pLCEvent->addCollection(pStartVertexCollection, m_settings.m_startVertexCollectionName.c_str());

    return pandora::STATUS_CODE_SUCCESS;
}

//-------------------------------------------------x-----------------------------------------------------------------------------------------

void DDPfoCreator::InitialiseSubDetectorNames(pandora::StringVector &subDetectorNames) const
{
    subDetectorNames.push_back("ecal");
    subDetectorNames.push_back("hcal");
    subDetectorNames.push_back("yoke");
    subDetectorNames.push_back("lcal");
    subDetectorNames.push_back("lhcal");
    subDetectorNames.push_back("bcal");
}

//------------------------------------------------------------------------------------------------------------------------------------------
void DDPfoCreator::SetClusterSubDetectorEnergies(const pandora::StringVector &subDetectorNames, IMPL::ClusterImpl *const pLcioCluster,
    const pandora::CaloHitList &pandoraCaloHitList, pandora::FloatVector &hitE, pandora::FloatVector &hitX, pandora::FloatVector &hitY,
    pandora::FloatVector &hitZ) const
{   bool HCAL_hit = false;
    bool ECAL_hit = false;
    for (pandora::CaloHitList::const_iterator hIter = pandoraCaloHitList.begin(), hIterEnd = pandoraCaloHitList.end(); hIter != hIterEnd; ++hIter)
    {
        const pandora::CaloHit *const pPandoraCaloHit(*hIter);
        EVENT::CalorimeterHit *const pCalorimeterHit = (EVENT::CalorimeterHit*)(pPandoraCaloHit->GetParentAddress());
        pLcioCluster->addHit(pCalorimeterHit, 1.f);

        double R_hit = sqrt(pCalorimeterHit->getPosition()[0]*pCalorimeterHit->getPosition()[0]+pCalorimeterHit->getPosition()[1]*pCalorimeterHit->getPosition()[1]);
        double z_hit = sqrt(pCalorimeterHit->getPosition()[2]*pCalorimeterHit->getPosition()[2]);
        const float caloHitEnergy(pCalorimeterHit->getEnergy());

        hitE.push_back(caloHitEnergy);
        hitX.push_back(pCalorimeterHit->getPosition()[0]);
        hitY.push_back(pCalorimeterHit->getPosition()[1]);
        hitZ.push_back(pCalorimeterHit->getPosition()[2]);
        
	std::cout<<"*********Calo HIt Energy: "<<caloHitEnergy<<"**************"<<std::endl;
    std::cout<<"*********subDetectorNames.size(): "<<subDetectorNames.size()<<"**************"<<std::endl;
        //std::cout<<"*********Position: R= "<<sqrt(pCalorimeterHit->getPosition()[0]*pCalorimeterHit->getPosition()[0]+pCalorimeterHit->getPosition()[1]*pCalorimeterHit->getPosition()[1])<<" z= "<<pCalorimeterHit->getPosition()[2]<<std::endl;
        
        std::vector<float> &subDetectorEnergies = pLcioCluster->subdetectorEnergies();
        std::cout<<"*********pLcioCluster->subdetectorEnergies(): "<<&subDetectorEnergies[0]<<"**************"<<std::endl;
        std::cout<<"*********after &subDetectorEnergies = pLcioCluster->subdetectorEnergies(); *************"<<std::endl;
        subDetectorEnergies.resize(subDetectorNames.size());
	std::cout<<"*********Calo HIt Energy: after resize : "<<subDetectorEnergies.size()<<"**************"<<std::endl;
  
    if((R_hit>= 1500 && R_hit<=1702 && z_hit < 2210.) || (R_hit>= 310 && R_hit<=1700 && z_hit>= 2307. && z_hit<= 2509.)){
        //std::cout<<"ECAL"<<std::endl;
        ECAL_hit = true;
    }
    if((R_hit>= 1740. && R_hit<=3330. && z_hit < 2210.) || (R_hit>= 300 && R_hit<=3246 && z_hit>= 2539. && z_hit<= 4129.)){
        //std::cout<<"HCAL"<<std::endl;
        HCAL_hit = true;
    }
        std::cout << "pCalorimeterHit->getType(): "<< pCalorimeterHit->getType() << std::endl;
        std::cout << "CHT(pCalorimeterHit->getType()).caloID()" << CHT(pCalorimeterHit->getType()).caloID()<< std::endl;
        switch (CHT(pCalorimeterHit->getType()).caloID())
        {
            case CHT::ecal:  subDetectorEnergies[ECAL_INDEX ] += caloHitEnergy; break;
            case CHT::hcal:  subDetectorEnergies[HCAL_INDEX ] += caloHitEnergy; break;
            case CHT::yoke:  subDetectorEnergies[YOKE_INDEX ] += caloHitEnergy; break;
            case CHT::lcal:  subDetectorEnergies[LCAL_INDEX ] += caloHitEnergy; break;
            case CHT::lhcal: subDetectorEnergies[LHCAL_INDEX] += caloHitEnergy; break;
            case CHT::bcal:  subDetectorEnergies[BCAL_INDEX ] += caloHitEnergy; break;
            default: streamlog_out(WARNING) << "DDPfoCreator::SetClusterSubDetectorEnergies: no subdetector found for hit with type: " << pCalorimeterHit->getType() << std::endl;
        }
    
    }
    if(HCAL_hit && ECAL_hit ) std::cout<<"HE~cluster~"<<std::endl;
    else if  (HCAL_hit && !ECAL_hit)std::cout<<"Honly~cluster~"<<std::endl;
    else if  (!HCAL_hit && ECAL_hit) std::cout<<"Eonly~cluster~"<<std::endl;
    else std::cout<<"UK~cluster~"<<std::endl;

    
}
//------------------------------------------------------------------------------------------------------------------------------------------


void DDPfoCreator::SetClusterEnergyAndError(const pandora::ParticleFlowObject *const pPandoraPfo, const pandora::Cluster *const pPandoraCluster, 
    IMPL::ClusterImpl *const pLcioCluster, float &clusterCorrectEnergy) const
{
    const bool isEmShower((pandora::PHOTON == pPandoraPfo->GetParticleId()) || (pandora::E_MINUS == std::abs(pPandoraPfo->GetParticleId())));
    clusterCorrectEnergy = (isEmShower ? pPandoraCluster->GetCorrectedElectromagneticEnergy(m_pandora) : pPandoraCluster->GetCorrectedHadronicEnergy(m_pandora));

    if (clusterCorrectEnergy < std::numeric_limits<float>::epsilon())
        throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);

    const float stochasticTerm(isEmShower ? m_settings.m_emStochasticTerm : m_settings.m_hadStochasticTerm); 
    const float constantTerm(isEmShower ? m_settings.m_emConstantTerm : m_settings.m_hadConstantTerm);
    const float energyError(std::sqrt(stochasticTerm * stochasticTerm / clusterCorrectEnergy + constantTerm * constantTerm) * clusterCorrectEnergy);

    pLcioCluster->setEnergy(clusterCorrectEnergy);
    pLcioCluster->setEnergyError(energyError);
    std::cout<<" clusterCorrectEnergy="<<clusterCorrectEnergy<<std::endl;
}

/*
void DDPfoCreator::SetClusterEnergyAndError(const pandora::ParticleFlowObject *const pPandoraPfo, const pandora::Cluster *const pPandoraCluster,
                                            IMPL::ClusterImpl *const pLcioCluster, float &clusterCorrectEnergy, const pandora::CaloHitList &pandoraCaloHitList) const
{
    const bool isEmShower((pandora::PHOTON == pPandoraPfo->GetParticleId()) || (pandora::E_MINUS == std::abs(pPandoraPfo->GetParticleId())));
    //    clusterCorrectEnergy = (isEmShower ? pPandoraCluster->GetCorrectedElectromagneticEnergy(m_pandora) : pPandoraCluster->GetCorrectedHadronicEnergy(m_pandora));
    std::cout<<"isEMShower="<<isEmShower<<std::endl;
    std::cout<<" pandora calo Hit list="<<pandoraCaloHitList.size()<<std::endl;
    float a=27.46;
    float b=1/0.8843;
    float c=-10.81/a;
    clusterCorrectEnergy=std::pow((pandoraCaloHitList.size()/a),b) ;   
    if (clusterCorrectEnergy < std::numeric_limits<float>::epsilon())
        throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);

    const float stochasticTerm(isEmShower ? m_settings.m_emStochasticTerm : m_settings.m_hadStochasticTerm);
    const float constantTerm(isEmShower ? m_settings.m_emConstantTerm : m_settings.m_hadConstantTerm);
    const float energyError(std::sqrt(stochasticTerm * stochasticTerm / clusterCorrectEnergy + constantTerm * constantTerm) * clusterCorrectEnergy);

    pLcioCluster->setEnergy(clusterCorrectEnergy);
    pLcioCluster->setEnergyError(energyError);
 
    //    std::cout<<" clusterCorrectEnergy="<<clusterCorrectEnergy<<" my digi calo energy="<<myClEn<<std::endl;
    std::cout<<" clusterCorrectEnergy="<<clusterCorrectEnergy<<std::endl;
    
}
*/

void DDPfoCreator::SetDigiCalClusterEnergyAndError(const pandora::ParticleFlowObject *const pPandoraPfo, const pandora::Cluster *const pPandoraCluster,
                                            IMPL::ClusterImpl *const pLcioCluster, float &clusterCorrectEnergy, const pandora::CaloHitList &pandoraCaloHitList) const
{
    const bool isEmShower((pandora::PHOTON == pPandoraPfo->GetParticleId()) || (pandora::E_MINUS == std::abs(pPandoraPfo->GetParticleId())));
    //    clusterCorrectEnergy = (isEmShower ? pPandoraCluster->GetCorrectedElectromagneticEnergy(m_pandora) : pPandoraCluster->GetCorrectedHadronicEnergy(m_pandora));     
    //std::cout<<"isEMShower="<<isEmShower<<std::endl;
    //std::cout<<" pandora calo Hit list size="<<pandoraCaloHitList.size()<<std::endl;

    //for Digi RO 
    float a_digi=27.46;
    float b_digi=1/0.8843;
    float c_digi=-10.81/a_digi;

    //for SemiDigi RO 
    int N1 = 0;
    int N2 = 0;
    int N3 = 0;
    int N; 
    int N_ECAL = 0;
    float a;
    float b;
    float c;

    // MPGD1x1_CRILIN (0.2-6-12) noES sample
    /*float ap0 =	0.0791714;
    float ap1 = 0.;
    float ap2 =  0.;

    float bp0 = 0.162359;
    float bp1 = 0.00011747;
    float bp2 = 0.; 

    float cp0 =  0.191908;
    float cp1 = 0.00115602; 
    float cp2 = -5.12013e-07;*/

    /*// MPGD1x1_CRILIN (0.2-2-12) noES sample
    float ap0 =	0.0791714;
    float ap1 = -5.15947e-06;
    float ap2 =  4.90034e-08;

    float bp0 = 0.0261397;
    float bp1 = 0.000203803;
    float bp2 = -5.40581e-07; 

    float cp0 = 0.29502;
    float cp1 = 0.000774884; 
    float cp2 = 6.51271e-07;*/

   // MPGD1x1_CRILIN (0.2-3.9-11.9) noES sample
    /* float ap0 =	0.071391;
    float ap1 = -3.73869e-05;
    float ap2 = 4.14876e-08;

    float bp0 = 0.0373246;
    float bp1 = 0.000707441;
    float bp2 = -1.15117e-06; 

    float cp0 = 0.245835;
    float cp1 = 0.000821071; 
    float cp2 = 1.04599e-07;*/

    //MPGD1x1_CRILIN (0.2-1-3) noES sample
    /*float ap0 = 0.0644425;
    float ap1 = 0.000455192;
    float ap2 = -7.29697e-07;
    float bp0 = 0.0542693;
    float bp1 =  -0.000458562;
    float bp2 = 3.86825e-07;
    float cp0 = 0.174251;
    float cp1 = -5.8699e-05;
    float cp2 = 8.82761e-07;*/

    // MPGD1x1_CRILIN (0.2-3) noES sample
    /*float ap0 = 0.0852425;
    float ap1 = 0.000161769;
    float ap2 = -4.64141e-07;
    
    float cp0 = 0.0696009;
    float cp1 = -6.08462e-05;
    float cp2 = 3.76417e-07;*/
    
    /*// Anna's thesis estimates (0.2-3.9-11.9)
    float ap0 = 4.13118e-02;
    float ap1 = 1.85454e-05;
    float ap2 = -1.33300e-08;
    float bp0 = 9.50907e-02;
    float bp1 = 6.50793e-06;
    float bp2 = 5.38691e-08;
    float cp0 = 8.92423e-02;
    float cp1 = 1.79733e-04;
    float cp2 = 1.70854e-07;
    bool not_all_HCAL = false;*/



    float total_ECAL_ene = 0;
    for (pandora::CaloHitList::const_iterator hIter = pandoraCaloHitList.begin(), hIterEnd = pandoraCaloHitList.end(); hIter != hIterEnd; ++hIter){
        
        const pandora::CaloHit *const pPandoraCaloHit(*hIter);
        EVENT::CalorimeterHit *const pCalorimeterHit = (EVENT::CalorimeterHit*)(pPandoraCaloHit->GetParentAddress());
        const float caloHitEnergy(pCalorimeterHit->getEnergy());
        if (CHT(pCalorimeterHit->getType()).caloID() == 2){ // hits from HCAL
            if(caloHitEnergy == 10) N1++; 
            if(caloHitEnergy == 20) N2++;
            if(caloHitEnergy == 30) N3++;
            //std::cout<<"N1: " << N1 << std::endl;
            //std::cout<<"N2: " << N2 << std::endl;
            //std::cout<<"N3: " << N3 << std::endl;
            
        }
        else if (CHT(pCalorimeterHit->getType()).caloID() == 1){ //hits form ECAL
            total_ECAL_ene += caloHitEnergy;
            std::cout<<"total_ECAL_ene: " << total_ECAL_ene << std::endl;
            N_ECAL += 1;
        }
        //std::cout<<"alpha:"<<ap0 + ap1*(N1+N2+N3) + ap2*(N1+N2+N3)*(N1+N2+N3) <<std::endl;
        //std::cout<<"beta:"<<bp0 + bp1*(N1+N2+N3) + bp2*(N1+N2+N3)*(N1+N2+N3) <<std::endl;
        //std::cout<<"gamma:"<<cp0 + cp1*(N1+N2+N3) + cp2*(N1+N2+N3)*(N1+N2+N3) <<std::endl;
    }
    
    std::cout<<"total_ECAL_ene: " << total_ECAL_ene << std::endl;
    std::cout<<"n_ECAL_hit: " << N_ECAL << std::endl;
    std::cout<<"n_HCAL_hit: " << N1+N2+N3 << std::endl;
    std::cout<<"n_Clu_hit="<<pandoraCaloHitList.size()<<std::endl;

    //THREE THRESHOLDS
    /*N = N1 + N2 + N3;
    a = ap0 + ap1*N + ap2*N*N;
    b = bp0 + bp1*N + bp2*N*N;
    c = cp0 + cp1*N + cp2*N*N;

    clusterCorrectEnergy= a*N1 + b*N2 + c*N3 + total_ECAL_ene;*/

    //TWO THERSHOLDS
    /*N = N1 + N2 + N3;
    a = ap0 + ap1*N + ap2*N*N;
    c = cp0 + cp1*N + cp2*N*N;

    clusterCorrectEnergy= a*(N1+N2) + c*N3 + total_ECAL_ene;*/

    //DIGITAL
    N = N1 + N2 + N3;

    clusterCorrectEnergy= std::pow((N/a_digi),b_digi) + total_ECAL_ene;

    std::cout << "clusterCorrectEnergy: " << clusterCorrectEnergy << std::endl;
    //std::cout << "std::numeric_limits<float>::epsilon()" <<std::numeric_limits<float>::epsilon()<< std::endl;
    if (clusterCorrectEnergy < std::numeric_limits<float>::epsilon())
        throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);
        std::cout<<"numeric_limits EXCEPTION"<<std::endl;
        
    const float stochasticTerm(isEmShower ? m_settings.m_emStochasticTerm : m_settings.m_hadStochasticTerm);
    const float constantTerm(isEmShower ? m_settings.m_emConstantTerm : m_settings.m_hadConstantTerm);
    const float energyError(std::sqrt(stochasticTerm * stochasticTerm / clusterCorrectEnergy + constantTerm * constantTerm) * clusterCorrectEnergy);

    pLcioCluster->setEnergy(clusterCorrectEnergy);
    pLcioCluster->setEnergyError(energyError);

    //    std::cout<<" clusterCorrectEnergy="<<clusterCorrectEnergy<<" my digi calo energy="<<myClEn<<std::endl;                                                            
    std::cout<<" clusterCorrectEnergy="<<clusterCorrectEnergy<<std::endl;

}



//------------------------------------------------------------------------------------------------------------------------------------------

void DDPfoCreator::SetClusterPositionAndError(const unsigned int nHitsInCluster, pandora::FloatVector &hitE, pandora::FloatVector &hitX, 
    pandora::FloatVector &hitY, pandora::FloatVector &hitZ, IMPL::ClusterImpl *const pLcioCluster, pandora::CartesianVector &clusterPositionVec) const
{
    ClusterShapes *const pClusterShapes(new ClusterShapes(nHitsInCluster, hitE.data(), hitX.data(), hitY.data(), hitZ.data()));

    try
    {
        pLcioCluster->setIPhi(std::atan2(pClusterShapes->getEigenVecInertia()[1], pClusterShapes->getEigenVecInertia()[0]));
        pLcioCluster->setITheta(std::acos(pClusterShapes->getEigenVecInertia()[2]));
        pLcioCluster->setPosition(pClusterShapes->getCentreOfGravity());
        //ATTN these two lines below would only compile with ilcsoft HEAD V2015-10-13 and above
        //pLcioCluster->setPositionError(pClusterShapes->getCenterOfGravityErrors());
        //pLcioCluster->setDirectionError(pClusterShapes->getEigenVecInertiaErrors());
        clusterPositionVec.SetValues(pClusterShapes->getCentreOfGravity()[0], pClusterShapes->getCentreOfGravity()[1], pClusterShapes->getCentreOfGravity()[2]);
    }
    catch (...)
    {
        streamlog_out(WARNING) << "DDPfoCreator::SetClusterPositionAndError: unidentified exception caught." << std::endl;
    }

    delete pClusterShapes;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode DDPfoCreator::CalculateTrackBasedReferencePoint(const pandora::ParticleFlowObject *const pPandoraPfo, pandora::CartesianVector &referencePoint) const
{
    const pandora::TrackList &trackList(pPandoraPfo->GetTrackList());

    float totalTrackMomentumAtDca(0.f), totalTrackMomentumAtStart(0.f);
    pandora::CartesianVector referencePointAtDCAWeighted(0.f, 0.f, 0.f), referencePointAtStartWeighted(0.f, 0.f, 0.f);

    bool hasSiblings(false);
    for (pandora::TrackList::const_iterator tIter = trackList.begin(), tIterEnd = trackList.end(); tIter != tIterEnd; ++tIter)
    {
        const pandora::Track *const pPandoraTrack(*tIter);

        if (!this->IsValidParentTrack(pPandoraTrack, trackList))
            continue;

        if (this->HasValidSiblingTrack(pPandoraTrack, trackList))
        {
            // Presence of sibling tracks typically represents a conversion
            const pandora::CartesianVector &trackStartPoint((pPandoraTrack->GetTrackStateAtStart()).GetPosition());
            const float trackStartMomentum(((pPandoraTrack->GetTrackStateAtStart()).GetMomentum()).GetMagnitude());
            referencePointAtStartWeighted += trackStartPoint * trackStartMomentum;
            totalTrackMomentumAtStart += trackStartMomentum;
            hasSiblings = true;
        }
        else
        {
            const EVENT::Track *const pLcioTrack = (EVENT::Track*)(pPandoraTrack->GetParentAddress());
            const float z0(pPandoraTrack->GetZ0());
            pandora::CartesianVector intersectionPoint(0.f, 0.f, 0.f);

            intersectionPoint.SetValues(pLcioTrack->getD0() * std::cos(pLcioTrack->getPhi()), pLcioTrack->getD0() * std::sin(pLcioTrack->getPhi()), z0);
            const float trackMomentumAtDca((pPandoraTrack->GetMomentumAtDca()).GetMagnitude());
            referencePointAtDCAWeighted += intersectionPoint * trackMomentumAtDca;
            totalTrackMomentumAtDca += trackMomentumAtDca;
        }
    }

    if (hasSiblings)
    {
        if (totalTrackMomentumAtStart < std::numeric_limits<float>::epsilon())
        {
            streamlog_out(WARNING) << "DDPfoCreator::CalculateTrackBasedReferencePoint: invalid track momentum " << totalTrackMomentumAtStart << std::endl;
            throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);
        }
        else
        {
            referencePoint = referencePointAtStartWeighted * (1.f / totalTrackMomentumAtStart);
        }
    }
    else
    {
        if (totalTrackMomentumAtDca < std::numeric_limits<float>::epsilon())
        {
            streamlog_out(WARNING) << "DDPfoCreator::CalculateTrackBasedReferencePoint: invalid track momentum " << totalTrackMomentumAtDca << std::endl;
            throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);
        }
        else
        {
            referencePoint = referencePointAtDCAWeighted * (1.f / totalTrackMomentumAtDca);
        }
    }

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DDPfoCreator::IsValidParentTrack(const pandora::Track *const pPandoraTrack, const pandora::TrackList &allTrackList) const
{
    const pandora::TrackList &parentTrackList(pPandoraTrack->GetParentList());

    for (pandora::TrackList::const_iterator iter = parentTrackList.begin(), iterEnd = parentTrackList.end(); iter != iterEnd; ++iter)
    {
        if (allTrackList.end() != std::find(allTrackList.begin(), allTrackList.end(), *iter))
            continue;

        // ATTN This track must have a parent not in the all track list; still use it if it is the closest to the ip
        streamlog_out(WARNING) << "DDPfoCreator::IsValidParentTrack: mismatch in track relationship information, use information as available " << std::endl;

        if (this->IsClosestTrackToIP(pPandoraTrack, allTrackList))
            return true;

        return false;
    }

    // Ideal case: All parents are associated to same pfo
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DDPfoCreator::HasValidSiblingTrack(const pandora::Track *const pPandoraTrack, const pandora::TrackList &allTrackList) const
{
    const pandora::TrackList &siblingTrackList(pPandoraTrack->GetSiblingList());

    for (pandora::TrackList::const_iterator iter = siblingTrackList.begin(), iterEnd = siblingTrackList.end(); iter != iterEnd; ++iter)
    {
        if (allTrackList.end() != std::find(allTrackList.begin(), allTrackList.end(), *iter))
            continue;

        // ATTN This track must have a sibling not in the all track list; still use it if it has a second sibling that is in the list
        streamlog_out(WARNING) << "DDPfoCreator::HasValidSiblingTrack: mismatch in track relationship information, use information as available " << std::endl;

        if (this->AreAnyOtherSiblingsInList(pPandoraTrack, allTrackList))
            return true;

        return false;
    }

    // Ideal case: All siblings associated to same pfo
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DDPfoCreator::IsClosestTrackToIP(const pandora::Track *const pPandoraTrack, const pandora::TrackList &allTrackList) const 
{
    const pandora::Track *pClosestTrack(NULL);
    float closestTrackDisplacement(std::numeric_limits<float>::max()); 

    for (pandora::TrackList::const_iterator iter = allTrackList.begin(), iterEnd = allTrackList.end(); iter != iterEnd; ++iter)
    {
        const pandora::Track *const pTrack(*iter);
        const float trialTrackDisplacement(pTrack->GetTrackStateAtStart().GetPosition().GetMagnitude());

        if (trialTrackDisplacement < closestTrackDisplacement)
        {
            closestTrackDisplacement = trialTrackDisplacement;
            pClosestTrack = pTrack;
        }
    }

    return (pPandoraTrack == pClosestTrack);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DDPfoCreator::AreAnyOtherSiblingsInList(const pandora::Track *const pPandoraTrack, const pandora::TrackList &allTrackList) const
{
    const pandora::TrackList &siblingTrackList(pPandoraTrack->GetSiblingList());

    for (pandora::TrackList::const_iterator iter = siblingTrackList.begin(), iterEnd = siblingTrackList.end(); iter != iterEnd; ++iter)
    {
        if (allTrackList.end() != std::find(allTrackList.begin(), allTrackList.end(), *iter))
            return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DDPfoCreator::SetRecoParticleReferencePoint(const pandora::CartesianVector &referencePoint, IMPL::ReconstructedParticleImpl *const pReconstructedParticle) const
{
    const float referencePointArray[3] = {referencePoint.GetX(), referencePoint.GetY(), referencePoint.GetZ()};
    pReconstructedParticle->setReferencePoint(referencePointArray);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DDPfoCreator::AddTracksToRecoParticle(const pandora::ParticleFlowObject *const pPandoraPfo, IMPL::ReconstructedParticleImpl *const pReconstructedParticle) const
{
    const pandora::TrackList &trackList(pPandoraPfo->GetTrackList());

    for (pandora::TrackList::const_iterator tIter = trackList.begin(), tIterEnd = trackList.end(); tIter != tIterEnd; ++tIter)
    {
        const pandora::Track *const pTrack(*tIter);
        pReconstructedParticle->addTrack((EVENT::Track*)(pTrack->GetParentAddress()));
        std::cout<<"pTrack->GetMomentum()"<<pTrack->GetMomentumAtDca() << std::endl;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DDPfoCreator::SetRecoParticlePropertiesFromPFO(const pandora::ParticleFlowObject *const pPandoraPfo, IMPL::ReconstructedParticleImpl *const pReconstructedParticle) const
{
    const float momentum[3] = {pPandoraPfo->GetMomentum().GetX(), pPandoraPfo->GetMomentum().GetY(), pPandoraPfo->GetMomentum().GetZ()};
    pReconstructedParticle->setMomentum(momentum);
    pReconstructedParticle->setEnergy(pPandoraPfo->GetEnergy());
    pReconstructedParticle->setMass(pPandoraPfo->GetMass());
    pReconstructedParticle->setCharge(pPandoraPfo->GetCharge());
    pReconstructedParticle->setType(pPandoraPfo->GetParticleId());
    std::cout<<" energy from set reco part prop="<<pPandoraPfo->GetEnergy()<<std::endl;
    std::cout<<" px="<<pPandoraPfo->GetMomentum().GetX()<<" py="<<pPandoraPfo->GetMomentum().GetY()<<" pz="<<pPandoraPfo->GetMomentum().GetZ()<<std::endl;
    std::cout<<" pt="<<std::pow((pPandoraPfo->GetMomentum().GetX()*pPandoraPfo->GetMomentum().GetX() + pPandoraPfo->GetMomentum().GetY()*pPandoraPfo->GetMomentum().GetY()),0.5);
}

//------------------------------------------------------------------------------------------------------------------------------------------

DDPfoCreator::Settings::Settings():
    m_clusterCollectionName(""),
    m_pfoCollectionName(""),
    m_startVertexCollectionName(""),
    m_startVertexAlgName(""),
    m_emStochasticTerm(0.17f),
    m_hadStochasticTerm(0.6f),
    m_emConstantTerm(0.01f),
    m_hadConstantTerm(0.03f),
    _digitalCalo(1)
{
}
