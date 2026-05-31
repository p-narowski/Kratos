//
// Authors:
// Miguel Angel Celigueta maceli@cimne.upc.edu
// Salvador Latorre latorre@cimne.upc.edu
// Miquel Santasusana msantasusana@cimne.upc.edu
// Guillermo Casas gcasas@cimne.upc.edu
// Chengshun Shang cshang@cimne.upc.edu
//
// Distribution or reproduction of any portion of this source code or
// documentation is not allowed without the explicit written consent of
// the CIMNE Structural Mechanics group.

// System includes
#include <string>
#include <iostream>
#include <cmath>
#include <algorithm>

// External includes

// Kratos includes
#include "includes/define.h"
#include "includes/kratos_flags.h"
#include "spheric_particle.h"
#include "custom_utilities/GeometryFunctions.h"
#include "custom_utilities/AuxiliaryFunctions.h"
#include "DEM_application_variables.h"
#include "custom_elements/Particle_Contact_Element.h"
#include "custom_constitutive/DEM_discontinuum_constitutive_law.h"
#include "utilities/quaternion.h"
#include "includes/variables.h"

namespace Kratos
{
// using namespace GeometryFunctions;

SphericParticle::SphericParticle()
    : DiscreteElement(), mRealMass(0), mRadius(0.0), mSearchControlVariable(0) {}

SphericParticle::SphericParticle(IndexType NewId, GeometryType::Pointer pGeometry)
    : DiscreteElement(NewId, pGeometry), mRealMass(0), mRadius(0.0), mSearchControlVariable(0) {}

SphericParticle::SphericParticle(IndexType NewId, NodesArrayType const& ThisNodes)
    : DiscreteElement(NewId, ThisNodes), mRealMass(0), mRadius(0.0), mSearchControlVariable(0) {}

SphericParticle::SphericParticle(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : DiscreteElement(NewId, pGeometry, pProperties), mRealMass(0), mRadius(0.0), mSearchControlVariable(0) {}

SphericParticle::SphericParticle(Element const& rOther, PropertiesType::Pointer pProperties)
    : DiscreteElement(rOther, pProperties), mRealMass(0), mRadius(0.0), mSearchControlVariable(0) {}

SphericParticle::~SphericParticle() {}

SphericParticle& SphericParticle::operator=(SphericParticle const& rOther) {
    DiscreteElement::operator=(rOther);
    return *this;
}

Element::Pointer SphericParticle::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const {
    return Element::Pointer(new SphericParticle(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

void SphericParticle::Initialize(const ProcessInfo& r_process_info) {
    SetValue(NEIGHBOUR_IDS, DenseVector<int>());
}

void SphericParticle::MemberDeclarationFirstStep(const ProcessInfo& r_process_info) {
    if (!r_process_info[IS_RESTARTED]) {
        mDiscontinuumConstitutiveLaw = DEMDiscontinuumConstitutiveLaw::Pointer(new DEMDiscontinuumConstitutiveLaw());
        if (GetProperties().Has(DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME)) {
            std::string constitutive_law_name = GetProperties()[DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME];
            Kratos::DEMDiscontinuumConstitutiveLaw::Pointer p_constitutive_law = KratosComponents<Kratos::DEMDiscontinuumConstitutiveLaw>::Get(constitutive_law_name).Clone();
            mDiscontinuumConstitutiveLaw = p_constitutive_law;
        }
    }
}

void SphericParticle::SetRadius(double radius) {
    mRadius = radius;
    if (GetGeometry().size() > 0) { // typically 1 node
        GetGeometry()[0].FastGetSolutionStepValue(RADIUS) = radius;
    }
}

void SphericParticle::SetRadius() {
    mRadius = GetGeometry()[0].FastGetSolutionStepValue(RADIUS);
}

double SphericParticle::GetRadius() {
    return mRadius;
}

double SphericParticle::CalculateVolume() {
    return (4.0 / 3.0) * Globals::Pi * mRadius * mRadius * mRadius;
}

double SphericParticle::GetInteractionRadius(const int radius_index) {
    return mRadius;
}

double SphericParticle::GetSearchRadius() {
    return this->GetValue(SEARCH_RADIUS);
}

void SphericParticle::SetSearchRadius(const double r) {
    this->SetValue(SEARCH_RADIUS, r);
}

void SphericParticle::Initialize(const ProcessInfo& r_process_info, const bool Flag) {}

void SphericParticle::CalculateRightHandSide(const ProcessInfo& r_process_info,
                                               double dt,
                                               const array_1d<double,3>& gravity) {
    DEM_COPY_SECOND_TO_FIRST_3(this->GetGeometry()[0].FastGetSolutionStepValue(TOTAL_FORCES), ZeroVector(3))
    DEM_COPY_SECOND_TO_FIRST_3(this->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_MOMENT), ZeroVector(3))
}

void SphericParticle::FirstCalculateRightHandSide(const ProcessInfo& r_process_info, double dt) {}
void SphericParticle::CollectCalculateRightHandSide(const ProcessInfo& r_process_info) {}
void SphericParticle::FinalCalculateRightHandSide(const ProcessInfo& r_process_info, double dt, const array_1d<double,3>& gravity) {}
void SphericParticle::InitializeForceComputation(const ProcessInfo& r_process_info) {}
void SphericParticle::FinalizeForceComputation(ParticleDataBuffer& rDataBuffer) {}

void SphericParticle::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& r_process_info) const {}
void SphericParticle::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& r_process_info) {}

void SphericParticle::CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& r_process_info) {
    rMassMatrix.resize(1, 1, false);
    rMassMatrix(0, 0) = GetMass();
}

void SphericParticle::GetDofList(DofsVectorType& ElementalDofList, const ProcessInfo& r_process_info) const {}

std::string SphericParticle::Info() const {
    std::stringstream buffer;
    buffer << "SphericParticle #" << Id();
    return buffer.str();
}

void SphericParticle::PrintInfo(std::ostream& rOStream) const { rOStream << "SphericParticle #" << Id(); }
void SphericParticle::PrintData(std::ostream& rOStream) const {}

void SphericParticle::ComputeNewNeighboursHistoricalData(DenseVector<int>& temp_neighbours_ids,
                                                          std::vector<array_1d<double, 3>>& temp_neighbour_elastic_contact_forces) {
    std::vector<array_1d<double, 3>> temp_old_forces(mNeighbourElements.size());
    DenseVector<int> temp_old_ids = GetValue(NEIGHBOUR_IDS);

    size_t old_size = mNeighbourElasticContactForces.size();
    for (size_t i = 0; i < mNeighbourElements.size(); ++i) {
        if (!mNeighbourElements[i]) {
            noalias(temp_old_forces[i]) = ZeroVector(3);
            continue;
        }
        int id = mNeighbourElements[i]->Id();
        bool found = false;
        for (unsigned int j = 0; j < temp_old_ids.size(); ++j) {
            if (temp_old_ids[j] == id && j < old_size) {
                noalias(temp_old_forces[i]) = mNeighbourElasticContactForces[j];
                found = true;
                break;
            }
        }
        if (!found) noalias(temp_old_forces[i]) = ZeroVector(3);
    }

    temp_neighbours_ids.resize(mNeighbourElements.size());
    for (size_t i = 0; i < mNeighbourElements.size(); ++i) {
        temp_neighbours_ids[i] = mNeighbourElements[i] ? mNeighbourElements[i]->Id() : -1;
    }
    temp_neighbour_elastic_contact_forces = temp_old_forces;
}

void SphericParticle::ComputeNewRigidFaceNeighboursHistoricalData() {
    array_1d<double, 3> zero = ZeroVector(3);
    size_t size = mNeighbourRigidFaces.size();
    if (mNeighbourRigidFacesTotalContactForce.size() != size) {
        mNeighbourRigidFacesTotalContactForce.resize(size);
        mNeighbourRigidFacesElasticContactForce.resize(size);
        for (size_t i = 0; i < size; ++i) {
            noalias(mNeighbourRigidFacesTotalContactForce[i]) = zero;
            noalias(mNeighbourRigidFacesElasticContactForce[i]) = zero;
        }
    }
}

void SphericParticle::ContactCalculations::RecordNewFEMContacts() {}
void SphericParticle::ContactCalculations::FinalizeForceComputation(ParticleDataBuffer& rDataBuffer) {}

ParticleDataBuffer* SphericParticle::CreateParticleDataBuffer(SphericParticle* p_this_particle) {
    return new ParticleDataBuffer(p_this_particle);
}

void SphericParticle::InitializeDataBuffer(ParticleDataBuffer& rDataBuffer) {
    rDataBuffer.mpOtherParticle         = nullptr;
    rDataBuffer.mMultiStageRHS          = false;
    rDataBuffer.mDeltTime               = 0.0;
    rDataBuffer.mDistance               = 0.0;
    rDataBuffer.mRadiusSum              = 0.0;
    rDataBuffer.mOtherRadius            = 0.0;
    rDataBuffer.mIndentation            = 0.0;
    noalias(rDataBuffer.mOtherToMeVector) = ZeroVector(3);
    noalias(rDataBuffer.mLocalCoordSystem[0]) = ZeroVector(3);
    noalias(rDataBuffer.mLocalCoordSystem[1]) = ZeroVector(3);
    noalias(rDataBuffer.mLocalCoordSystem[2]) = ZeroVector(3);
    noalias(rDataBuffer.mOldLocalCoordSystem[0]) = ZeroVector(3);
    noalias(rDataBuffer.mOldLocalCoordSystem[1]) = ZeroVector(3);
    noalias(rDataBuffer.mOldLocalCoordSystem[2]) = ZeroVector(3);
}

void SphericParticle::ComputeBallToBallContactForceAndMoment(ParticleDataBuffer& data_buffer,
                                                               const ProcessInfo& r_process_info,
                                                               array_1d<double, 3>& rElasticForce,
                                                               array_1d<double, 3>& rContactForce) {
    // Base implementation intentionally empty — derived classes override
}

void SphericParticle::ComputeBallToRigidFaceContactForceAndMoment(array_1d<double, 3>& rElasticForce,
                                                                    array_1d<double, 3>& rContactForce,
                                                                    array_1d<double, 3>& rigid_element_force,
                                                                    const ProcessInfo& r_process_info,
                                                                    double mTimeStep,
                                                                    bool has_mpi) {}

void SphericParticle::ComputeConditionBasedContactForces(array_1d<double, 3>& force) {}

void SphericParticle::InitializeSolutionStep(const ProcessInfo& r_process_info) {}
void SphericParticle::FinalizeSolutionStep(const ProcessInfo& r_process_info) {}

double SphericParticle::GetMass() { return mRealMass; }
void SphericParticle::SetMass(double real_mass) { mRealMass = real_mass; }

double SphericParticle::CalculateLocalMaxPeriod(const bool has_mpi, const ProcessInfo& r_process_info) { return 0.0; }

void SphericParticle::Move(const double delta_t, const bool rotation_option, const double force_reduction_factor, const int StepFlag) {}

void SphericParticle::SetIntegrationScheme(DEMIntegrationScheme::Pointer& translational_integration_scheme,
                                             DEMIntegrationScheme::Pointer& rotational_integration_scheme) {}

void SphericParticle::UpdatePositionAndOrientation(const ProcessInfo& r_process_info) {}
void SphericParticle::PrepareForPrinting(const ProcessInfo& r_process_info) {}
void SphericParticle::MemberDeclarationSecondStep(const ProcessInfo& r_process_info) {}

double SphericParticle::GetYoung() { return GetProperties()[YOUNG_MODULUS]; }
double SphericParticle::GetPoisson() { return GetProperties()[POISSON_RATIO]; }
double SphericParticle::GetDensity() { return GetProperties()[PARTICLE_DENSITY]; }
int SphericParticle::GetParticleMaterial() { return GetProperties()[PARTICLE_MATERIAL]; }
double SphericParticle::GetParticleRealYoungModulusRatio() { return 1.0; }

void SphericParticle::Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& r_process_info) {}
void SphericParticle::Calculate(const Variable<array_1d<double,3>>& rVariable, array_1d<double,3>& Output, const ProcessInfo& r_process_info) {}
void SphericParticle::Calculate(const Variable<Vector>& rVariable, Vector& Output, const ProcessInfo& r_process_info) {}
void SphericParticle::Calculate(const Variable<Matrix>& rVariable, Matrix& Output, const ProcessInfo& r_process_info) {}

void SphericParticle::AdditionalCalculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& r_process_info) {}

int SphericParticle::GetClusterId() { return GetValue(CLUSTER_ID); }
void SphericParticle::SetClusterId(int givenId) { SetValue(CLUSTER_ID, givenId); }

// The actual implementation lives in the .cpp built from spheric_particle.cpp
// (the full Kratos DEM source). The stubs above satisfy the linker for tests.
// When building the full Kratos DEM application the real implementation below
// is compiled instead.

} // namespace Kratos
