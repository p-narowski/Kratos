//  DPD-like simplified DEM constitutive law for small particles

#include "custom_constitutive/dem_dpd_sph_like_cl.h"
#include "DEM_application_variables.h"
#include "custom_elements/spheric_particle.h"
#include <cmath>

namespace Kratos {

DEMDiscontinuumConstitutiveLaw::Pointer DEM_DPD_SPH_LIKE::Clone() const {
    DEMDiscontinuumConstitutiveLaw::Pointer p_clone(new DEM_DPD_SPH_LIKE(*this));
    return p_clone;
}

std::unique_ptr<DEMDiscontinuumConstitutiveLaw> DEM_DPD_SPH_LIKE::CloneUnique() {
    return Kratos::make_unique<DEM_DPD_SPH_LIKE>(*this);
}

std::string DEM_DPD_SPH_LIKE::GetTypeOfLaw() {
    return "DPD_SPH_LIKE";
}

void DEM_DPD_SPH_LIKE::Check(Properties::Pointer pProp) const {
    DEM_D_Linear_viscous_Coulomb::Check(pProp);

    KRATOS_ERROR_IF_NOT(pProp->Has(DPD_SMOOTHING_LENGTH))
        << "Missing DPD_SMOOTHING_LENGTH in properties " << pProp->Id() << std::endl;

    KRATOS_ERROR_IF_NOT(pProp->Has(DPD_CUTOFF_RADIUS))
        << "Missing DPD_CUTOFF_RADIUS in properties " << pProp->Id() << std::endl;

    KRATOS_ERROR_IF_NOT(pProp->Has(DPD_CONSERVATIVE_COEFF))
        << "Missing DPD_CONSERVATIVE_COEFF in properties " << pProp->Id() << std::endl;

    KRATOS_ERROR_IF_NOT(pProp->Has(DPD_DISSIPATIVE_COEFF))
        << "Missing DPD_DISSIPATIVE_COEFF in properties " << pProp->Id() << std::endl;
}

void DEM_DPD_SPH_LIKE::InitializeContact(SphericParticle* const element1,
                                         SphericParticle* const element2,
                                         const double indentation) {
    mKn = 0.0;
    mKt = 0.0;
}

void DEM_DPD_SPH_LIKE::InitializeContactWithFEM(SphericParticle* const element,
                                                Condition* const wall,
                                                const double indentation,
                                                const double ini_delta) {
    mKn = 0.0;
    mKt = 0.0;
}

double DEM_DPD_SPH_LIKE::KernelWeight(const double r, const double h, const double rc) const {
    if (r >= rc) return 0.0;

    const double q = r / h;

    if (q >= 0.0 && q < 1.0) {
        return 1.0 - 1.5*q*q + 0.75*q*q*q;
    } else if (q >= 1.0 && q < 2.0) {
        const double a = 2.0 - q;
        return 0.25 * a * a * a;
    }

    return 0.0;
}

double DEM_DPD_SPH_LIKE::DissipativeWeight(const double r, const double h, const double rc) const {
    const double w = KernelWeight(r, h, rc);
    return w * w;
}

double DEM_DPD_SPH_LIKE::CalculateNormalForce(SphericParticle* const element1,
                                              SphericParticle* const element2,
                                              const double indentation,
                                              double LocalCoordSystem[3][3]) {
    return 0.0;
}

double DEM_DPD_SPH_LIKE::CalculateNormalForce(SphericParticle* const element,
                                              Condition* const wall,
                                              const double indentation) {
    return 0.0;
}

double DEM_DPD_SPH_LIKE::CalculateNormalForce(const double indentation) {
    return 0.0;
}

double DEM_DPD_SPH_LIKE::CalculateCohesiveNormalForce(SphericParticle* const element1,
                                                      SphericParticle* const element2,
                                                      const double indentation) {
    return 0.0;
}

double DEM_DPD_SPH_LIKE::CalculateCohesiveNormalForceWithFEM(SphericParticle* const element,
                                                             Condition* const wall,
                                                             const double indentation) {
    return 0.0;
}

void DEM_DPD_SPH_LIKE::CalculateForces(const ProcessInfo& r_process_info,
                                       const double OldLocalElasticContactForce[3],
                                       double LocalElasticContactForce[3],
                                       double LocalDeltDisp[3],
                                       double LocalRelVel[3],
                                       double indentation,
                                       double previous_indentation,
                                       double ViscoDampingLocalContactForce[3],
                                       double& cohesive_force,
                                       SphericParticle* element1,
                                       SphericParticle* element2,
                                       bool& sliding,
                                       double LocalCoordSystem[3][3]) {
    KRATOS_TRY

    LocalElasticContactForce[0] = 0.0;
    LocalElasticContactForce[1] = 0.0;
    LocalElasticContactForce[2] = 0.0;

    ViscoDampingLocalContactForce[0] = 0.0;
    ViscoDampingLocalContactForce[1] = 0.0;
    ViscoDampingLocalContactForce[2] = 0.0;

    cohesive_force = 0.0;
    sliding = false;

    Properties& properties_of_this_contact =
        element1->GetProperties().GetSubProperties(element2->GetProperties().Id());

    const double h     = properties_of_this_contact[DPD_SMOOTHING_LENGTH];
    const double rc    = properties_of_this_contact[DPD_CUTOFF_RADIUS];
    const double a     = properties_of_this_contact[DPD_CONSERVATIVE_COEFF];
    const double gamma = properties_of_this_contact[DPD_DISSIPATIVE_COEFF];

    const auto& x1 = element1->GetGeometry()[0].Coordinates();
    const auto& x2 = element2->GetGeometry()[0].Coordinates();

    const double dx = x1[0] - x2[0];
    const double dy = x1[1] - x2[1];
    const double dz = x1[2] - x2[2];
    const double r2 = dx*dx + dy*dy + dz*dz;
    const double r  = std::sqrt(r2);

    if (r < 1.0e-15 || r >= rc) {
        return;
    }

    const double wc = KernelWeight(r, h, rc);
    const double wd = DissipativeWeight(r, h, rc);

    // In your branch, local normal direction is component 2
    const double vrel_n = LocalRelVel[2];

    const double f_cons = a * wc;
    const double f_diss = -gamma * wd * vrel_n;

    LocalElasticContactForce[2]      = f_cons;
    ViscoDampingLocalContactForce[2] = f_diss;

    KRATOS_CATCH("")
}

void DEM_DPD_SPH_LIKE::CalculateForcesWithFEM(const ProcessInfo& r_process_info,
                                              const double OldLocalElasticContactForce[3],
                                              double LocalElasticContactForce[3],
                                              double LocalDeltDisp[3],
                                              double LocalRelVel[3],
                                              double indentation,
                                              double previous_indentation,
                                              double ViscoDampingLocalContactForce[3],
                                              double& cohesive_force,
                                              SphericParticle* const element,
                                              Condition* const wall,
                                              bool& sliding) {
    LocalElasticContactForce[0] = 0.0;
    LocalElasticContactForce[1] = 0.0;
    LocalElasticContactForce[2] = 0.0;

    ViscoDampingLocalContactForce[0] = 0.0;
    ViscoDampingLocalContactForce[1] = 0.0;
    ViscoDampingLocalContactForce[2] = 0.0;

    cohesive_force = 0.0;
    sliding = false;
}

} // namespace Kratos