// DPD-like simplified DEM constitutive law for small particles
// Based on DEM_D_Linear_viscous_Coulomb interface in this Kratos branch

#if !defined(DEM_DPD_SPH_LIKE_CL_H_INCLUDED)
#define DEM_DPD_SPH_LIKE_CL_H_INCLUDED

#include "includes/define.h"
#include "includes/serializer.h"
#include "custom_constitutive/DEM_D_Linear_viscous_Coulomb_CL.h"

namespace Kratos
{

    class SphericParticle;

    class KRATOS_API(DEM_APPLICATION) DEM_DPD_SPH_LIKE : public DEM_D_Linear_viscous_Coulomb
    {
    public:
        using DEMDiscontinuumConstitutiveLaw::CalculateNormalForce;

        KRATOS_CLASS_POINTER_DEFINITION(DEM_DPD_SPH_LIKE);

        DEM_DPD_SPH_LIKE() = default;
        ~DEM_DPD_SPH_LIKE() override = default;

        std::string GetTypeOfLaw() override;

        bool IsRangeForce() const override
        {
            return true;
        }

        void Check(Properties::Pointer pProp) const override;

        DEMDiscontinuumConstitutiveLaw::Pointer Clone() const override;
        std::unique_ptr<DEMDiscontinuumConstitutiveLaw> CloneUnique() override;

        void InitializeContact(SphericParticle *const element1,
                               SphericParticle *const element2,
                               const double indentation) override;

        void InitializeContactWithFEM(SphericParticle *const element,
                                      Condition *const wall,
                                      const double indentation,
                                      const double ini_delta = 0.0) override;

        void CalculateForces(const ProcessInfo &r_process_info,
                             const double OldLocalElasticContactForce[3],
                             double LocalElasticContactForce[3],
                             double LocalDeltDisp[3],
                             double LocalRelVel[3],
                             double indentation,
                             double previous_indentation,
                             double ViscoDampingLocalContactForce[3],
                             double &cohesive_force,
                             SphericParticle *element1,
                             SphericParticle *element2,
                             bool &sliding,
                             double LocalCoordSystem[3][3]) override;

        void CalculateForcesWithFEM(const ProcessInfo &r_process_info,
                                    const double OldLocalElasticContactForce[3],
                                    double LocalElasticContactForce[3],
                                    double LocalDeltDisp[3],
                                    double LocalRelVel[3],
                                    double indentation,
                                    double previous_indentation,
                                    double ViscoDampingLocalContactForce[3],
                                    double &cohesive_force,
                                    SphericParticle *const element,
                                    Condition *const wall,
                                    bool &sliding) override;

        double CalculateNormalForce(SphericParticle *const element1,
                                    SphericParticle *const element2,
                                    const double indentation,
                                    double LocalCoordSystem[3][3]) override;

        double CalculateNormalForce(SphericParticle *const element,
                                    Condition *const wall,
                                    const double indentation) override;

        double CalculateNormalForce(const double indentation) override;

        double CalculateCohesiveNormalForce(SphericParticle *const element1,
                                            SphericParticle *const element2,
                                            const double indentation) override;

        double CalculateCohesiveNormalForceWithFEM(SphericParticle *const element,
                                                   Condition *const wall,
                                                   const double indentation) override;

    private:
        double KernelWeight(const double r, const double h, const double rc) const;
        double DissipativeWeight(const double r, const double h, const double rc) const;

        friend class Serializer;

        void save(Serializer &rSerializer) const override
        {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, DEM_D_Linear_viscous_Coulomb);
        }

        void load(Serializer &rSerializer) override
        {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, DEM_D_Linear_viscous_Coulomb);
        }
    };

} // namespace Kratos

#endif // DEM_DPD_SPH_LIKE_CL_H_INCLUDED
